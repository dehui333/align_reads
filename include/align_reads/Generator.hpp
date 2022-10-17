#ifndef ALIGN_READS_GENERATOR_HPP_
#define ALIGN_READS_GENERATOR_HPP_

#include <iostream>

#include "align_reads/Aligner.hpp"
#include "align_reads/AlignCounter.hpp"
#include "align_reads/CountsConverter.hpp"
#include "align_reads/Converter.hpp"
#include "align_reads/IndexedSequences.hpp"
#include "align_reads/Inputs.hpp"
#include "align_reads/Overlapper.hpp"
#include "align_reads/Utilities.hpp"

/*
 * Issues:
 * 1. Add haplotype label to AlignmentSegment & target sequence
 *
 * **> need to keep on producing features and keeping
 *  an 'inventory' even when python level is doing work with a batch
 *
 *
 * !!!!!!!!!!!!!!!!! try to move width index stuff to MultiAlignment and add
 *  multi alignment aware printing functionality there
 *
 * >> The ground truth info can either be passed from python level
 * (actually can't, i don't know which sequence is coming up and passing everything will be too large)
 * or found in C++ level using ram....?
 * possible ways:
 * 1. find on the spot with ram
 * 2. put the truth reads in memory....
 * 3. get the truth reads from indexed fasta
 * 4. get from bam somehow (actually if we don't load whole bam into mem, as good
 *  as indexing fasta and getting from there?)
 */

/*
 * The 'Main' component that interacts with other components.
 */

namespace align_reads
{

    class Generator
    {
    public:
        // The arguments are vector of paths to the sequences // obsolete
        Generator(std::vector<std::string> &reads0, std::vector<std::string> &reads1,
                  std::vector<std::string> &haplotype0, std::vector<std::string> &haplotype1,
                  std::shared_ptr<thread_pool::ThreadPool> pool, std::uint16_t window_size);
        Generator(std::vector<std::string> &reads, std::shared_ptr<thread_pool::ThreadPool> pool, std::uint16_t window_length, bool debug_printing = false);
        // Data produce_data();
        //  The first vector<PyObject*> contains the truth matrices,
        //  the second the ..., etc.
        std::vector<std::vector<PyObject *>> produce_data();

        void index_truth(std::string &truth_path);

        ~Generator() {delete for_print;}

        void print_window(std::uint32_t width_idx, std::uint32_t len)
        {
            for_print->print_in_window(width_idx, len);
        }

    private:
        bool has_truth = false;
        std::uint32_t start_of_reads1 = 0; // id from which second set of reads start from
        std::uint32_t current_target = 0;  // may need atomic if parallelize?

        std::shared_ptr<thread_pool::ThreadPool> pool; // pool of workers
        Inputs inputs;                                 // Stores input sequences
        Overlapper overlapper;
        std::uint16_t window_length;

        IndexedSequences indexed_sequences;
        bool debug_printing;
        MultiAlignment *for_print = nullptr;

        inline void prepare_for_print(std::string &target_string, std::vector<clipped_alignment<EdlibAlignResult>> &align_results)
        {

            auto futures = Futures<AlignmentSegment>(pool, align_results.size());
            for (auto &result : align_results)
            {
                futures.add_inputs(get_alignment_segment2, result, target_string, false);
            }
            std::vector<AlignmentSegment> segments = futures.get();
            if (for_print != nullptr)
            {
                delete for_print;
            }
            for_print = new MultiAlignment(target_string, std::move(segments));
        }

        inline AlignCounter get_align_counter(std::unique_ptr<biosoup::NucleicAcid> &target, std::string &target_string)
        {
            auto overlaps = overlapper.find_overlaps(target, READS_GROUP);
            auto align_results = align_overlaps(overlaps, overlaps.size(), inputs, target_string, pool, READS_GROUP);
            if (debug_printing)
            {
                prepare_for_print(target_string, align_results);
            }
            return {target_string, align_results};
        }

        inline std::vector<PyObject *> get_ground_truth(std::unique_ptr<biosoup::NucleicAcid> &target, std::string &target_string, AlignCounter &align_counter, std::uint32_t &left_clip, std::uint32_t &right_clip, std::uint32_t &alignment_length, std::uint32_t &num_matrices)
        {
            // output vector
            std::vector<PyObject *> output;

            // get truth sequence and align to target
            std::string seq_name = target->name;
            std::string truth_string = indexed_sequences.get_sequence(seq_name);
            auto align_segment = get_alignment_segment(truth_string, 0, truth_string.size(),
                                                       target_string, 0, target_string.size(), EDLIB_MODE_GLOBAL, EDLIB_TASK_PATH);

            // find the length on the sides without truth alignment
            left_clip = align_segment.start_on_target;
            right_clip = target_string.size() - align_segment.end_on_target - 1;

            // calculate length with alignment and number of matrices
            alignment_length = align_counter.alignment_length - left_clip - right_clip;
            num_matrices = alignment_length % window_length == 0 ? alignment_length / window_length : alignment_length / window_length + 1;

            // prepare matrices
            npy_intp dims[2];
            dims[0] = 1;
            dims[1] = window_length;
            output.reserve(num_matrices);
            for (std::uint32_t i = 0; i < num_matrices; i++)
            {
                output.push_back(PyArray_SimpleNew(2, dims, NPY_UINT8));
            }

            // fill up matrices
            std::uint32_t matrix_idx = 0;
            std::uint32_t col_idx = 0;

            std::uint32_t target_idx = align_segment.start_on_target;
            std::uint16_t ins_idx = 0;

            auto &counts = align_counter.counts;
            auto end = counts.end() - right_clip;
            // The iteration below is just to achieve same number of positions
            // as the counts matrices. Don't actually need the counts.

            std::uint8_t *value_ptr;
            // each position on the target
            for (auto it = counts.begin() + left_clip; it < end; it++)
            {
                auto &counts_at_pos = *it;
                // aligned base (or del) plus any ins
                for (auto &counter : counts_at_pos)
                {
                    value_ptr = (uint8_t *)PyArray_GETPTR2(output[matrix_idx], 0, col_idx++);
                    *value_ptr = ENCODER[static_cast<std::uint8_t>(align_segment.get_at_target_pos(target_idx, ins_idx))];

                    if (col_idx == window_length)
                    {
                        matrix_idx++;
                        col_idx = 0;
                    }
                    ins_idx++;
                }
                target_idx++;
                ins_idx = 0;
            }

            for (; col_idx < window_length;)
            {
                value_ptr = (uint8_t *)PyArray_GETPTR2(output[matrix_idx], 0, col_idx++);
                *value_ptr = 5; // pad
            }

            return output;
        }
    };
} // namespace align_reads

#endif // ALIGN_READS_GENERATOR_HPP_
