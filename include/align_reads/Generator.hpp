#ifndef ALIGN_READS_GENERATOR_HPP_
#define ALIGN_READS_GENERATOR_HPP_


#include "align_reads/Aligner.hpp"
#include "align_reads/AlignCounter.hpp"
#include "align_reads/CountsConverter.hpp"
#include "align_reads/Converter.hpp"
#include "align_reads/IndexedSequences.hpp"
#include "align_reads/Inputs.hpp"
#include "align_reads/Overlapper.hpp"

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
        // The arguments are vector of paths to the sequences
        Generator(std::vector<std::string> &reads0, std::vector<std::string> &reads1,
                  std::vector<std::string> &haplotype0, std::vector<std::string> &haplotype1,
                  std::shared_ptr<thread_pool::ThreadPool> pool, std::uint16_t window_size);
        Generator(std::vector<std::string> &reads, std::shared_ptr<thread_pool::ThreadPool> pool, std::uint16_t window_length);
        // Data produce_data();
        //  The first vector<PyObject*> contains the truth matrices,
        //  the second the ..., etc.
        std::vector<std::vector<PyObject *>> produce_data();

        void index_truth(std::string &truth_path);

    private:
        bool has_truth = false;
        std::uint32_t start_of_reads1 = 0; // id from which second set of reads start from
        std::uint32_t current_target = 0;  // may need atomic if parallelize?

        std::shared_ptr<thread_pool::ThreadPool> pool; // pool of workers
        Inputs inputs;                                 // Stores input sequences
        Overlapper overlapper;
        std::uint16_t window_length;

        IndexedSequences indexed_sequences;

        inline AlignCounter get_align_counter(std::unique_ptr<biosoup::NucleicAcid> &target, std::string &target_string)
        {
            auto overlaps = overlapper.find_overlaps(target, READS_GROUP);
            auto align_results = align_overlaps(overlaps, overlaps.size(), inputs, target_string, pool, READS_GROUP);
            return {target_string, align_results};
        }

        inline std::vector<PyObject *> get_ground_truth(std::unique_ptr<biosoup::NucleicAcid> &target, std::string &target_string, AlignCounter &align_counter, std::uint32_t& left_clip, std::uint32_t& right_clip, std::uint32_t& alignment_length, std::uint32_t& num_matrices)
        {
            std::vector<PyObject *> output;
            
            npy_intp dims[2];
            dims[0] = 1;
            dims[1] = window_length;
            std::string seq_name = target->name;
            std::string truth_string = indexed_sequences.get_sequence(seq_name);
            auto align_segment = get_alignment_segment(truth_string, 0, truth_string.size(),
                                                       target_string, 0, target_string.size(), EDLIB_MODE_GLOBAL, EDLIB_TASK_PATH);
            left_clip = align_segment.start_on_target;
            right_clip = target_string.size() - align_segment.end_on_target - 1;
            
            alignment_length = align_counter.alignment_length - left_clip - right_clip;
            num_matrices = alignment_length % window_length == 0 ? alignment_length/window_length : alignment_length/window_length + 1;
            output.reserve(num_matrices);
            for (std::uint32_t i = 0; i < num_matrices; i++)
            {
                output.push_back(PyArray_SimpleNew(2, dims, NPY_UINT16));
            }
            // need to fill up + pass to python level
            return output;
        }
    };
} // namespace align_reads

#endif // ALIGN_READS_GENERATOR_HPP_
