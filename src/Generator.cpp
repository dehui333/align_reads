#include "align_reads/Aligner.hpp"
#include "align_reads/AlignCounter.hpp"
#include "align_reads/CountsConverter.hpp"
#include "align_reads/Converter.hpp"
#include "align_reads/Generator.hpp"
#include "align_reads/MultiAlignment.hpp"
#include "align_reads/Utilities.hpp"

#include <iostream>
namespace align_reads
{

    Generator::Generator(std::vector<std::string> &reads0, std::vector<std::string> &reads1,
                         std::vector<std::string> &haplotype0, std::vector<std::string> &haplotype1,
                         std::shared_ptr<thread_pool::ThreadPool> pool, std::uint16_t window_length)
        : has_truth(true), pool(pool), inputs(3), overlapper(3, pool), window_length(window_length)
    {
        inputs.append_to_group(READS_GROUP, reads0, pool);
        start_of_reads1 = inputs.get_group(0).size();
        inputs.append_to_group(READS_GROUP, reads1, pool);
        inputs.append_to_group(HAP0_GROUP, haplotype0, pool);
        inputs.append_to_group(HAP1_GROUP, haplotype1, pool);
        inputs.index_group(READS_GROUP);
        inputs.index_group(HAP0_GROUP);
        inputs.index_group(HAP1_GROUP);
        overlapper.index_sequences(inputs.get_group(READS_GROUP), READS_GROUP);
        overlapper.index_sequences(inputs.get_group(HAP0_GROUP), HAP0_GROUP);
        overlapper.index_sequences(inputs.get_group(HAP1_GROUP), HAP1_GROUP);
    }

    Generator::Generator(std::vector<std::string> &reads, std::shared_ptr<thread_pool::ThreadPool> pool, std::uint16_t window_length, bool debug_printing)
        : has_truth(false), pool(pool), inputs(1), overlapper(1, pool), window_length(window_length), debug_printing(debug_printing) 
    {
        inputs.append_to_group(READS_GROUP, reads, pool);
        inputs.index_group(READS_GROUP);
        overlapper.index_sequences(inputs.get_group(READS_GROUP), READS_GROUP);
    }

    void Generator::index_truth(std::string &truth_path)
    {
        has_truth = true;
        indexed_sequences.index_reads(truth_path);
    }

    std::vector<std::vector<PyObject *>> Generator::produce_data()
    {
        std::uint32_t left_clip = 0;
        std::uint32_t right_clip = 0;
        std::vector<std::vector<PyObject *>> output;
        auto &target = inputs.get_id_in_group(READS_GROUP, current_target++);
        std::string target_string = target->InflateData();
        AlignCounter align_counter = get_align_counter(target, target_string);
        std::uint32_t alignment_length = align_counter.alignment_length;
        std::uint32_t num_matrices = alignment_length % window_length == 0 ? alignment_length/window_length : alignment_length/window_length + 1;
        if (has_truth)
        {
            output.push_back(get_ground_truth(target, target_string, align_counter, left_clip, right_clip, alignment_length, num_matrices));
        }
        output.push_back(CountsConverter::get_counts_matrices(align_counter, window_length, left_clip, right_clip, num_matrices, alignment_length));
        return output;
    }

    /*
        Data Generator::produce_data()
        {
            // Get the next target read
            auto& target = inputs.get_id_in_group(READS_GROUP, current_target++);
            std::string target_string = target->InflateData();

            // Find its overlaps with other reads
            auto overlaps = overlapper.find_overlaps(target, READS_GROUP);
            // currently just align all
            // TODO:
            // 1. throw away badly aligned results
            // 2. align in batches, keep track distribution across windows, stop when enough
            // 3. align in smaller windows?

            // align the target with overlapping
            auto align_results = align_overlaps(overlaps, overlaps.size(), inputs, target_string, pool, READS_GROUP);
            // transform overlapping to targets align results into AlignmentSegment
            auto futures = Futures<AlignmentSegment>(pool, align_results.size());
            for (auto& result : align_results)
            {
                futures.add_inputs(get_alignment_segment2, result, target_string);
            }
            std::vector<AlignmentSegment> segments = futures.get();
            // put into MultiAlignment
            // currently no truth alignments used
            // TODO: add truth info if exist
            MultiAlignment multi_align {std::move(target_string), std::move(segments)};
            // Give MultiAlignment to AlignmentConverter
            // currently no additional info on aligned reads given
            // TODO: 1. add option to find them and input
            //       2. Add option for changing height and width of matrix
            AlignmentConverter converter {multi_align, 5, 10};
            // produce features
            //
            Data data = converter.produce_data(pool, false);
            return data;
        }
    */
} // namespace align_reads