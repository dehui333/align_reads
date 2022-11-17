#include "align_reads/Aligner.hpp"
#include "align_reads/AlignCounter.hpp"
#include "align_reads/CountsConverter.hpp"
#include "align_reads/Converter.hpp"
#include "align_reads/Generator.hpp"
#include "align_reads/MultiAlignment.hpp"
#include "align_reads/Utilities.hpp"

#include <iostream>
#include <chrono>

namespace align_reads
{
    std::mutex mtx;
    Generator::Generator(std::vector<std::string> &reads0, std::vector<std::string> &reads1,
                         std::vector<std::string> &haplotype0, std::vector<std::string> &haplotype1,
                         std::shared_ptr<thread_pool::ThreadPool> pool, std::uint16_t window_length)
        : has_truth(true), pool(pool), inputs(3), overlapper(3, pool), window_length(window_length)
    {
        current_target = 0;
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
        // 5808
        current_target = 0;
        inputs.append_to_group(READS_GROUP, reads, pool);
        inputs.index_group(READS_GROUP);
        overlapper.index_sequences(inputs.get_group(READS_GROUP), READS_GROUP);
    }

    void Generator::index_truth(std::string &truth_path)
    {
        has_truth = true;
        indexed_sequences.index_reads(truth_path);
    }

    bool Generator::has_next()
    {
        return current_target < inputs.get_group(READS_GROUP).size();
    }

    inline std::vector<PyObject*> get_target_info(std::string& target_name, std::uint32_t num_matrices) 
    {
        int int_num_matrices = static_cast<int>(num_matrices);
        mtx.lock();
        auto tuple = Py_BuildValue("si", target_name.c_str(), int_num_matrices);
        mtx.unlock();
        return {tuple};
    }

    std::vector<std::vector<PyObject *>> Generator::produce_data()
    {
        //auto p1 = std::chrono::steady_clock::now();
        std::uint32_t left_clip = 0;
        std::uint32_t right_clip = 0;
        std::vector<std::vector<PyObject *>> output;
        auto &target = inputs.get_id_in_group(READS_GROUP, current_target++);
        //auto &target = inputs.get_id_in_group(READS_GROUP, 3319);
        std::string target_string = target->InflateData();
        //auto p2 = std::chrono::steady_clock::now();
        AlignCounter align_counter = get_align_counter(target, target_string);
        std::uint32_t alignment_length = align_counter.alignment_length;
        std::uint32_t num_matrices = alignment_length % window_length == 0 ? alignment_length/window_length : alignment_length/window_length + 1;
        //std::cout << 1 << std::endl;
        //auto p3 = std::chrono::steady_clock::now();
        /*
        if (has_truth)
        {
            output.push_back(get_ground_truth(target, target_string, align_counter, left_clip, right_clip, alignment_length, num_matrices));
        }*/
        //auto p4 = std::chrono::steady_clock::now();
        //std::cout << 2 << std::endl;
        output.push_back(CountsConverter::get_counts_matrices(align_counter, window_length, left_clip, right_clip, num_matrices, alignment_length));
        //auto p5 = std::chrono::steady_clock::now();
        //std::cout << 3 << std::endl;
        output.push_back(CountsConverter::get_stats_matrices(align_counter, window_length, left_clip, right_clip, num_matrices, alignment_length));
        //auto p6 = std::chrono::steady_clock::now();
        //std::cout << 4 << std::endl;
        output.push_back(get_target_info(target->name, num_matrices));
        //std::cout << 5 << std::endl;
        //auto p7 = std::chrono::steady_clock::now();
        //std::cout << " 2- 1 " << std::chrono::duration_cast<std::chrono::microseconds>(p2 - p1).count() << std::endl;
        //std::cout << " 3- 2 " << std::chrono::duration_cast<std::chrono::microseconds>(p3 - p2).count() << std::endl;
        //std::cout << " 4- 3 " << std::chrono::duration_cast<std::chrono::microseconds>(p4 - p3).count() << std::endl;
        //std::cout << " 5- 4 " << std::chrono::duration_cast<std::chrono::microseconds>(p5 - p4).count() << std::endl;
        //std::cout << " 6- 5 " << std::chrono::duration_cast<std::chrono::microseconds>(p6 - p5).count() << std::endl;
        //std::cout << " 7- 6 " << std::chrono::duration_cast<std::chrono::microseconds>(p7 - p6).count() << std::endl;
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