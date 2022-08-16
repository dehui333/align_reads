#include "align_reads/Aligner.hpp"
#include "align_reads/Converter.hpp"
#include "align_reads/Generator.hpp"
#include "align_reads/MultiAlignment.hpp"
#include "align_reads/Utilities.hpp"

namespace align_reads
{

    Generator::Generator(std::vector<std::string> &reads0, std::vector<std::string> &reads1,
                         std::vector<std::string> &haplotype0, std::vector<std::string> &haplotype1,
                         std::shared_ptr<thread_pool::ThreadPool> pool)
        : has_haplotypes(true), pool(pool), inputs(3), overlapper(3, pool)
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

    Generator::Generator(std::vector<std::string> &reads, std::shared_ptr<thread_pool::ThreadPool> pool)
        : has_haplotypes(false), pool(pool), inputs(1), overlapper(1, pool)
    {
        inputs.append_to_group(READS_GROUP, reads, pool);
        inputs.index_group(READS_GROUP);
        overlapper.index_sequences(inputs.get_group(READS_GROUP), READS_GROUP);
    }

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

} // namespace align_reads