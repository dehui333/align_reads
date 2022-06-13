#include "align_reads/Generator.hpp"

namespace align_reads
{

    Generator::Generator(std::vector<std::string> &reads0, std::vector<std::string> &reads1,
                         std::vector<std::string> &haplotype0, std::vector<std::string> &haplotype1,
                         std::shared_ptr<thread_pool::ThreadPool> &pool)
        : has_haplotypes(true), pool(pool), inputs(3), overlapper(3, pool)
    {
        inputs.append_to_group(READS_INDEX, reads0, pool);
        start_of_reads1 = inputs.get_group(0).size();
        inputs.append_to_group(READS_INDEX, reads1, pool);
        inputs.append_to_group(HAP0_INDEX, haplotype0, pool);
        inputs.append_to_group(HAP1_INDEX, haplotype1, pool);
        inputs.index_group(READS_INDEX);
        inputs.index_group(HAP0_INDEX);
        inputs.index_group(HAP1_INDEX);
        overlapper.index_sequences(inputs.get_group(READS_INDEX), READS_INDEX);
        overlapper.index_sequences(inputs.get_group(HAP0_INDEX), HAP0_INDEX);
        overlapper.index_sequences(inputs.get_group(HAP1_INDEX), HAP1_INDEX);
    }

    Generator::Generator(std::vector<std::string> &reads, std::shared_ptr<thread_pool::ThreadPool> &pool)
        : has_haplotypes(false), pool(pool), inputs(1), overlapper(1, pool)
    {
        inputs.append_to_group(READS_INDEX, reads, pool);
        inputs.index_group(READS_INDEX);
        overlapper.index_sequences(inputs.get_group(READS_INDEX), READS_INDEX);
    }

} // namespace align_reads