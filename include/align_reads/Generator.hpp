#ifndef ALIGN_READS_GENERATOR_HPP_
#define ALIGN_READS_GENERATOR_HPP_

#include "align_reads/Inputs.hpp"
#include "align_reads/Overlapper.hpp"

#define READS_INDEX 0
#define HAP0_INDEX 1
#define HAP1_INDEX 2

/*
 * Issues:
 * 1. Need to ignore ins after the last target char as it can be of arbitary length for diff 
 *    aligned sequences. -> add a test to check this
 */


/*
 * The 'Main' component that interacts with other components. 
 */

namespace align_reads
{
    class Generator
    {
    public:
        Generator(std::vector<std::string> &reads0, std::vector<std::string> &reads1,
                  std::vector<std::string> &haplotype0, std::vector<std::string> &haplotype1,
                  std::shared_ptr<thread_pool::ThreadPool> &pool);
        Generator(std::vector<std::string> &reads, std::shared_ptr<thread_pool::ThreadPool> &pool);

    private:
        bool has_haplotypes = false; // is there haplotype info?
        std::uint32_t start_of_reads1 = 0; // id from which second set of reads start from

        std::shared_ptr<thread_pool::ThreadPool> pool = nullptr; // pool of workers
        Inputs inputs; // Stores input sequences
        Overlapper overlapper;
        

        
    };

} // namespace align_reads

#endif // ALIGN_READS_GENERATOR_HPP_
