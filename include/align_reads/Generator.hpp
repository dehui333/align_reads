#ifndef ALIGN_READS_GENERATOR_HPP_
#define ALIGN_READS_GENERATOR_HPP_

#include "align_reads/Converter.hpp"
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
                  std::shared_ptr<thread_pool::ThreadPool> pool);
        Generator(std::vector<std::string> &reads, std::shared_ptr<thread_pool::ThreadPool> pool);

        Data produce_data();


    private:
        bool has_haplotypes = false; // is there haplotype info?
        std::uint32_t start_of_reads1 = 0; // id from which second set of reads start from
        std::uint32_t current_target = 0; // may need atomic if parallelize?

        std::shared_ptr<thread_pool::ThreadPool> pool; // pool of workers
        Inputs inputs; // Stores input sequences
        Overlapper overlapper;
        

        
    };

} // namespace align_reads

#endif // ALIGN_READS_GENERATOR_HPP_
