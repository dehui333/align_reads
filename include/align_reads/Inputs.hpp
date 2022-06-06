#ifndef ALIGN_READS_INPUTS_HPP_
#define ALIGN_READS_INPUTS_HPP_
#include <memory> // unique_ptr, shared_ptr

#include "biosoup/nucleic_acid.hpp"
#include "thread_pool/thread_pool.hpp"

/*
 *  Manages reading in of and access to input sequences.
 *  Sequences are separated into distinct groups.
 *  
 */

namespace align_reads
{

    class Inputs
    {
        public:
            
        private:
            //std::vector<std::vector<std::unique_ptr<biosoup::NucleicAcid>>> groups_of_sequences;

       
    };

} // namespace align_reads
#endif // ALIGN_READS_INPUTS_HPP_
