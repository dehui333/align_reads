#ifndef ALIGN_READS_INPUTS_HPP_
#define ALIGN_READS_INPUTS_HPP_
#include <memory> // unique_ptr, shared_ptr

#include "biosoup/nucleic_acid.hpp"
#include "thread_pool/thread_pool.hpp"

#define READS_GROUP 0
#define HAP0_GROUP 1
#define HAP1_GROUP 2

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
        // split into this number of groups
        Inputs(std::uint8_t num_groups);

        // add sequences from paths to a group
        void append_to_group(std::uint32_t group_id, std::vector<std::string> &paths,
                             std::shared_ptr<thread_pool::ThreadPool> &pool);

        // return the sequences in a group. They may not be ordered
        std::vector<std::unique_ptr<biosoup::NucleicAcid>> &get_group(std::uint32_t group_id);

        // index a group
        void index_group(std::uint32_t group_id);

        // return reference to the sequence with the specified id in the specified group
        // requires group to be indexed.
        // If not indexed, the sequence_id-th sequence in the group will be returned, 
        // which will be correct if they are ordered by sequence id.
        std::unique_ptr<biosoup::NucleicAcid> &get_id_in_group(std::uint32_t group_id,
                                                               std::uint32_t sequence_id);
    private:
        std::vector<std::vector<std::unique_ptr<biosoup::NucleicAcid>>> groups_of_sequences;
        std::vector<std::vector<std::uint32_t>> indices_of_groups;
    };

} // namespace align_reads
#endif // ALIGN_READS_INPUTS_HPP_
