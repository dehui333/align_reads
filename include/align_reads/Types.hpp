#include <string>
#include <vector>

#include "edlib.h"

#ifndef ALIGN_READS_TYPES_HPP_
#define ALIGN_READS_TYPES_HPP_

/*
 * Types for representing various objects.
 */

namespace align_reads
{
    /* Stores alignment info of a query w.r.t. to a target.
     * Always exist w.r.t. to a target.
     */
    class alignment_segment
    {
    public:
        // q_start/t_start are offsets from the start of query/target
        // when aligning with edlib. Result should contain path.
        // Frees(consumes) the edlib result.
        alignment_segment(std::string& query, std::string& target,
         std::uint32_t q_start, std::uint32_t t_start, EdlibAlignResult& result);

        // Get the characters (possibily '_' ) aligned to the target segment.
        std::string& get_aligned_chars();

        // Get the insertion segment at the specified index of the target.
        // -1 for ins segment before 0.
        std::string& get_ins_segment_at(int index);

    private:
        std::uint32_t start_on_target; // The start index of the aligned segment on the target, 0-based
        std::uint32_t align_len_on_target; // Align with how many from start (inclusive).
        std::string aligned_chars;
        std::vector<std::string> ins_segments;
    };

} // namespace align_reads

#endif // ALIGN_READS_TYPES_HPP_