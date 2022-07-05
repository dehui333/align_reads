#ifndef ALIGN_READS_ALIGNMENTSEGMENT_HPP_
#define ALIGN_READS_ALIGNMENTSEGMENT_HPP_

#include "edlib.h"

namespace align_reads
{

    /* Stores alignment info of a query w.r.t. to a target.
     * Always exist w.r.t. to a target.
     */
    class AlignmentSegment
    {
    public:
        // q_start/t_start are offsets from the start of query/target
        // when aligning with edlib. Result should contain path.
        // Frees(consumes) the edlib result.
        AlignmentSegment(std::string &query, std::uint32_t q_start,
                         std::string &target, std::uint32_t t_start, EdlibAlignResult &result);

        /*AlignmentSegment(std::string &query, std::string &target,
                          std::uint32_t q_start, std::uint32_t t_start,
                          std::uint32_t q_len, std::uint32_t t_len,
                           EdlibAlignMode mode, EdlibAlignTask task);*/

        // Get the characters (possibily '_' ) aligned to the target segment.
        std::string &get_aligned_chars();

        // Get the insertion segment at the specified index of the target.
        // -1 for ins segment before 0.
        std::string &get_ins_segment_at(int index);

        // Get the aligned/inserted char w.r.t. to positions on the target
        // Assumes the position is within the alignment span of this segment excluding
        // ins before and after.
        // (>= (start_on_target, 0) and <= (end_on_target, 0) )
        // Returns gap for ins position in span but which this alignment has no value.
        char get_at_target_pos(std::uint32_t target_index, std::uint32_t ins_index);

        // Print this alignment segment (target string is required)
        void print(std::string &target);

    private:
        std::uint32_t start_on_target; // The start index of the aligned segment on the target, 0-based
        std::uint32_t end_on_target;   // The end index, 0-based, inclusive
        std::string aligned_chars;
        std::vector<std::string> ins_segments;

        friend class MultiAlignment;
        friend class AlignmentConverter;

        struct aligned_pos
        {
            char c;
            std::uint32_t target_index;
            std::uint32_t ins_index;
            aligned_pos(char c, std::uint32_t t_idx, std::uint32_t i_idx) : c(c), target_index(t_idx), ins_index(i_idx) {}
        };
        // untested
        class AlignmentIterator
        {
        public:
            bool has_next()
            {
                return current_target_index <= segment.end_on_target;
            }
            aligned_pos next()
            {
                if (current_ins_index == 0)
                {
                    return {segment.aligned_chars[current_target_index - segment.start_on_target], current_target_index, current_ins_index};
                }
                else
                {

                    aligned_pos r{segment.ins_segments[current_target_index - segment.start_on_target + 1][current_ins_index - 1], current_target_index, current_ins_index};
                    if (current_ins_index == current_max_ins)
                    {
                        current_target_index++;
                        current_ins_index = 0;
                    }
                    else
                    {
                        current_ins_index++;
                    }
                    return r;
                }
            }
            AlignmentIterator(AlignmentSegment &segment, std::uint32_t start_t_idx, std::uint32_t start_i_idx)
                : segment(segment), current_target_index(start_t_idx), current_ins_index(start_i_idx)
            {
                current_max_ins = segment.ins_segments[start_t_idx - segment.start_on_target + 1].size();
            }

        private:
            AlignmentSegment &segment;
            std::uint32_t current_target_index;
            std::uint32_t current_ins_index;
            std::uint32_t current_max_ins;
        };
    };

}

#endif // ALIGN_READS_ALIGNMENTSEGMENT_HPP_
