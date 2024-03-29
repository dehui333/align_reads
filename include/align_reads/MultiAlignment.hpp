#ifndef ALIGN_READS_MULTIALIGNMENT_HPP_
#define ALIGN_READS_MULTIALIGNMENT_HPP_

#include "align_reads/AlignmentSegment.hpp"

// ** Maybe can have iterator over each segment under the context of all the segments
// - will include positions where some segments have gaps and some have non gap chars.

// iterator over target string taking into account gaps due to ins over all aligned? (which will include pairwise)

namespace align_reads
{

    // Stores info on how a set of sequences align to a single target.
    class MultiAlignment
    {
    public:
        MultiAlignment(std::string &target, std::vector<AlignmentSegment> &segments);

        MultiAlignment(std::string &target, std::vector<AlignmentSegment> &&segments);
        // steals content of inputs
        MultiAlignment(std::string &&target, std::vector<AlignmentSegment> &&segments);

        MultiAlignment(std::string &&target, std::vector<AlignmentSegment> &&segments, std::vector<AlignmentSegment> &&truth);

        class MultiAlignmentIterator
        {
        public:
            bool has_next();
            char next();

            MultiAlignmentIterator(std::uint32_t start_width_idx,
                                   AlignmentSegment &segment,
                                   std::vector<std::pair<std::uint32_t, std::uint32_t>> &width_idx_to_pos_idx,
                                   std::uint32_t start_t_idx, std::uint32_t start_i_idx);

        private:
            std::uint32_t width_idx;
            std::vector<std::pair<std::uint32_t, std::uint32_t>> &width_idx_to_pos_idx;
            AlignmentSegment::AlignmentIterator iter;
            bool withholding = false;
            std::uint32_t waiting_for;
            char stored_char;
        };

        MultiAlignmentIterator iterator(std::uint32_t alignment_idx, std::uint32_t start_width_idx);
        MultiAlignmentIterator target_iterator(AlignmentSegment &target_segment, std::uint32_t start_width_idx);

        void print_in_window(std::uint32_t start_width_index, std::uint32_t length);

    private:
        std::string target; // The target sequence
        std::vector<AlignmentSegment> alignment_segments;
        std::vector<AlignmentSegment> truth_to_target;
        // std::uint32_t largest_truth_start;
        // std::uint32_t smallest_truth_end;
        //  Restricts the effective span to span with truth alignment
        std::uint32_t start_on_target;
        std::uint32_t end_on_target;

        // a mapping from 'width index' to 'pos index'
        // width index: the column indices along the width of the alignment
        // pos index: An index which takes into account which position is aligned
        // to an existing base on the target and which are insertions etc
        std::vector<std::pair<std::uint32_t, std::uint32_t>> width_idx_to_pos_idx;

        friend class AlignmentConverter;

        void initialize();
    };

} // namespace align_reads

#endif // ALIGN_READS_MULTIALIGNMENT_HPP_
