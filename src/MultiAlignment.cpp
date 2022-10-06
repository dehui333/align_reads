#include "align_reads/MultiAlignment.hpp"

#define GAP_CHAR '_'
#define PAD_CHAR '*'
namespace align_reads
{

    //----------------- MultiAlignment---------------------------

    MultiAlignment::MultiAlignment(std::string &target,
                                   std::vector<AlignmentSegment> &segments) : target(target), alignment_segments(segments)
    {
        initialize();
    }
    MultiAlignment::MultiAlignment(std::string &&target,
                                   std::vector<AlignmentSegment> &&segments) : target(std::move(target)), alignment_segments(std::move(segments))
    {
        initialize();
    }
    MultiAlignment::MultiAlignment(std::string &&target,
                                   std::vector<AlignmentSegment> &&segments,
                                   std::vector<AlignmentSegment> &&truth)
        : target(std::move(target)), alignment_segments(std::move(segments)), truth_to_target(std::move(truth))
    {
        initialize();
    }

    void MultiAlignment::initialize()
    {
        // The span of the alignment target (start/end) to consider may vary due
        // to unreliable alignments at the ends etc
        start_on_target = 0;               // start index on target from which to consider
        end_on_target = target.size() - 1; // inclusive end index

        if (!truth_to_target.empty())
        {
            for (auto &s : truth_to_target)
            {
                if (s.start_on_target > start_on_target)
                    start_on_target = s.start_on_target;
                if (s.end_on_target < end_on_target)
                    end_on_target = s.end_on_target;
            }
        }

        std::uint32_t i = 0;
        std::uint32_t j = 0;
        // record the highest number of ins at each target pos
        std::vector<std::uint32_t> max_ins_at_pos;
        max_ins_at_pos.resize(target.size(), 0); // * ignoring ins to the left of target

        for (auto &s : alignment_segments)
        {
            std::uint32_t index_on_target = s.start_on_target;
            if (start_on_target > index_on_target)
                index_on_target = start_on_target;

            for (i = index_on_target - s.start_on_target + 1; i < s.ins_segments.size(); i++) // * ignoring ins to the left
            {
                auto &ins_segment = s.ins_segments[i];
                if (ins_segment.size() > max_ins_at_pos[index_on_target])
                {
                    max_ins_at_pos[index_on_target] = ins_segment.size();
                }
                index_on_target++;
                if (index_on_target > end_on_target)
                    break;
            }
        }

        max_ins_at_pos[end_on_target] = 0; // * ignore ins to the right
        // Label each position with its window index

        width_idx_to_pos_idx.reserve(end_on_target - start_on_target + 1);
        for (i = start_on_target; i <= end_on_target; i++)
        {
            width_idx_to_pos_idx.emplace_back(i, 0);
            for (j = 0; j < max_ins_at_pos[i]; j++)
            {

                width_idx_to_pos_idx.emplace_back(i, j + 1);
            }
        }
    }
    MultiAlignment::MultiAlignmentIterator::MultiAlignmentIterator(std::uint32_t start_width_idx,
                                                                   AlignmentSegment &segment, std::vector<std::pair<std::uint32_t, std::uint32_t>> &width_idx_to_pos_idx)
        : width_idx(start_width_idx), width_idx_to_pos_idx(width_idx_to_pos_idx), iter(segment, segment.start_on_target, 0) {}

    char MultiAlignment::MultiAlignmentIterator::next()
    {
        if (withholding)
        {
            std::pair<std::uint32_t, std::uint32_t> &next_pos_idx = width_idx_to_pos_idx[width_idx++];
            if (next_pos_idx.first == waiting_for)
            {
                withholding = false;
                return stored_char;
            }
            else
            {
                return GAP_CHAR;
            }
        }
        else
        {
            if (iter.has_next())
            {
                std::pair<std::uint32_t, std::uint32_t> &next_pos_idx = width_idx_to_pos_idx[width_idx++];
                aligned_pos possible_next = iter.next();
                if (possible_next.target_index != next_pos_idx.first)
                {
                    withholding = true;
                    stored_char = possible_next.c;
                    waiting_for = possible_next.target_index;
                    return GAP_CHAR;
                }
                return possible_next.c;
            }
            else
            {
                width_idx++;
                return GAP_CHAR;
            }
        }
    }

    bool MultiAlignment::MultiAlignmentIterator::has_next()
    {
        return width_idx < width_idx_to_pos_idx.size();
    }

    MultiAlignment::MultiAlignmentIterator MultiAlignment::iterator(std::uint32_t alignment_idx, std::uint32_t start_width_idx)
    {
        return {start_width_idx, alignment_segments[alignment_idx], width_idx_to_pos_idx};
    }

} // namespace align_reads