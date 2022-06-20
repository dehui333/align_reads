#include <iostream>

#include "align_reads/Converter.hpp"

namespace align_reads
{
    AlignmentConverter::AlignmentConverter(MultiAlignment &alignment,
                                           std::uint32_t matrix_width, std::uint32_t matrix_height)
        : alignment_ptr(&alignment), matrix_width(matrix_width), matrix_height(matrix_height)
    {
        // record the highest number of ins at each target pos
        std::vector<std::uint32_t> max_ins_at_pos;
        max_ins_at_pos.resize(alignment_ptr->target.size(), 0);
        for (auto &s : alignment_ptr->alignment_segments)
        {
            std::uint32_t index_on_target = s.start_on_target;
            for (auto &ins_segment : s.ins_segments)
            {
                if (ins_segment.size() > max_ins_at_pos[index_on_target])
                {
                    max_ins_at_pos[index_on_target] = ins_segment.size();
                }
                index_on_target++;
            }
        }

        // Label each position with its window index
        std::vector<std::uint32_t> target_pos_window_index;
        target_pos_window_index.reserve(alignment_ptr->target.size());
        std::vector<std::vector<std::uint32_t>> ins_pos_window_index;
        ins_pos_window_index.resize(alignment_ptr->target.size());
        std::uint32_t i = 0;
        std::uint32_t j = 0;
        for (auto &s : ins_pos_window_index)
        {
            s.reserve(max_ins_at_pos[i++]);
        }
        std::uint32_t width_index = 0;
        for (i = 0; i < alignment_ptr->target.size(); i++)
        {
            target_pos_window_index.push_back((width_index++) / matrix_width);
            for (j = 0; j < max_ins_at_pos[i]; j++)
            {
                ins_pos_window_index[i].push_back((width_index++) / matrix_width);
            }
        }

        // record the segments that fall into each window
        segments_in_windows.resize(width_index / matrix_width + (width_index % matrix_width == 0 ? 0 : 1));
        i = 0;
        for (auto &segment : alignment_ptr->alignment_segments)
        {
            std::uint32_t first_window = target_pos_window_index[segment.start_on_target];
            std::uint32_t last_window = target_pos_window_index[segment.end_on_target];
            if (segment.ins_segments.back().size() > 0)
            {
                last_window = ins_pos_window_index[segment.end_on_target][segment.ins_segments.back().size() - 1];
            }
           
            for (j = first_window; j <= last_window; j++)
            {
                segments_in_windows[j].push_back(i);
            }

            i++;
        }
    }

} // namespace align_reads