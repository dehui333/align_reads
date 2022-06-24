#include <stdlib.h>  // srand, rand
#include <time.h>  // time

#include "align_reads/Converter.hpp"

#define PAD_CODE 5 
#define PAD_CHAR '*' // Also defined in Aligner.cpp
namespace align_reads
{

    constexpr static std::uint8_t ENCODER[] = {
        255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
        255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
        255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
        255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
        255, 255, 5, 255, 255, 255, 255, 255, 255, 255,
        255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
        255, 255, 255, 255, 255, 0, 255, 1, 255, 255,
        255, 2, 255, 255, 255, 255, 255, 255, 255, 255,
        255, 255, 255, 255, 3, 255, 255, 255, 255, 255,
        255, 255, 255, 255, 255, 4, 255, 255, 255, 255,
        255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
        255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
        255, 255, 255, 255, 255, 255, 255, 255};
    constexpr static char DECODER[] = {
        'A', 'C', 'G', 'T', '_', '*'};

    AlignmentConverter::AlignmentConverter(MultiAlignment &alignment,
                                           std::uint32_t matrix_height, std::uint32_t matrix_width)
        : alignment_ptr(&alignment), matrix_width(matrix_width), matrix_height(matrix_height)
    {
        std::uint32_t i = 0;
        std::uint32_t j = 0;
        // record the highest number of ins at each target pos
        std::vector<std::uint32_t> max_ins_at_pos;
        max_ins_at_pos.resize(alignment_ptr->target.size(), 0); // * ignoring ins to the left of target
        for (auto &s : alignment_ptr->alignment_segments)
        {
            std::uint32_t index_on_target = s.start_on_target;
            for (i = 1; i < s.ins_segments.size(); i++) // * ignoring ins to the left
            {
                auto &ins_segment = s.ins_segments[i];
                if (ins_segment.size() > max_ins_at_pos[index_on_target])
                {
                    max_ins_at_pos[index_on_target] = ins_segment.size();
                }
                index_on_target++;
            }
        }
        max_ins_at_pos.back() = 0; // * ignore ins to the right
        // Label each position with its window index
        
        std::vector<std::uint32_t> target_pos_window_index;
        target_pos_window_index.reserve(alignment_ptr->target.size());
        std::vector<std::vector<std::uint32_t>> ins_pos_window_index;
        ins_pos_window_index.resize(alignment_ptr->target.size());
        std::uint32_t total_width = alignment_ptr->target.size();
        i = 0;
        for (auto &s : ins_pos_window_index)
        {
            total_width += max_ins_at_pos[i];
            s.reserve(max_ins_at_pos[i++]);
        }
        width_idx_to_pos_idx.reserve(total_width);
        std::uint32_t width_index = 0;
        for (i = 0; i < alignment_ptr->target.size(); i++)
        {
            if (width_index % matrix_width == 0)
                start_of_windows.emplace_back(i, 0); // record start of window
            width_idx_to_pos_idx.emplace_back(i, 0);
            target_pos_window_index.push_back((width_index++) / matrix_width);
            for (j = 0; j < max_ins_at_pos[i]; j++)
            {
                if (width_index % matrix_width == 0)
                    start_of_windows.emplace_back(i, j + 1); // record start of window
                width_idx_to_pos_idx.emplace_back(i, j + 1);
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
            /*
            do not consider ins to right
            if (segment.ins_segments.back().size() > 0)
            {
                last_window = ins_pos_window_index[segment.end_on_target][segment.ins_segments.back().size() - 1];
            }*/

            for (j = first_window; j <= last_window; j++)
            {
                segments_in_windows[j].push_back(i);
            }

            i++;
        }
    }

    void AlignmentConverter::fill_row_from_alignment(PyObject *matrix, std::uint32_t window,
                                                     std::uint32_t row, std::uint32_t alignment_idx)
    {
        uint8_t *value_ptr;
        std::uint16_t col_index = 0;
        std::uint32_t width_index = window * matrix_width;
        AlignmentSegment& segment = alignment_ptr->alignment_segments[alignment_idx];
        auto current_index = width_idx_to_pos_idx[width_index++];
        std::uint32_t target_index = current_index.first;
        std::uint32_t ins_index = current_index.second;
        // pad on left
        while(target_index < segment.start_on_target)
        {
            value_ptr = (uint8_t*) PyArray_GETPTR2(matrix, row, col_index++);
            *value_ptr = PAD_CODE;
            if (col_index == matrix_width) return;
            current_index = width_idx_to_pos_idx[width_index++];
            target_index = current_index.first;
            ins_index = current_index.second;
        }
        // aligned segment 
        while (target_index <= segment.end_on_target)
        {
            value_ptr = (uint8_t*) PyArray_GETPTR2(matrix, row, col_index++);    
            *value_ptr = ENCODER[static_cast<std::uint8_t>(segment.get_at_target_pos(target_index, ins_index))];
            if (col_index == matrix_width) return;
            if (width_index == width_idx_to_pos_idx.size()) break;
            current_index = width_idx_to_pos_idx[width_index++];
            target_index = current_index.first;        
            ins_index = current_index.second;
        }
        while (col_index != matrix_width)
        {
            value_ptr = (uint8_t*) PyArray_GETPTR2(matrix, row, col_index++);
            *value_ptr = PAD_CODE;
        }
    }

    std::vector<std::vector<std::uint32_t>> AlignmentConverter::choose_segments(std::uint32_t num_reserved_for_target, bool sample_target)
    {
        std::vector<std::vector<std::uint32_t>> result;
        result.resize(segments_in_windows.size());
        for (std::uint32_t i = 0; i < segments_in_windows.size(); i++)
        {
            std::uint32_t num_in_window = segments_in_windows[i].size();
            if (num_in_window == 0) continue; // if no segment falls into this window
            result[i].reserve(matrix_height);
            result[i].resize(num_reserved_for_target, alignment_ptr->alignment_segments.size());
            std::uint32_t num_choices = sample_target ? alignment_ptr->alignment_segments.size() + 1 : alignment_ptr->alignment_segments.size();  
            for (std::uint32_t j = 0; j < matrix_height - num_reserved_for_target; j++)
            {
                result[i].push_back(rand() % num_choices);
            }
        }
        return result;
    }

} // namespace align_reads