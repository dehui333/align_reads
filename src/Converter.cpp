#include <algorithm>
#include <stdlib.h> // srand, rand
#include <time.h>   // time
#include <iostream>

#include "align_reads/Converter.hpp"
#include "align_reads/Utilities.hpp"

#define PAD_CODE 5
#define PAD_CHAR '*' // Also defined in Aligner.cpp

#define ROWS_FOR_TARGET 0  // Reserve this number of rows for the target seq in each matrix
#define SAMPLE_TARGET true // Do we sample the target seq too when filling matrix?
#define MIN_ALIGNED 1 // Min number of aligned segments for a window's matrix to be used
#define MAX_UINT32 static_cast<std::uint32_t>(-1)
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
        info_ptr = nullptr;
        
        // The span of the alignment target (start/end) to consider may vary due 
        // to unreliable alignments at the ends etc
        std::uint32_t start_on_target = alignment.start_on_target; // start index on target from which to consider
        std::uint32_t end_on_target = alignment.end_on_target; // inclusive end index
        
        std::vector<std::uint32_t> target_pos_window_index;
        target_pos_window_index.reserve(alignment.target.size());
        target_pos_window_index.resize(start_on_target, -1);
        std::uint32_t width = alignment.width_idx_to_pos_idx.size(); 
        std::uint32_t i;
        std::uint32_t j;
        for (i = 0; i < width; i++)
        {
            auto index_pair = alignment.width_idx_to_pos_idx[i];
            if (index_pair.second == 0)
            {
                target_pos_window_index.push_back(i/matrix_width);
            }
        }

        // record the segments that fall into each window
        segments_in_windows.resize(width / matrix_width + (width % matrix_width == 0 ? 0 : 1));
        i = 0;
        for (auto &segment : alignment.alignment_segments)
        {
            std::uint32_t first_window = target_pos_window_index[std::max(segment.start_on_target, start_on_target)];
            std::uint32_t last_window = target_pos_window_index[std::min(segment.end_on_target, end_on_target)];
            for (j = first_window; j <= last_window; j++)
            {
                segments_in_windows[j].push_back(i);
            }
            i++;
        }


        /*
        // The span of the alignment target (start/end) to consider may vary due 
        // to unreliable alignments at the ends etc
        std::uint32_t start_on_target = alignment.start_on_target; // start index on target from which to consider
        std::uint32_t end_on_target = alignment.end_on_target; // inclusive end index


        std::uint32_t i = 0;
        std::uint32_t j = 0;
        // record the highest number of ins at each target pos
        std::vector<std::uint32_t> max_ins_at_pos;
        max_ins_at_pos.resize(alignment_ptr->target.size(), 0); // * ignoring ins to the left of target
    
        for (auto &s : alignment_ptr->alignment_segments)
        {
            std::uint32_t index_on_target = s.start_on_target;
            if (start_on_target > index_on_target) index_on_target = start_on_target;
            
            for (i = index_on_target - s.start_on_target + 1; i < s.ins_segments.size(); i++) // * ignoring ins to the left
            {
                auto &ins_segment = s.ins_segments[i];
                if (ins_segment.size() > max_ins_at_pos[index_on_target])
                {
                    max_ins_at_pos[index_on_target] = ins_segment.size();
                }
                index_on_target++;
                if (index_on_target > end_on_target) break;
            }
        }
        
        max_ins_at_pos[end_on_target] = 0; // * ignore ins to the right
        // Label each position with its window index

        std::vector<std::uint32_t> target_pos_window_index;
        target_pos_window_index.reserve(alignment_ptr->target.size());
        target_pos_window_index.resize(start_on_target, -1);
        i = start_on_target;

        width_idx_to_pos_idx.reserve(end_on_target - start_on_target + 1);
        std::uint32_t width_index = 0;
        for (i = start_on_target; i <= end_on_target; i++)
        {
            width_idx_to_pos_idx.emplace_back(i, 0);
            target_pos_window_index.push_back((width_index++) / matrix_width);
            for (j = 0; j < max_ins_at_pos[i]; j++)
            {
                 
                width_idx_to_pos_idx.emplace_back(i, j + 1);
                width_index++;
            }
        }
        

        // record the segments that fall into each window
        segments_in_windows.resize(width_index / matrix_width + (width_index % matrix_width == 0 ? 0 : 1));
        i = 0;
        for (auto &segment : alignment_ptr->alignment_segments)
        {
            std::uint32_t first_window = target_pos_window_index[std::max(segment.start_on_target, start_on_target)];
            std::uint32_t last_window = target_pos_window_index[std::min(segment.end_on_target, end_on_target)];
            for (j = first_window; j <= last_window; j++)
            {
                segments_in_windows[j].push_back(i);
            }
            i++;
        }
        */
        
        
    }

    AlignmentConverter::AlignmentConverter(MultiAlignment &alignment, std::uint32_t matrix_height, std::uint32_t matrix_width, Info &info)
        : AlignmentConverter(alignment, matrix_height, matrix_width)
    {
        info_ptr = &info;
    }

    void AlignmentConverter::fill_row_from_alignment(PyObject *matrix, std::uint32_t window,
                                                     std::uint32_t row, AlignmentSegment &segment)
    {
        uint8_t *value_ptr;
        std::uint16_t col_index = 0;
        std::uint32_t width_index = window * matrix_width;
        auto& width_idx_to_pos_idx = alignment_ptr->width_idx_to_pos_idx;
        auto current_index = width_idx_to_pos_idx[width_index++];
        std::uint32_t target_index = current_index.first;
        std::uint32_t ins_index = current_index.second;
        
        // pad on left
        while (target_index < segment.start_on_target)
        {
            value_ptr = (uint8_t *)PyArray_GETPTR2(matrix, row, col_index++);
            *value_ptr = PAD_CODE;
            if (col_index == matrix_width)
                return;
            current_index = width_idx_to_pos_idx[width_index++];
            target_index = current_index.first;
            ins_index = current_index.second;
        }

        // aligned segment
        while (target_index <= segment.end_on_target)
        {
            value_ptr = (uint8_t *)PyArray_GETPTR2(matrix, row, col_index++);
            *value_ptr = ENCODER[static_cast<std::uint8_t>(segment.get_at_target_pos(target_index, ins_index))];
            if (col_index == matrix_width)
                return;
            if (width_index == width_idx_to_pos_idx.size())
                break;
            current_index = width_idx_to_pos_idx[width_index++];
            target_index = current_index.first;
            ins_index = current_index.second;
        }
        while (col_index != matrix_width)
        {
            value_ptr = (uint8_t *)PyArray_GETPTR2(matrix, row, col_index++);
            *value_ptr = PAD_CODE;
        }
    }

    void AlignmentConverter::fill_row_from_target(PyObject *matrix, std::uint32_t window,
                                                  std::uint32_t row)
    {
        uint8_t *value_ptr;
        std::uint16_t col_index = 0;
        std::uint32_t width_index = window * matrix_width;
        auto& width_idx_to_pos_idx = alignment_ptr->width_idx_to_pos_idx;
        auto current_index = width_idx_to_pos_idx[width_index++];
        std::uint32_t target_index = current_index.first;
        std::uint32_t ins_index = current_index.second;
        auto &target = alignment_ptr->target;

        while (col_index < matrix_width)
        {
            value_ptr = (uint8_t *)PyArray_GETPTR2(matrix, row, col_index++);
            if (ins_index == 0)
            {
                *value_ptr = ENCODER[static_cast<std::uint8_t>(target[target_index])];
            }
            else
            {
                *value_ptr = ENCODER[static_cast<std::uint8_t>('_')];
            }
            if (width_index == width_idx_to_pos_idx.size())
                break;
            current_index = width_idx_to_pos_idx[width_index++];
            target_index = current_index.first;
            ins_index = current_index.second;
        }
        while (col_index != matrix_width)
        {
            value_ptr = (uint8_t *)PyArray_GETPTR2(matrix, row, col_index++);
            *value_ptr = PAD_CODE;
        }
    }
    
    std::vector<std::vector<std::uint32_t>> AlignmentConverter::choose_segments(std::uint32_t num_reserved_for_target, bool sample_target)
    {
        std::vector<std::vector<std::uint32_t>> result;
        result.resize(segments_in_windows.size());
        
        for (std::uint32_t i = 0; i < segments_in_windows.size(); i++)
        {
            if (segments_in_windows[i].size() < MIN_ALIGNED)
                continue; // if no segment falls into this window
            std::uint32_t num_alignments_in_window = segments_in_windows[i].size();
            result[i].reserve(matrix_height);
            result[i].resize(num_reserved_for_target, num_alignments_in_window);
            std::uint32_t num_choices = sample_target ? num_alignments_in_window + 1 : num_alignments_in_window;
            for (std::uint32_t j = 0; j < matrix_height - num_reserved_for_target; j++)
            {
                std::uint16_t random_number = static_cast<std::uint16_t>(rand() % num_choices);
                if (random_number == num_alignments_in_window)
                {
                    result[i].push_back(-1);
                }
                else
                {
                    result[i].push_back(segments_in_windows[i][random_number    ]);
                }
                
            }
        }
        return result;
    }

    Data AlignmentConverter::produce_data(std::shared_ptr<thread_pool::ThreadPool> &pool, bool with_labels)
    {
        std::vector<std::vector<std::uint32_t>> chosen = choose_segments(ROWS_FOR_TARGET, SAMPLE_TARGET);
        Data data;
        /*
        std::cout << "chosen in window 0 " << std::endl;
        for (auto i: chosen[0])
        {
            std::cout << i << std::endl;
        }
        std::cout << "segments in window :" << std::endl;
        for (auto i : segments_in_windows[0])
        {
            std::cout << "segment " << i << std::endl;
            alignment_ptr->alignment_segments[i].print(alignment_ptr->target);
        }*/
        data.Xs = produce_alignment_matrices(chosen, pool);
        if (with_labels)
        {
            // TODO
            // attach label info to the struct
        }
        return data;
    }

    std::vector<PyObject *> AlignmentConverter::produce_alignment_matrices(std::vector<std::vector<std::uint32_t>> &chosen, std::shared_ptr<thread_pool::ThreadPool> &pool)
    {
        //std::uint32_t num_alignment = alignment_ptr->alignment_segments.size();
        npy_intp dims[2];
        dims[0] = matrix_height;
        dims[1] = matrix_width;
        std::uint32_t i = 0;
        std::uint32_t m_idx = 0;
        std::vector<PyObject *> Xs;
        // ---> I think PyArray_SimpleNew cannot have multiple called in parallel!!!
        Xs.reserve(chosen.size());
        for (auto &list : chosen)
        {
            if (!list.empty())
                Xs.push_back(PyArray_SimpleNew(2, dims, NPY_UINT8));
        }

        auto produce_matrix = [&](std::vector<std::uint32_t> &alignment_indices, std::uint32_t window_index, std::uint32_t matrix_idx)
        {
            uint32_t row_idx = 0;
            for (auto alignment_idx : alignment_indices)
            {
                if (alignment_idx == MAX_UINT32)
                {
                    this->fill_row_from_target(Xs[matrix_idx], window_index, row_idx++);
                }
                else
                {
                    this->fill_row_from_alignment(Xs[matrix_idx], window_index, row_idx++, alignment_ptr->alignment_segments[alignment_idx]);
                }
            }
        };

        if (pool == nullptr)
        {
            std::uint32_t row_idx = 0;
            for (; i < chosen.size(); i++)
            {
                if (chosen[i].empty())
                    continue;
                row_idx = 0;
                for (auto alignment_idx : chosen[i])
                {
                    if (alignment_idx == MAX_UINT32)
                    {
                        fill_row_from_target(Xs[m_idx], i, row_idx++);
                    }
                    else
                    {
                        fill_row_from_alignment(Xs[m_idx], i, row_idx++, alignment_ptr->alignment_segments[alignment_idx]);
                    }
                }
                m_idx++;
            }
        }
        else
        {
            Futures<void> futures(pool, chosen.size());

            for (; i < chosen.size(); i++)
            {
                if (chosen[i].empty())
                    continue;
                futures.add_inputs(produce_matrix, chosen[i], i, m_idx++);
            }
            futures.finish();
        }
        return Xs;
    }

    std::vector<PyObject *> AlignmentConverter::produce_truth_matrices(std::vector<std::vector<std::uint32_t>> &chosen, std::shared_ptr<thread_pool::ThreadPool> &pool)
    {
        //std::uint32_t num_alignment = alignment_ptr->alignment_segments.size();
        npy_intp dims[2];
        dims[0] = matrix_height;
        dims[1] = matrix_width;
        std::uint32_t i = 0;
        std::uint32_t m_idx = 0;
        std::vector<PyObject *> Ys;
        // ---> I think PyArray_SimpleNew cannot have multiple called in parallel!!!
        Ys.reserve(chosen.size());
        for (auto &list : chosen)
        {
            if (!list.empty())
                Ys.push_back(PyArray_SimpleNew(2, dims, NPY_UINT8));
        }

        auto produce_matrix = [&](std::vector<std::uint32_t> &alignment_indices, std::uint32_t window_index, std::uint32_t matrix_idx, Info* info_ptr)
        {
            uint32_t row_idx = 0;
            for (auto alignment_idx : alignment_indices)
            {
                if (alignment_idx == MAX_UINT32)
                {
                    this->fill_row_from_alignment(Ys[matrix_idx], window_index, row_idx++, alignment_ptr->truth_to_target[info_ptr->hap_of_target]);
                }
                else
                {
                    this->fill_row_from_alignment(Ys[matrix_idx], window_index, row_idx++, alignment_ptr->truth_to_target[info_ptr->hap_of_aligned[alignment_idx]]);
                }
            }
        };

        if (pool == nullptr)
        {
            std::uint32_t row_idx = 0;
            for (; i < chosen.size(); i++)
            {
                if (chosen[i].empty())
                    continue;
                row_idx = 0;
                for (auto alignment_idx : chosen[i])
                {
                    if (alignment_idx == MAX_UINT32)
                    {
                        this->fill_row_from_alignment(Ys[m_idx], i, row_idx++, alignment_ptr->truth_to_target[info_ptr->hap_of_target]);
                    }
                    else
                    {
                        this->fill_row_from_alignment(Ys[m_idx], i, row_idx++, alignment_ptr->truth_to_target[info_ptr->hap_of_aligned[alignment_idx]]);
                    }
                }
                m_idx++;
            }
        }
        else
        {
            Futures<void> futures(pool, chosen.size());

            for (; i < chosen.size(); i++)
            {
                if (chosen[i].empty())
                    continue;
                futures.add_inputs(produce_matrix, chosen[i], i, m_idx++, this->info_ptr);
            }
            futures.finish();
        }

        return Ys;
    }

    void AlignmentConverter::print_alignments_in_window(std::uint32_t window)
    {
        std::uint32_t width_idx = window * matrix_width;
        alignment_ptr->print_in_window(width_idx, matrix_width);
    }

} // namespace align_reads