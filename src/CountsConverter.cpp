#include <iostream>

#include "align_reads/CountsConverter.hpp"

#define NUM_COUNTS 5
#define NUM_STATS 5

namespace align_reads
{

    inline void fill_in_column(std::uint16_t column_idx, std::uint16_t num_rows, PyObject *matrix, std::vector<std::uint16_t> &counter)
    {
        uint16_t *value_ptr;
        for (std::uint16_t row_idx = 0; row_idx < num_rows; row_idx++)
        {
            value_ptr = (uint16_t *)PyArray_GETPTR2(matrix, row_idx, column_idx);
            *value_ptr = counter[row_idx];
        }
    }

    inline void fill_in_column_stats(std::uint16_t column_idx, std::uint16_t num_rows, PyObject *matrix, std::vector<float> &stat)
    {
        float *value_ptr;
        for (std::uint16_t row_idx = 0; row_idx < num_rows; row_idx++)
        {
            value_ptr = (float *)PyArray_GETPTR2(matrix, row_idx, column_idx);
            *value_ptr = stat[row_idx];
        }
    }

    // start padding in this matrix from col_idx
    inline void pad_zeros(PyObject *matrix, std::uint16_t col_idx, std::uint16_t window_length, std::uint16_t num_rows)
    {
        uint16_t *value_ptr;
        for (; col_idx < window_length; col_idx++)
        {
            for (std::uint16_t row_idx = 0; row_idx < num_rows; row_idx++)
            {
                value_ptr = (uint16_t *)PyArray_GETPTR2(matrix, row_idx, col_idx);
                *value_ptr = 0;
            }
        }
    }

    // start padding in this matrix from col_idx
    inline void pad_zeros_stats(PyObject *matrix, std::uint16_t col_idx, std::uint16_t window_length, std::uint16_t num_rows)
    {
        float *value_ptr;
        for (; col_idx < window_length; col_idx++)
        {
            for (std::uint16_t row_idx = 0; row_idx < num_rows; row_idx++)
            {
                value_ptr = (float *)PyArray_GETPTR2(matrix, row_idx, col_idx);
                *value_ptr = 0;
            }
        }
    }

    std::vector<PyObject *> CountsConverter::get_counts_matrices(AlignCounter &counter, std::uint16_t window_length, std::uint32_t left_clip, std::uint32_t right_clip, std::uint32_t num_matrices, std::uint32_t alignment_length)
    {
        npy_intp dims[2];
        dims[0] = NUM_COUNTS;
        dims[1] = window_length;
        std::vector<PyObject *> matrices;
        auto &counts = counter.counts;
        matrices.reserve(num_matrices);
        std::uint16_t matrix_idx = 0;
        std::uint16_t col_idx = 0;
        auto end = counts.end() - right_clip;
        for (auto it = counts.begin() + left_clip; it < end; it++)
        {
            auto& counts_at_pos = *it;
            for (auto &counter : counts_at_pos)
            {
                if (col_idx == 0)
                {
                    matrices.push_back(PyArray_SimpleNew(2, dims, NPY_UINT16));
                }
                fill_in_column(col_idx++, NUM_COUNTS, matrices[matrix_idx], counter);
                if (col_idx == window_length)
                {
                    matrix_idx++;
                    col_idx = 0;
                }
            }
        }
        std::uint16_t num_to_pad = num_matrices * window_length - alignment_length;
        pad_zeros(matrices.back(), window_length-num_to_pad, window_length, NUM_COUNTS); 

        return matrices;
    }

    std::vector<PyObject *> CountsConverter::get_stats_matrices(AlignCounter &counter, std::uint16_t window_length, std::uint32_t left_clip, std::uint32_t right_clip, std::uint32_t num_matrices, std::uint32_t alignment_length)
    {
        npy_intp dims[2];
        dims[0] = NUM_STATS;
        dims[1] = window_length;
        std::vector<PyObject *> matrices;
        auto &stats = counter.stats;
        matrices.reserve(num_matrices);
        std::uint16_t matrix_idx = 0;
        std::uint16_t col_idx = 0;
        auto end = stats.end() - right_clip;
        for (auto it = stats.begin() + left_clip; it < end; it++)
        {
            auto& stats_at_pos = *it;
            for (auto &stat : stats_at_pos)
            {
                if (col_idx == 0)
                {
                    matrices.push_back(PyArray_SimpleNew(2, dims, NPY_FLOAT32));
                }
                fill_in_column_stats(col_idx++, NUM_STATS, matrices[matrix_idx], stat);
                if (col_idx == window_length)
                {
                    matrix_idx++;
                    col_idx = 0;
                }
            }
        }
        std::uint16_t num_to_pad = num_matrices * window_length - alignment_length;
        pad_zeros_stats(matrices.back(), window_length-num_to_pad, window_length, NUM_STATS); 

        return matrices;
    }

} // namespace align_reads