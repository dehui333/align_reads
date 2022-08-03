#include <gtest/gtest.h>
#include <Python.h>


#define PY_ARRAY_UNIQUE_SYMBOL gen_ARRAY_API
#include "numpy/arrayobject.h"

#include "align_reads/Utilities.hpp"

#define private public
#include "../src/Converter.cpp"

#define PAD_CODE 5

using namespace align_reads;


/*
 *
 * Actually depends on both Aligner & 
 */


extern std::shared_ptr<thread_pool::ThreadPool> pool;

PyMODINIT_FUNC PyInit_align_reads_gen(void) {
    Py_Initialize();
    import_array();
    return 0;
}

TEST(Converter, Converter_all_match)
{
    
    std::string t  = "TAGGCATACCGG";
    std::string q1 = "TCGG";
    std::string q2 = "CACA";
    std::string q3 = "CCGGAAAAAAAAAAAAAAAAAAA";

    Futures<AlignmentSegment> futures(pool, 3);
    futures.add_inputs(get_alignment_segment, q1, 0, 4, t, 0, 4, EDLIB_MODE_INFIX, EDLIB_TASK_PATH);
    futures.add_inputs(get_alignment_segment, q2, 0, 4, t, 4, 4, EDLIB_MODE_INFIX, EDLIB_TASK_PATH);
    futures.add_inputs(get_alignment_segment, q3, 0, 4, t, 8, 4, EDLIB_MODE_INFIX, EDLIB_TASK_PATH);
    std::vector<AlignmentSegment> segments = futures.get();
    
    EXPECT_EQ(segments[0].aligned_chars, "TCGG");
    EXPECT_EQ(segments[1].aligned_chars, "CACA");
    EXPECT_EQ(segments[2].aligned_chars, "CCGG");
    MultiAlignment m_align {std::move(t), std::move(segments)};
    AlignmentConverter converter {m_align, 3, 4};

    EXPECT_EQ(converter.segments_in_windows.size(), 3);
    std::set<std::uint32_t> s0(converter.segments_in_windows[0].begin(), converter.segments_in_windows[0].end());
    std::set<std::uint32_t> a0 {0};
    std::set<std::uint32_t> s1(converter.segments_in_windows[1].begin(), converter.segments_in_windows[1].end());
    std::set<std::uint32_t> a1 {1};
    std::set<std::uint32_t> s2(converter.segments_in_windows[2].begin(), converter.segments_in_windows[2].end());
    std::set<std::uint32_t> a2 {2};
    
    EXPECT_EQ(s0, a0);
    EXPECT_EQ(s1, a1);
    EXPECT_EQ(s2, a2);

    t  = "TAGGCATACCGG";
    q1 = "AGGC";
    q2 = "GCAT";
    q3 = "ACCG";

    futures.add_inputs(get_alignment_segment, q1, 0, 4, t, 0, 6, EDLIB_MODE_INFIX, EDLIB_TASK_PATH);
    futures.add_inputs(get_alignment_segment, q2, 0, 4, t, 2, 5, EDLIB_MODE_INFIX, EDLIB_TASK_PATH);
    futures.add_inputs(get_alignment_segment, q3, 0, 4, t, 6, 5, EDLIB_MODE_INFIX, EDLIB_TASK_PATH);
    segments = futures.get();
    
    m_align = {std::move(t), std::move(segments)};
    converter = {m_align, 3, 4};
    EXPECT_EQ(converter.segments_in_windows.size(), 3);
    s0 = {converter.segments_in_windows[0].begin(), converter.segments_in_windows[0].end()};
    a0 = {0, 1};
    s1 = {converter.segments_in_windows[1].begin(), converter.segments_in_windows[1].end()};
    a1 = {0, 1, 2};
    s2 = {converter.segments_in_windows[2].begin(), converter.segments_in_windows[2].end()};
    a2 = {2};
    
    EXPECT_EQ(s0, a0);
    EXPECT_EQ(s1, a1);
    EXPECT_EQ(s2, a2);
    
}

TEST(Converter, Converter_del)
{
    
    std::string t  = "TAGGCATACAGG";
    std::string q1 = "TGGC";
    std::string q2 = "CAACA";
    std::string q3 = "TAAGG";

    Futures<AlignmentSegment> futures(pool, 3);
    futures.add_inputs(get_alignment_segment, q1, 0, 4, t, 0, 5, EDLIB_MODE_INFIX, EDLIB_TASK_PATH);
    futures.add_inputs(get_alignment_segment, q2, 0, 5, t, 4, 6, EDLIB_MODE_INFIX, EDLIB_TASK_PATH);
    futures.add_inputs(get_alignment_segment, q3, 0, 5, t, 6, 6, EDLIB_MODE_INFIX, EDLIB_TASK_PATH);
    std::vector<AlignmentSegment> segments = futures.get();
    
    EXPECT_EQ(segments[0].aligned_chars, "T_GGC");
    EXPECT_EQ(segments[1].aligned_chars, "CA_ACA");
    EXPECT_EQ(segments[2].aligned_chars, "TA_AGG");
    MultiAlignment m_align {std::move(t), std::move(segments)};
    AlignmentConverter converter {m_align, 3, 4};
    EXPECT_EQ(converter.segments_in_windows.size(), 3);
    std::set<std::uint32_t> s0(converter.segments_in_windows[0].begin(), converter.segments_in_windows[0].end());
    std::set<std::uint32_t> a0 {0};
    std::set<std::uint32_t> s1(converter.segments_in_windows[1].begin(), converter.segments_in_windows[1].end());
    std::set<std::uint32_t> a1 {0, 1, 2};
    std::set<std::uint32_t> s2(converter.segments_in_windows[2].begin(), converter.segments_in_windows[2].end());
    std::set<std::uint32_t> a2 {1, 2};
    
    EXPECT_EQ(s0, a0);
    EXPECT_EQ(s1, a1);
    EXPECT_EQ(s2, a2);  

}

TEST(Converter, Converter_ins)
{
    
    std::string t  = "TAGGCATACAGG";
    std::string q1 = "TAGTG";
    std::string q2 = "CAACATGG";
    std::string q3 = "TGGCATATCA";
    // TAG_|GCAT|A_CA|_GG

    Futures<AlignmentSegment> futures(pool, 3);
    futures.add_inputs(get_alignment_segment, q1, 0, 5, t, 0, 4, EDLIB_MODE_INFIX, EDLIB_TASK_PATH);
    futures.add_inputs(get_alignment_segment, q2, 0, 8, t, 4, 8, EDLIB_MODE_INFIX, EDLIB_TASK_PATH);
    futures.add_inputs(get_alignment_segment, q3, 0, 10, t, 0, 12, EDLIB_MODE_INFIX, EDLIB_TASK_PATH);
    std::vector<AlignmentSegment> segments = futures.get();
    
    EXPECT_EQ(segments[0].aligned_chars, "TAGG");
    EXPECT_EQ(segments[1].aligned_chars, "CA_ACAGG");
    EXPECT_EQ(segments[2].aligned_chars, "T_GGCATACA");
    EXPECT_EQ(segments[0].get_ins_segment_at(2), "T");
    EXPECT_EQ(segments[1].get_ins_segment_at(5), "T");
    EXPECT_EQ(segments[2].get_ins_segment_at(7), "T");
    MultiAlignment m_align {std::move(t), std::move(segments)};
    
    AlignmentConverter converter {m_align, 3, 4};
    EXPECT_EQ(converter.segments_in_windows.size(), 4);
    std::set<std::uint32_t> s0(converter.segments_in_windows[0].begin(), converter.segments_in_windows[0].end());
    std::set<std::uint32_t> a0 {0, 2};
    std::set<std::uint32_t> s1(converter.segments_in_windows[1].begin(), converter.segments_in_windows[1].end());
    std::set<std::uint32_t> a1 {0, 1, 2};
    std::set<std::uint32_t> s2(converter.segments_in_windows[2].begin(), converter.segments_in_windows[2].end());
    std::set<std::uint32_t> a2 {1, 2};
    std::set<std::uint32_t> s3(converter.segments_in_windows[3].begin(), converter.segments_in_windows[3].end());
    std::set<std::uint32_t> a3 {1};
    
    EXPECT_EQ(s0, a0);
    EXPECT_EQ(s1, a1);
    EXPECT_EQ(s2, a2);  
    EXPECT_EQ(s3, a3);  


    EXPECT_EQ(converter.width_idx_to_pos_idx[0].first, 0);
    EXPECT_EQ(converter.width_idx_to_pos_idx[0].second, 0);
    EXPECT_EQ(converter.width_idx_to_pos_idx[1].first, 1);
    EXPECT_EQ(converter.width_idx_to_pos_idx[1].second, 0);
    EXPECT_EQ(converter.width_idx_to_pos_idx[2].first, 2);
    EXPECT_EQ(converter.width_idx_to_pos_idx[2].second, 0);
    EXPECT_EQ(converter.width_idx_to_pos_idx[3].first, 2);
    EXPECT_EQ(converter.width_idx_to_pos_idx[3].second, 1);
    EXPECT_EQ(converter.width_idx_to_pos_idx[4].first, 3);
    EXPECT_EQ(converter.width_idx_to_pos_idx[4].second, 0);
    EXPECT_EQ(converter.width_idx_to_pos_idx[5].first, 4);
    EXPECT_EQ(converter.width_idx_to_pos_idx[5].second, 0);
    EXPECT_EQ(converter.width_idx_to_pos_idx[6].first, 5);
    EXPECT_EQ(converter.width_idx_to_pos_idx[6].second, 0);
    EXPECT_EQ(converter.width_idx_to_pos_idx[7].first, 6);
    EXPECT_EQ(converter.width_idx_to_pos_idx[7].second, 0);
    EXPECT_EQ(converter.width_idx_to_pos_idx[8].first, 7);
    EXPECT_EQ(converter.width_idx_to_pos_idx[8].second, 0);
    EXPECT_EQ(converter.width_idx_to_pos_idx[9].first, 7);
    EXPECT_EQ(converter.width_idx_to_pos_idx[9].second, 1);
    EXPECT_EQ(converter.width_idx_to_pos_idx[10].first, 8);
    EXPECT_EQ(converter.width_idx_to_pos_idx[10].second, 0);
    EXPECT_EQ(converter.width_idx_to_pos_idx[11].first, 9);
    EXPECT_EQ(converter.width_idx_to_pos_idx[11].second, 0);
    EXPECT_EQ(converter.width_idx_to_pos_idx[12].first, 9);
    EXPECT_EQ(converter.width_idx_to_pos_idx[12].second, 1);
    EXPECT_EQ(converter.width_idx_to_pos_idx[13].first, 10);
    EXPECT_EQ(converter.width_idx_to_pos_idx[13].second, 0);
    EXPECT_EQ(converter.width_idx_to_pos_idx[14].first, 11);
    EXPECT_EQ(converter.width_idx_to_pos_idx[14].second, 0);
}

TEST(Converter, numpy_array)
{
    PyInit_align_reads_gen();
    npy_intp dims[2];
    dims[0] = 2;
    dims[1] = 2;
    auto x = PyArray_SimpleNew(2, dims, NPY_UINT8);   
    uint8_t* value_ptr;
    value_ptr = (uint8_t*) PyArray_GETPTR2(x, 0, 0);
    *value_ptr = 0;
    value_ptr = (uint8_t*) PyArray_GETPTR2(x, 0, 1);
    *value_ptr = 1;
    value_ptr = (uint8_t*) PyArray_GETPTR2(x, 1, 0);
    *value_ptr = 10;
    value_ptr = (uint8_t*) PyArray_GETPTR2(x, 1, 1);
    *value_ptr = 11;
    value_ptr = (uint8_t*) PyArray_GETPTR2(x, 0, 0);
    EXPECT_EQ(*value_ptr, 0);
    value_ptr = (uint8_t*) PyArray_GETPTR2(x, 0, 1);
    EXPECT_EQ(*value_ptr, 1);
    value_ptr = (uint8_t*) PyArray_GETPTR2(x, 1, 0);
    EXPECT_EQ(*value_ptr, 10);
    value_ptr = (uint8_t*) PyArray_GETPTR2(x, 1, 1);
    EXPECT_EQ(*value_ptr, 11); 
}

TEST(Converter, fill_row_from_alignment)
{
    PyInit_align_reads_gen();  
    std::string t  = "TAGGCATACAGG";
    std::string q1 = "TAGTG";
    std::string q2 = "CAACATGGAAAAAAAAAAAAAAAA";
    std::string q3 = "TGGCATATCA";
    // TAG_|GCAT|A_CA|_GG

    Futures<AlignmentSegment> futures(pool, 3);
    futures.add_inputs(get_alignment_segment, q1, 0, 5, t, 0, 4, EDLIB_MODE_INFIX, EDLIB_TASK_PATH);
    futures.add_inputs(get_alignment_segment, q2, 0, 24, t, 4, 8, EDLIB_MODE_INFIX, EDLIB_TASK_PATH);
    futures.add_inputs(get_alignment_segment, q3, 0, 10, t, 0, 12, EDLIB_MODE_INFIX, EDLIB_TASK_PATH);
    futures.add_inputs(get_alignment_segment, t, 0, 12, t, 0, 12, EDLIB_MODE_INFIX, EDLIB_TASK_PATH);
    std::vector<AlignmentSegment> segments = futures.get();
    EXPECT_EQ(segments[1].get_aligned_chars(), "CA_ACAGG");
    EXPECT_EQ(segments[1].get_ins_segment_at(5), "T");
    EXPECT_EQ(segments[1].get_ins_segment_at(7), "AAAAAAAAAAAAAAAA");
    
    MultiAlignment m_align {std::move(t), std::move(segments)};
    
    AlignmentConverter converter {m_align, 4, 12};
    EXPECT_EQ(converter.segments_in_windows.size(), 2);;
    std::vector<AlignmentSegment>& ss = m_align.alignment_segments;
    uint8_t* value_ptr;
    npy_intp dims[2];
    dims[0] = 4;
    dims[1] = 12;
    auto x = PyArray_SimpleNew(2, dims, NPY_UINT8); 
    auto x2 = PyArray_SimpleNew(2, dims, NPY_UINT8);

    converter.fill_row_from_alignment(x, 0, 0, ss[0]);
    converter.fill_row_from_alignment(x, 0, 1, ss[1]);
    converter.fill_row_from_alignment(x, 0, 2, ss[2]);
    converter.fill_row_from_alignment(x, 0, 3, ss[3]);
  
    converter.fill_row_from_alignment(x2, 1, 1, ss[1]);
    converter.fill_row_from_alignment(x2, 1, 3, ss[3]);
    
    // check first matrix first row
    value_ptr = (uint8_t*) PyArray_GETPTR2(x, 0, 0);
    EXPECT_EQ(*value_ptr, ENCODER['T']);
    value_ptr = (uint8_t*) PyArray_GETPTR2(x, 0, 1);
    EXPECT_EQ(*value_ptr, ENCODER['A']);
    value_ptr = (uint8_t*) PyArray_GETPTR2(x, 0, 2);
    EXPECT_EQ(*value_ptr, ENCODER['G']);
    value_ptr = (uint8_t*) PyArray_GETPTR2(x, 0, 3);
    EXPECT_EQ(*value_ptr, ENCODER['T']);
    value_ptr = (uint8_t*) PyArray_GETPTR2(x, 0, 4);
    EXPECT_EQ(*value_ptr, ENCODER['G']);
    value_ptr = (uint8_t*) PyArray_GETPTR2(x, 0, 5);
    EXPECT_EQ(*value_ptr, PAD_CODE);
    value_ptr = (uint8_t*) PyArray_GETPTR2(x, 0, 6);
    EXPECT_EQ(*value_ptr, PAD_CODE);
    value_ptr = (uint8_t*) PyArray_GETPTR2(x, 0, 7);
    EXPECT_EQ(*value_ptr, PAD_CODE);
    value_ptr = (uint8_t*) PyArray_GETPTR2(x, 0, 8);
    EXPECT_EQ(*value_ptr, PAD_CODE);
    value_ptr = (uint8_t*) PyArray_GETPTR2(x, 0, 9);
    EXPECT_EQ(*value_ptr, PAD_CODE);
    value_ptr = (uint8_t*) PyArray_GETPTR2(x, 0, 10);
    EXPECT_EQ(*value_ptr, PAD_CODE);
    value_ptr = (uint8_t*) PyArray_GETPTR2(x, 0, 11);
    EXPECT_EQ(*value_ptr, PAD_CODE);
    

    // check first matrix second row
    value_ptr = (uint8_t*) PyArray_GETPTR2(x, 1, 0);
    EXPECT_EQ(*value_ptr, PAD_CODE);
    value_ptr = (uint8_t*) PyArray_GETPTR2(x, 1, 1);
    EXPECT_EQ(*value_ptr, PAD_CODE);
    value_ptr = (uint8_t*) PyArray_GETPTR2(x, 1, 2);
    EXPECT_EQ(*value_ptr, PAD_CODE);
    value_ptr = (uint8_t*) PyArray_GETPTR2(x, 1, 3);
    EXPECT_EQ(*value_ptr, PAD_CODE);
    value_ptr = (uint8_t*) PyArray_GETPTR2(x, 1, 4);
    EXPECT_EQ(*value_ptr, PAD_CODE);
    value_ptr = (uint8_t*) PyArray_GETPTR2(x, 1, 5);
    EXPECT_EQ(*value_ptr, ENCODER['C']);
    value_ptr = (uint8_t*) PyArray_GETPTR2(x, 1, 6);
    EXPECT_EQ(*value_ptr, ENCODER['A']);
    value_ptr = (uint8_t*) PyArray_GETPTR2(x, 1, 7);
    EXPECT_EQ(*value_ptr, ENCODER['_']);
    value_ptr = (uint8_t*) PyArray_GETPTR2(x, 1, 8);
    EXPECT_EQ(*value_ptr, ENCODER['A']);
    value_ptr = (uint8_t*) PyArray_GETPTR2(x, 1, 9);
    EXPECT_EQ(*value_ptr, ENCODER['_']);
    value_ptr = (uint8_t*) PyArray_GETPTR2(x, 1, 10);
    EXPECT_EQ(*value_ptr, ENCODER['C']);
    value_ptr = (uint8_t*) PyArray_GETPTR2(x, 1, 11);
    EXPECT_EQ(*value_ptr, ENCODER['A']); 

    // check first matrix fourth row
    value_ptr = (uint8_t*) PyArray_GETPTR2(x, 3, 0);
    EXPECT_EQ(*value_ptr, ENCODER['T']);
    value_ptr = (uint8_t*) PyArray_GETPTR2(x, 3, 1);
    EXPECT_EQ(*value_ptr, ENCODER['A']);
    value_ptr = (uint8_t*) PyArray_GETPTR2(x, 3, 2);
    EXPECT_EQ(*value_ptr, ENCODER['G']);
    value_ptr = (uint8_t*) PyArray_GETPTR2(x, 3, 3);
    EXPECT_EQ(*value_ptr, ENCODER['_']);
    value_ptr = (uint8_t*) PyArray_GETPTR2(x, 3, 4);
    EXPECT_EQ(*value_ptr, ENCODER['G']);
    value_ptr = (uint8_t*) PyArray_GETPTR2(x, 3, 5);
    EXPECT_EQ(*value_ptr, ENCODER['C']);
    value_ptr = (uint8_t*) PyArray_GETPTR2(x, 3, 6);
    EXPECT_EQ(*value_ptr, ENCODER['A']);
    value_ptr = (uint8_t*) PyArray_GETPTR2(x, 3, 7);
    EXPECT_EQ(*value_ptr, ENCODER['T']);
    value_ptr = (uint8_t*) PyArray_GETPTR2(x, 3, 8);
    EXPECT_EQ(*value_ptr, ENCODER['A']);
    value_ptr = (uint8_t*) PyArray_GETPTR2(x, 3, 9);
    EXPECT_EQ(*value_ptr, ENCODER['_']);
    value_ptr = (uint8_t*) PyArray_GETPTR2(x, 3, 10);
    EXPECT_EQ(*value_ptr, ENCODER['C']);
    value_ptr = (uint8_t*) PyArray_GETPTR2(x, 3, 11);
    EXPECT_EQ(*value_ptr, ENCODER['A']); 

    // check second matrix second row
    value_ptr = (uint8_t*) PyArray_GETPTR2(x2, 1, 0);
    EXPECT_EQ(*value_ptr, ENCODER['T']);
    value_ptr = (uint8_t*) PyArray_GETPTR2(x2, 1, 1);
    EXPECT_EQ(*value_ptr, ENCODER['G']);
    value_ptr = (uint8_t*) PyArray_GETPTR2(x2, 1, 2);
    EXPECT_EQ(*value_ptr, ENCODER['G']);
    value_ptr = (uint8_t*) PyArray_GETPTR2(x2, 1, 3);
    EXPECT_EQ(*value_ptr, PAD_CODE);
    value_ptr = (uint8_t*) PyArray_GETPTR2(x2, 1, 4);
    EXPECT_EQ(*value_ptr, PAD_CODE);
    value_ptr = (uint8_t*) PyArray_GETPTR2(x2, 1, 5);
    EXPECT_EQ(*value_ptr, PAD_CODE);
    value_ptr = (uint8_t*) PyArray_GETPTR2(x2, 1, 6);
    EXPECT_EQ(*value_ptr, PAD_CODE);
    value_ptr = (uint8_t*) PyArray_GETPTR2(x2, 1, 7);
    EXPECT_EQ(*value_ptr, PAD_CODE);
    value_ptr = (uint8_t*) PyArray_GETPTR2(x2, 1, 8);
    EXPECT_EQ(*value_ptr, PAD_CODE);
    value_ptr = (uint8_t*) PyArray_GETPTR2(x2, 1, 9);
    EXPECT_EQ(*value_ptr, PAD_CODE);
    value_ptr = (uint8_t*) PyArray_GETPTR2(x2, 1, 10);
    EXPECT_EQ(*value_ptr, PAD_CODE);
    value_ptr = (uint8_t*) PyArray_GETPTR2(x2, 1, 11);
    EXPECT_EQ(*value_ptr, PAD_CODE); 

    // check second matrix fourth row
    value_ptr = (uint8_t*) PyArray_GETPTR2(x2, 3, 0);
    EXPECT_EQ(*value_ptr, ENCODER['_']);
    value_ptr = (uint8_t*) PyArray_GETPTR2(x2, 3, 1);
    EXPECT_EQ(*value_ptr, ENCODER['G']);
    value_ptr = (uint8_t*) PyArray_GETPTR2(x2, 3, 2);
    EXPECT_EQ(*value_ptr, ENCODER['G']);
    value_ptr = (uint8_t*) PyArray_GETPTR2(x2, 3, 3);
    EXPECT_EQ(*value_ptr, PAD_CODE);
    value_ptr = (uint8_t*) PyArray_GETPTR2(x2, 3, 4);
    EXPECT_EQ(*value_ptr, PAD_CODE);
    value_ptr = (uint8_t*) PyArray_GETPTR2(x2, 3, 5);
    EXPECT_EQ(*value_ptr, PAD_CODE);
    value_ptr = (uint8_t*) PyArray_GETPTR2(x2, 3, 6);
    EXPECT_EQ(*value_ptr, PAD_CODE);
    value_ptr = (uint8_t*) PyArray_GETPTR2(x2, 3, 7);
    EXPECT_EQ(*value_ptr, PAD_CODE);
    value_ptr = (uint8_t*) PyArray_GETPTR2(x2, 3, 8);
    EXPECT_EQ(*value_ptr, PAD_CODE);
    value_ptr = (uint8_t*) PyArray_GETPTR2(x2, 3, 9);
    EXPECT_EQ(*value_ptr, PAD_CODE);
    value_ptr = (uint8_t*) PyArray_GETPTR2(x2, 3, 10);
    EXPECT_EQ(*value_ptr, PAD_CODE);
    value_ptr = (uint8_t*) PyArray_GETPTR2(x2, 3, 11);
    EXPECT_EQ(*value_ptr, PAD_CODE); 
}

TEST(Converter, choose_segments)
{
    srand(time(NULL));
    PyInit_align_reads_gen();  
    std::string t  = "TAGGCATACAGG";
    std::string q1 = "TAGTG";

    Futures<AlignmentSegment> futures(pool, 1);
    futures.add_inputs(get_alignment_segment, q1, 0, 5, t, 0, 4, EDLIB_MODE_INFIX, EDLIB_TASK_PATH);
    std::vector<AlignmentSegment> segments = futures.get();
    
    MultiAlignment m_align {std::move(t), std::move(segments)};
    
    AlignmentConverter converter {m_align, 1000, 12};
    EXPECT_EQ(converter.segments_in_windows.size(), 2);

    auto result = converter.choose_segments(1000, false);
    EXPECT_EQ(result.size(), 2);
    EXPECT_EQ(result[0].size(), 1000);
    EXPECT_EQ(result[1].size(), 0);
    for (auto i : result[0])
    {
        EXPECT_EQ(i, 1);
    }
    result = converter.choose_segments(1000, true);
    EXPECT_EQ(result.size(), 2);
    EXPECT_EQ(result[0].size(), 1000);
    EXPECT_EQ(result[1].size(), 0);
    for (auto i : result[0])
    {
        EXPECT_EQ(i, 1);
    }
    result = converter.choose_segments(0, false);
    EXPECT_EQ(result.size(), 2);
    EXPECT_EQ(result[0].size(), 1000);
    EXPECT_EQ(result[1].size(), 0);
    for (auto i : result[0])
    {
        EXPECT_EQ(i, 0);
    }
    
    result = converter.choose_segments(123, false);
    EXPECT_EQ(result.size(), 2);
    EXPECT_EQ(result[0].size(), 1000);
    EXPECT_EQ(result[1].size(), 0);
    auto j = 0;
    for (auto i : result[0])
    {
        if (j < 123) 
        {
            EXPECT_EQ(i, 1);
        } else 
        {
            EXPECT_EQ(i, 0);
        }
        j++;
    }

}


TEST(Converter, fill_row_from_target)
{
    PyInit_align_reads_gen();  
    std::string t  = "TAGGCATACAGG";
    std::string q1 = "TAGTG";
    std::string q2 = "CAACATGGAAAAAAAAAAAAAAAA";
    std::string q3 = "TGGCATATCA";
    // TAG_|GCAT|A_CA|_GG

    Futures<AlignmentSegment> futures(pool, 3);
    futures.add_inputs(get_alignment_segment, q1, 0, 5, t, 0, 4, EDLIB_MODE_INFIX, EDLIB_TASK_PATH);
    futures.add_inputs(get_alignment_segment, q2, 0, 24, t, 4, 8, EDLIB_MODE_INFIX, EDLIB_TASK_PATH);
    futures.add_inputs(get_alignment_segment, q3, 0, 10, t, 0, 12, EDLIB_MODE_INFIX, EDLIB_TASK_PATH);
    futures.add_inputs(get_alignment_segment, t, 0, 12, t, 0, 12, EDLIB_MODE_INFIX, EDLIB_TASK_PATH);
    std::vector<AlignmentSegment> segments = futures.get();
    EXPECT_EQ(segments[1].get_aligned_chars(), "CA_ACAGG");
    EXPECT_EQ(segments[1].get_ins_segment_at(5), "T");
    EXPECT_EQ(segments[1].get_ins_segment_at(7), "AAAAAAAAAAAAAAAA");
    
    MultiAlignment m_align {std::move(t), std::move(segments)};
    
    AlignmentConverter converter {m_align, 4, 12};
    EXPECT_EQ(converter.segments_in_windows.size(), 2);;
     
    uint8_t* value_ptr;
    npy_intp dims[2];
    dims[0] = 4;
    dims[1] = 12;
    auto x = PyArray_SimpleNew(2, dims, NPY_UINT8); 
    auto x2 = PyArray_SimpleNew(2, dims, NPY_UINT8);
    
    converter.fill_row_from_target(x, 0, 0);

  
    converter.fill_row_from_target(x2, 1, 0);
    
    // check first matrix first row
    value_ptr = (uint8_t*) PyArray_GETPTR2(x, 0, 0);
    EXPECT_EQ(*value_ptr, ENCODER['T']);
    value_ptr = (uint8_t*) PyArray_GETPTR2(x, 0, 1);
    EXPECT_EQ(*value_ptr, ENCODER['A']);
    value_ptr = (uint8_t*) PyArray_GETPTR2(x, 0, 2);
    EXPECT_EQ(*value_ptr, ENCODER['G']);
    value_ptr = (uint8_t*) PyArray_GETPTR2(x, 0, 3);
    EXPECT_EQ(*value_ptr, ENCODER['_']);
    value_ptr = (uint8_t*) PyArray_GETPTR2(x, 0, 4);
    EXPECT_EQ(*value_ptr, ENCODER['G']);
    value_ptr = (uint8_t*) PyArray_GETPTR2(x, 0, 5);
    EXPECT_EQ(*value_ptr, ENCODER['C']);
    value_ptr = (uint8_t*) PyArray_GETPTR2(x, 0, 6);
    EXPECT_EQ(*value_ptr, ENCODER['A']);
    value_ptr = (uint8_t*) PyArray_GETPTR2(x, 0, 7);
    EXPECT_EQ(*value_ptr, ENCODER['T']);
    value_ptr = (uint8_t*) PyArray_GETPTR2(x, 0, 8);
    EXPECT_EQ(*value_ptr, ENCODER['A']);
    value_ptr = (uint8_t*) PyArray_GETPTR2(x, 0, 9);
    EXPECT_EQ(*value_ptr, ENCODER['_']);
    value_ptr = (uint8_t*) PyArray_GETPTR2(x, 0, 10);
    EXPECT_EQ(*value_ptr, ENCODER['C']);
    value_ptr = (uint8_t*) PyArray_GETPTR2(x, 0, 11);
    EXPECT_EQ(*value_ptr, ENCODER['A']);
    
    // check second matrix second row
    value_ptr = (uint8_t*) PyArray_GETPTR2(x2, 0, 0);
    EXPECT_EQ(*value_ptr, ENCODER['_']);
    value_ptr = (uint8_t*) PyArray_GETPTR2(x2, 0, 1);
    EXPECT_EQ(*value_ptr, ENCODER['G']);
    value_ptr = (uint8_t*) PyArray_GETPTR2(x2, 0, 2);
    EXPECT_EQ(*value_ptr, ENCODER['G']);
    value_ptr = (uint8_t*) PyArray_GETPTR2(x2, 0, 3);
    EXPECT_EQ(*value_ptr, PAD_CODE);
    value_ptr = (uint8_t*) PyArray_GETPTR2(x2, 0, 4);
    EXPECT_EQ(*value_ptr, PAD_CODE);
    value_ptr = (uint8_t*) PyArray_GETPTR2(x2, 0, 5);
    EXPECT_EQ(*value_ptr, PAD_CODE);
    value_ptr = (uint8_t*) PyArray_GETPTR2(x2, 0, 6);
    EXPECT_EQ(*value_ptr, PAD_CODE);
    value_ptr = (uint8_t*) PyArray_GETPTR2(x2, 0, 7);
    EXPECT_EQ(*value_ptr, PAD_CODE);
    value_ptr = (uint8_t*) PyArray_GETPTR2(x2, 0, 8);
    EXPECT_EQ(*value_ptr, PAD_CODE);
    value_ptr = (uint8_t*) PyArray_GETPTR2(x2, 0, 9);
    EXPECT_EQ(*value_ptr, PAD_CODE);
    value_ptr = (uint8_t*) PyArray_GETPTR2(x2, 0, 10);
    EXPECT_EQ(*value_ptr, PAD_CODE);
    value_ptr = (uint8_t*) PyArray_GETPTR2(x2, 0, 11);
    EXPECT_EQ(*value_ptr, PAD_CODE);
}

TEST(Converter, produce_alignment_matrices)
{
    PyInit_align_reads_gen();  
    std::string t  = "TAGGCATACAGG";
    std::string q1 = "TAGTG";
    std::string q2 = "CAACATGGAAAAAAAAAAAAAAAA";
    std::string q3 = "TGGCATATCA";
    // TAG_|GCAT|A_CA|_GG

    Futures<AlignmentSegment> futures(pool, 3);
    futures.add_inputs(get_alignment_segment, q1, 0, 5, t, 0, 4, EDLIB_MODE_INFIX, EDLIB_TASK_PATH);
    futures.add_inputs(get_alignment_segment, q2, 0, 24, t, 4, 8, EDLIB_MODE_INFIX, EDLIB_TASK_PATH);
    futures.add_inputs(get_alignment_segment, q3, 0, 10, t, 0, 12, EDLIB_MODE_INFIX, EDLIB_TASK_PATH);
    futures.add_inputs(get_alignment_segment, t, 0, 12, t, 0, 12, EDLIB_MODE_INFIX, EDLIB_TASK_PATH);
    std::vector<AlignmentSegment> segments = futures.get();
    EXPECT_EQ(segments[1].get_aligned_chars(), "CA_ACAGG");
    EXPECT_EQ(segments[1].get_ins_segment_at(5), "T");
    EXPECT_EQ(segments[1].get_ins_segment_at(7), "AAAAAAAAAAAAAAAA");
    
    MultiAlignment m_align {std::move(t), std::move(segments)};
    
    AlignmentConverter converter {m_align, 4, 12};
    EXPECT_EQ(converter.segments_in_windows.size(), 2);;
     
    std::vector<std::vector<std::uint32_t>> chosen;
    chosen.resize(2);
    chosen[0] = {0, 1, 2, 3};
    chosen[1] = {1, 4, 4, 3};
    std::shared_ptr<thread_pool::ThreadPool> null_pool = nullptr;
    auto Xs = converter.produce_alignment_matrices(chosen, null_pool);
    EXPECT_EQ(Xs.size(), 2);
    auto x = Xs[0];
    auto x2 = Xs[1];
    uint8_t* value_ptr;
    // check first matrix first row
    value_ptr = (uint8_t*) PyArray_GETPTR2(x, 0, 0);
    EXPECT_EQ(*value_ptr, ENCODER['T']);
    value_ptr = (uint8_t*) PyArray_GETPTR2(x, 0, 1);
    EXPECT_EQ(*value_ptr, ENCODER['A']);
    value_ptr = (uint8_t*) PyArray_GETPTR2(x, 0, 2);
    EXPECT_EQ(*value_ptr, ENCODER['G']);
    value_ptr = (uint8_t*) PyArray_GETPTR2(x, 0, 3);
    EXPECT_EQ(*value_ptr, ENCODER['T']);
    value_ptr = (uint8_t*) PyArray_GETPTR2(x, 0, 4);
    EXPECT_EQ(*value_ptr, ENCODER['G']);
    value_ptr = (uint8_t*) PyArray_GETPTR2(x, 0, 5);
    EXPECT_EQ(*value_ptr, PAD_CODE);
    value_ptr = (uint8_t*) PyArray_GETPTR2(x, 0, 6);
    EXPECT_EQ(*value_ptr, PAD_CODE);
    value_ptr = (uint8_t*) PyArray_GETPTR2(x, 0, 7);
    EXPECT_EQ(*value_ptr, PAD_CODE);
    value_ptr = (uint8_t*) PyArray_GETPTR2(x, 0, 8);
    EXPECT_EQ(*value_ptr, PAD_CODE);
    value_ptr = (uint8_t*) PyArray_GETPTR2(x, 0, 9);
    EXPECT_EQ(*value_ptr, PAD_CODE);
    value_ptr = (uint8_t*) PyArray_GETPTR2(x, 0, 10);
    EXPECT_EQ(*value_ptr, PAD_CODE);
    value_ptr = (uint8_t*) PyArray_GETPTR2(x, 0, 11);
    EXPECT_EQ(*value_ptr, PAD_CODE);

    // check first matrix fourth row
    value_ptr = (uint8_t*) PyArray_GETPTR2(x, 3, 0);
    EXPECT_EQ(*value_ptr, ENCODER['T']);
    value_ptr = (uint8_t*) PyArray_GETPTR2(x, 3, 1);
    EXPECT_EQ(*value_ptr, ENCODER['A']);
    value_ptr = (uint8_t*) PyArray_GETPTR2(x, 3, 2);
    EXPECT_EQ(*value_ptr, ENCODER['G']);
    value_ptr = (uint8_t*) PyArray_GETPTR2(x, 3, 3);
    EXPECT_EQ(*value_ptr, ENCODER['_']);
    value_ptr = (uint8_t*) PyArray_GETPTR2(x, 3, 4);
    EXPECT_EQ(*value_ptr, ENCODER['G']);
    value_ptr = (uint8_t*) PyArray_GETPTR2(x, 3, 5);
    EXPECT_EQ(*value_ptr, ENCODER['C']);
    value_ptr = (uint8_t*) PyArray_GETPTR2(x, 3, 6);
    EXPECT_EQ(*value_ptr, ENCODER['A']);
    value_ptr = (uint8_t*) PyArray_GETPTR2(x, 3, 7);
    EXPECT_EQ(*value_ptr, ENCODER['T']);
    value_ptr = (uint8_t*) PyArray_GETPTR2(x, 3, 8);
    EXPECT_EQ(*value_ptr, ENCODER['A']);
    value_ptr = (uint8_t*) PyArray_GETPTR2(x, 3, 9);
    EXPECT_EQ(*value_ptr, ENCODER['_']);
    value_ptr = (uint8_t*) PyArray_GETPTR2(x, 3, 10);
    EXPECT_EQ(*value_ptr, ENCODER['C']);
    value_ptr = (uint8_t*) PyArray_GETPTR2(x, 3, 11);
    EXPECT_EQ(*value_ptr, ENCODER['A']);
    
    // check second matrix second row
    value_ptr = (uint8_t*) PyArray_GETPTR2(x2, 1, 0);
    EXPECT_EQ(*value_ptr, ENCODER['_']);
    value_ptr = (uint8_t*) PyArray_GETPTR2(x2, 1, 1);
    EXPECT_EQ(*value_ptr, ENCODER['G']);
    value_ptr = (uint8_t*) PyArray_GETPTR2(x2, 1, 2);
    EXPECT_EQ(*value_ptr, ENCODER['G']);
    value_ptr = (uint8_t*) PyArray_GETPTR2(x2, 1, 3);
    EXPECT_EQ(*value_ptr, PAD_CODE);
    value_ptr = (uint8_t*) PyArray_GETPTR2(x2, 1, 4);
    EXPECT_EQ(*value_ptr, PAD_CODE);
    value_ptr = (uint8_t*) PyArray_GETPTR2(x2, 1, 5);
    EXPECT_EQ(*value_ptr, PAD_CODE);
    value_ptr = (uint8_t*) PyArray_GETPTR2(x2, 1, 6);
    EXPECT_EQ(*value_ptr, PAD_CODE);
    value_ptr = (uint8_t*) PyArray_GETPTR2(x2, 1, 7);
    EXPECT_EQ(*value_ptr, PAD_CODE);
    value_ptr = (uint8_t*) PyArray_GETPTR2(x2, 1, 8);
    EXPECT_EQ(*value_ptr, PAD_CODE);
    value_ptr = (uint8_t*) PyArray_GETPTR2(x2, 1, 9);
    EXPECT_EQ(*value_ptr, PAD_CODE);
    value_ptr = (uint8_t*) PyArray_GETPTR2(x2, 1, 10);
    EXPECT_EQ(*value_ptr, PAD_CODE);
    value_ptr = (uint8_t*) PyArray_GETPTR2(x2, 1, 11);
    EXPECT_EQ(*value_ptr, PAD_CODE); 
}

TEST(Converter, produce_alignment_matrices_parallel)
{
    PyInit_align_reads_gen();  
    std::string t  = "TAGGCATACAGG";
    std::string q1 = "TAGTG";
    std::string q2 = "CAACATGGAAAAAAAAAAAAAAAA";
    std::string q3 = "TGGCATATCA";
    // TAG_|GCAT|A_CA|_GG

    Futures<AlignmentSegment> futures(pool, 3);
    futures.add_inputs(get_alignment_segment, q1, 0, 5, t, 0, 4, EDLIB_MODE_INFIX, EDLIB_TASK_PATH);
    futures.add_inputs(get_alignment_segment, q2, 0, 24, t, 4, 8, EDLIB_MODE_INFIX, EDLIB_TASK_PATH);
    futures.add_inputs(get_alignment_segment, q3, 0, 10, t, 0, 12, EDLIB_MODE_INFIX, EDLIB_TASK_PATH);
    futures.add_inputs(get_alignment_segment, t, 0, 12, t, 0, 12, EDLIB_MODE_INFIX, EDLIB_TASK_PATH);
    std::vector<AlignmentSegment> segments = futures.get();
    EXPECT_EQ(segments[1].get_aligned_chars(), "CA_ACAGG");
    EXPECT_EQ(segments[1].get_ins_segment_at(5), "T");
    EXPECT_EQ(segments[1].get_ins_segment_at(7), "AAAAAAAAAAAAAAAA");
    
    MultiAlignment m_align {std::move(t), std::move(segments)};
    
    AlignmentConverter converter {m_align, 4, 12};
    EXPECT_EQ(converter.segments_in_windows.size(), 2);;
     
    std::vector<std::vector<std::uint32_t>> chosen;
    chosen.resize(2);
    chosen[0] = {0, 1, 2, 3};
    chosen[1] = {1, 4, 4, 3};
    //std::shared_ptr<thread_pool::ThreadPool> n_pool = nullptr;
    auto Xs = converter.produce_alignment_matrices(chosen, pool);
    EXPECT_EQ(Xs.size(), 2);
    auto x = Xs[0];
    auto x2 = Xs[1];
    uint8_t* value_ptr;
    // check first matrix first row
    value_ptr = (uint8_t*) PyArray_GETPTR2(x, 0, 0);
    EXPECT_EQ(*value_ptr, ENCODER['T']);
    value_ptr = (uint8_t*) PyArray_GETPTR2(x, 0, 1);
    EXPECT_EQ(*value_ptr, ENCODER['A']);
    value_ptr = (uint8_t*) PyArray_GETPTR2(x, 0, 2);
    EXPECT_EQ(*value_ptr, ENCODER['G']);
    value_ptr = (uint8_t*) PyArray_GETPTR2(x, 0, 3);
    EXPECT_EQ(*value_ptr, ENCODER['T']);
    value_ptr = (uint8_t*) PyArray_GETPTR2(x, 0, 4);
    EXPECT_EQ(*value_ptr, ENCODER['G']);
    value_ptr = (uint8_t*) PyArray_GETPTR2(x, 0, 5);
    EXPECT_EQ(*value_ptr, PAD_CODE);
    value_ptr = (uint8_t*) PyArray_GETPTR2(x, 0, 6);
    EXPECT_EQ(*value_ptr, PAD_CODE);
    value_ptr = (uint8_t*) PyArray_GETPTR2(x, 0, 7);
    EXPECT_EQ(*value_ptr, PAD_CODE);
    value_ptr = (uint8_t*) PyArray_GETPTR2(x, 0, 8);
    EXPECT_EQ(*value_ptr, PAD_CODE);
    value_ptr = (uint8_t*) PyArray_GETPTR2(x, 0, 9);
    EXPECT_EQ(*value_ptr, PAD_CODE);
    value_ptr = (uint8_t*) PyArray_GETPTR2(x, 0, 10);
    EXPECT_EQ(*value_ptr, PAD_CODE);
    value_ptr = (uint8_t*) PyArray_GETPTR2(x, 0, 11);
    EXPECT_EQ(*value_ptr, PAD_CODE);
    
    // check first matrix fourth row
    value_ptr = (uint8_t*) PyArray_GETPTR2(x, 3, 0);
    EXPECT_EQ(*value_ptr, ENCODER['T']);
    value_ptr = (uint8_t*) PyArray_GETPTR2(x, 3, 1);
    EXPECT_EQ(*value_ptr, ENCODER['A']);
    value_ptr = (uint8_t*) PyArray_GETPTR2(x, 3, 2);
    EXPECT_EQ(*value_ptr, ENCODER['G']);
    value_ptr = (uint8_t*) PyArray_GETPTR2(x, 3, 3);
    EXPECT_EQ(*value_ptr, ENCODER['_']);
    value_ptr = (uint8_t*) PyArray_GETPTR2(x, 3, 4);
    EXPECT_EQ(*value_ptr, ENCODER['G']);
    value_ptr = (uint8_t*) PyArray_GETPTR2(x, 3, 5);
    EXPECT_EQ(*value_ptr, ENCODER['C']);
    value_ptr = (uint8_t*) PyArray_GETPTR2(x, 3, 6);
    EXPECT_EQ(*value_ptr, ENCODER['A']);
    value_ptr = (uint8_t*) PyArray_GETPTR2(x, 3, 7);
    EXPECT_EQ(*value_ptr, ENCODER['T']);
    value_ptr = (uint8_t*) PyArray_GETPTR2(x, 3, 8);
    EXPECT_EQ(*value_ptr, ENCODER['A']);
    value_ptr = (uint8_t*) PyArray_GETPTR2(x, 3, 9);
    EXPECT_EQ(*value_ptr, ENCODER['_']);
    value_ptr = (uint8_t*) PyArray_GETPTR2(x, 3, 10);
    EXPECT_EQ(*value_ptr, ENCODER['C']);
    value_ptr = (uint8_t*) PyArray_GETPTR2(x, 3, 11);
    EXPECT_EQ(*value_ptr, ENCODER['A']);
    
    // check second matrix second row
    value_ptr = (uint8_t*) PyArray_GETPTR2(x2, 1, 0);
    EXPECT_EQ(*value_ptr, ENCODER['_']);
    value_ptr = (uint8_t*) PyArray_GETPTR2(x2, 1, 1);
    EXPECT_EQ(*value_ptr, ENCODER['G']);
    value_ptr = (uint8_t*) PyArray_GETPTR2(x2, 1, 2);
    EXPECT_EQ(*value_ptr, ENCODER['G']);
    value_ptr = (uint8_t*) PyArray_GETPTR2(x2, 1, 3);
    EXPECT_EQ(*value_ptr, PAD_CODE);
    value_ptr = (uint8_t*) PyArray_GETPTR2(x2, 1, 4);
    EXPECT_EQ(*value_ptr, PAD_CODE);
    value_ptr = (uint8_t*) PyArray_GETPTR2(x2, 1, 5);
    EXPECT_EQ(*value_ptr, PAD_CODE);
    value_ptr = (uint8_t*) PyArray_GETPTR2(x2, 1, 6);
    EXPECT_EQ(*value_ptr, PAD_CODE);
    value_ptr = (uint8_t*) PyArray_GETPTR2(x2, 1, 7);
    EXPECT_EQ(*value_ptr, PAD_CODE);
    value_ptr = (uint8_t*) PyArray_GETPTR2(x2, 1, 8);
    EXPECT_EQ(*value_ptr, PAD_CODE);
    value_ptr = (uint8_t*) PyArray_GETPTR2(x2, 1, 9);
    EXPECT_EQ(*value_ptr, PAD_CODE);
    value_ptr = (uint8_t*) PyArray_GETPTR2(x2, 1, 10);
    EXPECT_EQ(*value_ptr, PAD_CODE);
    value_ptr = (uint8_t*) PyArray_GETPTR2(x2, 1, 11);
    EXPECT_EQ(*value_ptr, PAD_CODE); 
}

TEST(Converter, trim_window)
{
    PyInit_align_reads_gen();  
    std::string t  = "TAGGCATACAGG";
    std::string q1 = "TAGTG";
    std::string q2 = "CAACATGGAAAAAAAAAAAAAAAA";
    std::string q3 = "TGGCATATCA";
    // TAG_|GCAT|A_CA|_GG

    Futures<AlignmentSegment> futures(pool, 3);
    futures.add_inputs(get_alignment_segment, q1, 0, 5, t, 0, 4, EDLIB_MODE_INFIX, EDLIB_TASK_PATH);
    futures.add_inputs(get_alignment_segment, q2, 0, 24, t, 4, 8, EDLIB_MODE_INFIX, EDLIB_TASK_PATH);
    futures.add_inputs(get_alignment_segment, q3, 0, 10, t, 0, 12, EDLIB_MODE_INFIX, EDLIB_TASK_PATH);
    futures.add_inputs(get_alignment_segment, t, 0, 12, t, 0, 12, EDLIB_MODE_INFIX, EDLIB_TASK_PATH);
    std::vector<AlignmentSegment> segments = futures.get();
    EXPECT_EQ(segments[1].get_aligned_chars(), "CA_ACAGG");
    EXPECT_EQ(segments[1].get_ins_segment_at(5), "T");
    EXPECT_EQ(segments[1].get_ins_segment_at(7), "AAAAAAAAAAAAAAAA");
    EXPECT_EQ(segments[0].get_ins_segment_at(2), "T");

    std::string tr = "GGCATACA";
    futures.add_inputs(get_alignment_segment, tr, 0, 8, t, 2, 8, EDLIB_MODE_INFIX, EDLIB_TASK_PATH);
    auto truth = futures.get();

    
    MultiAlignment m_align {std::move(t), std::move(segments), std::move(truth)};
    
    AlignmentConverter converter {m_align, 4, 12};
    EXPECT_EQ(converter.segments_in_windows.size(), 1);;
    
    std::vector<std::vector<std::uint32_t>> chosen;
    chosen.resize(1);
    chosen[0] = {0, 1, 2, 3};
    std::shared_ptr<thread_pool::ThreadPool> n_pool = nullptr;
    auto Xs = converter.produce_alignment_matrices(chosen, n_pool);
    EXPECT_EQ(Xs.size(), 1);
    
    auto x = Xs[0];
    uint8_t* value_ptr;
    // check first matrix first row
    value_ptr = (uint8_t*) PyArray_GETPTR2(x, 0, 0);
    EXPECT_EQ(*value_ptr, ENCODER['G']);
    
    value_ptr = (uint8_t*) PyArray_GETPTR2(x, 0, 1);
    EXPECT_EQ(*value_ptr, ENCODER['T']);

    value_ptr = (uint8_t*) PyArray_GETPTR2(x, 0, 2);
    EXPECT_EQ(*value_ptr, ENCODER['G']);
    value_ptr = (uint8_t*) PyArray_GETPTR2(x, 0, 3);
    EXPECT_EQ(*value_ptr, PAD_CODE);
    value_ptr = (uint8_t*) PyArray_GETPTR2(x, 0, 4);
    EXPECT_EQ(*value_ptr, PAD_CODE);
    value_ptr = (uint8_t*) PyArray_GETPTR2(x, 0, 5);
    EXPECT_EQ(*value_ptr, PAD_CODE);
    value_ptr = (uint8_t*) PyArray_GETPTR2(x, 0, 6);
    EXPECT_EQ(*value_ptr, PAD_CODE);
    value_ptr = (uint8_t*) PyArray_GETPTR2(x, 0, 7);
    EXPECT_EQ(*value_ptr, PAD_CODE);
    value_ptr = (uint8_t*) PyArray_GETPTR2(x, 0, 8);
    EXPECT_EQ(*value_ptr, PAD_CODE);
    value_ptr = (uint8_t*) PyArray_GETPTR2(x, 0, 9);
    EXPECT_EQ(*value_ptr, PAD_CODE);
    value_ptr = (uint8_t*) PyArray_GETPTR2(x, 0, 10);
    EXPECT_EQ(*value_ptr, PAD_CODE);
    value_ptr = (uint8_t*) PyArray_GETPTR2(x, 0, 11);
    EXPECT_EQ(*value_ptr, PAD_CODE);
    
    // check first matrix fourth row
    value_ptr = (uint8_t*) PyArray_GETPTR2(x, 3, 0);
    EXPECT_EQ(*value_ptr, ENCODER['G']);
    value_ptr = (uint8_t*) PyArray_GETPTR2(x, 3, 1);
    EXPECT_EQ(*value_ptr, ENCODER['_']);
    value_ptr = (uint8_t*) PyArray_GETPTR2(x, 3, 2);
    EXPECT_EQ(*value_ptr, ENCODER['G']);
    value_ptr = (uint8_t*) PyArray_GETPTR2(x, 3, 3);
    EXPECT_EQ(*value_ptr, ENCODER['C']);
    value_ptr = (uint8_t*) PyArray_GETPTR2(x, 3, 4);
    EXPECT_EQ(*value_ptr, ENCODER['A']);
    value_ptr = (uint8_t*) PyArray_GETPTR2(x, 3, 5);
    EXPECT_EQ(*value_ptr, ENCODER['T']);
    value_ptr = (uint8_t*) PyArray_GETPTR2(x, 3, 6);
    EXPECT_EQ(*value_ptr, ENCODER['A']);
    value_ptr = (uint8_t*) PyArray_GETPTR2(x, 3, 7);
    EXPECT_EQ(*value_ptr, ENCODER['_']);
    value_ptr = (uint8_t*) PyArray_GETPTR2(x, 3, 8);
    EXPECT_EQ(*value_ptr, ENCODER['C']);
    value_ptr = (uint8_t*) PyArray_GETPTR2(x, 3, 9);
    EXPECT_EQ(*value_ptr, ENCODER['A']);
    value_ptr = (uint8_t*) PyArray_GETPTR2(x, 3, 10);
    EXPECT_EQ(*value_ptr, PAD_CODE);
    value_ptr = (uint8_t*) PyArray_GETPTR2(x, 3, 11);
    EXPECT_EQ(*value_ptr, PAD_CODE);
    
    // check first matrix second row
    value_ptr = (uint8_t*) PyArray_GETPTR2(x, 1, 0);
    EXPECT_EQ(*value_ptr, PAD_CODE);
    value_ptr = (uint8_t*) PyArray_GETPTR2(x, 1, 1);
    EXPECT_EQ(*value_ptr, PAD_CODE);
    value_ptr = (uint8_t*) PyArray_GETPTR2(x, 1, 2);
    EXPECT_EQ(*value_ptr, PAD_CODE);
    value_ptr = (uint8_t*) PyArray_GETPTR2(x, 1, 3);
    EXPECT_EQ(*value_ptr, ENCODER['C']);
    value_ptr = (uint8_t*) PyArray_GETPTR2(x, 1, 4);
    EXPECT_EQ(*value_ptr, ENCODER['A']);
    value_ptr = (uint8_t*) PyArray_GETPTR2(x, 1, 5);
    EXPECT_EQ(*value_ptr, ENCODER['_']);
    value_ptr = (uint8_t*) PyArray_GETPTR2(x, 1, 6);
    EXPECT_EQ(*value_ptr, ENCODER['A']);
    value_ptr = (uint8_t*) PyArray_GETPTR2(x, 1, 7);
    EXPECT_EQ(*value_ptr, ENCODER['_']);
    value_ptr = (uint8_t*) PyArray_GETPTR2(x, 1, 8);
    EXPECT_EQ(*value_ptr, ENCODER['C']);
    value_ptr = (uint8_t*) PyArray_GETPTR2(x, 1, 9);
    EXPECT_EQ(*value_ptr, ENCODER['A']);
    value_ptr = (uint8_t*) PyArray_GETPTR2(x, 1, 10);
    EXPECT_EQ(*value_ptr, PAD_CODE);
    value_ptr = (uint8_t*) PyArray_GETPTR2(x, 1, 11);
    EXPECT_EQ(*value_ptr, PAD_CODE); 
}

TEST(Converter, produce_truth_matrices)
{
    PyInit_align_reads_gen();  
    std::string t  = "TAGGCATACAGG";
    std::string q1 = "TAGTG";
    std::string q2 = "CAACATGGAAAAAAAAAAAAAAAA";
    std::string q3 = "TGGCATATCA";
    // TAG_|GCAT|A_CA|_GG

    Futures<AlignmentSegment> futures(pool, 3);
    futures.add_inputs(get_alignment_segment, q1, 0, 5, t, 0, 4, EDLIB_MODE_INFIX, EDLIB_TASK_PATH);
    futures.add_inputs(get_alignment_segment, q2, 0, 24, t, 4, 8, EDLIB_MODE_INFIX, EDLIB_TASK_PATH);
    futures.add_inputs(get_alignment_segment, q3, 0, 10, t, 0, 12, EDLIB_MODE_INFIX, EDLIB_TASK_PATH);
    futures.add_inputs(get_alignment_segment, t, 0, 12, t, 0, 12, EDLIB_MODE_INFIX, EDLIB_TASK_PATH);
    std::vector<AlignmentSegment> segments = futures.get();
    EXPECT_EQ(segments[1].get_aligned_chars(), "CA_ACAGG");
    EXPECT_EQ(segments[1].get_ins_segment_at(5), "T");
    EXPECT_EQ(segments[1].get_ins_segment_at(7), "AAAAAAAAAAAAAAAA");
    EXPECT_EQ(segments[0].get_ins_segment_at(2), "T");

    std::string tr = "GGCATACA";
    futures.add_inputs(get_alignment_segment, tr, 0, 8, t, 2, 8, EDLIB_MODE_INFIX, EDLIB_TASK_PATH);
    auto truth = futures.get();

    
    MultiAlignment m_align {std::move(t), std::move(segments), std::move(truth)};
    Info info;
    info.hap_of_target = 0;
    info.hap_of_aligned.resize(3, 0);
    AlignmentConverter converter {m_align, 4, 12, info};
    EXPECT_EQ(converter.segments_in_windows.size(), 1);;
    
    std::vector<std::vector<std::uint32_t>> chosen;
    chosen.resize(1);
    chosen[0] = {0, 1, 2, 3};
    std::shared_ptr<thread_pool::ThreadPool> n_pool = nullptr;
    auto Xs = converter.produce_truth_matrices(chosen, n_pool);
    EXPECT_EQ(Xs.size(), 1);
    
    auto x = Xs[0];
    uint8_t* value_ptr;
    
    // check first matrix first row
    value_ptr = (uint8_t*) PyArray_GETPTR2(x, 0, 0);
    EXPECT_EQ(*value_ptr, ENCODER['G']);
    
    value_ptr = (uint8_t*) PyArray_GETPTR2(x, 0, 1);
    EXPECT_EQ(*value_ptr, ENCODER['_']);

    value_ptr = (uint8_t*) PyArray_GETPTR2(x, 0, 2);
    EXPECT_EQ(*value_ptr, ENCODER['G']);
    value_ptr = (uint8_t*) PyArray_GETPTR2(x, 0, 3);
    EXPECT_EQ(*value_ptr, ENCODER['C']);
    value_ptr = (uint8_t*) PyArray_GETPTR2(x, 0, 4);
    EXPECT_EQ(*value_ptr, ENCODER['A']);
    value_ptr = (uint8_t*) PyArray_GETPTR2(x, 0, 5);
    EXPECT_EQ(*value_ptr, ENCODER['T']);
    
    value_ptr = (uint8_t*) PyArray_GETPTR2(x, 0, 6);
    EXPECT_EQ(*value_ptr, ENCODER['A']);
    value_ptr = (uint8_t*) PyArray_GETPTR2(x, 0, 7);
    EXPECT_EQ(*value_ptr, ENCODER['_']);
    value_ptr = (uint8_t*) PyArray_GETPTR2(x, 0, 8);
    EXPECT_EQ(*value_ptr, ENCODER['C']);
    value_ptr = (uint8_t*) PyArray_GETPTR2(x, 0, 9);
    EXPECT_EQ(*value_ptr, ENCODER['A']);
    value_ptr = (uint8_t*) PyArray_GETPTR2(x, 0, 10);
    EXPECT_EQ(*value_ptr, PAD_CODE);
    value_ptr = (uint8_t*) PyArray_GETPTR2(x, 0, 11);
    EXPECT_EQ(*value_ptr, PAD_CODE);
    
    // check first matrix fourth row
    value_ptr = (uint8_t*) PyArray_GETPTR2(x, 3, 0);
    EXPECT_EQ(*value_ptr, ENCODER['G']);
    
    value_ptr = (uint8_t*) PyArray_GETPTR2(x, 3, 1);
    EXPECT_EQ(*value_ptr, ENCODER['_']);

    value_ptr = (uint8_t*) PyArray_GETPTR2(x, 3, 2);
    EXPECT_EQ(*value_ptr, ENCODER['G']);
    value_ptr = (uint8_t*) PyArray_GETPTR2(x, 3, 3);
    EXPECT_EQ(*value_ptr, ENCODER['C']);
    value_ptr = (uint8_t*) PyArray_GETPTR2(x, 3, 4);
    EXPECT_EQ(*value_ptr, ENCODER['A']);
    value_ptr = (uint8_t*) PyArray_GETPTR2(x, 3, 5);
    EXPECT_EQ(*value_ptr, ENCODER['T']);
    
    value_ptr = (uint8_t*) PyArray_GETPTR2(x, 3, 6);
    EXPECT_EQ(*value_ptr, ENCODER['A']);
    value_ptr = (uint8_t*) PyArray_GETPTR2(x, 3, 7);
    EXPECT_EQ(*value_ptr, ENCODER['_']);
    value_ptr = (uint8_t*) PyArray_GETPTR2(x, 3, 8);
    EXPECT_EQ(*value_ptr, ENCODER['C']);
    value_ptr = (uint8_t*) PyArray_GETPTR2(x, 3, 9);
    EXPECT_EQ(*value_ptr, ENCODER['A']);
    value_ptr = (uint8_t*) PyArray_GETPTR2(x, 3, 10);
    EXPECT_EQ(*value_ptr, PAD_CODE);
    value_ptr = (uint8_t*) PyArray_GETPTR2(x, 3, 11);
    EXPECT_EQ(*value_ptr, PAD_CODE);
    
    // check first matrix second row
    value_ptr = (uint8_t*) PyArray_GETPTR2(x, 1, 0);
    EXPECT_EQ(*value_ptr, ENCODER['G']);
    
    value_ptr = (uint8_t*) PyArray_GETPTR2(x, 1, 1);
    EXPECT_EQ(*value_ptr, ENCODER['_']);

    value_ptr = (uint8_t*) PyArray_GETPTR2(x, 1, 2);
    EXPECT_EQ(*value_ptr, ENCODER['G']);
    value_ptr = (uint8_t*) PyArray_GETPTR2(x, 1, 3);
    EXPECT_EQ(*value_ptr, ENCODER['C']);
    value_ptr = (uint8_t*) PyArray_GETPTR2(x, 1, 4);
    EXPECT_EQ(*value_ptr, ENCODER['A']);
    value_ptr = (uint8_t*) PyArray_GETPTR2(x, 1, 5);
    EXPECT_EQ(*value_ptr, ENCODER['T']);
    
    value_ptr = (uint8_t*) PyArray_GETPTR2(x, 1, 6);
    EXPECT_EQ(*value_ptr, ENCODER['A']);
    value_ptr = (uint8_t*) PyArray_GETPTR2(x, 1, 7);
    EXPECT_EQ(*value_ptr, ENCODER['_']);
    value_ptr = (uint8_t*) PyArray_GETPTR2(x, 1, 8);
    EXPECT_EQ(*value_ptr, ENCODER['C']);
    value_ptr = (uint8_t*) PyArray_GETPTR2(x, 1, 9);
    EXPECT_EQ(*value_ptr, ENCODER['A']);
    value_ptr = (uint8_t*) PyArray_GETPTR2(x, 1, 10);
    EXPECT_EQ(*value_ptr, PAD_CODE);
    value_ptr = (uint8_t*) PyArray_GETPTR2(x, 1, 11);
    EXPECT_EQ(*value_ptr, PAD_CODE);
     
    
}