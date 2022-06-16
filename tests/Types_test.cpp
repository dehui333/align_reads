#include <gtest/gtest.h>

#include "align_reads/Aligner.hpp"
#define private public
#include "../src/Types.cpp"

TEST(Types_alignment_segment, global_all_match)
{
    std::string s1 = "ACGT";
    auto result = align_reads::get_edlib_result(s1.c_str(), s1.c_str(), 4, 4, EDLIB_MODE_GLOBAL, EDLIB_TASK_PATH);
    auto align_segment = align_reads::alignment_segment(s1, s1, 0, 0, result);
    EXPECT_EQ(align_segment.start_on_target, 0);
    EXPECT_EQ(align_segment.align_len_on_target, 4);
    EXPECT_EQ(align_segment.get_aligned_chars().compare(s1), 0);
    for (int i = -1; i < 4; i++)
    {
        EXPECT_EQ(align_segment.get_ins_segment_at(i).size(), 0);
    }
}

TEST(Types_alignment_segment, global_ins)
{
    std::string t = "GGGG";
    std::string q = "AGGGG";
    auto result = align_reads::get_edlib_result(q.c_str(), t.c_str(), q.size(), t.size(), EDLIB_MODE_GLOBAL, EDLIB_TASK_PATH);
    auto align_segment = align_reads::alignment_segment(q, t, 0, 0, result);
    EXPECT_EQ(align_segment.start_on_target, 0);
    EXPECT_EQ(align_segment.align_len_on_target, 4);
    EXPECT_EQ(align_segment.get_aligned_chars().compare(t), 0);

    EXPECT_EQ(align_segment.get_ins_segment_at(-1).size(), 1);
    int i = 0;
    for (; i < 4; i++)
    {
        EXPECT_EQ(align_segment.get_ins_segment_at(i).size(), 0);
    }
    EXPECT_EQ(align_segment.get_ins_segment_at(-1), "A");

    t = "GGGG";
    q = "AGAGGGAA";
    result = align_reads::get_edlib_result(q.c_str(), t.c_str(), q.size(), t.size(), EDLIB_MODE_GLOBAL, EDLIB_TASK_PATH);
    align_segment = align_reads::alignment_segment(q, t, 0, 0, result);
    EXPECT_EQ(align_segment.start_on_target, 0);
    EXPECT_EQ(align_segment.align_len_on_target, 4);
    EXPECT_EQ(align_segment.get_aligned_chars().compare(t), 0);

    EXPECT_EQ(align_segment.get_ins_segment_at(-1).size(), 1);
    EXPECT_EQ(align_segment.get_ins_segment_at(0).size(), 1);
    EXPECT_EQ(align_segment.get_ins_segment_at(3).size(), 2);
}

TEST(Types_alignment_segment, global_del)
{
    std::string t = "ACGT";
    std::string q = "AGT";
    auto result = align_reads::get_edlib_result(q.c_str(), t.c_str(), q.size(), t.size(), EDLIB_MODE_GLOBAL, EDLIB_TASK_PATH);
    auto align_segment = align_reads::alignment_segment(q, t, 0, 0, result);
    EXPECT_EQ(align_segment.start_on_target, 0);
    EXPECT_EQ(align_segment.align_len_on_target, 4);
    EXPECT_EQ(align_segment.get_aligned_chars(), "A_GT");

    int i = 0;
    for (i = -1; i < 4; i++)
    {
        EXPECT_EQ(align_segment.get_ins_segment_at(i).size(), 0);
    }

    t = "ACGT";
    q = "AT";
    result = align_reads::get_edlib_result(q.c_str(), t.c_str(), q.size(), t.size(), EDLIB_MODE_GLOBAL, EDLIB_TASK_PATH);
    align_segment = align_reads::alignment_segment(q, t, 0, 0, result);
    EXPECT_EQ(align_segment.start_on_target, 0);
    EXPECT_EQ(align_segment.align_len_on_target, 4);
    EXPECT_EQ(align_segment.get_aligned_chars(), "A__T");

    for (i = -1; i < 4; i++)
    {
        EXPECT_EQ(align_segment.get_ins_segment_at(i).size(), 0);
    }
}

TEST(Types_alignment_segment, infix)
{
    std::string t = "ACGT";
    std::string q = "CG";
    auto result = align_reads::get_edlib_result(q.c_str(), t.c_str(), q.size(), t.size(), EDLIB_MODE_INFIX, EDLIB_TASK_PATH);
    auto align_segment = align_reads::alignment_segment(q, t, 0, 0, result);
    EXPECT_EQ(align_segment.start_on_target, 1);
    EXPECT_EQ(align_segment.align_len_on_target, 3);
    EXPECT_EQ(align_segment.get_aligned_chars(), "CG");

    
    int i = 0;
    for (i = -1; i < 1; i++)
    {
        EXPECT_EQ(align_segment.get_ins_segment_at(i).size(), 0);
    }
    
    t = "ACGT";
    q ="TATG";
    result = align_reads::get_edlib_result(q.c_str(), t.c_str(), q.size(), t.size(), EDLIB_MODE_INFIX, EDLIB_TASK_PATH);
    align_segment = align_reads::alignment_segment(q, t, 0, 0, result);
    EXPECT_EQ(align_segment.start_on_target, 0);
    EXPECT_EQ(align_segment.align_len_on_target, 3);
    EXPECT_EQ(align_segment.get_aligned_chars(), "ATG");
    EXPECT_EQ(align_segment.get_ins_segment_at(-1), "T");   
}

TEST(Types_alignment_segment, prefix)
{
    std::string t = "ACGT";
    std::string q = "CG";
    auto result = align_reads::get_edlib_result(q.c_str(), t.c_str(), q.size(), t.size(), EDLIB_MODE_PREFIX, EDLIB_TASK_PATH);
    auto align_segment = align_reads::alignment_segment(q, t, 0, 0, result);
    EXPECT_EQ(align_segment.start_on_target, 0);
    EXPECT_EQ(align_segment.align_len_on_target, 3);
    EXPECT_EQ(align_segment.get_aligned_chars(), "_CG");

    
    int i = 0;
    for (i = -1; i < 3; i++)
    {
        EXPECT_EQ(align_segment.get_ins_segment_at(i).size(), 0);
    }
    
    t = "ACGT";
    q ="TA";
    result = align_reads::get_edlib_result(q.c_str(), t.c_str(), q.size(), t.size(), EDLIB_MODE_PREFIX, EDLIB_TASK_PATH);
    align_segment = align_reads::alignment_segment(q, t, 0, 0, result);
    EXPECT_EQ(align_segment.start_on_target, 0);
    EXPECT_EQ(align_segment.align_len_on_target, 1);
    EXPECT_EQ(align_segment.get_aligned_chars(), "A");
    EXPECT_EQ(align_segment.get_ins_segment_at(-1), "T");   
}