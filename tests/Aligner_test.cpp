#include <gtest/gtest.h>

#include "align_reads/Utilities.hpp"
#define private public
#include "../src/Aligner.cpp"

extern std::shared_ptr<thread_pool::ThreadPool> pool;

using namespace align_reads;

TEST(Aligner, get_edlib_results)
{
    std::string s1 = "ACGT";
    std::string s2 = "AC";
    std::vector<align_reads::EdlibTask> tasks;
    tasks.emplace_back(s1.c_str(), 4, s1.c_str(), 4, EDLIB_MODE_GLOBAL, EDLIB_TASK_DISTANCE);
    tasks.emplace_back(s1.c_str(), 4, s2.c_str(), 2, EDLIB_MODE_GLOBAL, EDLIB_TASK_DISTANCE);
    tasks.emplace_back(s2.c_str(), 2, s1.c_str(), 4, EDLIB_MODE_PREFIX, EDLIB_TASK_DISTANCE);

    auto results = align_reads::get_edlib_results(tasks, pool);
    EXPECT_EQ(results[0].editDistance, 0);
    EXPECT_EQ(results[1].editDistance, 2);
    EXPECT_EQ(results[2].editDistance, 0);

    auto result = align_reads::get_edlib_result(s1.c_str(), 4, s1.c_str(), 4, EDLIB_MODE_GLOBAL, EDLIB_TASK_DISTANCE);
    EXPECT_EQ(result.editDistance, 0);
}

TEST(Aligner, Futures)
{
    std::string s1 = "ACGT";
    std::string s2 = "AC";

    auto futures = align_reads::Futures<EdlibAlignResult>(pool, 3);
    futures.add_inputs(align_reads::get_edlib_result, s1.c_str(), 4, s1.c_str(), 4, EDLIB_MODE_GLOBAL, EDLIB_TASK_DISTANCE);
    futures.add_inputs(align_reads::get_edlib_result, s1.c_str(), 4, s2.c_str(), 2, EDLIB_MODE_GLOBAL, EDLIB_TASK_DISTANCE);
    futures.add_inputs(align_reads::get_edlib_result, s2.c_str(), 2, s1.c_str(), 4, EDLIB_MODE_PREFIX, EDLIB_TASK_DISTANCE);

    auto results = futures.get();

    EXPECT_EQ(results[0].editDistance, 0);
    EXPECT_EQ(results[1].editDistance, 2);
    EXPECT_EQ(results[2].editDistance, 0);
}

TEST(Aligner_AlignmentSegment, global_all_match)
{
    std::string s1 = "ACGT";
    auto result = align_reads::get_edlib_result(s1.c_str(), 4, s1.c_str(), 4, EDLIB_MODE_GLOBAL, EDLIB_TASK_PATH);
    auto align_segment = align_reads::AlignmentSegment(s1, 0, s1, 0, result);
    EXPECT_EQ(align_segment.start_on_target, 0);
    EXPECT_EQ(align_segment.end_on_target, 3);
    EXPECT_EQ(align_segment.get_aligned_chars().compare(s1), 0);
    for (int i = -1; i < 4; i++)
    {
        EXPECT_EQ(align_segment.get_ins_segment_at(i).size(), 0);
    }

    align_segment = align_reads::get_alignment_segment(s1, 0, 4, s1, 0, 4, EDLIB_MODE_GLOBAL, EDLIB_TASK_PATH);
    EXPECT_EQ(align_segment.start_on_target, 0);
    EXPECT_EQ(align_segment.end_on_target, 3);
    EXPECT_EQ(align_segment.get_aligned_chars().compare(s1), 0);
    for (int i = -1; i < 4; i++)
    {
        EXPECT_EQ(align_segment.get_ins_segment_at(i).size(), 0);
    }


}

TEST(Aligner_AlignmentSegment, global_ins)
{
    std::string t = "GGGG";
    std::string q = "AGGGG";
    auto result = align_reads::get_edlib_result(q.c_str(), q.size(), t.c_str(), t.size(), EDLIB_MODE_GLOBAL, EDLIB_TASK_PATH);
    auto align_segment = align_reads::AlignmentSegment(q, 0, t, 0, result);
    EXPECT_EQ(align_segment.start_on_target, 0);
    EXPECT_EQ(align_segment.end_on_target, 3);
    EXPECT_EQ(align_segment.get_aligned_chars().compare(t), 0);

    EXPECT_EQ(align_segment.get_ins_segment_at(-1).size(), 1);
    int i = 0;
    for (; i < 4; i++)
    {
        EXPECT_EQ(align_segment.get_ins_segment_at(i).size(), 0);
    }
    EXPECT_EQ(align_segment.get_ins_segment_at(-1), "A");

    align_segment = align_reads::get_alignment_segment(q, 0, 5, t, 0, 4, EDLIB_MODE_GLOBAL, EDLIB_TASK_PATH);
    EXPECT_EQ(align_segment.start_on_target, 0);
    EXPECT_EQ(align_segment.end_on_target, 3);
    EXPECT_EQ(align_segment.get_aligned_chars().compare(t), 0);

    EXPECT_EQ(align_segment.get_ins_segment_at(-1).size(), 1);

    for (i=0; i < 4; i++)
    {
        EXPECT_EQ(align_segment.get_ins_segment_at(i).size(), 0);
    }
    EXPECT_EQ(align_segment.get_ins_segment_at(-1), "A");


    t = "GGGG";
    q = "AGAGGGAA";
    result = align_reads::get_edlib_result(q.c_str(), q.size(), t.c_str(), t.size(), EDLIB_MODE_GLOBAL, EDLIB_TASK_PATH);
    align_segment = align_reads::AlignmentSegment(q, 0, t, 0, result);
    EXPECT_EQ(align_segment.start_on_target, 0);
    EXPECT_EQ(align_segment.end_on_target, 3);
    EXPECT_EQ(align_segment.get_aligned_chars().compare(t), 0);

    EXPECT_EQ(align_segment.get_ins_segment_at(-1).size(), 1);
    EXPECT_EQ(align_segment.get_ins_segment_at(0).size(), 1);
    EXPECT_EQ(align_segment.get_ins_segment_at(3).size(), 2);

    align_segment = align_reads::get_alignment_segment(q, 0, 8, t, 0, 4, EDLIB_MODE_GLOBAL, EDLIB_TASK_PATH);
    EXPECT_EQ(align_segment.start_on_target, 0);
    EXPECT_EQ(align_segment.end_on_target, 3);
    EXPECT_EQ(align_segment.get_aligned_chars().compare(t), 0);

    EXPECT_EQ(align_segment.get_ins_segment_at(-1).size(), 1);
    EXPECT_EQ(align_segment.get_ins_segment_at(0).size(), 1);
    EXPECT_EQ(align_segment.get_ins_segment_at(3).size(), 2);

}

TEST(Aligner_AlignmentSegment, global_del)
{
    std::string t = "ACGT";
    std::string q = "AGT";
    auto result = align_reads::get_edlib_result(q.c_str(), q.size(), t.c_str(), t.size(), EDLIB_MODE_GLOBAL, EDLIB_TASK_PATH);
    auto align_segment = align_reads::AlignmentSegment(q, 0, t, 0, result);
    EXPECT_EQ(align_segment.start_on_target, 0);
    EXPECT_EQ(align_segment.end_on_target, 3);
    EXPECT_EQ(align_segment.get_aligned_chars(), "A_GT");

    int i = 0;
    for (i = -1; i < 4; i++)
    {
        EXPECT_EQ(align_segment.get_ins_segment_at(i).size(), 0);
    }

    t = "ACGT";
    q = "AT";
    result = align_reads::get_edlib_result(q.c_str(), q.size(), t.c_str(), t.size(), EDLIB_MODE_GLOBAL, EDLIB_TASK_PATH);
    align_segment = align_reads::AlignmentSegment(q, 0, t, 0, result);
    EXPECT_EQ(align_segment.start_on_target, 0);
    EXPECT_EQ(align_segment.end_on_target, 3);
    EXPECT_EQ(align_segment.get_aligned_chars(), "A__T");

    for (i = -1; i < 4; i++)
    {
        EXPECT_EQ(align_segment.get_ins_segment_at(i).size(), 0);
    }
}

TEST(Aligner_AlignmentSegment, infix)
{
    std::string t = "ACGT";
    std::string q = "CG";
    auto result = align_reads::get_edlib_result(q.c_str(), q.size(), t.c_str(), t.size(), EDLIB_MODE_INFIX, EDLIB_TASK_PATH);
    auto align_segment = align_reads::AlignmentSegment(q, 0, t, 0, result);
    EXPECT_EQ(align_segment.start_on_target, 1);
    EXPECT_EQ(align_segment.end_on_target, 2);
    EXPECT_EQ(align_segment.get_aligned_chars(), "CG");

    int i = 0;
    for (i = -1; i < 2; i++)
    {
        EXPECT_EQ(align_segment.get_ins_segment_at(i).size(), 0);
    }

    align_segment = align_reads::get_alignment_segment(q, 0, 2, t, 1, 3, EDLIB_MODE_INFIX, EDLIB_TASK_PATH);
    EXPECT_EQ(align_segment.start_on_target, 1);
    EXPECT_EQ(align_segment.end_on_target, 2);
    EXPECT_EQ(align_segment.get_aligned_chars(), "CG");

    for (i = -1; i < 2; i++)
    {
        EXPECT_EQ(align_segment.get_ins_segment_at(i).size(), 0);
    }


    t = "ACGT";
    q = "TATG";
    result = align_reads::get_edlib_result(q.c_str(), q.size(), t.c_str(), t.size(), EDLIB_MODE_INFIX, EDLIB_TASK_PATH);
    align_segment = align_reads::AlignmentSegment(q, 0, t, 0, result);
    EXPECT_EQ(align_segment.start_on_target, 0);
    EXPECT_EQ(align_segment.end_on_target, 2);
    EXPECT_EQ(align_segment.get_aligned_chars(), "ATG");
    EXPECT_EQ(align_segment.get_ins_segment_at(-1), "T");
}

TEST(Aligner_AlignmentSegment, prefix)
{
    std::string t = "ACGT";
    std::string q = "CG";
    auto result = align_reads::get_edlib_result(q.c_str(), q.size(), t.c_str(), t.size(), EDLIB_MODE_PREFIX, EDLIB_TASK_PATH);
    auto align_segment = align_reads::AlignmentSegment(q, 0, t, 0, result);
    EXPECT_EQ(align_segment.start_on_target, 0);
    EXPECT_EQ(align_segment.end_on_target, 2);
    EXPECT_EQ(align_segment.get_aligned_chars(), "_CG");

    int i = 0;
    for (i = -1; i < 3; i++)
    {
        EXPECT_EQ(align_segment.get_ins_segment_at(i).size(), 0);
    }

    align_segment = align_reads::get_alignment_segment(q, 0, 2, t, 1, 3, EDLIB_MODE_PREFIX, EDLIB_TASK_PATH);
    EXPECT_EQ(align_segment.start_on_target, 1);
    EXPECT_EQ(align_segment.end_on_target, 2);
    EXPECT_EQ(align_segment.get_aligned_chars(), "CG");

    for (i = -1; i < 2; i++)
    {
        EXPECT_EQ(align_segment.get_ins_segment_at(i).size(), 0);
    }


    t = "ACGT";
    q = "TA";
    result = align_reads::get_edlib_result(q.c_str(), q.size(), t.c_str(), t.size(), EDLIB_MODE_PREFIX, EDLIB_TASK_PATH);
    align_segment = align_reads::AlignmentSegment(q, 0, t, 0, result);
    EXPECT_EQ(align_segment.start_on_target, 0);
    EXPECT_EQ(align_segment.end_on_target, 0);
    EXPECT_EQ(align_segment.get_aligned_chars(), "A");
    EXPECT_EQ(align_segment.get_ins_segment_at(-1), "T");
}

TEST(Aligner_AlignmentSegment, get_at_target_pos)
{
    std::string t  = "TAGGCATACAGG";
    std::string q1 = "TAGTG";
    std::string q2 = "CAACATGGAAAA";
    std::string q3 = "TGGCATATCA";
    // TAG_|GCAT|A_CA|_GG

    Futures<AlignmentSegment> futures(pool, 3);
    futures.add_inputs(get_alignment_segment, q1, 0, 5, t, 0, 4, EDLIB_MODE_INFIX, EDLIB_TASK_PATH);
    futures.add_inputs(get_alignment_segment, q2, 0, 8, t, 4, 8, EDLIB_MODE_INFIX, EDLIB_TASK_PATH);
    futures.add_inputs(get_alignment_segment, q3, 0, 10, t, 0, 12, EDLIB_MODE_INFIX, EDLIB_TASK_PATH);
    std::vector<AlignmentSegment> segments = futures.get();

    EXPECT_EQ(segments[0].get_at_target_pos(0, 0), 'T');
    EXPECT_EQ(segments[0].get_at_target_pos(1, 0), 'A');
    EXPECT_EQ(segments[0].get_at_target_pos(2, 0), 'G');
    EXPECT_EQ(segments[0].get_at_target_pos(2, 1), 'T');
    EXPECT_EQ(segments[0].get_at_target_pos(3, 0), 'G');

    EXPECT_EQ(segments[1].get_at_target_pos(4, 0), 'C');
    EXPECT_EQ(segments[1].get_at_target_pos(5, 0), 'A');
    EXPECT_EQ(segments[1].get_at_target_pos(6, 0), '_');
    EXPECT_EQ(segments[1].get_at_target_pos(7, 0), 'A');
    EXPECT_EQ(segments[1].get_at_target_pos(8, 0), 'C');
    EXPECT_EQ(segments[1].get_at_target_pos(9, 0), 'A');
    EXPECT_EQ(segments[1].get_at_target_pos(9, 1), 'T');
    EXPECT_EQ(segments[1].get_at_target_pos(10, 0), 'G');
    EXPECT_EQ(segments[1].get_at_target_pos(11, 0), 'G');
}
