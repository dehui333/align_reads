#include <gtest/gtest.h>

#include "align_reads/Utilities.hpp"

#define private public
#include "../src/Converter.cpp"

using namespace align_reads;

extern std::shared_ptr<thread_pool::ThreadPool> pool;

TEST(Converter, Converter_all_match)
{
    
    std::string t  = "TAGGCATACCGG";
    std::string q1 = "TAGG";
    std::string q2 = "CATA";
    std::string q3 = "CCGG";

    Futures<AlignmentSegment> futures(pool, 3);
    futures.add_inputs(get_alignment_segment, q1, 0, 4, t, 0, 4, EDLIB_MODE_INFIX, EDLIB_TASK_PATH);
    futures.add_inputs(get_alignment_segment, q2, 0, 4, t, 4, 4, EDLIB_MODE_INFIX, EDLIB_TASK_PATH);
    futures.add_inputs(get_alignment_segment, q3, 0, 4, t, 8, 4, EDLIB_MODE_INFIX, EDLIB_TASK_PATH);
    std::vector<AlignmentSegment> segments = futures.get();
    
    EXPECT_EQ(segments[0].aligned_chars, "TAGG");
    EXPECT_EQ(segments[1].aligned_chars, "CATA");
    EXPECT_EQ(segments[2].aligned_chars, "CCGG");
    MultiAlignment m_align {std::move(t), std::move(segments)};
    AlignmentConverter converter {m_align, 4, 3};

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
    converter = {m_align, 4, 3};

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
    AlignmentConverter converter {m_align, 4, 3};

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