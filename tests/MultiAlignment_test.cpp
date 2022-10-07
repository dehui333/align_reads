#include <gtest/gtest.h>
#include <iostream>
#include "align_reads/Utilities.hpp"
#define private public
#include "../src/MultiAlignment.cpp"
#include "align_reads/Aligner.hpp"

using namespace align_reads;
extern std::shared_ptr<thread_pool::ThreadPool> pool;

#define GAP_CHAR '_'
TEST(MultiAlignment, basic_test)
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
    MultiAlignment ma{std::move(t), std::move(segments)};

    /*    
    for (std::uint32_t i = 0; i < 3; i++)
    {
        auto it = ma.iterator(i, 0);
        std::string s;
        while (it.has_next())
        {
            s.push_back(it.next());
        }
        std::cout << s << std::endl;
    }*/
    

    auto it = ma.iterator(0, 0);
    EXPECT_TRUE(it.has_next());
    EXPECT_EQ(it.next(), 'T');
    
    EXPECT_TRUE(it.has_next());
    EXPECT_EQ(it.next(), 'A');
    
    EXPECT_TRUE(it.has_next());
    EXPECT_EQ(it.next(), 'G');
    
    EXPECT_TRUE(it.has_next());
    EXPECT_EQ(it.next(), 'T');
    
    EXPECT_TRUE(it.has_next());
    EXPECT_EQ(it.next(), 'G');
    
    EXPECT_TRUE(it.has_next());
    EXPECT_EQ(it.next(), GAP_CHAR);

    EXPECT_TRUE(it.has_next());
    EXPECT_EQ(it.next(), GAP_CHAR);
    
    EXPECT_TRUE(it.has_next());
    EXPECT_EQ(it.next(), GAP_CHAR);

    EXPECT_TRUE(it.has_next());
    EXPECT_EQ(it.next(), GAP_CHAR);

    EXPECT_TRUE(it.has_next());
    EXPECT_EQ(it.next(), GAP_CHAR);

    EXPECT_TRUE(it.has_next());
    EXPECT_EQ(it.next(), GAP_CHAR);

    EXPECT_TRUE(it.has_next());
    EXPECT_EQ(it.next(), GAP_CHAR);

    EXPECT_TRUE(it.has_next());
    EXPECT_EQ(it.next(), GAP_CHAR);

    EXPECT_TRUE(it.has_next());
    EXPECT_EQ(it.next(), GAP_CHAR);

    EXPECT_TRUE(it.has_next());
    EXPECT_EQ(it.next(), GAP_CHAR);

    EXPECT_FALSE(it.has_next());
    
    // move assignment deleted, maybe because of the reference variable member
    auto it2 = ma.iterator(1, 0);
    EXPECT_TRUE(it2.has_next());
    EXPECT_EQ(it2.next(), GAP_CHAR);
    
    EXPECT_TRUE(it2.has_next());
    EXPECT_EQ(it2.next(), GAP_CHAR);
    
    EXPECT_TRUE(it2.has_next());
    EXPECT_EQ(it2.next(), GAP_CHAR);
    
    EXPECT_TRUE(it2.has_next());
    EXPECT_EQ(it2.next(), GAP_CHAR);
    
    EXPECT_TRUE(it2.has_next());
    EXPECT_EQ(it2.next(), GAP_CHAR);
    
    EXPECT_TRUE(it2.has_next());
    EXPECT_EQ(it2.next(), 'C');

    EXPECT_TRUE(it2.has_next());
    EXPECT_EQ(it2.next(), 'A');
    
    EXPECT_TRUE(it2.has_next());
    EXPECT_EQ(it2.next(), GAP_CHAR);

    EXPECT_TRUE(it2.has_next());
    EXPECT_EQ(it2.next(), 'A');

    EXPECT_TRUE(it2.has_next());
    EXPECT_EQ(it2.next(), GAP_CHAR);

    EXPECT_TRUE(it2.has_next());
    EXPECT_EQ(it2.next(), 'C');

    EXPECT_TRUE(it2.has_next());
    EXPECT_EQ(it2.next(), 'A');

    EXPECT_TRUE(it2.has_next());
    EXPECT_EQ(it2.next(), 'T');

    EXPECT_TRUE(it2.has_next());
    EXPECT_EQ(it2.next(), 'G');

    EXPECT_TRUE(it2.has_next());
    EXPECT_EQ(it2.next(), 'G');

    EXPECT_FALSE(it2.has_next());

    auto it3 = ma.iterator(2, 0);
    EXPECT_TRUE(it3.has_next());
    EXPECT_EQ(it3.next(), 'T');
    
    EXPECT_TRUE(it3.has_next());
    EXPECT_EQ(it3.next(), GAP_CHAR);
    
    EXPECT_TRUE(it3.has_next());
    EXPECT_EQ(it3.next(), 'G');
    
    EXPECT_TRUE(it3.has_next());
    EXPECT_EQ(it3.next(), GAP_CHAR);
    
    EXPECT_TRUE(it3.has_next());
    EXPECT_EQ(it3.next(), 'G');
    
    EXPECT_TRUE(it3.has_next());
    EXPECT_EQ(it3.next(), 'C');

    EXPECT_TRUE(it3.has_next());
    EXPECT_EQ(it3.next(), 'A');
    
    EXPECT_TRUE(it3.has_next());
    EXPECT_EQ(it3.next(), 'T');

    EXPECT_TRUE(it3.has_next());
    EXPECT_EQ(it3.next(), 'A');

    EXPECT_TRUE(it3.has_next());
    EXPECT_EQ(it3.next(), 'T');

    EXPECT_TRUE(it3.has_next());
    EXPECT_EQ(it3.next(), 'C');

    EXPECT_TRUE(it3.has_next());
    EXPECT_EQ(it3.next(), 'A');

    EXPECT_TRUE(it3.has_next());
    EXPECT_EQ(it3.next(), GAP_CHAR);

    EXPECT_TRUE(it3.has_next());
    EXPECT_EQ(it3.next(), GAP_CHAR);

    EXPECT_TRUE(it3.has_next());
    EXPECT_EQ(it3.next(), GAP_CHAR);

    EXPECT_FALSE(it3.has_next());   
}

TEST(MultiAlignment, print)
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
    MultiAlignment ma{std::move(t), std::move(segments)};
    
  
    ma.print_in_window(0, 5);
}
