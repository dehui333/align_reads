#include <gtest/gtest.h>

#include "align_reads/Utilities.hpp"
#define private public
#include "align_reads/AlignmentSegment.hpp"
#include "../src/Aligner.cpp"
#include "align_reads/MultiAlignment.hpp"


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
