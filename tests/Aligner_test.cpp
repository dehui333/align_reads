#include <gtest/gtest.h>

#include "../src/Aligner.cpp"

extern std::shared_ptr<thread_pool::ThreadPool> pool;

TEST(Aligner, get_edlib_results) {
    std::string s1 = "ACGT";
    std::string s2 = "AC";
    std::vector<align_reads::edlib_task> tasks;
    tasks.emplace_back(s1.c_str(), s1.c_str(), 4, 4, EDLIB_MODE_GLOBAL, EDLIB_TASK_DISTANCE);
    tasks.emplace_back(s1.c_str(), s2.c_str(), 4, 2, EDLIB_MODE_GLOBAL, EDLIB_TASK_DISTANCE);
    tasks.emplace_back(s2.c_str(), s1.c_str(), 2, 4, EDLIB_MODE_PREFIX, EDLIB_TASK_DISTANCE);

    auto results = align_reads::get_edlib_results(tasks, pool);
    EXPECT_EQ(results[0].editDistance, 0);
    EXPECT_EQ(results[1].editDistance, 2);
    EXPECT_EQ(results[2].editDistance, 0);

    auto result = align_reads::get_edlib_result(s1.c_str(), s1.c_str(), 4, 4, EDLIB_MODE_GLOBAL, EDLIB_TASK_DISTANCE);
    EXPECT_EQ(result.editDistance, 0);

}