#include <gtest/gtest.h>

#include "align_reads/Utilities.hpp"

extern std::shared_ptr<thread_pool::ThreadPool> pool;

TEST(Utilities, class_Futures)
{
    auto f = [](int i, int j)
    {
        return i + j;
    };
    auto futures = align_reads::Futures<int>(pool, 10);
    for (auto i = 0; i < 10; i++)
    {
        futures.add_inputs(f, i, i * i);
    }
    auto r = futures.get();
    for (auto i = 0; i < 10; i++)
    {
        EXPECT_EQ(r[i], i + i * i);
    }
}