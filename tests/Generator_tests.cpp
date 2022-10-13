#include <gtest/gtest.h>

#define private public
#include "../src/Generator.cpp"


extern const char* fastq_path;
extern const char* fasta_path;
extern std::shared_ptr<thread_pool::ThreadPool> pool;
TEST(Generator, class_Generator) {
    std::vector<std::string> paths = {fastq_path};
    align_reads::Generator gen {paths, pool, 10};
    EXPECT_FALSE(gen.has_truth);
    EXPECT_EQ(gen.inputs.get_group(0).size(), 3165);
    EXPECT_EQ(gen.start_of_reads1, 0);

    align_reads::Generator gen2 {paths, paths, paths, paths, pool, 10};
    EXPECT_TRUE(gen2.has_truth);
    EXPECT_EQ(gen2.inputs.get_group(0).size(), 3165 * 2);
    EXPECT_EQ(gen2.start_of_reads1, 3165);

    auto& seq = gen2.inputs.get_id_in_group(0, 0);
    auto num1 = gen2.overlapper.find_overlaps(seq, 0).size();
    auto num2 = gen2.overlapper.find_overlaps(seq, 1).size();
    EXPECT_GT(num1, num2);

}