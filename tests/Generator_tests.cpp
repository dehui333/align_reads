#include <gtest/gtest.h>

#define private public
#include "../src/Generator.cpp"

extern const char* fastq_path;
extern const char* fasta_path;

TEST(BasicTests, class_Generator) {
    std::vector<std::string> paths = {fastq_path};
    std::shared_ptr<thread_pool::ThreadPool> pool = std::make_shared<thread_pool::ThreadPool>(10);
    align_reads::Generator gen {paths, pool};
    EXPECT_FALSE(gen.has_haplotypes);
    EXPECT_EQ(gen.inputs.get_group(0).size(), 3165);
    EXPECT_EQ(gen.start_of_reads1, 0);

    align_reads::Generator gen2 {paths, paths, paths, paths, pool};
    EXPECT_TRUE(gen2.has_haplotypes);
    EXPECT_EQ(gen2.inputs.get_group(0).size(), 3165 * 2);
    EXPECT_EQ(gen2.start_of_reads1, 3165);

    auto& seq = gen2.inputs.get_id_in_group(0, 0);
    auto num1 = gen2.overlapper.find_overlaps(seq, 0).size();
    auto num2 = gen2.overlapper.find_overlaps(seq, 1).size();
    EXPECT_GT(num1, num2);

}