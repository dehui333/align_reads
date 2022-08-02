#include <gtest/gtest.h>

#include "align_reads/Inputs.hpp"

#include "../src/Overlapper.cpp"

extern const char *fastq_path;
extern const char *fasta_path;
extern std::shared_ptr<thread_pool::ThreadPool> pool;

TEST(Overlapper, class_Overlapper)
{
    std::vector<std::string> paths = {fastq_path, fastq_path};
    align_reads::Inputs inputs(2);
    inputs.append_to_group(0, paths, pool);
    paths.push_back(fastq_path);
    inputs.append_to_group(1, paths, pool);
    inputs.index_group(0);
    inputs.index_group(1);

    align_reads::Overlapper finder{2, pool};
    finder.index_sequences(inputs.get_group(0), 0);
    finder.index_sequences(inputs.get_group(1), 1);

    auto &seq = inputs.get_id_in_group(0, 0);
    auto ov1 = finder.find_overlaps(seq, 0);
    auto ov2 = finder.find_overlaps(seq, 1);
    EXPECT_GT(ov2.size(), ov1.size());
}
/*
TEST(Overlapper, find_from_real)
{
    std::vector<std::string> paths = {"../ignored/reads/drosophila-f1-100k.fastq"};
    align_reads::Inputs inputs(1);
    inputs.append_to_group(0, paths, pool);

    align_reads::Overlapper finder{1, pool};
    finder.index_sequences(inputs.get_group(0), 0);

    auto &seq = inputs.get_id_in_group(0, 0);
    auto ov = finder.find_overlaps(seq, 0);

    std::cout << "num " << ov.size() << std::endl;

}*/
