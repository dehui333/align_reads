#include <gtest/gtest.h>

#include "align_reads/Utilities.hpp"
#define private public
#include "align_reads/AlignmentSegment.hpp"
#include "../src/Aligner.cpp"
#include "align_reads/Inputs.hpp"
#include "align_reads/MultiAlignment.hpp"
#include "align_reads/Overlapper.hpp"

extern std::shared_ptr<thread_pool::ThreadPool> pool;
const char *reads0_path = "../test_data/fake_reads0.fasta";
const char *align_reverse_path = "../test_data/align_reverse.fasta";

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

TEST(Aligner, align_overlap)
{
    
    biosoup::NucleicAcid::num_objects = 0;
    std::vector<std::string> paths = {reads0_path};
    align_reads::Inputs inputs(1);
    inputs.append_to_group(0, paths, pool);
    Overlapper ovl {1, pool};
    ovl.index_sequences(inputs.get_group(0), 0);
    auto& target = inputs.get_id_in_group(0, 0);
    auto overlaps = ovl.find_overlaps(target, 0);
    std::string target_string = target->InflateData();
    auto result = align_overlap(overlaps[0], &inputs, target_string.c_str(), target_string.size());
    EXPECT_TRUE(overlaps[0].strand);
    auto alignment = AlignmentSegment(result.clipped_query, 0, target_string, result.t_start, result.result);
    alignment.print(target_string);

}

TEST(Aligner, align_overlap_reverse)
{
    biosoup::NucleicAcid::num_objects = 0;
    std::vector<std::string> paths = {align_reverse_path};
    align_reads::Inputs inputs(1);
    inputs.append_to_group(0, paths, pool);
    Overlapper ovl {1, pool};
    ovl.index_sequences(inputs.get_group(0), 0);
    auto& target = inputs.get_id_in_group(0, 0);
    auto overlaps = ovl.find_overlaps(target, 0);
    std::string target_string = target->InflateData();
    auto result = align_overlap(overlaps[0], &inputs, target_string.c_str(), target_string.size());
    EXPECT_FALSE(overlaps[0].strand);
    auto alignment = AlignmentSegment(result.clipped_query, 0, target_string, result.t_start, result.result);
    alignment.print(target_string);

}

