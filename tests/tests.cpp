#include <gtest/gtest.h>

#include <Python.h>

#define PY_ARRAY_UNIQUE_SYMBOL gen_ARRAY_API
#include "numpy/arrayobject.h"

#define private public
#include "../include/align_reads/aligner.hpp"
#include "bioparser/fasta_parser.hpp"
#include "bioparser/fastq_parser.hpp"
#include "edlib.h"

const char* reads_fastq_path = "../test_data/reads.fastq";
const char* reads_fasta_path = "../test_data/reads.fasta";
const char* overlap_paf_path = "../test_data/overlap.paf";


// Demonstrate some basic assertions.
TEST(TrivialTests, BasicAssertions) {
    //init();
    // Expect two strings not to be equal.
    EXPECT_STRNE("hello", "world");
    // Expect equality.
    EXPECT_EQ(7 * 6, 42);
}

TEST(BasicTests, Parse_Fasta) {
    auto p = bioparser::Parser<biosoup::NucleicAcid>::Create<bioparser::FastaParser>(reads_fasta_path);
    auto s = p->Parse(-1);
    EXPECT_EQ(2, s.size());
    uint32_t last = s.size() - 1;
    EXPECT_EQ(s[0]->name, "read1");
    EXPECT_EQ(s[0]->InflateData(), "ACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTG");
    EXPECT_EQ(s[last]->name, "read2");
    EXPECT_EQ(s[last]->InflateData(), "ACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTG");
}

TEST(BasicTests, Parse_Fastq) {
    // Quality string ignored for now 
    auto p = bioparser::Parser<biosoup::NucleicAcid>::Create<bioparser::FastqParser>(reads_fastq_path);
    auto s = p->Parse(-1);
    EXPECT_EQ(3165, s.size());
    uint32_t last = s.size() - 1;
    EXPECT_EQ(s[0]->name, "49e202ae-bf1a-4a18-8590-3a3620b1b257");
    EXPECT_EQ(s[0]->InflateData(), "GTAGCCGCTTCGTTCAGTTGCGTATTGCTGCCACCGGCAGCACTTCAACTTGGCCCACCAGCACAAACAGCACAAAGCCCGGCAGCAACATTGCGGCACAGCAGAACAGCCCTCAAGTTCAAGGCGGCTATCGTTAATCAGATGTTTTATATCCCCGCACCAGTTCAACTTGTTCCTGCTGCTTCTCCTGTTGGCCAGCAGAGTCCGTGCTGCCAGCAAGCTGGCCACTCGGAGCAACAAAACATCGGGCTCATCGGCATCGGGAGTCGGAGCCGGAACAGCAGCAGCTTCGGCAGCGGCGGCGGCAGCGTCAGGAACACCCATTTTGGGATCATGCCATCGATTTGAATGACGATTTTCGCATACTATGGCAAATCGGTAATCAGGATATAACATTTGAAATACAGGCACGTACCCTGGGCTACGTCGGGTTCGGATTCTCGCCCGATGGCAATCTGGCTGATGCCGATATGGCCATCGGCTGGGTGAACAAGGGTCAAACCCTATTTTCAGGCCGAAGTGAAGCGCAAAATAAAACCCAGCACAAAGAAAGTTTTTGAGGGTGGGGTTTGGGCTGAATACCCGTAGATGGGGCTTCTGCTGGGCGCTAACGATTTTCATTCCCCGCAGAACGTTTCATTTTTCCAATGCGTGTGGACGAAAATTAAGGGCGGTCAGGCCGGCCGGCCAACTAAATGAAAAGTTGTCGCCCTAATATCGTAATTTATTCTGCTCATGGGCCTGCCGGAGTGGCCCAAGGAGCAGCGACCTCAGTGCCGCCGGAGCAATTGGTGCGTACATGGTGCTGGTACGCTGCTGGCAATCTGGAGATATTCAGCGAAAACTGCGTCAGAATGAACCTGGATTTGCCCGACTGGACACTGGAGCTGGTGGCAACCAACGATACGAGACGATGTGATTTGGTGGGTAAAAGTCAAGGATATGGTTTATAAGTTGCAAAAATAATGCGCATGTTCTTAACAGTTCCTGCAGAAGCACGCGATCCTGTTGGATGCAGCTGCTCTGGCAGCCATACTGACACTTCTCCCTGCTGCCCTTTCCGGAAATGCTCTACATGCAGGCCAGATAGCCAATCTCTCGGTGAGTACTTTCTTATGCTTGTTTTCCCAGCACCAGTAACAACTGAAATCCTGTCCAGACCACGGTGTTCCCGCAAGCCTTTCAGGAGGGACGGCGTCTATTCGCCTCCTCACACCAAACAGCGGAGGATGGACCTGCGCCGCGGCAGTTCCGTGAGATGTACAAAAAGAAACTGCAGATCAAGCGGAGCAATACGTA");
    EXPECT_EQ(s[last]->name, "1b34a65d-b51b-4595-9132-bf9090b546d7");
    EXPECT_EQ(s[last]->InflateData(), "GTACATGCTTCGTTCAGTTACGTATTGCTGAACGTTGCAGAAGTGCCAATGCCGCAAAGCAAAAGAAGCGAGACTTACAAATAAATAAATATACAATTGCAAAAAGTAAATCGTTTACCCTAAATTATGCTGGCTAACTTACACGTTTGTATCACTTTTCGCAAATAATTCGTTTGATTAAATTAATTTGAGGCCTTTGCACGAATACAAATAATGGTAATAGTTACAATAACTAAATCACGTTTTTGGTTGAAATCAATTTAGGGTTCCCCATTTCCAACTTCCATATGTTCTCCTAAAAGATATTGATCTGGAACTCGCTGACGATCCCCTCAAGGAAGTTCTTCAGGCTGTGGTACGAGTGGCTGCGCCAGCCAAATGTTGTTGTTCACCTGCTCCAGCGCCAACTCGATGGCCGATTCCACGGCATGGAGTCCCAGTTTACGGCAGGTAGCGGCCAAGTCCAGAAGCTGGTGCTTATGGTAGTCCTTGTTTAATAGATGGTCAGCGATTTGATCATCTCGGACAGCGTGGAGAAGCCATCGCCATAGCTATGGTGCAATTGAACATAACTAAATTAATATACATTCAACCAATCACACAAAGCATCACTTACTATTCGGCAATCTCCTTTATGTTGCCCTGCAGAAAATCAAAGGCTATCTCATGACCAATGGCATTGGAGGCCACAGCACGGAAGGCCAAGGCGCCATCCTGTTTTAGTACTTACAACGATGTTGGATTGATGGTCATGTTGAGGTACCTACAAAGTTGTTGCACATCAGTCTTTGGTCTTTGATGGATCAAACTATTGTAGTTACTCACTTGGACAGCAGCCAGGGTTTGGTGGTGCAGCCTTGAGGTCGAGTATCTCCTCCTTCTCGGAAGCACTCGCGGTTGACTTGTACTGTTTGTAGGCGAAATACCATTCCGGTGACGAGCCCTCCGCCAAGGAGGTGCAGTAGATCACAGACTTGGGGTTTGGCTTAATGCTGTTGATTTATGATTGCATCGTTTTTAATAAGTTTTTGTATTGTATTTCAGTTACAAGCTTACGGATTGTTTTTGGGATCACGCATCCACTCGCGGAACTTCATCTGCGCCTTTTGGGTACAGCAGTCGTAGTTGAACTTGCAGGCAAAGTAGGCCACCAAGGCACGGTGCTTCAATTGCAAGTGGGACTCGTGTCCGGCTCATGCAGGCCATAATGATCAAAGGCAGGACGTACGATGAATTTCATGAAGGCCTGGTTAAACAAAATATATTATAAAACACTCTAAAAGTTACTAGGTTTACTGGGTTCCACTACCCTGAAATCTCATAGGCAGGCTCCCTCTTCCAGGTTGTAGATCAGATAGTTGAGACCAGGTTTGGCGGCAATCCACAGCAGCTCATCATCCACAGCATCGAACAGCTCCATGAGAGTCCAATCTAAAGGCAAAAGTATTTTATGTTAATAATTTCATTTTACTTGTCTTACTTCTCAAACTTACGGTATGTCGTAGGTAAGATATTCCGCTTGCGACAGATGCAGTGCATCATCCAGCAACTGGGCCCTTGTGATCTGGGCAATGTGCTAAAAGTTCTTCTTGAGCGCCAGCCAGGAGGTCATATCGTAGTTGACCCTATAGTAACCCTGCCAGTTGAGATTCAGATAGATCACGTTATCGCTGTTGCTGCTATGCGCAAAGACATTGCCCACGATGAGCTCCTCTTCGTCCTGCTTCTCATCGGTGGGTGGGTATGTTGTCGCCCTTGCGCAACTCATCCCGTCTCGAGTGATGGGTATAAACCAGGTGCTCTGATCCGCAGTGTTCTTGGGAGGCAGCAGATAGCGTTCCTGGCGCAGCACGAGATCAGCACCACGACGCTCCATTGACCACCGGATAACCGGGCTGTGTGATCCACGAGTCCATGATCTGCTTGACACTCAGATCCTTGGGCAGAGTACCCTGTTCGTGACCATGGCGCGTGAGCATGGCCCACAGATCATCCTTGTCCATGTTTCCATAGGCGAACTTCTTTAGAAGATCGCGAGTGGCCGACAAGGCTACATCACCCACGATCGAATTGAGCATGCGCAGCAAGATAGTGCCCTTTGAGTAGCTGATGGGATCGAAAATCCGCCTGACATCGTTGGTGGAGCGCACATCAAAGGGGAAATGGCATGCGAGGTGTTGTCCGCATCGTGCTCCATCGACTCCTTAAACTCCAGCATGGTCAGTGTCCGCTCTGGAACTCCGGATGGGCGCTCCAGTGCCTTGTAGCTCATGTAGCAGGCGAAGCCCTCCTTCAGCCAGAGATCATCCCACCACTTCAGGGTCACTAGATTGCCGAGCCACTGATGGGCCAACTCGTGTGCAATGATTCCGGCCACCACCTGCATGTGTTCCGATGAGGACGCCAGCTGCAGATCCTCGGGCACCAGTAGCGCCGGATCGCGGAACGTTATAGCGGGTCCCCCAGTTTTCCATGGCAGCGAATCCAAAGTCGGCACGGACACCGTCAATTTGGGCAGCTTGTTCTTAATACCGAAGAGTCCTCGTAGTAGGGCAAGAATTTCGCACCATCTTGTACGCATAGTGAGTCATACCCACAAACTGGAGTCGCGTCCAGATCTCCACTCGCGGCGTCAACCCACTGTCACAGCTGGCAAGCCGCGAATCCACCATGTTGGACACGATGAAAGCCACGGGTAAGTGGGCATCTTCGGCGTGGTCTCGAAATCGTCTCTTATGAAACCACGGCGGGAAAACGCTTGCCCGACTTGGGCATGTTGGAAAGGGCCATCTTGAACTGCATGGGTCTGGCTTGATGCTGATCGAGAAGTTAGCTTTCATGTCCCGGACGGTGAAGCAGGGAAAGGCGCGACAGGCATCGACGGGCGAGAACTAGTGCTTATCATCCATCTAAAAGTAGGGAACTTTAGAAGGGGAACTTTCGTGGAGATGGTTCTGGACTTACTCCGCTTCATTCTTGGTGTCCGGATTGGTGTAGCTGGTCTTGTAGATGCGCTGCAGTGTATCCGTTACCTGGCTGACGAAATCAGACTTAGCAGCACTCTCAGTTGAGTCTCCACCGCCAAAGTCTTGCTCAAATTGATCACGAACGTGGCATTATCCTCCCCGTAGTCGCTGTAGAAATCCAAGTGCTTGCTCCTCGCTGGCATTGCTGGCGCCATCCGCAAGGGCACGGATCACCCGGGCATTGGAGATGCTCACGTTGTGCACGTCGAGCACGATGGGCTCCCAGCTGGTCACCTTGGACACATCCCGTTCGATCTCGATGGTCAGGCTGCGTTGCTGCTAGATGTTGCCACACTTGGCTCAATCAAAGGCTGGAAGATATCATATAGTTGAGAAGGTATCATATAAATAGATTCAAATAGGGAAATGCTCTGCATTTGCAGGTGAAATGCAAATGTTGGCCAAAAACGCCAGGTCGCTGCCATAAAAAGTGCCAATCAATCAATCAATGA");
}

TEST(BasicTests, Construct_Generator_Fasta) {
    const char* p[2] = {reads_fasta_path, nullptr};
    align_reads::Aligner gen_fasta {p, 1, 15, 5, 0.001};
    EXPECT_EQ(gen_fasta.sequences.size(), 2);
    
    uint32_t last1 = gen_fasta.sequences.size() - 1;
    EXPECT_EQ(gen_fasta.sequences[0]->name, "read1");
    EXPECT_EQ(gen_fasta.sequences[0]->InflateData(), "ACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTG");
    EXPECT_EQ(gen_fasta.sequences[last1]->name, "read2");
    EXPECT_EQ(gen_fasta.sequences[last1]->InflateData(), "ACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTG");    
}

TEST(BasicTests, Construct_Generator_Fastq) {
    const char* p[2] = {reads_fastq_path, nullptr};
    align_reads::Aligner gen_fastq {p, 3, 15, 10, 0.0002};
    EXPECT_EQ(gen_fastq.sequences.size(), 3165);
    std::uint32_t last_len = 1000000000;
    for (auto& s: gen_fastq.sequences) {
        EXPECT_LE(s->inflated_len, last_len);
        last_len = s->inflated_len;
    }
    std::uint32_t expected_pos_123 = gen_fastq.id_to_pos_index[123];
    EXPECT_EQ(gen_fastq.sequences[expected_pos_123]->id, 123);
}

TEST(ComponentTests, aligning) {
    std::string target = "GACTAGGAC";
    std::vector<std::string> queries = {"AAGACTTCACGGTACAA", "GACTCAGGGAC"};
    std::vector<std::pair<std::uint32_t, std::uint32_t>> v;
    auto result = align_reads::Aligner::pseudoMSA(queries, target, v);    
    result.print();
}

TEST(ComponentTests, align_overlapping) {
    const char* p[2] = {"../test_data/reads_align.fasta", nullptr};
    align_reads::Aligner gen_fasta {p, 3, 15, 5, 0.001};
    auto result = gen_fasta.align_overlapping(gen_fasta.sequences[0]);
    if (result.valid) result.alignment.print();
    
    const char* p2[2] = {"../test_data/fake_haplotype0.fasta", "../test_data/fake_haplotype1.fasta"};
    align_reads::Aligner gen_fasta2 {p2, 3, 15, 5, 0.001};
    auto result2 = gen_fasta2.align_overlapping(gen_fasta2.sequences[0]);
    if (result2.valid) result2.alignment.print();
}

TEST(ComponentTests, all_inputs) {
    const char* reads_path[2] = {"../test_data/fake_reads0.fasta", "../test_data/fake_reads1.fasta"};
    const char* haplotypes_path[2] = {"../test_data/fake_haplotype0.fasta", "../test_data/fake_haplotype1.fasta"};
    align_reads::Aligner gen {reads_path, 3, 15, 5, 0.001, haplotypes_path};
    EXPECT_EQ(gen.start_of_other_phase, 4);
    EXPECT_EQ(gen.id_to_pos_index[gen.haplotypes_sequences[0][0]->id], 0);
    EXPECT_EQ(gen.id_to_pos_index[gen.haplotypes_sequences[1][0]->id], 0);
    EXPECT_EQ(gen.haplotypes_sequences[0].size(), 1);
    EXPECT_EQ(gen.haplotypes_sequences[1].size(), 1);
    
    auto result = gen.align_overlapping(gen.sequences[0]);
    if (result.valid) result.alignment.print();
 

}

TEST(ComponentTests, align_hap) {
    const char* reads_path[2] = {"../test_data/fake_reads0.fasta", "../test_data/fake_reads1.fasta"};
    const char* haplotypes_path[2] = {"../test_data/fake_haplotype0.fasta", "../test_data/fake_haplotype1.fasta"};
    align_reads::Aligner gen {reads_path, 3, 15, 10, 0.0002, haplotypes_path};
    
    for (int i = 0; i < 7; i++) {
        auto result1 = gen.align_overlapping_plus_haplotypes(gen.sequences[i]);
        if (result1.valid) {
            result1.alignment.print();
        } else {
            std::cout << i << " is invalid " << std::endl;
        }
    }
    
}