#include <gtest/gtest.h>
#include <iostream>
#define private public
#include "../src/AlignCounter.cpp"
#include "align_reads/Aligner.hpp"

using namespace align_reads;

TEST(AlignCounter, construct)
{
  std::string t = "TAGGCATACAGG";
  std::string q1 = "TCGGTCTTAACA";
  std::string q2 = "AGGTTCATACAGG";
  std::string q3 = "CTAGGCAACAGGCC";
  /*
  TAGG  CATA CAGG
  TCGGT CTTAACA
   AGGTTCATA CAGG
 CTAGG  CA A CAGGCC
  */

  auto result1 = get_edlib_result(q1.c_str(), 12, t.c_str(), 12, EDLIB_MODE_INFIX, EDLIB_TASK_PATH);
  auto result2 = get_edlib_result(q2.c_str(), 13, t.c_str(), 12, EDLIB_MODE_INFIX, EDLIB_TASK_PATH);
  auto result3 = get_edlib_result(q3.c_str(), 14, t.c_str(), 12, EDLIB_MODE_INFIX, EDLIB_TASK_PATH);

  clipped_alignment<EdlibAlignResult> ca1;
  ca1.clipped_query = q1;
  ca1.q_start = 0;
  ca1.q_end = 11;
  ca1.t_start = 0;
  ca1.t_end = 11;
  ca1.result = result1;
  clipped_alignment<EdlibAlignResult> ca2;
  ca2.clipped_query = q2;
  ca2.q_start = 0;
  ca2.q_end = 12;
  ca2.t_start = 0;
  ca2.t_end = 11;
  ca2.result = result2;
  clipped_alignment<EdlibAlignResult> ca3;
  ca3.clipped_query = q3;
  ca3.q_start = 0;
  ca3.q_end = 13;
  ca3.t_start = 0;
  ca3.t_end = 11;
  ca3.result = result3;

  std::vector<clipped_alignment<EdlibAlignResult>> v{ca1, ca2, ca3};
  


  AlignCounter ac{t, v};
  /*
  for (auto i = 0; i < v.size(); i++)
  {
    std::cout << i << " iden " << 1 - (float) v[i].result.editDistance / v[i].clipped_query.size() << std::endl;
  }

  for (auto i = 0; i < ac.counts.size(); i++)
  {
    std::cout << "target pos " << i << std::endl;
    for (auto j = 0; j < ac.counts[i].size(); j++)
    {
      std::cout << "ins pos " << j << std::endl;
      for (auto k = 0; k < 5; k++)
      {
        std::cout << "char idx " << k << std::endl;
        std::cout << ac.stats[i][j][k] << std::endl;
      }
    }
  } */

  EXPECT_EQ(ac.max_ins_at(0), 0);
  EXPECT_EQ(ac.base_count_at(0, 0, 0), 0);
  EXPECT_EQ(ac.base_count_at(0, 0, 1), 0);
  EXPECT_EQ(ac.base_count_at(0, 0, 2), 0);
  EXPECT_EQ(ac.base_count_at(0, 0, 3), 3);
  EXPECT_EQ(ac.base_count_at(0, 0, 4), 0);

  EXPECT_EQ(ac.max_ins_at(1), 0);
  EXPECT_EQ(ac.base_count_at(1, 0, 0), 3);
  EXPECT_EQ(ac.base_count_at(1, 0, 1), 1);
  EXPECT_EQ(ac.base_count_at(1, 0, 2), 0);
  EXPECT_EQ(ac.base_count_at(1, 0, 3), 0);
  EXPECT_EQ(ac.base_count_at(1, 0, 4), 0);

  EXPECT_EQ(ac.max_ins_at(2), 0);
  EXPECT_EQ(ac.base_count_at(2, 0, 0), 0);
  EXPECT_EQ(ac.base_count_at(2, 0, 1), 0);
  EXPECT_EQ(ac.base_count_at(2, 0, 2), 4);
  EXPECT_EQ(ac.base_count_at(2, 0, 3), 0);
  EXPECT_EQ(ac.base_count_at(2, 0, 4), 0);

  EXPECT_EQ(ac.max_ins_at(3), 2);

  EXPECT_EQ(ac.base_count_at(3, 0, 0), 0);
  EXPECT_EQ(ac.base_count_at(3, 0, 1), 0);
  EXPECT_EQ(ac.base_count_at(3, 0, 2), 4);
  EXPECT_EQ(ac.base_count_at(3, 0, 3), 0);
  EXPECT_EQ(ac.base_count_at(3, 0, 4), 0);

  EXPECT_EQ(ac.base_count_at(3, 1, 0), 0);
  EXPECT_EQ(ac.base_count_at(3, 1, 1), 0);
  EXPECT_EQ(ac.base_count_at(3, 1, 2), 0);
  EXPECT_EQ(ac.base_count_at(3, 1, 3), 2);
  EXPECT_EQ(ac.base_count_at(3, 1, 4), 2);

  EXPECT_EQ(ac.base_count_at(3, 2, 0), 0);
  EXPECT_EQ(ac.base_count_at(3, 2, 1), 0);
  EXPECT_EQ(ac.base_count_at(3, 2, 2), 0);
  EXPECT_EQ(ac.base_count_at(3, 2, 3), 1);
  EXPECT_EQ(ac.base_count_at(3, 2, 4), 3);

  EXPECT_EQ(ac.max_ins_at(4), 0);
  EXPECT_EQ(ac.base_count_at(4, 0, 0), 0);
  EXPECT_EQ(ac.base_count_at(4, 0, 1), 4);
  EXPECT_EQ(ac.base_count_at(4, 0, 2), 0);
  EXPECT_EQ(ac.base_count_at(4, 0, 3), 0);
  EXPECT_EQ(ac.base_count_at(4, 0, 4), 0);

  EXPECT_EQ(ac.max_ins_at(5), 0);
  EXPECT_EQ(ac.base_count_at(5, 0, 0), 3);
  EXPECT_EQ(ac.base_count_at(5, 0, 1), 0);
  EXPECT_EQ(ac.base_count_at(5, 0, 2), 0);
  EXPECT_EQ(ac.base_count_at(5, 0, 3), 1);
  EXPECT_EQ(ac.base_count_at(5, 0, 4), 0);

  EXPECT_EQ(ac.max_ins_at(6), 0);
  EXPECT_EQ(ac.base_count_at(6, 0, 0), 0);
  EXPECT_EQ(ac.base_count_at(6, 0, 1), 0);
  EXPECT_EQ(ac.base_count_at(6, 0, 2), 0);
  EXPECT_EQ(ac.base_count_at(6, 0, 3), 3);
  EXPECT_EQ(ac.base_count_at(6, 0, 4), 1);

  EXPECT_EQ(ac.max_ins_at(7), 1);
  EXPECT_EQ(ac.base_count_at(7, 0, 0), 4);
  EXPECT_EQ(ac.base_count_at(7, 0, 1), 0);
  EXPECT_EQ(ac.base_count_at(7, 0, 2), 0);
  EXPECT_EQ(ac.base_count_at(7, 0, 3), 0);
  EXPECT_EQ(ac.base_count_at(7, 0, 4), 0);

  EXPECT_EQ(ac.base_count_at(7, 1, 0), 1);
  EXPECT_EQ(ac.base_count_at(7, 1, 1), 0);
  EXPECT_EQ(ac.base_count_at(7, 1, 2), 0);
  EXPECT_EQ(ac.base_count_at(7, 1, 3), 0);
  EXPECT_EQ(ac.base_count_at(7, 1, 4), 3);

  EXPECT_EQ(ac.max_ins_at(8), 0);
  EXPECT_EQ(ac.base_count_at(8, 0, 0), 0);
  EXPECT_EQ(ac.base_count_at(8, 0, 1), 4);
  EXPECT_EQ(ac.base_count_at(8, 0, 2), 0);
  EXPECT_EQ(ac.base_count_at(8, 0, 3), 0);
  EXPECT_EQ(ac.base_count_at(8, 0, 4), 0);

  EXPECT_EQ(ac.max_ins_at(9), 0);
  EXPECT_EQ(ac.base_count_at(9, 0, 0), 4);
  EXPECT_EQ(ac.base_count_at(9, 0, 1), 0);
  EXPECT_EQ(ac.base_count_at(9, 0, 2), 0);
  EXPECT_EQ(ac.base_count_at(9, 0, 3), 0);
  EXPECT_EQ(ac.base_count_at(9, 0, 4), 0);

  EXPECT_EQ(ac.max_ins_at(10), 0);
  EXPECT_EQ(ac.base_count_at(10, 0, 0), 0);
  EXPECT_EQ(ac.base_count_at(10, 0, 1), 0);
  EXPECT_EQ(ac.base_count_at(10, 0, 2), 3);
  EXPECT_EQ(ac.base_count_at(10, 0, 3), 0);
  EXPECT_EQ(ac.base_count_at(10, 0, 4), 0);

  EXPECT_EQ(ac.max_ins_at(11), 0);
  EXPECT_EQ(ac.base_count_at(11, 0, 0), 0);
  EXPECT_EQ(ac.base_count_at(11, 0, 1), 0);
  EXPECT_EQ(ac.base_count_at(11, 0, 2), 3);
  EXPECT_EQ(ac.base_count_at(11, 0, 3), 0);
  EXPECT_EQ(ac.base_count_at(11, 0, 4), 0);
}