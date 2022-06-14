#include <assert.h>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <iterator>
#include <math.h>
#include <stdexcept>
#include <stdio.h>  /* printf, scanf, puts, NULL */
#include <stdlib.h> /* srand, rand */
#include <string>
#include <time.h> /* time */
#include <unordered_set>
#include <utility>
#include <chrono>
#include <set>

#include "align_reads/aligner.hpp"
#include "bioparser/fasta_parser.hpp"
#include "bioparser/fastq_parser.hpp"
#include "biosoup/nucleic_acid.hpp"
#include "edlib.h"
#include "ram/minimizer_engine.hpp"
#include "thread_pool/thread_pool.hpp"

#define EXTRACT nanosim_extract
#define RANDOM_SEED 422
#define OVLP_THRES 0
#define TARGET_NUM 1000
#define MINHASH_BOOL false

#define TYPE_PREFIX 0
#define TYPE_SUFFIX 1
#define TYPE_OTHERS 2

std::atomic<std::uint32_t> biosoup::NucleicAcid::num_objects{0};

namespace align_reads
{

    // ------------------------------------

    struct seq_info
    {
        std::uint32_t id;
        bool forward;
        bool aligned;
        std::uint32_t start;
        std::uint32_t end;
        std::uint32_t left_extend;
        std::uint32_t right_extend;
        std::string contig;
        std::string name;
        std::uint32_t overlap_len;

        std::uint32_t idx_in_sequences;

        void print()
        {
            std::cout << "id: " << id << std::endl;
            std::cout << "forward?: " << forward << std::endl;
            std::cout << "aligned?: " << aligned << std::endl;
            std::cout << "start: " << start << std::endl;
            std::cout << "end: " << end << std::endl;
            std::cout << "left extend: " << left_extend << std::endl;
            std::cout << "right extend: " << right_extend << std::endl;
            std::cout << "contig: " << contig << std::endl;
        }
    };

    struct overlap_info
    {
        seq_info query_info;
        seq_info target_info;

        std::uint32_t q_clip_left;
        std::uint32_t q_clip_right;
        std::uint32_t t_clip_left;
        std::uint32_t t_clip_right;
        bool is_reverse = false;
        std::uint32_t overlap_len;
        std::uint32_t score = 0;

        std::uint32_t q_o_begin;
        std::uint32_t q_o_end;
        std::uint32_t t_o_begin;
        std::uint32_t t_o_end;

        overlap_info(seq_info q, seq_info t) : query_info(q), target_info(t) {}
    };

    struct comp
    {
        bool operator()(const overlap_info &lhs, const overlap_info &rhs)
        {
            if (lhs.query_info.id == rhs.query_info.id)
            {
                return lhs.target_info.id < rhs.target_info.id;
            }
            else if (lhs.query_info.id == rhs.target_info.id)
            {
                return lhs.target_info.id < rhs.query_info.id;
            }
            else
            {
                return lhs.query_info.id < rhs.query_info.id;
            }
        }
    };

    /*
    struct comp {
        bool operator() (const overlap_info& lhs, const overlap_info& rhs) {
            if (lhs.query_info.id != rhs.query_info.id) return lhs.query_info.id < rhs.query_info.id;
            return lhs.target_info.id < rhs.target_info.id;
        }
    };
    */

    struct id_pair
    {
        std::uint32_t id1 = 0;
        std::uint32_t id2 = 0;
        std::uint32_t overlap_len = 0;
        id_pair(std::uint32_t id1, std::uint32_t id2, std::uint32_t overlap_len) : id1(id1), id2(id2), overlap_len(overlap_len) {}
    };

    struct comp2
    {
        bool operator()(const id_pair &lhs, const id_pair &rhs)
        {
            if (lhs.id1 == rhs.id1)
            {
                return lhs.id2 < rhs.id2;
            }
            else if (lhs.id1 == rhs.id2)
            {
                return lhs.id2 < rhs.id1;
            }
            else
            {
                return lhs.id1 < rhs.id1;
            }
        }
    };

    std::set<overlap_info, comp> all_true_overlaps;
    std::set<overlap_info, comp> tps;
    std::set<overlap_info, comp> fps;
    std::set<overlap_info, comp> fns;
    std::set<overlap_info, comp> filtered_t;

    std::set<id_pair, comp2> accepts;
    std::set<id_pair, comp2> rejects;

    std::unordered_map<std::uint32_t, std::pair<std::uint16_t, std::uint16_t>> true_has;   
    std::unordered_map<std::uint32_t, std::pair<std::uint16_t, std::uint16_t>> true_lack;
    std::unordered_map<std::uint32_t, std::pair<std::uint16_t, std::uint16_t>> found_has;
    std::unordered_map<std::uint32_t, std::pair<std::uint16_t, std::uint16_t>> found_lack;

    
    std::unordered_map<std::uint32_t, std::pair<std::uint16_t, std::uint16_t>> found_lack_false;

    std::unordered_map<std::uint32_t, std::pair<std::uint16_t, std::uint16_t>> found_has_false;
    std::uint32_t num_true_prefix_o = 0;
    std::uint32_t num_true_suffix_o = 0;

    std::uint8_t overlap_type(std::uint32_t t_start, std::uint32_t t_end,
                              std::uint32_t q_start, std::uint32_t q_end)
    {
        if (q_start < t_start)
        {
            if (q_end > t_start && q_end < t_end)
                return TYPE_PREFIX;
        }
        else if (q_start < t_end)
        {
            if (q_end > t_end)
                return TYPE_SUFFIX;
        }
        return TYPE_OTHERS;
    }

    std::uint8_t overlap_type(const seq_info &target_info, const seq_info &query_info, bool t_forward)
    {
        std::uint8_t type = overlap_type(target_info.start, target_info.end, query_info.start, query_info.end);
        if (!t_forward)
        {
            if (type == TYPE_PREFIX)
            {
                return TYPE_SUFFIX;
            }
            if (type == TYPE_SUFFIX)
            {
                return TYPE_PREFIX;
            }
        }
        return type;
    }

    static seq_info nanosim_extract(std::string name)
    {
        //>NC-004354_15504699;aligned_0_R_9_1481_32
        //>NC-004354_8681245_unaligned_17149_F_0_1489_0

        seq_info info;

        std::uint32_t idx = 0;
        std::string temp;
        // extract contig name
        while (name[idx] != '_')
        {
            temp += name[idx++];
        }
        info.contig = temp;
        temp.clear();
        idx++;

        // extract start
        while (name[idx] != ';' && name[idx] != '_')
        {
            temp += name[idx++];
        }
        info.start = std::stoi(temp);
        temp.clear();
        idx++;
        if (name[idx] == 'a')
        {
            info.aligned = true;
        }
        else
        {
            info.aligned = false;
        }
        while (name[idx] != '_')
        {
            idx++;
        }
        idx++;
        // extract id
        while (name[idx] != '_')
        {
            temp += name[idx++];
        }
        info.id = std::stoi(temp);
        temp.clear();
        idx++;

        // extract forward/reverse
        if (name[idx] == 'F')
        {
            info.forward = true;
        }
        else
        {
            info.forward = false;
        }
        idx += 2;

        while (name[idx] != '_')
        {
            temp += name[idx++];
        }
        info.left_extend = stoi(temp);
        temp.clear();
        idx++;

        while (name[idx] != '_')
        {
            temp += name[idx++];
        }
        info.end = info.start + std::stoi(temp);
        temp.clear();
        idx++;
        while (idx != name.size())
        {
            temp += name[idx++];
        }
        info.right_extend = std::stoi(temp);
        info.name = std::move(name);
        return info;
    }

    static seq_info seqreq_extract(std::string name)
    {
        //>read=1,forward,position=19104036-19118431,length=14395,NT_033777.3
        seq_info info;

        std::uint32_t idx = 5;
        std::string temp;
        // extract read id string
        while (name[idx] != ',')
        {
            temp += name[idx++];
        }
        info.id = std::stoi(temp);
        temp.clear();
        idx++;

        // extract forward or reverse
        while (name[idx] != ',')
        {
            temp += name[idx++];
        }
        if (temp.compare("forward") == 0)
        {
            info.forward = true;
        }
        else
        {
            info.forward = false;
        }
        temp.clear();
        idx += 10;

        // extract start
        while (name[idx] != '-')
        {
            temp += name[idx++];
        }
        info.start = std::stoi(temp);
        temp.clear();
        idx++;

        // extract end
        while (name[idx] != ',')
        {
            temp += name[idx++];
        }
        info.end = std::stoi(temp);
        temp.clear();
        idx++;
        while (name[idx] != ',')
        {
            idx++;
        }
        idx++;

        // extract contig
        while (idx < name.size())
        {
            temp += name[idx++];
        }
        info.contig = std::move(temp);

        return info;
    }

    std::uint32_t overlap_len(seq_info query_info, seq_info target_info)
    {
        if (query_info.contig != target_info.contig)
            return 0;
        if (query_info.start >= target_info.start && query_info.start < target_info.end)
        {
            return query_info.end < target_info.end ? query_info.end - query_info.start : target_info.end - query_info.start;
        }
        else if (query_info.end <= target_info.end && query_info.end > target_info.start)
        {
            return query_info.start < target_info.start ? query_info.end - target_info.start : query_info.end - query_info.start;
        }
        else if (target_info.start >= query_info.start && target_info.end <= query_info.end)
        {
            return target_info.end - target_info.start;
        }
        return 0;
    }

    bool locality_bad(EdlibAlignResult &result)
    {
        std::uint16_t window_size = 200;
        std::uint16_t thres = 0.3 * window_size;
        std::uint16_t count;

        for (std::uint32_t i = 0; i < result.alignmentLength; i++)
        {
            if (i % window_size == 0)
                count = 0;
            if (result.alignment[i] != 0)
                count++;
            if (count >= thres)
                return true;
        }

        return false;
    }

    //------------------------------------------------

    void Aligner::filter_overlaps(std::vector<biosoup::Overlap> &overlaps)
    {
        std::uint16_t overhang_thres = 500;
        double dist_thres = 0.2;
        std::vector<biosoup::Overlap> filtered;
        filtered.reserve(overlaps.size());
        /*std::uint16_t num_prefix_pass = 0;
        std::uint16_t num_suffix_pass = 0;
        std::vector<biosoup::Overlap> failed_prefix_o;
        std::vector<biosoup::Overlap> failed_suffix_o;
        failed_prefix_o.reserve(overlaps.size() / 10);
        failed_suffix_o.reserve(overlaps.size() / 10);
        std::uint16_t thres = 1;*/
        std::set<overlap_info, comp> non_filtered_t;
        for (auto &o : overlaps)
        {

            auto &query_sequence = sequences[id_to_pos_index[o.rhs_id]];
            auto &target_sequence = sequences[id_to_pos_index[o.lhs_id]];
            
            seq_info query_info = EXTRACT(query_sequence->name);
            seq_info target_info = EXTRACT(target_sequence->name);
            target_info.idx_in_sequences = id_to_pos_index[o.lhs_id];
            query_info.idx_in_sequences = id_to_pos_index[o.rhs_id];
            std::uint32_t ovlp_len = overlap_len(query_info, target_info);
            overlap_info o_info {query_info, target_info};
            
            std::uint32_t t_o_begin = o.lhs_begin;
            std::uint32_t t_o_end = o.lhs_end;
            std::uint32_t q_o_begin = o.rhs_begin;
            
            std::uint32_t q_o_end = o.rhs_end;

            if (!o.strand)
            {
                o_info.is_reverse = true;
                query_sequence->ReverseAndComplement();
                q_o_end = query_sequence->inflated_len - o.rhs_begin;
                q_o_begin = query_sequence->inflated_len - o.rhs_end;
            }
            int protrude_left = q_o_begin - t_o_begin;
            int protrude_right = (query_sequence->inflated_len - q_o_end) - (target_sequence->inflated_len - t_o_end);

            std::uint32_t q_clip_left;
            std::uint32_t q_clip_right;
            std::uint32_t t_clip_left;
            std::uint32_t t_clip_right;
            if (protrude_left > 0)
            {
                q_clip_left = protrude_left;
                t_clip_left = 0;
            }
            else
            {
                q_clip_left = 0;
                t_clip_left = -protrude_left;
            }
            if (protrude_right > 0)
            {
                q_clip_right = protrude_right;
                t_clip_right = 0;
            }
            else
            {
                q_clip_right = 0;
                t_clip_right = -protrude_right;
            }

            /*
            if (protrude_left > 0) {
                if (protrude_right < 0) num_prefix++;
            } else if (protrude_right > 0) {
                num_suffix++;
            }
            */

            o_info.q_clip_left = q_clip_left;
            o_info.q_clip_right = q_clip_right;
            o_info.t_clip_left = t_clip_left;
            o_info.t_clip_right = t_clip_right;

            o_info.q_o_begin = q_o_begin;
            o_info.q_o_end = q_o_end;
            o_info.t_o_begin = t_o_begin;
            o_info.t_o_end = t_o_end;
            o_info.score = o.score;

            std::uint32_t left_overhang = q_o_begin - q_clip_left;
            std::uint32_t right_overhang = query_sequence->inflated_len - q_clip_right - q_o_end;

            std::uint32_t larger = std::max(left_overhang, right_overhang);

            std::uint32_t query_string_len = query_sequence->inflated_len - q_clip_left - q_clip_right;
            std::uint32_t target_string_len = target_sequence->inflated_len - t_clip_left - t_clip_right;

            // double overhang_ratio = (double)(left_overhang + right_overhang) / std::min(query_string_len, target_string_len);
            double norm_score = (double)o.score / std::max(o.rhs_end - o.rhs_begin, o.lhs_end - o.lhs_begin);
            std::string target_string = target_sequence->InflateData(t_clip_left, target_string_len);
            std::string query_string = query_sequence->InflateData(q_clip_left, query_string_len);

            bool declared_bad = false;

            if (false /*norm_score < 0.08*/)
            {
                declared_bad = true;
            }
            else
            {
                if (larger > overhang_thres)
                {
                    declared_bad = true;
                }

                if (declared_bad)
                {
                    bool left_ok = true;
                    bool right_ok = true;
                    if (left_overhang > overhang_thres)
                    {
                        left_ok = false;
                        EdlibAlignResult result = edlibAlign(query_string.c_str() + 100, 100,
                                                             target_string.c_str(), 300,
                                                             edlibNewAlignConfig(-1, EDLIB_MODE_HW, EDLIB_TASK_DISTANCE, NULL, 0));
                        if ((double)result.editDistance / 100 < dist_thres)
                            left_ok = true;
                        edlibFreeAlignResult(result);
                    }

                    if (right_overhang > overhang_thres)
                    {
                        right_ok = false;
                        EdlibAlignResult result = edlibAlign(query_string.c_str() + query_string.size() - 200, 100,
                                                             target_string.c_str() + target_string.size() - 300, 300,
                                                             edlibNewAlignConfig(-1, EDLIB_MODE_HW, EDLIB_TASK_DISTANCE, NULL, 0));
                        if ((double)result.editDistance / 100 < dist_thres)
                            right_ok = true;
                        edlibFreeAlignResult(result);
                    }
                    if (left_ok && right_ok)
                        declared_bad = false;
                }
            }

            if (!declared_bad)
            {
                filtered.push_back(o);
                /*if (is_prefix) {
                    num_prefix_pass++;
                } else if (is_suffix) {
                    num_suffix_pass++;
                }*/
            } /*else {
                if (is_prefix) {
                    failed_prefix_o.push_back(o);
                } else if(is_suffix) {
                    failed_suffix_o.push_back(o);
                }
            }*/
            if (!declared_bad && ovlp_len > OVLP_THRES) non_filtered_t.insert(o_info);
            if (declared_bad && ovlp_len > OVLP_THRES) filtered_t.insert(o_info); 
            if (!o.strand)
                query_sequence->ReverseAndComplement();
        }
        for (auto& i : non_filtered_t) {
            filtered_t.erase(i);
        }

        /*if (num_prefix_pass < thres) {
            filtered.insert(filtered.end(), failed_prefix_o.begin(), failed_prefix_o.begin() + std::min((uint16_t) failed_prefix_o.size(), thres - num_prefix_pass));
        }
        if (num_suffix_pass < thres) {

            filtered.insert(filtered.end(), failed_suffix_o.begin(), failed_suffix_o.begin() + std::min((uint16_t) failed_suffix_o.size(), thres - num_suffix_pass));
        }*/
        filtered.shrink_to_fit();
        overlaps = std::move(filtered);
    }

    void account_prefix_suffix()
    {
        std::uint32_t num_prefix = 0;
        std::uint32_t num_suffix = 0;
        for (auto& o_info: tps) {
            auto type = overlap_type(o_info.target_info, o_info.query_info, o_info.target_info.forward);
            if (type == TYPE_PREFIX) {
                num_prefix++;
            } else if (type == TYPE_SUFFIX) {
                num_suffix++;
            }
        }
        std::cout << "num true prefix: " << num_true_prefix_o << std::endl;
        std::cout << "num true suffix: " << num_true_suffix_o << std::endl;
        std::cout << "num true prefix found:" << num_prefix << std::endl;
        std::cout << "num true suffix found:" << num_suffix << std::endl;
        std::cout << "recall of prefix-suffix: " << (double) (num_prefix + num_suffix) / (num_true_prefix_o + num_true_suffix_o) << std::endl; 
    }

    void Aligner::compare_break() {
        std::cout << "true lack: " << std::endl;
        std::cout << "-------------------------" << std::endl;
        for (auto& item: true_lack) {
            auto info = EXTRACT(sequences[item.first]->name);
            auto len = info.end - info.start;
            std::cout << "num p: " << item.second.first << "num s: " << item.second.second << " len: " << len << std::endl;
        }

        std::cout << "------------------------" << std::endl;
        std::cout << "extra lack in found: " << std::endl;
        std::cout << "--------------------------" << std::endl;

        std::uint32_t num_500 = 0;
        std::uint32_t num_1000 = 0;
        std::uint32_t num_2000 = 0;
        std::uint32_t num_more = 0;
        for (auto& item : found_lack) {
            if (true_lack.find(item.first) != true_lack.end()) {
                continue;
            }
            std::cout << "id " << item.first << std::endl;
            auto info = EXTRACT(sequences[item.first]->name);
            auto len = info.end - info.start;
            if (len <= 500) {
                std::cout << "cat 500" << std::endl;
                num_500++;
                
            } else if (len <= 1000) {
                std::cout << "cat 1000" << std::endl;
                num_1000++;
            } else if (len <= 2000) {
                std::cout << "cat 2000" << std::endl;
                num_2000++;
            } else {
                std::cout << "cat more" << std::endl;
                num_more++;
            }
            std::cout << sequences[item.first]->name << std::endl;
            std::cout << "num p: " << item.second.first << " num s: " << item.second.second << " len: " << len << std::endl;
            auto& item_before = true_has[item.first];
            std::cout << "before: num p: " << item_before.first << " num s: " << item_before.second  << std::endl;
            auto& item_false = found_lack_false[item.first];
            std::cout << "false: num p: " << item_false.first << " num s: " << item_false.second  << std::endl;
 

        }


        std::cout << "num (0, 500]: " << num_500 << std::endl; 
        std::cout << "num (500, 1000]: " << num_1000 << std::endl; 
        std::cout << "num (1000, 2000]: " << num_2000 << std::endl; 
        std::cout << "num (2000, ?]: " << num_more << std::endl; 
        std::cout << "-------------FOUND HAS --------------" << std::endl;
        for (auto& item : found_has) {
            std::cout << "has id " << item.first << std::endl;
            auto info = EXTRACT(sequences[item.first]->name);
            auto len = info.end - info.start;
            std::cout << sequences[item.first]->name << std::endl;
            std::cout << "num p: " << item.second.first << " num s: " << item.second.second << " len: " << len << std::endl;
            auto& item_before = true_has[item.first];
            std::cout << "before: num p: " << item_before.first << " num s: " << item_before.second  << std::endl;
            auto& item_false = found_has_false[item.first];
            std::cout << "false: num p: " << item_false.first << " num s: " << item_false.second  << std::endl;
 

        }


    }


    void Aligner::run()
    {
        // RAM_overlaps_simulated_reads();
        find_true_overlaps();
        // within_each();
        //find_RAM_overlaps(false);
        //compare_break();
        //find_RAM_overlaps(true);
        //filtered_trues();
        //account_prefix_suffix();
        // find_RAM_overlaps_real(true);
        //  true_positives_align_part();
        //  false_positives_align_part();
        //  true_positives();
        // false_positives();
        //false_negatives();
    }

    void Aligner::find_true_overlaps()
    {   
        std::uint32_t min_num_prefix = 10000;
        std::uint32_t min_num_suffix = 10000;
        std::uint16_t num_lacking_any = 0;
        std::uint64_t total_len_lacking = 0;
        std::uint32_t longest_len_lacking = 0;

        std::uint16_t num_200 = 0;
        std::uint16_t num_500 = 0;
        std::uint16_t num_1000 = 0;
        std::uint16_t num_2000 = 0;
        std::uint16_t num_longer = 0;

        std::uint16_t num_o_200 = 0;
        std::uint16_t num_o_500 = 0;
        std::uint16_t num_o_1000 = 0;
        std::uint16_t num_o_2000 = 0;
        std::uint16_t num_o_longer = 0;

        std::uint32_t total_o_len_200 = 0;
        std::uint32_t total_o_len_500 = 0;
        std::uint32_t total_o_len_1000 = 0;
        std::uint32_t total_o_len_2000 = 0;
        std::uint32_t total_o_len_longer = 0;

        for (auto i = 0; i < TARGET_NUM; i++)
        {
            std::uint32_t num_prefix = 0;
            std::uint32_t num_suffix = 0;

            auto target_info = EXTRACT(sequences[i]->name);
            auto s_len = target_info.end - target_info.start;
            
            //
            /*if (i == 432 || i == 538) {
                std::cout << "-----------------" << std::endl;
                std::cout << sequences[i]->name << std::endl;
            }*/
            if (s_len <= 200) {
                num_200++;
            } else if (s_len <= 500) {
                num_500++;
            } else if (s_len <= 1000) {
                num_1000++;
            } else if (s_len <= 2000) {
                num_2000++;
            } else {
                num_longer++;
            }


            for (auto j = 0; j < sequences.size(); j++)
            {
                if (j == i)
                    continue;
                auto query_info = EXTRACT(sequences[j]->name);
                std::uint32_t len = overlap_len(query_info, target_info);

                if (len > OVLP_THRES)
                {
                    query_info.idx_in_sequences = j;
                    target_info.idx_in_sequences = i;
                    overlap_info o_info{query_info, target_info};
                    o_info.overlap_len = len;
                    all_true_overlaps.insert(o_info);

                    auto o_type = overlap_type(target_info, query_info, target_info.forward);
                    if (o_type == TYPE_PREFIX)
                    {
                        num_true_prefix_o++;
                        num_prefix++;
                            
                        //
                        /*if (i == 432 || i == 538) {
                            std::cout << "prefix ovlp len " << len << std::endl;
                            std::cout << query_info.name << std::endl;
                        }*/
                        if (s_len <= 200) {
                            num_o_200++;
                            total_o_len_200 += len;

                        } else if (s_len <= 500) {
                            num_o_500++;
                            total_o_len_500 += len;

                        } else if (s_len <= 1000) {
                            num_o_1000++;
                            total_o_len_1000 += len;

                        } else if (s_len <= 2000) {
                            num_o_2000++;
                            total_o_len_2000 += len;    
                        } else {
                            num_o_longer++;
                            total_o_len_longer += len;    
                        }
                    }
                    else if (o_type == TYPE_SUFFIX)
                    {
                        num_true_suffix_o++;
                        num_suffix++;

                        if (s_len <= 200) {
                            num_o_200++;
                            total_o_len_200 += len;

                        } else if (s_len <= 500) {
                            num_o_500++;
                            total_o_len_500 += len;

                        } else if (s_len <= 1000) {
                            num_o_1000++;
                            total_o_len_1000 += len;

                        } else if (s_len <= 2000) {
                            num_o_2000++;
                            total_o_len_2000 += len;    
                        } else {
                            num_o_longer++;
                            total_o_len_longer += len;    
                        }
                        /*if (i == 432 || i == 538) {
                            //
                            std::cout << "suffix ovlp len " << len << std::endl;
                            std::cout << query_info.name << std::endl;
                        }*/


                    }
                }
            }
            if (num_prefix < min_num_prefix) min_num_prefix = num_prefix;
            if (num_suffix < min_num_suffix) min_num_suffix = num_suffix;
            if (num_prefix == 0 || num_suffix == 0) {
                num_lacking_any++;
                std::uint32_t len = (target_info.end - target_info.start);
                total_len_lacking += len;
                if (len > longest_len_lacking) longest_len_lacking = len;
                true_lack.emplace(i, std::make_pair(num_prefix, num_suffix));
            } else {

                true_has.emplace(i, std::make_pair(num_prefix, num_suffix));
            }
        }
        std::cout << "min num prefix: " << min_num_prefix << std::endl;
        std::cout << "min num suffix: " << min_num_suffix << std::endl;
        std::cout << "avg num prefix: " << (double) num_true_prefix_o / TARGET_NUM << std::endl;
        std::cout << "avg num suffix: " << (double) num_true_suffix_o / TARGET_NUM << std::endl;
        std::cout << "num lacking any in true: "  << num_lacking_any << std::endl;
        std::cout << "avg len lacking in true: " << (double) total_len_lacking / num_lacking_any << std::endl;
        std::cout << "longest len lacking: " << longest_len_lacking << std::endl;
        
        std::cout << "num segment (0, 200]: " << num_200 << std::endl;
        std::cout << "num segment (200, 500]: " << num_500 << std::endl;
        std::cout << "num segment (500, 1000]: " << num_1000 << std::endl;
        std::cout << "num segment (1000, 2000]: " << num_2000 << std::endl;
        std::cout << "num segment (2000, ?]: " << num_longer << std::endl;

        std::cout << "avg num prefix/suffix for (0, 200]: " << (double) num_o_200/num_200  << std::endl;
        std::cout << "avg num prefix/suffix for (200, 500]: " << (double) num_o_500/num_500 << std::endl;
        std::cout << "avg num prefix/suffix for (500, 1000]: " << (double) num_o_1000/num_1000 << std::endl;
        std::cout << "avg num prefix/suffix for (1000, 2000]: " << (double) num_o_2000/num_2000 << std::endl;
        std::cout << "avg num prefix/suffix for (2000, ?]: " << (double) num_o_longer/num_longer << std::endl;
    
        std::cout << "avg prefix/suffix overlap len for (0, 200]: " << (double) total_o_len_200 / num_o_200 << std::endl;
        std::cout << "avg prefix/suffix overlap len for (200, 500]: " << (double) total_o_len_500 / num_o_500 << std::endl;
        std::cout << "avg prefix/suffix overlap len for (500, 1000]: " << (double) total_o_len_1000 / num_o_1000 << std::endl;
        std::cout << "avg prefix/suffix overlap len for (1000, 2000]: " << (double) total_o_len_2000 / num_o_2000 << std::endl;
        std::cout << "avg prefix/suffix overlap len for (2000, ?]: " << (double) total_o_len_longer / num_o_longer << std::endl;
    }

    void Aligner::within_each()
    {

        for (auto i = 0; i < TARGET_NUM; i++)
        {
            auto &target = sequences[i];
            std::uint32_t target_id = target->id;
            auto target_info = EXTRACT(target->name);
            target_info.idx_in_sequences = id_to_pos_index[target_id];
            std::vector<biosoup::Overlap> overlaps = minimizer_engine.Map(target, true, false, true);

            std::uint32_t num_tp = 0;
            std::uint32_t num_fp = 0;

            double total_norm_dist_tp = 0;
            double total_norm_score_whole_tp = 0;
            double total_norm_score_part_tp = 0;

            double total_norm_dist_fp = 0;
            double total_norm_score_whole_fp = 0;
            double total_norm_score_part_fp = 0;

            for (auto &o : overlaps)
            {
                auto &s = sequences[id_to_pos_index[o.rhs_id]];
                auto query_info = EXTRACT(s->name);
                query_info.idx_in_sequences = id_to_pos_index[o.rhs_id];
                overlap_info o_info{query_info, target_info};

                std::uint32_t t_begin = o.lhs_begin;
                std::uint32_t t_end = o.lhs_end;
                std::uint32_t q_begin = o.rhs_begin;
                std::uint32_t q_end = o.rhs_end;

                if (!o.strand)
                {
                    s->ReverseAndComplement();
                    o_info.is_reverse = true;
                    q_end = s->inflated_len - o.rhs_begin;
                    q_begin = s->inflated_len - o.rhs_end;
                }
                int protrude_left = q_begin - t_begin;
                int protrude_right = (s->inflated_len - q_end) - (target->inflated_len - t_end);

                std::uint32_t q_clip_left;
                std::uint32_t q_clip_right;
                std::uint32_t t_clip_left;
                std::uint32_t t_clip_right;
                if (protrude_left > 0)
                {
                    q_clip_left = protrude_left;
                    t_clip_left = 0;
                }
                else
                {
                    q_clip_left = 0;
                    t_clip_left = -protrude_left;
                }
                if (protrude_right > 0)
                {
                    q_clip_right = protrude_right;
                    t_clip_right = 0;
                }
                else
                {
                    q_clip_right = 0;
                    t_clip_right = -protrude_right;
                }

                o_info.q_clip_left = q_clip_left;
                o_info.q_clip_right = q_clip_right;
                o_info.t_clip_left = t_clip_left;
                o_info.t_clip_right = t_clip_right;

                o_info.q_o_begin = q_begin;
                o_info.q_o_end = q_end;
                o_info.t_o_begin = t_begin;
                o_info.t_o_end = t_end;
                o_info.score = o.score;

                std::uint32_t query_string_len = s->inflated_len - q_clip_left - q_clip_right;
                std::uint32_t target_string_len = target->inflated_len - t_clip_left - t_clip_right;
                std::string query_string = s->InflateData(q_clip_left, query_string_len);
                std::string target_string = target->InflateData(t_clip_left, target_string_len);
                EdlibAlignResult result = edlibAlign(query_string.c_str(), query_string.size(),
                                                     target_string.c_str(), target_string.size(),
                                                     edlibNewAlignConfig(-1, EDLIB_MODE_HW, EDLIB_TASK_PATH, NULL, 0));

                if (!o.strand)
                {
                    s->ReverseAndComplement();
                }

                if (overlap_len(query_info, target_info) > OVLP_THRES)
                {
                    tps.insert(o_info);
                    total_norm_dist_tp += (double)result.editDistance / query_string_len;
                    total_norm_score_whole_tp += (double)o.score / std::max(query_string_len, target_string_len);
                    total_norm_score_part_tp += (double)o.score / std::max(o.rhs_end - o.rhs_begin, o.lhs_end - o.rhs_begin);
                    num_tp++;
                }
                else
                {
                    fps.insert(o_info);
                    total_norm_dist_fp += (double)result.editDistance / query_string_len;
                    total_norm_score_whole_fp += (double)o.score / std::max(query_string_len, target_string_len);
                    total_norm_score_part_fp += (double)o.score / std::max(o.rhs_end - o.rhs_begin, o.lhs_end - o.rhs_begin);
                    num_fp++;
                }

                edlibFreeAlignResult(result);
            }
            std::cout << "target_name: " << target_info.name << std::endl;
            std::cout << "-------------------------------------------" << std::endl;
            std::cout << "num tp: " << num_tp << std::endl;
            std::cout << "avg norm dist tp : " << total_norm_dist_tp / num_tp << std::endl;
            std::cout << "avg norm score whole tp: " << total_norm_score_whole_tp / num_tp << std::endl;
            std::cout << "avg norm score part tp: " << total_norm_score_part_tp / num_tp << std::endl;

            std::cout << "-----------------------" << std::endl;
            std::cout << "num fp: " << num_fp << std::endl;
            std::cout << "avg norm dist fp : " << total_norm_dist_fp / num_fp << std::endl;
            std::cout << "avg norm score whole fp: " << total_norm_score_whole_fp / num_fp << std::endl;
            std::cout << "avg norm score part fp: " << total_norm_score_part_fp / num_fp << std::endl;
            std::cout << "-------------------------------------------" << std::endl;
        }

        for (auto &o_info : all_true_overlaps)
        {
            if (tps.find(o_info) == tps.end())
            {
                fns.insert(o_info);
            }
        }

        std::cout << "precision: " << ((double)tps.size() / (tps.size() + fps.size())) << std::endl;
        std::cout << "recall: " << ((double)tps.size() / all_true_overlaps.size()) << std::endl;
    }

    void Aligner::RAM_overlaps_simulated_reads()
    {

        std::uint32_t num_good = 0;
        std::uint32_t num_bad = 0;
        double bad_thres = 0.2;
        double overhang_thres = 500;
        std::uint32_t num_declared_bad = 0;
        std::uint32_t num_declared_correctly = 0;
        std::uint32_t num_declared_wrongly = 0;
        std::uint32_t num_not_declared_correctly = 0;
        std::uint32_t num_not_declared_wrongly = 0;

        double total_good_overhang_ratio = 0;
        double total_bad_overhang_ratio = 0;

        std::uint32_t num_tp = 0;
        std::uint32_t num_fp = 0;
        std::uint32_t num_rejected_tp = 0;
        std::uint32_t num_rejected_fp = 0;
        double total_tp_norm_dist = 0;
        double total_fp_norm_dist = 0;
        std::uint32_t num_good_tp = 0;
        std::uint32_t num_bad_tp = 0;

        std::uint32_t num_locality_good = 0;
        std::uint32_t num_locality_bad = 0;
        std::uint32_t locality_good_reject = 0;
        std::uint32_t locality_bad_reject = 0;
        std::uint32_t locality_good_tp = 0;
        std::uint32_t locality_bad_tp = 0;

        auto sort_by_qid = [](biosoup::Overlap &o1, biosoup::Overlap &o2)
        {
            return o1.rhs_id < o2.rhs_id;
        };

        for (auto i = 0; i < TARGET_NUM; i++)
        {
            auto &target_sequence = sequences[i];
            std::vector<biosoup::Overlap> overlaps = minimizer_engine.Map(target_sequence, true, false, true);
            std::sort(overlaps.begin(), overlaps.end(), sort_by_qid);
            std::uint32_t j = 0;
            auto target_info = EXTRACT(target_sequence->name);
            while (j < overlaps.size())
            {
                std::vector<biosoup::Overlap> for_same_seq;

                while (true)
                {
                    for_same_seq.push_back(overlaps[j]);
                    j++;
                    if (j >= overlaps.size() || overlaps[j].rhs_id != overlaps[j - 1].rhs_id)
                        break;
                }
                for (auto &o : for_same_seq)
                {
                    std::cout << "---------------------------" << std::endl;
                    auto &query_sequence = sequences[id_to_pos_index[o.rhs_id]];
                    auto query_info = EXTRACT(query_sequence->name);
                    // std::cout << "q_o_begin: " << o.rhs_begin << " q_o_end: " << o.rhs_end << std::endl;
                    // std::cout << "t_o_begin: " << o.lhs_begin << " t_o_end: " << o.lhs_end << std::endl;

                    std::uint32_t t_o_begin = o.lhs_begin;
                    std::uint32_t t_o_end = o.lhs_end;
                    std::uint32_t q_o_begin = o.rhs_begin;
                    std::uint32_t q_o_end = o.rhs_end;

                    std::cout << "q_full_len: " << query_sequence->inflated_len << std::endl;
                    std::cout << "t_full_len: " << target_sequence->inflated_len << std::endl;
                    // std::cout << "t_begin: " << o.lhs_begin << " t_end: " << o.lhs_end << std::endl;
                    // std::cout << "q_begin: " << o.rhs_begin << " q_end: " << o.rhs_end << std::endl;

                    if (!o.strand)
                    {
                        query_sequence->ReverseAndComplement();
                        q_o_end = query_sequence->inflated_len - o.rhs_begin;
                        q_o_begin = query_sequence->inflated_len - o.rhs_end;
                    }
                    int protrude_left = q_o_begin - t_o_begin;
                    int protrude_right = (query_sequence->inflated_len - q_o_end) - (target_sequence->inflated_len - t_o_end);

                    std::uint32_t q_clip_left;
                    std::uint32_t q_clip_right;
                    std::uint32_t t_clip_left;
                    std::uint32_t t_clip_right;
                    if (protrude_left > 0)
                    {
                        q_clip_left = protrude_left;
                        t_clip_left = 0;
                    }
                    else
                    {
                        q_clip_left = 0;
                        t_clip_left = -protrude_left;
                    }
                    if (protrude_right > 0)
                    {
                        q_clip_right = protrude_right;
                        t_clip_right = 0;
                    }
                    else
                    {
                        q_clip_right = 0;
                        t_clip_right = -protrude_right;
                    }
                    std::uint32_t left_overhang = q_o_begin - q_clip_left;
                    std::uint32_t right_overhang = query_sequence->inflated_len - q_clip_right - q_o_end;

                    std::uint32_t larger = std::max(left_overhang, right_overhang);

                    std::uint32_t query_string_len = query_sequence->inflated_len - q_clip_left - q_clip_right;
                    std::uint32_t target_string_len = target_sequence->inflated_len - t_clip_left - t_clip_right;

                    double overhang_ratio = (double)(left_overhang + right_overhang) / std::min(query_string_len, target_string_len);

                    std::string target_string = target_sequence->InflateData(t_clip_left, target_string_len);
                    std::string query_string = query_sequence->InflateData(q_clip_left, query_string_len);

                    EdlibAlignResult result = edlibAlign(query_string.c_str(), query_string.size(),
                                                         target_string.c_str(), target_string.size(),
                                                         edlibNewAlignConfig(-1, EDLIB_MODE_HW, EDLIB_TASK_PATH, NULL, 0));

                    double norm_dist = (double)result.editDistance / query_string_len;

                    bool is_locality_bad = locality_bad(result);
                    std::uint32_t ovlp_len = overlap_len(query_info, target_info);
                    bool is_tp = ovlp_len > OVLP_THRES;
                    if (is_tp)
                    {
                        num_tp++;
                        total_tp_norm_dist += norm_dist;
                    }
                    else
                    {
                        num_fp++;
                        total_fp_norm_dist += norm_dist;
                    }

                    bool declared_bad = false;
                    bool bad = false;
                    if (norm_dist > bad_thres)
                    {
                        bad = true;
                        if (is_tp)
                        {
                            num_bad_tp++;
                            std::cout << "BAD TP" << std::endl;
                        }
                    }
                    else
                    {
                        if (is_tp)
                        {
                            num_good_tp++;
                            std::cout << "GOOD TP" << std::endl;
                        }
                    }

                    if (larger > overhang_thres)
                    {
                        declared_bad = true;
                    }

                    if (declared_bad)
                    {
                        bool left_ok = true;
                        bool right_ok = true;
                        if (left_overhang > overhang_thres)
                        {
                            left_ok = false;
                            EdlibAlignResult result = edlibAlign(query_string.c_str() + 100, 100,
                                                                 target_string.c_str(), 300,
                                                                 edlibNewAlignConfig(-1, EDLIB_MODE_HW, EDLIB_TASK_PATH, NULL, 0));
                            if ((double)result.editDistance / 100 < 0.2)
                                left_ok = true;
                            std::cout << "----left ok: " << (double)result.editDistance / 100 << std::endl;
                            edlibFreeAlignResult(result);
                        }

                        if (right_overhang > overhang_thres)
                        {
                            right_ok = false;
                            EdlibAlignResult result = edlibAlign(query_string.c_str() + query_string.size() - 200, 100,
                                                                 target_string.c_str() + target_string.size() - 300, 300,
                                                                 edlibNewAlignConfig(-1, EDLIB_MODE_HW, EDLIB_TASK_PATH, NULL, 0));
                            if ((double)result.editDistance / 100 < 0.2)
                                right_ok = true;
                            std::cout << "----right ok: " << (double)result.editDistance / 100 << std::endl;
                            edlibFreeAlignResult(result);
                        }
                        if (left_ok && right_ok)
                        {
                            declared_bad = false;
                        }
                        else
                        {
                            num_declared_bad++;
                        }
                    }
                    // ---------------
                    // declared_bad = false;
                    //-------------
                    if (declared_bad)
                    {
                        if (is_tp)
                        {
                            num_rejected_tp++;
                            std::cout << "REJECT TP" << std::endl;
                        }
                        else
                        {
                            num_rejected_fp++;
                        }
                    }

                    if (declared_bad)
                    {
                        rejects.emplace(query_info.id, target_info.id, ovlp_len);
                    }
                    else
                    {

                        accepts.emplace(query_info.id, target_info.id, ovlp_len);
                    }

                    if (is_locality_bad)
                    {
                        num_locality_bad++;
                        std::cout << "LOCALITY BAD ";
                        if (declared_bad)
                        {
                            locality_bad_reject++;
                            std::cout << "REJECT";
                        }
                        else
                        {
                            std::cout << "ACCEPT";
                        }
                        if (is_tp)
                        {
                            locality_bad_tp++;
                            std::cout << "LOCALITY BAD TP" << std::endl;
                        }
                    }
                    else
                    {
                        num_locality_good++;
                        std::cout << "LOCALITY GOOD ";
                        if (declared_bad)
                        {
                            locality_good_reject++;
                            std::cout << "REJECT";
                        }
                        else
                        {
                            std::cout << "ACCEPT";
                        }
                        if (is_tp)
                        {
                            locality_good_tp++;
                            std::cout << "LOCALITY GOOD TP" << std::endl;
                        }
                    }
                    std::cout << std::endl;

                    if (bad)
                    {
                        num_bad++;
                        total_bad_overhang_ratio += overhang_ratio;
                        if (declared_bad)
                        {
                            num_declared_correctly++;
                            std::cout << "CORRECTLY DECLARED" << std::endl;
                        }
                        else
                        {
                            num_not_declared_wrongly++;
                            std::cout << "WRONGLY NOT DECLARED" << std::endl;
                        }
                    }
                    else
                    {
                        num_good++;
                        total_good_overhang_ratio += overhang_ratio;
                        if (declared_bad)
                        {
                            num_declared_wrongly++;
                            std::cout << "WRONGLY DECLARED" << std::endl;
                        }
                        else
                        {
                            num_not_declared_correctly++;
                            std::cout << "CORRECTLY NOT DECLARED" << std::endl;
                        }
                    }
                    if (is_tp)
                    {
                        std::cout << "TP" << std::endl;
                    }
                    else
                    {
                        std::cout << "FP" << std::endl;
                    }
                    std::cout << "qid: " << query_sequence->id << " tid: " << target_sequence->id << std::endl;
                    std::cout << "norm dist " << norm_dist << std::endl;
                    std::cout << "reverse?: " << !o.strand << std::endl;
                    std::cout << "q_clip_left: " << q_clip_left << " query_clip_right: " << q_clip_right << std::endl;
                    std::cout << "t_clip_left: " << t_clip_left << " target_clip_right: " << t_clip_right << std::endl;
                    std::cout << "q_o_begin: " << q_o_begin << " q_o_end: " << q_o_end << std::endl;
                    std::cout << "t_o_begin: " << t_o_begin << " t_o_end: " << t_o_end << std::endl;
                    std::cout << "clipped q len: " << query_string.size() << " clipped t len: " << target_string.size() << std::endl;
                    std::cout << "left overhang: " << left_overhang << std::endl;
                    std::cout << "right overhang: " << right_overhang << std::endl;
                    std::cout << "overlap len: " << overlap_len(query_info, target_info) << std::endl;

                    std::vector<std::string> overlapping_reads;
                    overlapping_reads.push_back(std::move(query_string));
                    std::vector<std::pair<std::uint32_t, std::uint32_t>> clips;
                    clips.emplace_back(0, 0);
                    std::vector<EdlibAlignResult> rs = {result};
                    auto alignment = multi_align(overlapping_reads, target_string, clips, rs, false);
                    alignment.print();
                    // edlibFreeAlignResult(result);
                    if (!o.strand)
                        query_sequence->ReverseAndComplement();
                }
            }
        }
        for (auto &e : accepts)
        {
            rejects.erase(e);
        }
        std::uint32_t pair_tp = 0;
        std::uint32_t pair_fp = 0;
        std::uint32_t pair_tn = 0;
        std::uint32_t pair_fn = 0;
        for (auto &e : accepts)
        {
            if (e.overlap_len > OVLP_THRES)
            {
                pair_tp++;
            }
            else
            {
                pair_fp++;
            }
        }

        for (auto &e : rejects)
        {
            if (e.overlap_len > OVLP_THRES)
            {
                pair_fn++;
            }
            else
            {
                pair_tn++;
            }
        }
        std::cout << "pair tp : " << pair_tp << " pair fp: " << pair_fp << std::endl;
        std::cout << "pair tn: " << pair_tn << " pair fn: " << pair_fn << std::endl;
        std::cout << "pair prec: " << (double)pair_tp / (pair_tp + pair_fp) << std::endl;
        std::cout << "pair recall: " << (double)pair_tp / (pair_tp + pair_fn) << std::endl;

        std::cout << "---------------" << std::endl;
        std::cout << "num good: " << num_good << std::endl;
        std::cout << "num bad: " << num_bad << std::endl;
        std::cout << "avg good overhang ratio: " << total_good_overhang_ratio / num_good << std::endl;
        std::cout << "avg bad overhang ratio: " << total_bad_overhang_ratio / num_bad << std::endl;
        std::cout << "num declared bad: " << num_declared_bad << std::endl;
        std::cout << "num declared correctly: " << num_declared_correctly << std::endl;
        std::cout << "num declared wrongly: " << num_declared_wrongly << std::endl;
        std::cout << "num not declared correctly: " << num_not_declared_correctly << std::endl;
        std::cout << "num not declared wrongly: " << num_not_declared_wrongly << std::endl;
        std::cout << "recall of good: " << (double)num_not_declared_correctly / (num_not_declared_correctly + num_declared_wrongly) << std::endl;
        std::cout << "precision of good: " << (double)num_not_declared_correctly / (num_not_declared_correctly + num_not_declared_wrongly) << std::endl;
        std::cout << "-----------------" << std::endl;
        std::cout << "num tp: " << num_tp << " num fp: " << num_fp << std::endl;
        std::cout << "avg norm dist tp: " << (double)total_tp_norm_dist / num_tp << " avg norm dist fp: " << (double)total_fp_norm_dist / num_fp << std::endl;
        std::cout << "num rejected tp: " << num_rejected_tp << " num rejected fp: " << num_rejected_fp << std::endl;
        std::cout << "num good tp: " << num_good_tp << " num bad tp:" << num_bad_tp << std::endl;
        std::cout << "---------------" << std::endl;
        std::cout << "num locality good: " << num_locality_good << std::endl;
        std::cout << "num locality bad: " << num_locality_bad << std::endl;
        std::cout << "num locality good reject: " << locality_good_reject << std::endl;
        std::cout << "num locality bad reject: " << locality_bad_reject << std::endl;
        std::cout << "num locality good tp: " << locality_good_tp << std::endl;
        std::cout << "num locality bad tp: " << locality_bad_tp << std::endl;
    }

    void Aligner::RAM_overlaps_true_reads()
    {

        std::uint32_t num_good = 0;
        std::uint32_t num_bad = 0;
        double bad_thres = 0.2;
        double overhang_thres = 500;
        std::uint32_t num_declared_bad = 0;
        std::uint32_t num_declared_correctly = 0;
        std::uint32_t num_declared_wrongly = 0;
        std::uint32_t num_not_declared_correctly = 0;
        std::uint32_t num_not_declared_wrongly = 0;

        std::uint32_t num_locality_good = 0;
        std::uint32_t num_locality_bad = 0;
        std::uint32_t locality_good_reject = 0;
        std::uint32_t locality_bad_reject = 0;

        double total_good_overhang_ratio = 0;
        double total_bad_overhang_ratio = 0;

        auto sort_by_qid = [](biosoup::Overlap &o1, biosoup::Overlap &o2)
        {
            return o1.rhs_id < o2.rhs_id;
        };

        for (auto i = 0; i < TARGET_NUM; i++)
        {
            auto &target_sequence = sequences[i];
            std::vector<biosoup::Overlap> overlaps = minimizer_engine.Map(target_sequence, true, false, true);
            std::sort(overlaps.begin(), overlaps.end(), sort_by_qid);
            std::uint32_t j = 0;
            while (j < overlaps.size())
            {
                std::vector<biosoup::Overlap> for_same_seq;

                while (true)
                {
                    for_same_seq.push_back(overlaps[j]);
                    j++;
                    if (j >= overlaps.size() || overlaps[j].rhs_id != overlaps[j - 1].rhs_id)
                        break;
                }
                for (auto &o : for_same_seq)
                {
                    std::cout << "---------------------------" << std::endl;
                    auto &query_sequence = sequences[id_to_pos_index[o.rhs_id]];
                    // std::cout << "q_o_begin: " << o.rhs_begin << " q_o_end: " << o.rhs_end << std::endl;
                    // std::cout << "t_o_begin: " << o.lhs_begin << " t_o_end: " << o.lhs_end << std::endl;

                    std::uint32_t t_o_begin = o.lhs_begin;
                    std::uint32_t t_o_end = o.lhs_end;
                    std::uint32_t q_o_begin = o.rhs_begin;
                    std::uint32_t q_o_end = o.rhs_end;

                    std::cout << "q_full_len: " << query_sequence->inflated_len << std::endl;
                    std::cout << "t_full_len: " << target_sequence->inflated_len << std::endl;
                    // std::cout << "t_begin: " << o.lhs_begin << " t_end: " << o.lhs_end << std::endl;
                    // std::cout << "q_begin: " << o.rhs_begin << " q_end: " << o.rhs_end << std::endl;

                    if (!o.strand)
                    {
                        query_sequence->ReverseAndComplement();
                        q_o_end = query_sequence->inflated_len - o.rhs_begin;
                        q_o_begin = query_sequence->inflated_len - o.rhs_end;
                    }
                    int protrude_left = q_o_begin - t_o_begin;
                    int protrude_right = (query_sequence->inflated_len - q_o_end) - (target_sequence->inflated_len - t_o_end);

                    std::uint32_t q_clip_left;
                    std::uint32_t q_clip_right;
                    std::uint32_t t_clip_left;
                    std::uint32_t t_clip_right;
                    if (protrude_left > 0)
                    {
                        q_clip_left = protrude_left;
                        t_clip_left = 0;
                    }
                    else
                    {
                        q_clip_left = 0;
                        t_clip_left = -protrude_left;
                    }
                    if (protrude_right > 0)
                    {
                        q_clip_right = protrude_right;
                        t_clip_right = 0;
                    }
                    else
                    {
                        q_clip_right = 0;
                        t_clip_right = -protrude_right;
                    }
                    std::uint32_t left_overhang = q_o_begin - q_clip_left;
                    std::uint32_t right_overhang = query_sequence->inflated_len - q_clip_right - q_o_end;

                    std::uint32_t larger = std::max(left_overhang, right_overhang);

                    std::uint32_t query_string_len = query_sequence->inflated_len - q_clip_left - q_clip_right;
                    std::uint32_t target_string_len = target_sequence->inflated_len - t_clip_left - t_clip_right;

                    double overhang_ratio = (double)(left_overhang + right_overhang) / std::min(query_string_len, target_string_len);

                    std::string target_string = target_sequence->InflateData(t_clip_left, target_string_len);
                    std::string query_string = query_sequence->InflateData(q_clip_left, query_string_len);

                    EdlibAlignResult result = edlibAlign(query_string.c_str(), query_string.size(),
                                                         target_string.c_str(), target_string.size(),
                                                         edlibNewAlignConfig(-1, EDLIB_MODE_HW, EDLIB_TASK_PATH, NULL, 0));

                    double norm_dist = (double)result.editDistance / query_string_len;
                    bool is_locality_bad = locality_bad(result);

                    bool declared_bad = false;
                    bool bad = false;
                    if (norm_dist > bad_thres)
                        bad = true;

                    if (larger > overhang_thres)
                    {
                        declared_bad = true;
                    }

                    if (declared_bad)
                    {
                        bool left_ok = true;
                        bool right_ok = true;
                        if (left_overhang > overhang_thres)
                        {
                            left_ok = false;
                            EdlibAlignResult result = edlibAlign(query_string.c_str() + 100, 100,
                                                                 target_string.c_str(), 300,
                                                                 edlibNewAlignConfig(-1, EDLIB_MODE_HW, EDLIB_TASK_PATH, NULL, 0));
                            if ((double)result.editDistance / 100 < 0.2)
                                left_ok = true;
                            std::cout << "----left ok: " << (double)result.editDistance / 100 << std::endl;
                            edlibFreeAlignResult(result);
                        }

                        if (right_overhang > overhang_thres)
                        {
                            right_ok = false;
                            EdlibAlignResult result = edlibAlign(query_string.c_str() + query_string.size() - 200, 100,
                                                                 target_string.c_str() + target_string.size() - 300, 300,
                                                                 edlibNewAlignConfig(-1, EDLIB_MODE_HW, EDLIB_TASK_PATH, NULL, 0));
                            if ((double)result.editDistance / 100 < 0.2)
                                right_ok = true;
                            std::cout << "----right ok: " << (double)result.editDistance / 100 << std::endl;
                            edlibFreeAlignResult(result);
                        }
                        if (left_ok && right_ok)
                        {
                            declared_bad = false;
                        }
                        else
                        {
                            num_declared_bad++;
                        }
                    }

                    if (is_locality_bad)
                    {
                        num_locality_bad++;
                        std::cout << "LOCALITY BAD ";
                        if (declared_bad)
                        {
                            locality_bad_reject++;
                            std::cout << "REJECT";
                        }
                        else
                        {
                            std::cout << "ACCEPT";
                        }
                    }
                    else
                    {
                        num_locality_good++;
                        std::cout << "LOCALITY GOOD ";
                        if (declared_bad)
                        {
                            locality_good_reject++;
                            std::cout << "REJECT";
                        }
                        else
                        {
                            std::cout << "ACCEPT";
                        }
                    }
                    std::cout << std::endl;

                    if (bad)
                    {
                        num_bad++;
                        total_bad_overhang_ratio += overhang_ratio;
                        if (declared_bad)
                        {
                            num_declared_correctly++;
                            std::cout << "CORRECTLY DECLARED" << std::endl;
                        }
                        else
                        {
                            num_not_declared_wrongly++;
                            std::cout << "WRONGLY NOT DECLARED" << std::endl;
                        }
                    }
                    else
                    {
                        num_good++;
                        total_good_overhang_ratio += overhang_ratio;
                        if (declared_bad)
                        {
                            num_declared_wrongly++;
                            std::cout << "WRONGLY DECLARED" << std::endl;
                        }
                        else
                        {
                            num_not_declared_correctly++;
                            std::cout << "CORRECTLY NOT DECLARED" << std::endl;
                        }
                    }

                    std::cout << "qid: " << query_sequence->id << " tid: " << target_sequence->id << std::endl;
                    std::cout << "norm dist " << norm_dist << std::endl;
                    std::cout << "reverse?: " << !o.strand << std::endl;
                    std::cout << "q_clip_left: " << q_clip_left << " query_clip_right: " << q_clip_right << std::endl;
                    std::cout << "t_clip_left: " << t_clip_left << " target_clip_right: " << t_clip_right << std::endl;
                    std::cout << "q_o_begin: " << q_o_begin << " q_o_end: " << q_o_end << std::endl;
                    std::cout << "t_o_begin: " << t_o_begin << " t_o_end: " << t_o_end << std::endl;
                    std::cout << "clipped q len: " << query_string.size() << " clipped t len: " << target_string.size() << std::endl;
                    std::cout << "left overhang: " << left_overhang << std::endl;
                    std::cout << "right overhang: " << right_overhang << std::endl;

                    std::vector<std::string> overlapping_reads;
                    overlapping_reads.push_back(std::move(query_string));
                    std::vector<std::pair<std::uint32_t, std::uint32_t>> clips;
                    clips.emplace_back(0, 0);
                    std::vector<EdlibAlignResult> rs = {result};
                    auto alignment = multi_align(overlapping_reads, target_string, clips, rs, false);
                    alignment.print();
                    // edlibFreeAlignResult(result);
                    if (!o.strand)
                        query_sequence->ReverseAndComplement();
                }
            }
        }
        std::cout << "num good: " << num_good << std::endl;
        std::cout << "num bad: " << num_bad << std::endl;
        std::cout << "avg good overhang ratio: " << total_good_overhang_ratio / num_good << std::endl;
        std::cout << "avg bad overhang ratio: " << total_bad_overhang_ratio / num_bad << std::endl;
        std::cout << "num declared bad: " << num_declared_bad << std::endl;
        std::cout << "num declared correctly: " << num_declared_correctly << std::endl;
        std::cout << "num declared wrongly: " << num_declared_wrongly << std::endl;
        std::cout << "num not declared correctly: " << num_not_declared_correctly << std::endl;
        std::cout << "num not declared wrongly: " << num_not_declared_wrongly << std::endl;
        std::cout << "recall of good: " << (double)num_not_declared_correctly / (num_not_declared_correctly + num_declared_wrongly) << std::endl;
        std::cout << "precision of good: " << (double)num_not_declared_correctly / (num_not_declared_correctly + num_not_declared_wrongly) << std::endl;
        std::cout << "num locality good: " << num_locality_good << std::endl;
        std::cout << "num locality bad: " << num_locality_bad << std::endl;
        std::cout << "num locality good reject: " << locality_good_reject << std::endl;
        std::cout << "num locality bad reject: " << locality_bad_reject << std::endl;
    }
    void Aligner::find_RAM_overlaps_real(bool filter)
    {
        std::uint32_t num_lacking = 0;

        if (filter)
        {
            std::cout << "filter" << std::endl;
        }
        else
        {
            std::cout << "no filter" << std::endl;
        }
        for (auto i = 0; i < TARGET_NUM; i++)
        {
            std::uint32_t num_prefix = 0;
            std::uint32_t num_suffix = 0;

            auto &target = sequences[i];
            std::uint32_t target_id = target->id;
            std::vector<biosoup::Overlap> overlaps = minimizer_engine.Map(target, true, false, MINHASH_BOOL);
            if (filter)
            {
                filter_overlaps(overlaps);
            }
            for (auto &o : overlaps)
            {
                auto &s = sequences[id_to_pos_index[o.rhs_id]];

                std::uint32_t t_begin = o.lhs_begin;
                std::uint32_t t_end = o.lhs_end;
                std::uint32_t q_begin = o.rhs_begin;
                std::uint32_t q_end = o.rhs_end;

                if (!o.strand)
                {
                    q_end = s->inflated_len - o.rhs_begin;
                    q_begin = s->inflated_len - o.rhs_end;
                }
                int protrude_left = q_begin - t_begin;
                int protrude_right = (s->inflated_len - q_end) - (target->inflated_len - t_end);

                std::uint32_t q_clip_left;
                std::uint32_t q_clip_right;
                std::uint32_t t_clip_left;
                std::uint32_t t_clip_right;
                if (protrude_left > 0)
                {
                    q_clip_left = protrude_left;
                    t_clip_left = 0;
                }
                else
                {
                    q_clip_left = 0;
                    t_clip_left = -protrude_left;
                }
                if (protrude_right > 0)
                {
                    q_clip_right = protrude_right;
                    t_clip_right = 0;
                }
                else
                {
                    q_clip_right = 0;
                    t_clip_right = -protrude_right;
                }

                /* if (protrude_left == 0 && protrude_right ==0) {
                     num_contain++;

                 } else if (protrude_left > 0 && protrude_right > 0) {
                     num_contained++;
                 } else if (protrude_left > 0) {
                     num_prefix++;
                 } else {
                     num_suffix++;
                 }*/

                if (protrude_left != protrude_right)
                {
                    if (protrude_left > 0)
                    {
                        num_prefix++;
                    }
                    else
                    {
                        num_suffix++;
                    }
                }
            }
            if (num_prefix == 0 || num_suffix == 0)
            {
                num_lacking++;
            }
        }
        std::cout << "num lacking: " << num_lacking << std::endl;
    }

    void Aligner::find_RAM_overlaps(bool filter)
    {
        tps.clear();
        fps.clear();
        fns.clear();
        std::uint32_t num_lacking = 0;
        std::uint32_t total_num_prefix = 0;
        std::uint32_t total_num_suffix = 0;

        std::uint32_t total_num_false_prefix = 0;
        
        std::uint32_t total_num_false_suffix = 0;

        if (filter)
        {
            std::cout << "filter" << std::endl;
        }
        else
        {
            std::cout << "no filter" << std::endl;
        }
        for (auto i = 0; i < TARGET_NUM; i++)
        {
            std::uint32_t num_prefix = 0;
            std::uint32_t num_suffix = 0;
            std::uint32_t num_false_prefix = 0;
            std::uint32_t num_false_suffix = 0;
            auto &target = sequences[i];
            std::uint32_t target_id = target->id;
            auto target_info = EXTRACT(target->name);
            target_info.idx_in_sequences = id_to_pos_index[target_id];
            std::vector<biosoup::Overlap> overlaps = minimizer_engine.Map(target, true, false, MINHASH_BOOL);
            if (filter)
            {
                filter_overlaps(overlaps);
            }
            for (auto &o : overlaps)
            {
                auto &s = sequences[id_to_pos_index[o.rhs_id]];
                auto query_info = EXTRACT(s->name);
                query_info.idx_in_sequences = id_to_pos_index[o.rhs_id];
                overlap_info o_info{query_info, target_info};

                std::uint32_t t_begin = o.lhs_begin;
                std::uint32_t t_end = o.lhs_end;
                std::uint32_t q_begin = o.rhs_begin;
                std::uint32_t q_end = o.rhs_end;

                if (!o.strand)
                {
                    o_info.is_reverse = true;
                    q_end = s->inflated_len - o.rhs_begin;
                    q_begin = s->inflated_len - o.rhs_end;
                }
                int protrude_left = q_begin - t_begin;
                int protrude_right = (s->inflated_len - q_end) - (target->inflated_len - t_end);

                std::uint32_t q_clip_left;
                std::uint32_t q_clip_right;
                std::uint32_t t_clip_left;
                std::uint32_t t_clip_right;
                if (protrude_left > 0)
                {
                    q_clip_left = protrude_left;
                    t_clip_left = 0;
                }
                else
                {
                    q_clip_left = 0;
                    t_clip_left = -protrude_left;
                }
                if (protrude_right > 0)
                {
                    q_clip_right = protrude_right;
                    t_clip_right = 0;
                }
                else
                {
                    q_clip_right = 0;
                    t_clip_right = -protrude_right;
                }

                o_info.q_clip_left = q_clip_left;
                o_info.q_clip_right = q_clip_right;
                o_info.t_clip_left = t_clip_left;
                o_info.t_clip_right = t_clip_right;

                o_info.q_o_begin = q_begin;
                o_info.q_o_end = q_end;
                o_info.t_o_begin = t_begin;
                o_info.t_o_end = t_end;
                o_info.score = o.score;
                /*
                if (protrude_left == 0 && protrude_right == 0)
                {
                    num_contain++;
                }
                else if (protrude_left > 0 && protrude_right > 0)
                {
                    num_contained++;
                }
                else if (protrude_left > 0)
                {
                    num_prefix++;
                }
                else
                {
                    num_suffix++;
                }*/

                
                if (overlap_len(query_info, target_info) > OVLP_THRES)
                {
                    tps.insert(o_info);
                    auto o_type = overlap_type(target_info, query_info, target_info.forward);
                    if (o_type == TYPE_PREFIX)
                    {
                   
                        num_prefix++;
                    }
                    else if (o_type == TYPE_SUFFIX)
                    {
                       
                        num_suffix++;
                    }

                }
                else
                {
                    fps.insert(o_info);
                    if (protrude_left > 0) {
                         if (protrude_right < 0) {
                             total_num_false_prefix++;
                             num_false_prefix++;
                         }
                    } else if (protrude_right > 0) {
                         total_num_false_suffix++;
                         num_false_suffix++;
                    }

                }
            }
            if (num_prefix == 0 || num_suffix == 0)
            {
                num_lacking++;
                found_lack.emplace(i, std::make_pair(num_prefix, num_suffix));
                found_lack_false.emplace(i, std::make_pair(num_false_prefix, num_false_suffix));
                 
            } else {
                
                found_has.emplace(i, std::make_pair(num_prefix, num_suffix)); 
                
                found_has_false.emplace(i, std::make_pair(num_false_prefix, num_false_suffix));
            }
            total_num_prefix += num_prefix;
            total_num_suffix += num_suffix;
        }

        for (auto &o_info : all_true_overlaps)
        {
            if (tps.find(o_info) == tps.end())
            {
                fns.insert(o_info);
            }
        }
        std::cout << "num prefix found: " << total_num_prefix << std::endl;
        std::cout << "num suffix found: " << total_num_suffix << std::endl;
        std::cout << "avg num recovered prefix: " << (double) total_num_prefix / TARGET_NUM << std::endl;
        std::cout << "avg num recovered suffix: " << (double) total_num_suffix / TARGET_NUM << std::endl;
        std::cout << "avg num false prefix: " << (double) total_num_false_prefix / TARGET_NUM << std::endl;
        std::cout << "avg num false suffix: " << (double) total_num_false_suffix / TARGET_NUM << std::endl;
        std::cout << "num lacking: " << num_lacking << std::endl;
        std::cout << "precision: " << ((double)tps.size() / (tps.size() + fps.size())) << std::endl;
        std::cout << "recall: " << ((double)tps.size() / all_true_overlaps.size()) << std::endl;
    }

    void Aligner::false_negatives()
    {
        std::uint32_t num_100 = 0;
        std::uint32_t num_200 = 0;        
        std::uint32_t num_300 = 0;
        std::uint32_t num_400 = 0;
        std::uint32_t num_500 = 0;
        std::uint32_t num_1000 = 0;
        std::uint32_t num_2000 = 0;
        std::uint32_t num_4000 = 0;
        std::uint32_t num_10000 = 0;
        std::uint32_t more = 0;
        std::cout << "----false negatives----" << std::endl;
        for (auto &o_info : fns)
        {
            if (!o_info.query_info.aligned || !o_info.target_info.aligned)
            {
                continue;
            }

            auto &query_sequence = sequences[o_info.query_info.idx_in_sequences];
            auto &target_sequence = sequences[o_info.target_info.idx_in_sequences];

            if (!o_info.query_info.forward)
                query_sequence->ReverseAndComplement();
            if (!o_info.target_info.forward)
                target_sequence->ReverseAndComplement();
            std::uint32_t q_start = o_info.query_info.start;
            std::uint32_t q_end = o_info.query_info.end;
            std::uint32_t t_start = o_info.target_info.start;
            std::uint32_t t_end = o_info.target_info.end;

            int q_clip_left = t_start - q_start;
            int q_clip_right = q_end - t_end;
            int t_clip_left;
            int t_clip_right;
            if (q_clip_left >= 0)
            {
                t_clip_left = 0;
            }
            else
            {
                t_clip_left = -q_clip_left;
                q_clip_left = 0;
            }
            if (q_clip_right >= 0)
            {
                t_clip_right = 0;
            }
            else
            {
                t_clip_right = -q_clip_right;
                q_clip_right = 0;
            }

            q_clip_left += o_info.query_info.left_extend;
            q_clip_right += o_info.query_info.right_extend;
            t_clip_left += o_info.target_info.left_extend;
            t_clip_right += o_info.target_info.right_extend;

            std::uint32_t ovlp_len = overlap_len(o_info.query_info, o_info.target_info);

            int query_string_len = query_sequence->inflated_len - q_clip_left - q_clip_right;
            int target_string_len = target_sequence->inflated_len - t_clip_left - t_clip_right;
            if (query_string_len == 0 || target_string_len == 0 || ovlp_len < 0.1 * std::max(query_sequence->inflated_len, target_sequence->inflated_len))
            {
                continue;
            }

            std::cout << "ovlp " << ovlp_len << std::endl;
            if (ovlp_len <= 100) {
                num_100++;         
                std::cout << "cat 100" << std::endl;

            } else if (ovlp_len <= 200) {
                num_200++;

                std::cout << "cat 200" << std::endl;
            } else if (ovlp_len <= 300) {
                num_300++;         
                std::cout << "cat 300" << std::endl;
            } else if (ovlp_len <= 400) {
                num_400++;         
                std::cout << "cat 400" << std::endl;
            } else if (ovlp_len <= 500) {
                num_500++;
                std::cout << "cat 500" << std::endl;
            } else if (ovlp_len <= 1000) {
                num_1000++;
                std::cout << "cat 1000" << std::endl;
            } else if (ovlp_len <= 2000) {
                num_2000++;

                std::cout << "cat 2000" << std::endl;
            } else if (ovlp_len <= 4000) {
                num_4000++;

                std::cout << "cat 4000" << std::endl;
            } else if (ovlp_len <= 10000) {
                num_10000++;

                std::cout << "cat 10000" << std::endl;
            } else {
                more++;

                std::cout << "cat more" << std::endl;
            }






            std::cout << "query name: " << o_info.query_info.name << std::endl;
            std::cout << "target name: " << o_info.target_info.name << std::endl;

            std::cout << "pre-clip q-len: " << query_sequence->inflated_len << " pre-clip t-len: " << target_sequence->inflated_len << std::endl;
            std::cout << "q_start: " << o_info.query_info.start << " query_end: " << o_info.query_info.end << std::endl;
            std::cout << "query_clip_left: " << q_clip_left << " query_clip_right: " << q_clip_right << std::endl;
            std::cout << "t_start: " << o_info.target_info.start << " target_end: " << o_info.target_info.end << std::endl;
            std::cout << "t_clip_left: " << t_clip_left << " target_clip_right: " << t_clip_right << std::endl;
            std::cout << "overlap len: " << ovlp_len << std::endl;

            std::string query_string = query_sequence->InflateData(q_clip_left, query_string_len);
            std::string target_string = target_sequence->InflateData(t_clip_left, target_string_len);

            if (!o_info.query_info.forward)
                query_sequence->ReverseAndComplement();
            if (!o_info.target_info.forward)
                target_sequence->ReverseAndComplement();

            std::vector<std::string> overlapping_reads;

            EdlibAlignResult result = edlibAlign(query_string.c_str(), query_string.size(),
                                                 target_string.c_str(), target_string.size(),
                                                 edlibNewAlignConfig(-1, EDLIB_MODE_HW, EDLIB_TASK_PATH, NULL, 0));

            overlapping_reads.push_back(std::move(query_string));
            std::vector<std::pair<std::uint32_t, std::uint32_t>> clips;
            clips.emplace_back(0, 0);

            std::vector<EdlibAlignResult> rs = {result};

            auto alignment = multi_align(overlapping_reads, target_string, clips, rs, false);

            alignment.print();
        }
        std::cout << "(0, 100]: " << num_100 << std::endl;
        std::cout << "(100, 200]: " << num_200 << std::endl;
        std::cout << "(200, 300]: " << num_300 << std::endl;        
        std::cout << "(300, 400]: " << num_400 << std::endl;
        std::cout << "(400, 500]: " << num_500 << std::endl;
        std::cout << "(500, 1000]: " << num_1000 << std::endl;
        std::cout << "(1000, 2000]: " << num_2000 << std::endl;
        std::cout << "(2000, 4000]: " << num_4000 << std::endl;
        std::cout << "(4000, 10000]: " << num_10000 << std::endl;
        std::cout << "(10000, ?]: " << more << std::endl;
    }
    void Aligner::false_positives_align_part()
    {

        std::cout << "----false positives----" << std::endl;

        double total_norm_dist = 0;
        double total_norm_score_whole = 0;
        double total_norm_score_part = 0;
        std::uint32_t count_0001 = 0;
        std::uint32_t count_0102 = 0;
        std::uint32_t count_0203 = 0;
        std::uint32_t count_0304 = 0;
        std::uint32_t count_greater = 0;

        for (auto &o_info : fps)
        {

            auto &query_sequence = sequences[o_info.query_info.idx_in_sequences];
            auto &target_sequence = sequences[o_info.target_info.idx_in_sequences];

            if (o_info.is_reverse)
                query_sequence->ReverseAndComplement();

            std::uint32_t t_o_begin = o_info.t_o_begin;
            std::uint32_t t_o_end = o_info.t_o_end;
            std::uint32_t q_o_begin = o_info.q_o_begin;
            std::uint32_t q_o_end = o_info.q_o_end;
            std::uint32_t query_string_len = q_o_end - q_o_begin;
            std::uint32_t target_string_len = t_o_end - t_o_begin;
            std::string target_string = target_sequence->InflateData(t_o_begin, target_string_len);
            std::string query_string = query_sequence->InflateData(q_o_begin, query_string_len);

            EdlibAlignResult result = edlibAlign(query_string.c_str(), query_string.size(),
                                                 target_string.c_str(), target_string.size(),
                                                 edlibNewAlignConfig(-1, EDLIB_MODE_HW, EDLIB_TASK_PATH, NULL, 0));

            double norm_dist = (double)result.editDistance / query_string_len;
            total_norm_dist += norm_dist;
            total_norm_score_whole = (double)o_info.score / std::max(query_string_len, target_string_len);
            total_norm_score_part = (double)o_info.score / std::max(o_info.t_o_end - o_info.t_o_begin, o_info.q_o_end - o_info.q_o_begin);

            std::vector<std::string> overlapping_reads;
            overlapping_reads.push_back(std::move(query_string));
            std::vector<std::pair<std::uint32_t, std::uint32_t>> clips;
            clips.emplace_back(0, 0);
            std::vector<EdlibAlignResult> rs = {result};
            auto alignment = multi_align(overlapping_reads, target_string, clips, rs, false);
            std::cout << "query name: " << o_info.query_info.name << std::endl;
            std::cout << "target name: " << o_info.target_info.name << std::endl;
            std::cout << "query_start: " << o_info.query_info.start << " query_end: " << o_info.query_info.end << std::endl;
            std::cout << "q_o_begin: " << o_info.q_o_begin << " q_o_end: " << o_info.q_o_end << std::endl;
            std::cout << "t_start: " << o_info.target_info.start << " target_end: " << o_info.target_info.end << std::endl;
            std::cout << "t_o_begin: " << o_info.t_o_begin << " t_o_end: " << o_info.t_o_end << std::endl;
            std::uint32_t ovlp_len = overlap_len(o_info.query_info, o_info.target_info);
            std::cout << "overlap len: " << ovlp_len << std::endl;
            std::cout << "norm dist: " << norm_dist << std::endl;
            std::cout << "norm score whole: " << (double)o_info.score / std::max(query_string_len, target_string_len) << std::endl;
            std::cout << "norm score part: " << (double)o_info.score / std::max(o_info.t_o_end - o_info.t_o_begin, o_info.q_o_end - o_info.q_o_begin) << std::endl;
            // std::cout << "q_o_left: " << o_info.q_o_left << " query_o_right: " << o_info.q_o_right << std::endl;
            // std::cout << "t_o_left: " << o_info.t_o_left << " target_o_right: " << o_info.t_o_right << std::endl;
            if (norm_dist < 0.1)
            {
                count_0001++;
            }
            else if (norm_dist < 0.2)
            {
                count_0102++;
            }
            else if (norm_dist < 0.3)
            {
                count_0203++;
            }
            else if (norm_dist < 0.4)
            {
                count_0304++;
            }
            else
            {
                count_greater++;
            }

            alignment.print();

            if (o_info.is_reverse)
                query_sequence->ReverseAndComplement();
        }

        std::cout << "avg norm dist: " << total_norm_dist / tps.size() << std::endl;
        std::cout << "avg norm score whole: " << total_norm_score_whole / tps.size() << std::endl;
        std::cout << "avg norm score part: " << total_norm_score_part / tps.size() << std::endl;
        std::cout << "count [0, 0.1): " << count_0001 << std::endl;
        std::cout << "count [0.1, 0.2): " << count_0102 << std::endl;
        std::cout << "count [0.2, 0.3): " << count_0203 << std::endl;
        std::cout << "count [0.3, 0.4): " << count_0304 << std::endl;
        std::cout << "count >= 0.4: " << count_greater << std::endl;

        std::cout << "----false positives----" << std::endl;
    }

    void Aligner::true_positives_align_part()
    {

        std::cout << "----true positives----" << std::endl;

        double total_norm_dist = 0;
        double total_norm_score_whole = 0;
        double total_norm_score_part = 0;
        std::uint32_t count_0001 = 0;
        std::uint32_t count_0102 = 0;
        std::uint32_t count_0203 = 0;
        std::uint32_t count_0304 = 0;
        std::uint32_t count_greater = 0;

        // true positives
        for (auto &o_info : tps)
        {

            auto &query_sequence = sequences[o_info.query_info.idx_in_sequences];
            auto &target_sequence = sequences[o_info.target_info.idx_in_sequences];

            if (o_info.is_reverse)
                query_sequence->ReverseAndComplement();

            std::uint32_t t_o_begin = o_info.t_o_begin;
            std::uint32_t t_o_end = o_info.t_o_end;
            std::uint32_t q_o_begin = o_info.q_o_begin;
            std::uint32_t q_o_end = o_info.q_o_end;

            std::uint32_t query_string_len = q_o_end - q_o_begin;
            std::uint32_t target_string_len = t_o_end - t_o_begin;
            std::string target_string = target_sequence->InflateData(t_o_begin, target_string_len);
            std::string query_string = query_sequence->InflateData(q_o_begin, query_string_len);

            EdlibAlignResult result = edlibAlign(query_string.c_str(), query_string.size(),
                                                 target_string.c_str(), target_string.size(),
                                                 edlibNewAlignConfig(-1, EDLIB_MODE_HW, EDLIB_TASK_PATH, NULL, 0));

            double norm_dist = (double)result.editDistance / query_string_len;
            total_norm_dist += norm_dist;
            total_norm_score_whole = (double)o_info.score / std::max(query_string_len, target_string_len);
            total_norm_score_part = (double)o_info.score / std::max(o_info.t_o_end - o_info.t_o_begin, o_info.q_o_end - o_info.q_o_begin);

            std::vector<std::string> overlapping_reads;
            overlapping_reads.push_back(std::move(query_string));
            std::vector<std::pair<std::uint32_t, std::uint32_t>> clips;
            clips.emplace_back(0, 0);
            std::vector<EdlibAlignResult> rs = {result};
            auto alignment = multi_align(overlapping_reads, target_string, clips, rs, false);
            std::cout << "query name: " << o_info.query_info.name << std::endl;
            std::cout << "target name: " << o_info.target_info.name << std::endl;
            std::cout << "query_start: " << o_info.query_info.start << " query_end: " << o_info.query_info.end << std::endl;
            std::cout << "q_o_begin: " << o_info.q_o_begin << " q_o_end: " << o_info.q_o_end << std::endl;
            std::cout << "t_start: " << o_info.target_info.start << " target_end: " << o_info.target_info.end << std::endl;
            std::cout << "t_o_begin: " << o_info.t_o_begin << " t_o_end: " << o_info.t_o_end << std::endl;
            std::uint32_t ovlp_len = overlap_len(o_info.query_info, o_info.target_info);
            std::cout << "overlap len: " << ovlp_len << std::endl;
            std::cout << "norm dist: " << norm_dist << std::endl;
            std::cout << "norm score whole: " << (double)o_info.score / std::max(query_string_len, target_string_len) << std::endl;
            std::cout << "norm score part: " << (double)o_info.score / std::max(o_info.t_o_end - o_info.t_o_begin, o_info.q_o_end - o_info.q_o_begin) << std::endl;
            // std::cout << "q_o_left: " << o_info.q_o_left << " query_o_right: " << o_info.q_o_right << std::endl;
            // std::cout << "t_o_left: " << o_info.t_o_left << " target_o_right: " << o_info.t_o_right << std::endl;

            if (norm_dist < 0.1)
            {
                count_0001++;
            }
            else if (norm_dist < 0.2)
            {
                count_0102++;
            }
            else if (norm_dist < 0.3)
            {
                count_0203++;
            }
            else if (norm_dist < 0.4)
            {
                count_0304++;
            }
            else
            {
                count_greater++;
            }

            alignment.print();
            if (o_info.is_reverse)
                query_sequence->ReverseAndComplement();
        }

        std::cout << "avg norm dist: " << total_norm_dist / tps.size() << std::endl;
        std::cout << "avg norm score whole: " << total_norm_score_whole / tps.size() << std::endl;
        std::cout << "avg norm score part: " << total_norm_score_part / tps.size() << std::endl;
        std::cout << "count [0, 0.1): " << count_0001 << std::endl;
        std::cout << "count [0.1, 0.2): " << count_0102 << std::endl;
        std::cout << "count [0.2, 0.3): " << count_0203 << std::endl;
        std::cout << "count [0.3, 0.4): " << count_0304 << std::endl;
        std::cout << "count >= 0.4: " << count_greater << std::endl;

        std::cout << "----true positives----" << std::endl;
    }

    void Aligner::true_positives()
    {

        std::cout << "----true positives----" << std::endl;

        double total_norm_dist = 0;
        double total_norm_score_whole = 0;
        double total_norm_score_part = 0;
        double total_hang_ratio = 0;
        std::uint32_t count_0001 = 0;
        std::uint32_t count_0102 = 0;
        std::uint32_t count_0203 = 0;
        std::uint32_t count_0304 = 0;
        std::uint32_t count_greater = 0;

        std::uint32_t count0to100 = 0;
        std::uint32_t count100to500 = 0;
        std::uint32_t count500to2000 = 0;
        std::uint32_t count_greater2000 = 0;

        std::uint32_t count_h_0001 = 0;
        std::uint32_t count_h_0102 = 0;
        std::uint32_t count_h_0203 = 0;
        std::uint32_t count_h_0304 = 0;
        std::uint32_t count_h_greater = 0;

        // true positives
        for (auto &o_info : tps)
        {

            auto &query_sequence = sequences[o_info.query_info.idx_in_sequences];
            auto &target_sequence = sequences[o_info.target_info.idx_in_sequences];

            if (o_info.is_reverse)
                query_sequence->ReverseAndComplement();

            std::uint32_t left_overhang = o_info.q_o_begin - o_info.q_clip_left;
            std::uint32_t right_overhang = query_sequence->inflated_len - o_info.q_clip_right - o_info.q_o_end;
            // double bad_left_ratio = (double) left_overhang / (query_sequence->inflated_len - o_info.q_clip_left - o_info.q_clip_right);
            // double bad_right_ratio = (double) right_overhang / (query_sequence->inflated_len - o_info.q_clip_left - o_info.q_clip_right);
            std::uint32_t larger = std::max(left_overhang, right_overhang);
            if (larger < 100)
            {
                count0to100++;
            }
            else if (larger < 500)
            {
                count100to500++;
            }
            else if (larger < 2000)
            {
                count500to2000++;
            }
            else
            {
                count_greater2000++;
            }

            std::uint32_t query_string_len = query_sequence->inflated_len - o_info.q_clip_left - o_info.q_clip_right;
            double hang_ratio = (double)(left_overhang + right_overhang) / query_string_len;
            total_hang_ratio += hang_ratio;
            if (hang_ratio < 0.1)
            {
                count_h_0001++;
            }
            else if (hang_ratio < 0.2)
            {
                count_h_0102++;
            }
            else if (hang_ratio < 0.3)
            {
                count_h_0203++;
            }
            else if (hang_ratio < 0.4)
            {
                count_h_0304++;
            }
            else
            {
                count_h_greater++;
            }

            std::uint32_t target_string_len = target_sequence->inflated_len - o_info.t_clip_left - o_info.t_clip_right;
            std::string target_string = target_sequence->InflateData(o_info.t_clip_left, target_string_len);
            std::string query_string = query_sequence->InflateData(o_info.q_clip_left, query_string_len);

            EdlibAlignResult result = edlibAlign(query_string.c_str(), query_string.size(),
                                                 target_string.c_str(), target_string.size(),
                                                 edlibNewAlignConfig(-1, EDLIB_MODE_HW, EDLIB_TASK_PATH, NULL, 0));

            double norm_dist = (double)result.editDistance / query_string_len;
            total_norm_dist += norm_dist;
            total_norm_score_whole = (double)o_info.score / std::max(query_string_len, target_string_len);
            total_norm_score_part = (double)o_info.score / std::max(o_info.t_o_end - o_info.t_o_begin, o_info.q_o_end - o_info.q_o_begin);

            std::vector<std::string> overlapping_reads;
            overlapping_reads.push_back(std::move(query_string));
            std::vector<std::pair<std::uint32_t, std::uint32_t>> clips;
            clips.emplace_back(0, 0);
            std::vector<EdlibAlignResult> rs = {result};
            auto alignment = multi_align(overlapping_reads, target_string, clips, rs, false);
            std::cout << "query name: " << o_info.query_info.name << std::endl;
            std::cout << "target name: " << o_info.target_info.name << std::endl;
            std::cout << "query_start: " << o_info.query_info.start << " query_end: " << o_info.query_info.end << std::endl;
            std::cout << "q_clip_left: " << o_info.q_clip_left << " query_clip_right: " << o_info.q_clip_right << std::endl;
            std::cout << "t_start: " << o_info.target_info.start << " target_end: " << o_info.target_info.end << std::endl;
            std::cout << "t_clip_left: " << o_info.t_clip_left << " target_clip_right: " << o_info.t_clip_right << std::endl;
            std::uint32_t ovlp_len = overlap_len(o_info.query_info, o_info.target_info);
            std::cout << "overlap len: " << ovlp_len << std::endl;
            std::cout << "norm dist: " << norm_dist << std::endl;
            std::cout << "norm score whole: " << (double)o_info.score / std::max(query_string_len, target_string_len) << std::endl;
            std::cout << "norm score part: " << (double)o_info.score / std::max(o_info.t_o_end - o_info.t_o_begin, o_info.q_o_end - o_info.q_o_begin) << std::endl;
            std::cout << "left overhang: " << left_overhang << std::endl;
            std::cout << "right overhang: " << right_overhang << std::endl;
            // std::cout << "q_o_left: " << o_info.q_o_left << " query_o_right: " << o_info.q_o_right << std::endl;
            // std::cout << "t_o_left: " << o_info.t_o_left << " target_o_right: " << o_info.t_o_right << std::endl;
            if (norm_dist < 0.1)
            {
                count_0001++;
            }
            else if (norm_dist < 0.2)
            {
                count_0102++;
            }
            else if (norm_dist < 0.3)
            {
                count_0203++;
            }
            else if (norm_dist < 0.4)
            {
                count_0304++;
            }
            else
            {
                count_greater++;
            }

            alignment.print();
            if (o_info.is_reverse)
                query_sequence->ReverseAndComplement();
        }

        std::cout << "avg norm dist: " << total_norm_dist / tps.size() << std::endl;
        std::cout << "avg norm score whole: " << total_norm_score_whole / tps.size() << std::endl;
        std::cout << "avg norm score part: " << total_norm_score_part / tps.size() << std::endl;
        std::cout << "avg hang ratio: " << total_hang_ratio / tps.size() << std::endl;
        std::cout << "count [0, 0.1): " << count_0001 << std::endl;
        std::cout << "count [0.1, 0.2): " << count_0102 << std::endl;
        std::cout << "count [0.2, 0.3): " << count_0203 << std::endl;
        std::cout << "count [0.3, 0.4): " << count_0304 << std::endl;
        std::cout << "count >= 0.4: " << count_greater << std::endl;

        std::cout << "count h[0, 0.1): " << count_h_0001 << std::endl;
        std::cout << "count h[0.1, 0.2): " << count_h_0102 << std::endl;
        std::cout << "count h[0.2, 0.3): " << count_h_0203 << std::endl;
        std::cout << "count h[0.3, 0.4): " << count_h_0304 << std::endl;
        std::cout << "count h>= 0.4: " << count_h_greater << std::endl;

        std::cout << "count [0, 100): " << count0to100 << std::endl;
        std::cout << "count [100, 500): " << count100to500 << std::endl;
        std::cout << "count [500, 2000): " << count500to2000 << std::endl;
        std::cout << "count >= 2000: " << count_greater2000 << std::endl;

        std::cout << "----true positives----" << std::endl;
    }

    void Aligner::filtered_trues()
    {

        std::cout << "----filtered true positives----" << std::endl;

        double total_norm_dist = 0;
        double total_hang_ratio = 0;

        std::uint32_t count0to100 = 0;
        std::uint32_t count100to500 = 0;
        std::uint32_t count500to2000 = 0;
        std::uint32_t count_greater2000 = 0;
        std::uint32_t min_large = 10000;
        double min_norm_dist = 1;


        // true positives
        for (auto &o_info : filtered_t)
        {
            auto &query_sequence = sequences[o_info.query_info.idx_in_sequences];
            auto &target_sequence = sequences[o_info.target_info.idx_in_sequences];

            if (o_info.is_reverse)
                query_sequence->ReverseAndComplement();


            
            std::uint32_t large_extend_q = std::max(o_info.query_info.left_extend, o_info.query_info.right_extend);
            std::uint32_t large_extend_t = std::max(o_info.target_info.left_extend, o_info.target_info.right_extend);
            std::uint32_t large_extend = std::max(large_extend_q, large_extend_t);
            if (large_extend < 100) {
                std::cout << "< 100" << std::endl;       
            } else if (large_extend < 200) {
                std::cout << "< 200" << std::endl;
            } else if (large_extend < 300) {
                std::cout << "< 300" << std::endl;
            } else if (large_extend < 400) {
                std::cout << "< 400" << std::endl;
            } else if (large_extend < 500) {
                std::cout << "< 500" << std::endl;
            } 
            if (large_extend < min_large) min_large = large_extend;
            std::uint32_t left_overhang = o_info.q_o_begin - o_info.q_clip_left;
            std::uint32_t right_overhang = query_sequence->inflated_len - o_info.q_clip_right - o_info.q_o_end;
            // double bad_left_ratio = (double) left_overhang / (query_sequence->inflated_len - o_info.q_clip_left - o_info.q_clip_right);
            // double bad_right_ratio = (double) right_overhang / (query_sequence->inflated_len - o_info.q_clip_left - o_info.q_clip_right);
            std::uint32_t larger = std::max(left_overhang, right_overhang);
            if (larger < 100)
            {
                count0to100++;
            }
            else if (larger < 500)
            {
                count100to500++;
            }
            else if (larger < 2000)
            {
                count500to2000++;
            }
            else
            {
                count_greater2000++;
            }

            std::uint32_t query_string_len = query_sequence->inflated_len - o_info.q_clip_left - o_info.q_clip_right;
            double hang_ratio = (double)(left_overhang + right_overhang) / query_string_len;
            total_hang_ratio += hang_ratio;
            std::uint32_t target_string_len = target_sequence->inflated_len - o_info.t_clip_left - o_info.t_clip_right;
            std::string target_string = target_sequence->InflateData(o_info.t_clip_left, target_string_len);
            std::string query_string = query_sequence->InflateData(o_info.q_clip_left, query_string_len); 
            EdlibAlignResult result = edlibAlign(query_string.c_str(), query_string.size(),
                                                 target_string.c_str(), target_string.size(),
                                                 edlibNewAlignConfig(-1, EDLIB_MODE_HW, EDLIB_TASK_PATH, NULL, 0));
            double norm_dist = (double)result.editDistance / query_string_len;
            if (norm_dist < min_norm_dist) min_norm_dist = norm_dist;
            total_norm_dist += norm_dist;

            std::vector<std::string> overlapping_reads;
            overlapping_reads.push_back(std::move(query_string));
            std::vector<std::pair<std::uint32_t, std::uint32_t>> clips;
            clips.emplace_back(0, 0);
            std::vector<EdlibAlignResult> rs = {result};
            auto alignment = multi_align(overlapping_reads, target_string, clips, rs, false);
            std::cout << "query name: " << o_info.query_info.name << std::endl;
            std::cout << "target name: " << o_info.target_info.name << std::endl;
            std::cout << "qob: " << o_info.q_o_begin << " qoe: " << o_info.q_o_end << std::endl;
            std::cout << "tob: " << o_info.t_o_begin << " toe: " << o_info.t_o_end << std::endl;
            std::cout << "query_start: " << o_info.query_info.start << " query_end: " << o_info.query_info.end << std::endl;
            std::cout << "q_clip_left: " << o_info.q_clip_left << " query_clip_right: " << o_info.q_clip_right << std::endl;
            std::cout << "t_start: " << o_info.target_info.start << " target_end: " << o_info.target_info.end << std::endl;
            std::cout << "t_clip_left: " << o_info.t_clip_left << " target_clip_right: " << o_info.t_clip_right << std::endl;
            std::uint32_t ovlp_len = overlap_len(o_info.query_info, o_info.target_info);
            std::cout << "overlap len: " << ovlp_len << std::endl;
            std::cout << "norm dist: " << norm_dist << std::endl;
            std::cout << "norm score whole: " << (double)o_info.score / std::max(query_string_len, target_string_len) << std::endl;
            std::cout << "norm score part: " << (double)o_info.score / std::max(o_info.t_o_end - o_info.t_o_begin, o_info.q_o_end - o_info.q_o_begin) << std::endl;
            std::cout << "left overhang: " << left_overhang << std::endl;
            std::cout << "right overhang: " << right_overhang << std::endl;
            // std::cout << "q_o_left: " << o_info.q_o_left << " query_o_right: " << o_info.q_o_right << std::endl;
            // std::cout << "t_o_left: " << o_info.t_o_left << " target_o_right: " << o_info.t_o_right << std::endl;
            
            alignment.print();
            
            if (o_info.is_reverse)
                query_sequence->ReverseAndComplement();
        }

        std::cout << "avg norm dist: " << total_norm_dist / filtered_t.size() << std::endl;
        std::cout << "avg hang ratio: " << total_hang_ratio / filtered_t.size() << std::endl;

        std::cout << "count [0, 100): " << count0to100 << std::endl;
        std::cout << "count [100, 500): " << count100to500 << std::endl;
        std::cout << "count [500, 2000): " << count500to2000 << std::endl;
        std::cout << "count >= 2000: " << count_greater2000 << std::endl;
        std::cout << "min large : " << min_large << std::endl;
        std::cout << "min norm dist: " << min_norm_dist << std::endl;
        std::cout << "----filtered true positives----" << std::endl;
    }

    void Aligner::false_positives()
    {

        std::cout << "----false positives----" << std::endl;

        double total_norm_dist = 0;
        double total_norm_score_whole = 0;
        double total_norm_score_part = 0;
        double total_hang_ratio = 0;
        std::uint32_t count_h_0001 = 0;
        std::uint32_t count_h_0102 = 0;
        std::uint32_t count_h_0203 = 0;
        std::uint32_t count_h_0304 = 0;
        std::uint32_t count_h_greater = 0;

        std::uint32_t count_0001 = 0;
        std::uint32_t count_0102 = 0;
        std::uint32_t count_0203 = 0;
        std::uint32_t count_0304 = 0;
        std::uint32_t count_greater = 0;

        std::uint32_t count0to100 = 0;
        std::uint32_t count100to500 = 0;
        std::uint32_t count500to2000 = 0;
        std::uint32_t count_greater2000 = 0;

        // true positives
        for (auto &o_info : fps)
        {

            auto &query_sequence = sequences[o_info.query_info.idx_in_sequences];
            auto &target_sequence = sequences[o_info.target_info.idx_in_sequences];

            if (o_info.is_reverse)
                query_sequence->ReverseAndComplement();

            std::uint32_t left_overhang = o_info.q_o_begin - o_info.q_clip_left;
            std::uint32_t right_overhang = query_sequence->inflated_len - o_info.q_clip_right - o_info.q_o_end;
            std::uint32_t larger = std::max(left_overhang, right_overhang);
            if (larger < 100)
            {
                count0to100++;
            }
            else if (larger < 500)
            {
                count100to500++;
            }
            else if (larger < 2000)
            {
                count500to2000++;
            }
            else
            {
                count_greater2000++;
            }

            std::uint32_t query_string_len = query_sequence->inflated_len - o_info.q_clip_left - o_info.q_clip_right;
            double hang_ratio = (double)(left_overhang + right_overhang) / query_string_len;
            total_hang_ratio += hang_ratio;
            if (hang_ratio < 0.1)
            {
                count_h_0001++;
            }
            else if (hang_ratio < 0.2)
            {
                count_h_0102++;
            }
            else if (hang_ratio < 0.3)
            {
                count_h_0203++;
            }
            else if (hang_ratio < 0.4)
            {
                count_h_0304++;
            }
            else
            {
                count_h_greater++;
            }

            std::uint32_t target_string_len = target_sequence->inflated_len - o_info.t_clip_left - o_info.t_clip_right;
            std::string target_string = target_sequence->InflateData(o_info.t_clip_left, target_string_len);
            std::string query_string = query_sequence->InflateData(o_info.q_clip_left, query_string_len);
            EdlibAlignResult result = edlibAlign(query_string.c_str(), query_string.size(),
                                                 target_string.c_str(), target_string.size(),
                                                 edlibNewAlignConfig(-1, EDLIB_MODE_HW, EDLIB_TASK_PATH, NULL, 0));

            double norm_dist = (double)result.editDistance / query_string_len;
            total_norm_dist += norm_dist;
            total_norm_score_whole = (double)o_info.score / std::max(query_string_len, target_string_len);
            total_norm_score_part = (double)o_info.score / std::max(o_info.t_o_end - o_info.t_o_begin, o_info.q_o_end - o_info.q_o_begin);

            std::vector<std::string> overlapping_reads;
            overlapping_reads.push_back(std::move(query_string));
            std::vector<std::pair<std::uint32_t, std::uint32_t>> clips;
            clips.emplace_back(0, 0);
            std::vector<EdlibAlignResult> rs = {result};
            auto alignment = multi_align(overlapping_reads, target_string, clips, rs, false);
            std::cout << "query name: " << o_info.query_info.name << std::endl;
            std::cout << "target name: " << o_info.target_info.name << std::endl;
            std::cout << "query_start: " << o_info.query_info.start << " query_end: " << o_info.query_info.end << std::endl;
            std::cout << "q_clip_left: " << o_info.q_clip_left << " query_clip_right: " << o_info.q_clip_right << std::endl;
            std::cout << "t_start: " << o_info.target_info.start << " target_end: " << o_info.target_info.end << std::endl;
            std::cout << "t_clip_left: " << o_info.t_clip_left << " target_clip_right: " << o_info.t_clip_right << std::endl;
            std::uint32_t ovlp_len = overlap_len(o_info.query_info, o_info.target_info);
            std::cout << "overlap len: " << ovlp_len << std::endl;
            std::cout << "norm dist: " << norm_dist << std::endl;
            std::cout << "norm score whole: " << (double)o_info.score / std::max(query_string_len, target_string_len) << std::endl;
            std::cout << "norm score part: " << (double)o_info.score / std::max(o_info.t_o_end - o_info.t_o_begin, o_info.q_o_end - o_info.q_o_begin) << std::endl;
            std::cout << "left overhang: " << left_overhang << std::endl;
            std::cout << "right overhang: " << right_overhang << std::endl;

            if (norm_dist < 0.1)
            {
                count_0001++;
            }
            else if (norm_dist < 0.2)
            {
                count_0102++;
            }
            else if (norm_dist < 0.3)
            {
                count_0203++;
            }
            else if (norm_dist < 0.4)
            {
                count_0304++;
            }
            else
            {
                count_greater++;
            }

            // std::cout << "q_o_left: " << o_info.q_o_left << " query_o_right: " << o_info.q_o_right << std::endl;
            // std::cout << "t_o_left: " << o_info.t_o_left << " target_o_right: " << o_info.t_o_right << std::endl;
            alignment.print();
            if (o_info.is_reverse)
                query_sequence->ReverseAndComplement();
        }

        std::cout << "avg norm dist: " << total_norm_dist / fps.size() << std::endl;
        std::cout << "avg norm score whole: " << total_norm_score_whole / fps.size() << std::endl;
        std::cout << "avg norm score part: " << total_norm_score_part / fps.size() << std::endl;
        std::cout << "avg hang ratio: " << total_hang_ratio / fps.size() << std::endl;
        std::cout << "count [0, 0.1): " << count_0001 << std::endl;
        std::cout << "count [0.1, 0.2): " << count_0102 << std::endl;
        std::cout << "count [0.2, 0.3): " << count_0203 << std::endl;
        std::cout << "count [0.3, 0.4): " << count_0304 << std::endl;
        std::cout << "count >= 0.4: " << count_greater << std::endl;
        std::cout << "count h[0, 0.1): " << count_h_0001 << std::endl;
        std::cout << "count h[0.1, 0.2): " << count_h_0102 << std::endl;
        std::cout << "count h[0.2, 0.3): " << count_h_0203 << std::endl;
        std::cout << "count h[0.3, 0.4): " << count_h_0304 << std::endl;
        std::cout << "count h>= 0.4: " << count_h_greater << std::endl;

        std::cout << "count [0, 100): " << count0to100 << std::endl;
        std::cout << "count [100, 500): " << count100to500 << std::endl;
        std::cout << "count [500, 2000): " << count500to2000 << std::endl;
        std::cout << "count >= 2000: " << count_greater2000 << std::endl;
        std::cout << "----false positives----" << std::endl;
    }

    Aligner::Aligner(const char **sequences_paths, std::shared_ptr<thread_pool::ThreadPool> &pool, std::uint8_t kmer_len,
                     std::uint8_t window_len, double freq)
        : pool(pool), minimizer_engine(pool, kmer_len, window_len, 500, 4, 100, 10000)
    {

        // srand (time(NULL));
        srand(RANDOM_SEED);
        auto is_suffix = [](const std::string &s, const std::string &suff)
        {
            return s.size() < suff.size() ? false : s.compare(s.size() - suff.size(), suff.size(), suff) == 0;
        };
        for (std::uint8_t i = 0; i < 2; i++)
        {
            const char *path = sequences_paths[i];

            if (path == nullptr)
                break;
            std::string seq_path{path};
            if (is_suffix(seq_path, ".fasta"))
            {
                try
                {
                    auto p = bioparser::Parser<biosoup::NucleicAcid>::Create<bioparser::FastaParser>(seq_path);
                    auto s = p->Parse(-1);
                    sequences.insert(sequences.end(), std::make_move_iterator(s.begin()), std::make_move_iterator(s.end()));
                }
                catch (const std::invalid_argument &exception)
                {
                    std::cerr << exception.what() << std::endl;
                }
            }
            else if (is_suffix(seq_path, ".fastq"))
            {
                try
                {
                    auto p = bioparser::Parser<biosoup::NucleicAcid>::Create<bioparser::FastqParser>(seq_path);
                    auto s = p->Parse(-1);
                    sequences.insert(sequences.end(), std::make_move_iterator(s.begin()), std::make_move_iterator(s.end()));
                }
                catch (const std::invalid_argument &exception)
                {
                    std::cerr << exception.what() << std::endl;
                }
            }
            else
            {
                throw std::invalid_argument("[align_reads::Aligner::Aligner] Error: Invalid sequences file format.");
            }
        }
        auto long_seq_first = [](const std::unique_ptr<biosoup::NucleicAcid> &s1, const std::unique_ptr<biosoup::NucleicAcid> &s2)
        {
            return s1->inflated_len > s2->inflated_len;
        };


        std::random_shuffle(sequences.begin(), sequences.end());

        id_to_pos_index.resize(sequences.size());
        std::uint32_t pos_index = 0;

        for (auto &s : sequences)
        {
            id_to_pos_index[s->id] = pos_index++;
        }

        minimizer_engine.Minimize(sequences.begin(), sequences.end(), MINHASH_BOOL);
        minimizer_engine.Filter(freq);
    }

    Aligner::align_result Aligner::align_to_target_clip(std::vector<std::string> &queries, std::string &target,
                                                        std::vector<std::pair<std::uint32_t, std::uint32_t>> &clips, std::vector<std::vector<std::uint32_t>> &inserters,
                                                        std::vector<std::uint32_t> &ins_at_least2, bool has_hap, std::vector<EdlibAlignResult> &edlib_results)
    {
        align_result result;
        // Fill in the bases from the target into the first row

        result.target_columns.resize(target.size());
        for (std::uint32_t i = 0; i < result.target_columns.size(); i++)
        { // for each column
            auto &column = result.target_columns[i];
            column.resize(queries.size() + 1, '_');
            column[0] = target[i]; // the first item of each row - makes up the first row
        }
        // Allowance for insertion columns after each target position
        result.ins_columns.resize(target.size());

        // Initialize width of the alignment to the length of the target
        result.width = target.size(); // will be incremented along with adding insertion columns

        // Initialize the record of who is inserting at each position

        // Aligning queries to target and filling up/adding columns
        for (std::uint32_t k = 0; k < queries.size(); k++)
        {
            auto &clip = clips[k];
            std::uint32_t left_clip = clip.first;
            std::uint32_t right_clip = clip.second;
            auto &query = queries[k];
            EdlibAlignResult &align = edlib_results[k];
            /*edlibAlign(query.c_str(), query.size(),
                target.c_str() + left_clip, target.size() - left_clip - right_clip,
                edlibNewAlignConfig(-1, EDLIB_MODE_HW, EDLIB_TASK_PATH, NULL, 0));*/

            std::uint32_t next_query_index = 0;                                    // upcoming query position
            std::uint32_t next_target_index = align.startLocations[0] + left_clip; // upcoming target position
            std::uint32_t next_ins_index = 0;                                      // upcoming ins index for the current segment
            std::uint32_t align_start = next_target_index;
            std::uint32_t align_end = static_cast<std::uint32_t>(align.endLocations[0] + left_clip);

            for (int i = 0; i < align.alignmentLength; i++)
            {

                switch (align.alignment[i])
                {
                case 0:
                { // match

                    // place matching query base from the query into column

                    result.target_columns[next_target_index++][k + 1] = query[next_query_index++];

                    next_ins_index = 0;
                    break;
                }
                case 1:
                { // ins
                    std::uint32_t ins_columns_index = next_target_index - 1;

                    // insertions in the query to the sides of target is ignored
                    if (next_target_index == 0)
                    {
                        next_query_index++;
                        continue;
                    }
                    else if (next_target_index == target.size())
                    {
                        i = align.alignmentLength; // if at ins of the query to the right of target, terminate loop early
                        continue;
                    }

                    // check if we have enough columns ready for this target position
                    auto &ins_columns = result.ins_columns[ins_columns_index];
                    if (ins_columns.size() < next_ins_index + 1)
                    {
                        // ignore extra ins column due to hap sequence
                        if (has_hap && k >= queries.size() - 2)
                        {
                            next_query_index++;
                            continue;
                        }
                        // if not enough, create new column
                        ins_columns.resize(next_ins_index + 1);
                        ins_columns[next_ins_index].resize(queries.size() + 1, '_');
                        result.width++;
                    }

                    // Record existence of ins at certain positions
                    if (next_ins_index == 0)
                    {
                        auto &v = inserters[ins_columns_index];
                        v.push_back(k);
                        if (v.size() == 2)
                        {
                            ins_at_least2.push_back(ins_columns_index);
                        }
                    }

                    ins_columns[next_ins_index++][k + 1] = query[next_query_index++];
                    break;
                }
                case 2:
                { // del
                    next_target_index++;
                    next_ins_index = 0;
                    break;
                }
                case 3:
                { // mismatch
                    // place mismatching query base into column
                    result.target_columns[next_target_index++][k + 1] = query[next_query_index++];
                    next_ins_index = 0;
                    break;
                }
                default:
                {
                    std::cerr << "Unknown alignment result by edlib!" << std::endl;
                }
                }
            }

            edlibFreeAlignResult(align);
        }

        return result;
    }

    Aligner::align_result Aligner::align_to_target_no_clip(std::vector<std::string> &queries, std::string &target, bool has_hap)
    {

        align_result result;
        // Fill in the bases from the target into the first row
        result.target_columns.resize(target.size());
        for (std::uint32_t i = 0; i < result.target_columns.size(); i++)
        { // for each column
            auto &column = result.target_columns[i];
            column.resize(queries.size() + 1, '_');
            column[0] = target[i]; // the first item of each row - makes up the first row
        }

        // Allowance for insertion columns after each target position
        result.ins_columns.reserve(target.size() + 1);
        result.ins_columns.resize(target.size());

        // Initialize width of the alignment to the length of the target
        result.width = target.size(); // will be incremented along with adding insertion columns

        // Aligning queries to target and filling up/adding columns
        for (std::uint32_t k = 0; k < queries.size(); k++)
        {

            auto &query = queries[k];
            EdlibAlignResult align = edlibAlign(query.c_str(), query.size(),
                                                target.c_str(), target.size(), edlibNewAlignConfig(-1, EDLIB_MODE_NW, EDLIB_TASK_PATH, NULL, 0));
            std::uint32_t next_query_index = 0;                        // upcoming query position
            std::uint32_t next_target_index = align.startLocations[0]; // upcoming target position
            std::uint32_t next_ins_index = 0;                          // upcoming ins index for the current segment
            std::uint32_t align_start = next_target_index;
            std::uint32_t align_end = static_cast<std::uint32_t>(align.endLocations[0]);
            bool contained = true;

            for (int i = 0; i < align.alignmentLength; i++)
            {

                switch (align.alignment[i])
                {
                case 0:
                { // match

                    // place matching query base from the query into column

                    result.target_columns[next_target_index++][k + 1] = query[next_query_index++];

                    next_ins_index = 0;
                    break;
                }
                case 1:
                { // ins

                    std::uint32_t ins_columns_index = next_target_index - 1;

                    // Ins to left of target.
                    if (next_target_index == 0)
                    {
                        contained = false;
                        if (has_hap && k >= queries.size() - 2)
                        {
                            next_query_index++;
                            continue;
                        }
                        result.ins_columns.resize(target.size() + 1);
                        ins_columns_index = target.size();
                    }

                    // check if we have enough columns ready for this target position
                    auto &ins_columns = result.ins_columns[ins_columns_index];
                    if (ins_columns.size() < next_ins_index + 1)
                    {
                        if (has_hap && k >= queries.size() - 2)
                        {
                            next_query_index++;
                            continue;
                        }
                        // if not enough, create new column
                        ins_columns.resize(next_ins_index + 1);
                        ins_columns[next_ins_index].resize(queries.size() + 1, '_');
                        result.width++;
                    }
                    ins_columns[next_ins_index++][k + 1] = query[next_query_index++];
                    break;
                }
                case 2:
                { // del

                    next_target_index++;
                    next_ins_index = 0;
                    break;
                }
                case 3:
                { // mismatch

                    // place mismatching query base into column
                    result.target_columns[next_target_index++][k + 1] = query[next_query_index++];
                    next_ins_index = 0;
                    break;
                }
                default:
                {
                    std::cerr << "Unknown alignment result by edlib!" << std::endl;
                }
                }
            }

            edlibFreeAlignResult(align);
        }

        return result;
    }

    Aligner::align_result Aligner::multi_align(std::vector<std::string> &queries, std::string &target,
                                               std::vector<std::pair<std::uint32_t, std::uint32_t>> &clips, std::vector<EdlibAlignResult> &edlib_results, bool has_hap)
    {
        std::vector<std::vector<std::uint32_t>> inserters; // for each position, contains positional index of queries that has ins there
        inserters.resize(target.size());
        std::vector<std::uint32_t> ins_at_least2; // target positions where there are at least 2 reads with ins
        align_result result = align_to_target_clip(queries, target, clips, inserters, ins_at_least2, has_hap, edlib_results);
        for (auto &ins_pos : ins_at_least2)
        { // for all positions that have at least two queries with insertions,
            // extract the ins segments, and record the longest one
            std::uint32_t position_index_longest = 0; // index in ins_segments
            std::uint32_t longest_ins_len = 0;

            std::vector<std::string> ins_segments;
            auto inserters_at_pos = inserters[ins_pos]; // all the inserters here
            ins_segments.resize(inserters_at_pos.size());
            auto &ins_columns = result.ins_columns[ins_pos];
            for (auto &column : ins_columns)
            { // for each column, extract the non gap bases for the queries with ins there
                for (std::uint32_t i = 0; i < inserters_at_pos.size(); i++)
                {
                    char &base = column[inserters_at_pos[i] + 1];
                    auto &segment = ins_segments[i];
                    if (base != '_')
                    {
                        segment.push_back(base);
                        if (segment.size() > longest_ins_len && (!has_hap || inserters_at_pos[i] < queries.size() - 2))
                        {
                            position_index_longest = i;
                            longest_ins_len = segment.size();
                        }
                    }
                }
            }

            // align the rest to the longest segment
            std::uint32_t query_position_longest_ins = inserters_at_pos[position_index_longest];
            inserters_at_pos.erase(inserters_at_pos.begin() + position_index_longest); // exclude longest from list of queries
            std::string longest_ins_segment = std::move(ins_segments[position_index_longest]);
            ins_segments.erase(ins_segments.begin() + position_index_longest);
            auto sub_result = align_to_target_no_clip(ins_segments, longest_ins_segment, has_hap);

            std::uint32_t t = ins_columns.size();

            // get more ins columns for a target position if needed
            ins_columns.resize(sub_result.width);
            result.width += sub_result.width - t;
            for (; t < ins_columns.size(); t++)
            {
                ins_columns[t].resize(queries.size() + 1, '_');
            }

            std::uint32_t to_column_id = 0;

            // if there is ins to the left of target...
            if (sub_result.ins_columns.size() > sub_result.target_columns.size())
            {
                auto &sub_ins_columns = sub_result.ins_columns[sub_result.target_columns.size()];
                for (std::uint32_t j = 0; j < sub_ins_columns.size(); j++)
                {
                    auto &to_column = ins_columns[to_column_id++];
                    auto &from_column = sub_ins_columns[j];
                    to_column[query_position_longest_ins + 1] = from_column[0];
                    for (std::uint32_t j = 0; j < inserters_at_pos.size(); j++)
                    {
                        to_column[inserters_at_pos[j] + 1] = from_column[j + 1];
                    }
                }
            }

            for (std::uint32_t i = 0; i < sub_result.target_columns.size(); i++)
            {
                auto &to_column = ins_columns[to_column_id++];
                auto &from_column = sub_result.target_columns[i];

                // traversing the from columns
                to_column[query_position_longest_ins + 1] = from_column[0]; // the target is the top row of the alignment

                for (std::uint32_t j = 0; j < inserters_at_pos.size(); j++)
                {
                    to_column[inserters_at_pos[j] + 1] = from_column[j + 1];
                }

                auto &sub_ins_columns = sub_result.ins_columns[i];

                for (std::uint32_t j = 0; j < sub_ins_columns.size(); j++)
                {
                    auto &to_column = ins_columns[to_column_id++];

                    auto &from_column = sub_ins_columns[j];

                    to_column[query_position_longest_ins + 1] = from_column[0];

                    for (std::uint32_t j = 0; j < inserters_at_pos.size(); j++)
                    {
                        to_column[inserters_at_pos[j] + 1] = from_column[j + 1];
                    }
                }
            }
        }

        return result;
    }

    void Aligner::align_result::print()
    {
        std::uint32_t block_width = 100;
        std::vector<std::vector<std::string>> output_blocks;
        std::uint32_t num_rows = this->target_columns[0].size() * 2 - 1;
        std::uint32_t num_blocks = this->width / block_width;
        num_blocks = this->width % block_width == 0 ? num_blocks : num_blocks + 1;
        output_blocks.resize(num_blocks);
        for (auto &block : output_blocks)
        {
            block.resize(num_rows);
            for (auto &row : block)
            {
                row.reserve(width);
            }
        }
        std::vector<std::uint32_t> start_of_blocks;
        start_of_blocks.reserve(num_blocks);
        start_of_blocks.push_back(0);
        std::uint32_t col_index = 0;
        std::uint32_t current_block_index = 0;
        if (this->ins_columns.size() > this->target_columns.size())
        {

            auto &ins_cols = this->ins_columns[this->target_columns.size()];
            for (std::uint32_t j = 0; j < ins_cols.size(); j++)
            { // for each ins column here
                auto &ins_column = ins_cols[j];
                std::uint32_t block_index = col_index++ / block_width;
                auto &block = output_blocks[block_index];
                if (block_index > current_block_index)
                {
                    current_block_index++;
                    start_of_blocks.push_back(0);
                }
                for (std::uint32_t k = 0; k < num_rows; k++)
                { // move down the column
                    if (k % 2 == 0)
                    {
                        block[k].push_back(ins_column[k / 2]);
                        if (k != 0)
                        {
                            if (block[k].back() != '_' && block[k - 2].back() != '_' && block[k].back() != block[k - 2].back())
                            {
                                block[k - 1].push_back('!');
                            }
                            else
                            {
                                block[k - 1].push_back(' ');
                            }
                        }
                    }
                }
            }
        }
        for (std::uint32_t i = 0; i < this->target_columns.size(); i++)
        { // for each target position
            auto &column = this->target_columns[i];

            std::uint32_t block_index = col_index++ / block_width;
            auto &block = output_blocks[block_index];
            if (block_index > current_block_index)
            {
                current_block_index++;
                start_of_blocks.push_back(i);
            }

            for (std::uint32_t k = 0; k < num_rows; k++)
            { // move down the column

                if (k % 2 == 0)
                {
                    block[k].push_back(column[k / 2]);

                    if (k != 0)
                    {

                        if (block[k].back() != '_' && block[k - 2].back() != '_' && block[k].back() != block[k - 2].back())
                        {

                            block[k - 1].push_back('!');
                        }
                        else
                        {

                            block[k - 1].push_back(' ');
                        }
                    }
                }
            }

            auto &ins_cols = this->ins_columns[i];

            for (std::uint32_t j = 0; j < ins_cols.size(); j++)
            { // for each ins column here
                std::uint32_t block_index = col_index++ / block_width;
                auto &block = output_blocks[block_index];
                if (block_index > current_block_index)
                {
                    current_block_index++;
                    start_of_blocks.push_back(i);
                }

                auto &ins_column = ins_cols[j];
                for (std::uint32_t k = 0; k < num_rows; k++)
                { // move down the column
                    if (k % 2 == 0)
                    {
                        block[k].push_back(ins_column[k / 2]);
                        if (k != 0)
                        {
                            if (block[k].back() != '_' && block[k - 2].back() != '_' && block[k].back() != block[k - 2].back())
                            {
                                block[k - 1].push_back('!');
                            }
                            else
                            {
                                block[k - 1].push_back(' ');
                            }
                        }
                    }
                }
            }
        }

        std::uint32_t counter = 0;
        // printing
        for (auto &block : output_blocks)
        {
            std::cout << start_of_blocks[counter++] << std::endl;
            for (auto &row : block)
            {
                std::cout << row << std::endl;
            }
            std::cout << std::endl;
        }
    }

} // namespace align_reads
