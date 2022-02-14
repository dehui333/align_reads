#include <iostream>

#include "bioparser/fasta_parser.hpp"
#include "bioparser/fastq_parser.hpp"
#include "bioparser/paf_parser.hpp"


struct FastqSequence {
    std::string name;
    std::string seq;
    std::string qual;
  
    FastqSequence(
      const char* n, std::uint32_t ln,
      const char* s, std::uint32_t ls,
      const char* q, std::uint32_t lq
    ) :name(n, ln), seq(s, ls), qual(q, lq) {}
};

struct FastaSequence {
    std::string name;
    std::string seq;
  
    FastaSequence(
      const char* n, std::uint32_t ln,
      const char* s, std::uint32_t ls
    ) :name(n, ln), seq(s, ls){}
};


struct Overlap {
    std::string q_name;
    std::uint32_t q_len;
    std::uint32_t q_start;
    std::uint32_t q_end;
    char dir;
    std::string t_name;
    std::uint32_t t_len;
    std::uint32_t t_start;
    std::uint32_t t_end;
    std::uint32_t num_match;
    std::uint32_t overlap_len;
    std::uint32_t mapq;
    
    Overlap(  
        const char* qn, std::uint32_t qnl,
        std::uint32_t ql, 
        std::uint32_t qs, std::uint32_t qe,
        char d,
        const char* tn, std::uint32_t tnl,
        std::uint32_t tl,
        std::uint32_t ts, std::uint32_t te,
        std::uint32_t nm,
        std::uint32_t ol,
        std::uint32_t q
    ) : q_name(qn, qnl), q_len(ql), q_start(qs), q_end(qe), dir(d), t_name(tn, tnl),
        t_len(tl), t_start(ts), t_end(te), num_match(nm), overlap_len(ol), mapq(q) {}
};