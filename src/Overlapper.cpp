#include "ram/minimizer_engine.hpp"

#include "align_reads/Overlapper.hpp"

#define MINHASH_BOOL true
#define MINIMIZER_K 15
#define MINIMIZER_W 5
#define MINIMIZER_B 500
#define MINIMIZER_N 4
#define MINIMIZER_M 100
#define MINIMIZER_G 10000

namespace align_reads
{
    Overlapper::Overlapper(std::uint8_t num_groups,
                           std::shared_ptr<thread_pool::ThreadPool> &pool)
    {
        indices.reserve(num_groups);
        for (auto i = 0; i < num_groups; i++)
        {
            indices.emplace_back(pool, MINIMIZER_K, MINIMIZER_W, MINIMIZER_B, MINIMIZER_N,
                                 MINIMIZER_M, MINIMIZER_G);
        }
    }
    std::vector<biosoup::Overlap> Overlapper::find_overlaps(std::unique_ptr<biosoup::NucleicAcid> &sequence,
                                                            std::uint32_t index_id)
    {
        std::vector<biosoup::Overlap> overlaps = indices[index_id].Map(sequence, true, false, MINHASH_BOOL);
        return overlaps;
    }

    void Overlapper::index_sequences(std::vector<std::unique_ptr<biosoup::NucleicAcid>> &sequences,
                                     std::uint32_t index_id)
    {
        indices[index_id].Minimize(sequences.begin(), sequences.end(), MINHASH_BOOL);
    }


} // namespace align_reads