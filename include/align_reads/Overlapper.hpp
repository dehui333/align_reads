#ifndef ALIGN_READS_OVERLAPPER_HPP_
#define ALIGN_READS_OVERLAPPER_HPP_

#include "ram/minimizer_engine.hpp"

#include "align_reads/Overlapper.hpp"

namespace align_reads
{
    class Overlapper
    {
    public:
        Overlapper(std::uint8_t num_groups, std::shared_ptr<thread_pool::ThreadPool> &pool);

        // Index the sequences so that other sequences can be matched with them
        void index_sequences(std::vector<std::unique_ptr<biosoup::NucleicAcid>> &sequences,
                             std::uint32_t index_id);

        // Find overlaps to the input sequence with sequences indexed
        // -> maybe can filter overlaps
        std::vector<biosoup::Overlap> find_overlaps(std::unique_ptr<biosoup::NucleicAcid> &sequence,
                                                    std::uint32_t index_id);

    private:
        std::vector<ram::MinimizerEngine> indices;
    };

} // namespace align_reads

#endif // ALIGN_READS_OVERLAPPER_HPP_