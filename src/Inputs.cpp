#include <iostream> // std::cerr

#include "bioparser/fasta_parser.hpp"
#include "bioparser/fastq_parser.hpp"
#include "biosoup/nucleic_acid.hpp"

#include "align_reads/Inputs.hpp"

/*
 *   -> Beware of counter variable biosoup::NucleicAcid::num_objects
 *
 */

namespace align_reads
{
    auto is_suffix = [](const std::string &s, const std::string &suff)
    {
        return s.size() < suff.size() ? false : s.compare(s.size() - suff.size(), suff.size(), suff) == 0;
    };

    std::vector<std::unique_ptr<biosoup::NucleicAcid>> read_fasta(const char *path)
    {
        try
        {
            auto parser = bioparser::Parser<biosoup::NucleicAcid>::Create<bioparser::FastaParser>(path);
            return parser->Parse(-1);
        }
        catch (const std::invalid_argument &exception)
        {
            std::cerr << exception.what() << std::endl;
        }
        return {};
    }

    std::vector<std::unique_ptr<biosoup::NucleicAcid>> read_fastq(const char *path)
    {
        try
        {
            auto parser = bioparser::Parser<biosoup::NucleicAcid>::Create<bioparser::FastqParser>(path);
            return parser->Parse(-1);
        }
        catch (const std::invalid_argument &exception)
        {
            std::cerr << exception.what() << std::endl;
        }
        return {};
    }

    std::vector<std::unique_ptr<biosoup::NucleicAcid>> read_fasta_or_fastq(const char *path)
    {
        if (is_suffix(path, "fasta"))
        {
            return read_fasta(path);
        }
        else if (is_suffix(path, "fastq"))
        {
            return read_fastq(path);
        }
        else
        {
            throw std::invalid_argument("Error: Invalid sequences file format.");
        }
    }

    /* Maps id of sequences to their position in the vector.
     * Meaningful only when indices of sequences are out of order(parsing in parallel?).
     * Input: sequences with ids in {0, 1, ... sequences.size() - 1}
     * Output: vector with nth element being the position index of the sequence with id n
     * in 'sequences'.
     */
    std::vector<std::uint32_t> construct_index(std::vector<std::unique_ptr<biosoup::NucleicAcid>> &sequences)
    {
        std::vector<std::uint32_t> id_to_pos_index(sequences.size());
        std::uint32_t pos_index = 0;
        for (auto &s : sequences)
        {
            id_to_pos_index[s->id] = pos_index++;
        }
        return id_to_pos_index;
    }

    /*
     * Read a batch of sequences into a destination vector, possibly from multiple paths (in parallel).
     * The sequences will occupy contiguous ids starting from 'start_id' but may be out of order
     * if a pool is used.
     * Only 1 instance of this function should run at a time.
     *
     */
    void read_sequences(std::vector<std::unique_ptr<biosoup::NucleicAcid>> &dest, std::vector<std::string> &paths, std::uint32_t start_id, std::shared_ptr<thread_pool::ThreadPool> &pool)
    {
        if (paths.empty())
            return;
        biosoup::NucleicAcid::num_objects = start_id;
        if (paths.size() == 1)
        {
            auto &path = paths[0];
            auto seqs = read_fasta_or_fastq(path.c_str());
            dest.insert(dest.end(), std::make_move_iterator(seqs.begin()), std::make_move_iterator(seqs.end()));
        }
        else if (pool != nullptr)
        {
            // read in parallel
            struct item
            {
                std::vector<std::unique_ptr<biosoup::NucleicAcid>> sequences;
                item(std::vector<std::unique_ptr<biosoup::NucleicAcid>> &sequences) : sequences(std::move(sequences)) {}
            };

            auto thread_task = [](const char *path)
            {
                auto sequences = read_fasta_or_fastq(path);
                item it{sequences};
                return it;
            };           
            std::vector<std::future<item>> futures;
            for (auto &path : paths)
            {
                futures.emplace_back(pool->Submit(thread_task, path.c_str()));
            }

            for (auto &f : futures)
            {
                auto part = f.get().sequences;
                dest.insert(dest.end(), std::make_move_iterator(part.begin()), std::make_move_iterator(part.end()));
            }
        }
        else
        {
            // read sequentially
            for (auto &path : paths)
            {
                auto part = read_fasta_or_fastq(path.c_str());
                dest.insert(dest.end(), std::make_move_iterator(part.begin()), std::make_move_iterator(part.end()));
            }
        }
    }

} // namespace align_reads
