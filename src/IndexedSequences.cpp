

#include "align_reads/IndexedSequences.hpp"

namespace align_reads
{

    IndexedSequences::IndexedSequences(std::string &sequence_file_name) : seq_file_stream(sequence_file_name)
    {
        std::ifstream index_file_stream(sequence_file_name + ".fai");
        std::string seq_name;
        std::uint32_t seq_len;
        std::uint64_t offset;
        std::uint32_t line_bases;
        std::uint32_t line_width;
        while (index_file_stream >> seq_name >> seq_len >> offset >> line_bases >> line_width)
        {
            index.emplace(std::piecewise_construct, std::forward_as_tuple(std::move(seq_name)),
                          std::forward_as_tuple(seq_len, offset, line_bases));
        }

        index_file_stream.close();
    }

    std::string IndexedSequences::get_sequence(std::string &seq_name)
    {
        auto &info = index.at(seq_name);
        std::string output;
        std::uint64_t offset = info.offset;
        std::uint32_t total_length = info.length;
        std::uint32_t line_bases = info.line_bases;

        output.reserve(total_length);
        std::string buffer;
        std::uint32_t length_read = 0;
        seq_file_stream.seekg(offset);
        while (length_read < total_length)
        {
            seq_file_stream >> buffer;
            output += buffer;
            length_read += line_bases;
        }
        return output;
    }

} // namespace align_reads