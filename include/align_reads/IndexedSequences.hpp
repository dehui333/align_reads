#ifndef INDEXED_SEQUENCES_HPP_
#define INDEXED_SEQUENCES_HPP_

#include <fstream>
#include <string>
#include <unordered_map>

namespace align_reads
{

    struct indexed_info
    {
        std::uint32_t length;
        std::uint64_t offset;
        std::uint32_t line_bases; 
        // std::uint32_t line_width;
        //std::uint64_t qual_offset;

        indexed_info(std::uint32_t len, std::uint64_t offset, std::uint32_t line_bases)
            : length(len), offset(offset), line_bases(line_bases) {}
    };

    class IndexedSequences
    {
    public:

       
        // Path to fasta file, index file name is path_to_file.fai
        // only supports fasta
        IndexedSequences(std::string &path_to_file);
        IndexedSequences() = default;
        ~IndexedSequences()
        {
            seq_file_stream.close();
        }
        std::string get_sequence(std::string& seq_name);

        void index_reads(std::string& sequence_file_path);
    private:


        std::unordered_map<std::string, indexed_info> index; // seq name to info
        std::ifstream seq_file_stream;
        // unordered map or sth mapping from name to indexed info
        // file handle to fasta/q
    };

}

#endif // INDEXED_SEQUENCES_HPP_