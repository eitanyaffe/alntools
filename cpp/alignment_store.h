#pragma once

#include "aln_types.h"
#include <vector>
#include <string>
#include <memory>
#include <unordered_map>
#include <cstdint>
#include <fstream>

using std::string;
using std::unordered_map;
using std::vector;

// Enum for mutation position integer type
enum class MutationIntType
{
    UINT16,
    UINT32
};

class AlignmentStore
{
private:
    std::vector<Contig> contigs_;
    std::vector<Read> reads_;
    std::vector<Alignment> alignments_;
    unordered_map<string, size_t> read_id_to_index;
    unordered_map<string, size_t> contig_id_to_index;

    // Helper functions for saving/loading positions
    static void save_position(std::ofstream &file, uint32_t position, MutationIntType type);
    static uint32_t load_position(std::ifstream &file, MutationIntType type);

public:
    // Add methods
    void add_contig(const Contig &contig) { contigs_.push_back(contig); }
    void add_read(const Read &read) { reads_.push_back(read); }
    void add_alignment(const Alignment &alignment) { alignments_.push_back(alignment); }

    // Getter methods
    const std::vector<Contig> &get_contigs() const { return contigs_; }
    const std::vector<Read> &get_reads() const { return reads_; }
    const std::vector<Alignment> &get_alignments() const { return alignments_; }

    // Non-const getters for modification
    std::vector<Contig> &get_contigs() { return contigs_; }
    std::vector<Read> &get_reads() { return reads_; }
    std::vector<Alignment> &get_alignments() { return alignments_; }

    void export_tab_delimited(const string &prefix);

    // Add missing methods
    void save(const string &filename);
    void load(const string &filename);
    size_t get_alignment_count() const { return alignments_.size(); }
    size_t get_read_count() const { return reads_.size(); }

    size_t add_or_get_read_index(const string &read_id, uint32_t length);
    size_t add_or_get_contig_index(const string &contig_id, uint32_t length);

    size_t get_read_index(const string &read_id);
    size_t get_contig_index(const string &contig_id);

    // Get id by index
    const string &get_read_id(size_t read_index) const;
    const string &get_contig_id(size_t contig_index) const;
};
