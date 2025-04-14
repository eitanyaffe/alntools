#pragma once

#include "aln_types.h"
#include <vector>
#include <string>
#include <memory>

class AlignmentStore 
{
private:
    std::vector<Contig> contigs_;
    std::vector<Read> reads_;
    std::vector<Alignment> alignments_;

public:
    // Add methods
    void add_contig(const Contig& contig) { contigs_.push_back(contig); }
    void add_read(const Read& read) { reads_.push_back(read); }
    void add_alignment(const Alignment& alignment) { alignments_.push_back(alignment); }

    // Getter methods
    const std::vector<Contig>& get_contigs() const { return contigs_; }
    const std::vector<Read>& get_reads() const { return reads_; }
    const std::vector<Alignment>& get_alignments() const { return alignments_; }

    // Non-const getters for modification
    std::vector<Contig>& get_contigs() { return contigs_; }
    std::vector<Read>& get_reads() { return reads_; }
    std::vector<Alignment>& get_alignments() { return alignments_; }
    
    void export_tab_delimited(const string& prefix);
    
    // Add missing methods
    void save(const string& filename);
    void load(const string& filename);
    size_t get_alignment_count() const { return alignments_.size(); }
    size_t get_read_count() const { return reads_.size(); }
}; 