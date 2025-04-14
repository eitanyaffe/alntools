#pragma once

#include <string>
#include <vector>
#include "aln_types.h"

class AlignmentStore {
private:
    std::vector<Contig> contigs;
    std::vector<Read> reads;
    std::vector<Alignment> alignments;

public:
    void add_contig(const Contig& contig) {
        contigs.emplace_back(contig);
    }

    void add_read(const Read& read) {
        reads.emplace_back(read);
    }

    void add_alignment(const Alignment& alignment) {
        alignments.emplace_back(alignment);
    }

    const std::vector<Contig>& get_contigs() const {
        return contigs;
    }

    std::vector<Read>& get_reads() {
        return reads;
    }

    const std::vector<Alignment>& get_alignments() const {
        return alignments;
    }
    
    void export_tab_delimited(const std::string& prefix);
    
    // Add missing methods
    void save(const std::string& filename);
    void load(const std::string& filename);
    size_t get_alignment_count() const { return alignments.size(); }
    size_t get_read_count() const { return reads.size(); }
}; 