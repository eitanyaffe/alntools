#ifndef VERIFY_REARRANGE_H
#define VERIFY_REARRANGE_H

#include "alignment_store.h"
#include "aln_types.h"
#include "Rearrange.h"
#include <string>
#include <unordered_map>

class VerifyRearrange {
private:
    std::unordered_map<std::string, std::string> contigs;
    std::unordered_map<std::string, std::string> reads;
    bool sequences_loaded;
    
public:
    VerifyRearrange();
    
    // Load sequences from FASTA files
    void load_sequences(const std::string& ifn_contigs, const std::string& ifn_reads);
    
    // Load sequences from maps (for R interface)
    void load_sequences_from_maps(const std::unordered_map<std::string, std::string>& contig_map,
                                 const std::unordered_map<std::string, std::string>& read_map);
    
    // Verify a single rearrangement event
    bool verify_event(const ReadEvent& event, const Alignment& A, const Alignment& B, 
                     const Alignment* X, const AlignmentStore& store);
    
private:
    // Extract and apply mutations to a contig segment
    std::string get_mutated_contig_segment(const Alignment& alignment, const AlignmentStore& store);
    
    // Execute the rearrangement operation (insert/delete/invert)
    std::string execute_rearrangement(const ReadEvent& event, const std::string& segment_A,
                                    const std::string& segment_B, const std::string& segment_X);
    
    // Get the expected read segment for comparison
    std::string get_read_segment(const ReadEvent& event, const AlignmentStore& store);
    
    // Compare sequences and report mismatches
    bool compare_sequences(const std::string& expected, const std::string& actual, 
                          const std::string& event_desc);
};

#endif // VERIFY_REARRANGE_H
