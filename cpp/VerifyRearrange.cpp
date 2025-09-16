#include "VerifyRearrange.h"
#include "utils.h"
#include <iostream>
#include <fstream>

using namespace std;

VerifyRearrange::VerifyRearrange() : sequences_loaded(false)
{
}

void VerifyRearrange::load_sequences(const std::string& ifn_contigs, const std::string& ifn_reads)
{
    cout << "loading sequences for rearrangement verification..." << endl;
    
    // Load contigs
    unordered_set<string> all_contig_ids; // empty set means load all
    read_fasta(ifn_contigs, all_contig_ids, contigs);
    cout << "loaded " << contigs.size() << " contigs" << endl;
    
    // Load reads  
    unordered_set<string> all_read_ids; // empty set means load all
    read_fastq(ifn_reads, all_read_ids, reads);
    cout << "loaded " << reads.size() << " reads" << endl;
    
    sequences_loaded = true;
}

void VerifyRearrange::load_sequences_from_maps(const std::unordered_map<std::string, std::string>& contig_map,
                                              const std::unordered_map<std::string, std::string>& read_map)
{
    cout << "loading sequences from maps for rearrangement verification..." << endl;
    
    // Copy sequences from maps
    contigs = contig_map;
    reads = read_map;
    
    cout << "loaded " << contigs.size() << " contigs and " << reads.size() << " reads from maps" << endl;
    
    sequences_loaded = true;
}

bool VerifyRearrange::verify_event(const ReadEvent& event, const Alignment& A, const Alignment& B,
                                  const Alignment* X, const AlignmentStore& store)
{
    if (!sequences_loaded) {
        cerr << "error: sequences not loaded for verification" << endl;
        return false;
    }
    
    // Get read and contig IDs
    string read_id = store.get_read_id(A.read_index);
    string contig_id = event.contig_id;
    
    // Check that sequences exist
    if (contigs.find(contig_id) == contigs.end()) {
        cerr << "error: contig " << contig_id << " not found in loaded sequences" << endl;
        return false;
    }
    if (reads.find(read_id) == reads.end()) {
        cerr << "error: read " << read_id << " not found in loaded sequences" << endl;
        return false;
    }
    
    // Extract mutated contig segments
    string segment_A = get_mutated_contig_segment(A, store);
    string segment_B = get_mutated_contig_segment(B, store);
    string segment_X = "";
    
    if (X && event.type != RearrangementType::LARGE_DELETE) {
        segment_X = get_mutated_contig_segment(*X, store);
    }
    
    // Execute the rearrangement to get expected sequence
    string expected_sequence = execute_rearrangement(event, segment_A, segment_B, segment_X);
    
    // Get actual read segment
    string actual_sequence = get_read_segment(event, store);
    
    // Compare sequences
    string event_desc = "event " + string(1, static_cast<char>(event.type)) + " in read " + read_id;
    return compare_sequences(expected_sequence, actual_sequence, event_desc);
}

std::string VerifyRearrange::get_mutated_contig_segment(const Alignment& alignment, const AlignmentStore& store)
{
    string contig_id = store.get_contig_id(alignment.contig_index);
    string read_id = store.get_read_id(alignment.read_index);
    
    // Extract contig fragment
    string contig_fragment = contigs[contig_id].substr(alignment.contig_start,
                                                      alignment.contig_end - alignment.contig_start);
    
    // Apply mutations using existing utility function
    string mutated_segment = apply_mutations(contig_fragment, alignment.mutations, store, 
                                           alignment, read_id, contig_id);
    
    return mutated_segment;
}

std::string VerifyRearrange::execute_rearrangement(const ReadEvent& event, const std::string& segment_A,
                                                  const std::string& segment_B, const std::string& segment_X)
{
    string result;
    
    switch (event.type) {
    case RearrangementType::LARGE_INSERT:
        // Insert element X between A and B
        result = segment_A + segment_X + segment_B;
        break;
        
    case RearrangementType::LARGE_INVERT:
        // Insert reverse complement of element X between A and B
        result = segment_A + reverse_complement(segment_X) + segment_B;
        break;
        
    case RearrangementType::LARGE_DELETE:
        // Concatenate A and B directly (no element)
        result = segment_A + segment_B;
        break;
        
    default:
        cerr << "error: unknown rearrangement type" << endl;
        return "";
    }
    
    return result;
}

std::string VerifyRearrange::get_read_segment(const ReadEvent& event, const AlignmentStore& /* store */)
{
    string read_id = event.read_id;
    
    // Extract read segment from read_clip_out to read_clip_in
    // Note: read coordinates are already in the correct orientation from the event
    uint32_t start = event.read_clip_out;
    uint32_t end = event.read_clip_in;
    
    // Ensure start <= end
    if (start > end) {
        swap(start, end);
    }
    
    string read_segment = reads[read_id].substr(start, end - start);
    
    // Apply reverse complement if the event is on reverse strand
    // Following paf_reader pattern: reverse complement the read segment for reverse alignments
    if (event.strand == "-") {
        read_segment = reverse_complement(read_segment);
    }
    
    return read_segment;
}

bool VerifyRearrange::compare_sequences(const std::string& expected, const std::string& actual,
                                       const std::string& event_desc)
{
    // Check length first
    if (expected.size() != actual.size()) {
        cout << "verification failed for " << event_desc << ": length mismatch" << endl;
        cout << "expected length: " << expected.size() << ", actual length: " << actual.size() << endl;
        return false;
    }
    
    // Check sequence content
    for (size_t i = 0; i < expected.size(); ++i) {
        if (expected[i] != actual[i]) {
            // Calculate start and end indices for context (like paf_reader)
            size_t start = (i >= 8) ? i - 8 : 0;
            size_t end = (i + 8 < expected.size()) ? i + 8 : expected.size() - 1;
            
            cout << "verification failed for " << event_desc << ": sequence mismatch at position " << i << endl;
            cout << "expected: " << expected.substr(start, end - start + 1) << endl;
            cout << "actual  : " << actual.substr(start, end - start + 1) << endl;
            cout << "          " << string(i - start, ' ') << "^" << endl;
            return false;
        }
    }
    
    cout << "verification passed for " << event_desc << endl;
    return true;
}
