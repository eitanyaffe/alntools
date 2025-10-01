#ifndef REARRANGE_READ_EVENT_H
#define REARRANGE_READ_EVENT_H

#include "alignment_store.h"
#include "utils.h"
#include "Sequences.h"
#include <string>
#include <vector>
#include <unordered_map>
#include <set>
#include <map>

// Forward declaration
class RearrangeVerify;

// Event test results for rejection tracking
enum class EventTestResult {
    FOUND_EVENT,
    REJECT_ALIGNMENT_OVERLAP,
    REJECT_INSERT_ELEMENT_REINSERTED,
    REJECT_INSERT_ELEMENT_OVERLAPS,
    REJECT_ELEMENT_TOO_SHORT,
    REJECT_CONTAINED_ELEMENT_MARGIN_TOO_LARGE,
    REJECT_READ_GAP_TOO_LARGE,
    REJECT_INSERT_GAP_TOO_LARGE
};

// convert EventTestResult to string for rejection counting
std::string event_test_result_to_string(EventTestResult result);

enum class RearrangementType {
    LARGE_INSERT,
    LARGE_INVERT,
    LARGE_DELETE
};

// add operator<< for RearrangementType
inline std::ostream& operator<<(std::ostream& os, const RearrangementType& type)
{
    switch (type) {
    case RearrangementType::LARGE_INSERT:
        os << "large_insert";
        break;
    case RearrangementType::LARGE_INVERT:
        os << "large_invert";
        break;
    case RearrangementType::LARGE_DELETE:
        os << "large_delete";
        break;
    default:
        os << "unknown";
        break;
    }
    return os;
}

struct ReadEvent {
    std::string read_id;
    RearrangementType type;
    std::string contig_id;        // contig of anchors A,B
    std::string read_strand;      // "+" or "-" for read orientation
    uint32_t out_clip;           // contig coordinates
    uint32_t in_clip;            // contig coordinates  
    uint32_t read_clip_out;      // read coordinates
    uint32_t read_clip_in;       // read coordinates
    uint32_t span_start;         // contig coordinates - start of A
    uint32_t span_end;           // contig coordinates - end of B
    uint32_t read_span_start;    // read coordinates - start of A
    uint32_t read_span_end;      // read coordinates - end of B
    std::string element_contig;   // contig of element X (empty for deletions)
    std::string element_strand;  // "+" or "-" for element X (empty for deletions)
    uint32_t element_start;      // contig coordinates (0 for deletions)
    uint32_t element_end;        // contig coordinates (0 for deletions)
    std::string read_seams;      // seams in read coordinates and sequences (populated by resolve_seams)
    std::string assembly_seams;  // seams in assembly coordinates and sequences
    std::vector<SeamInterval> read_seam_intervals;  // intervals to be resolved
    
    // alignment copies for verification (A, B, and optionally X)
    Alignment alignment_A;
    Alignment alignment_B;
    bool has_alignment_X;
    Alignment alignment_X;  // only valid if has_alignment_X is true
    
    ReadEvent(const std::string& read_id = "", RearrangementType type = RearrangementType::LARGE_DELETE,
              const std::string& contig_id = "", const std::string& read_strand = "+",
              uint32_t out_clip = 0, uint32_t in_clip = 0,
              uint32_t read_clip_out = 0, uint32_t read_clip_in = 0,
              uint32_t span_start = 0, uint32_t span_end = 0,
              uint32_t read_span_start = 0, uint32_t read_span_end = 0,
              const std::string& element_contig = "", const std::string& element_strand = "",
              uint32_t element_start = 0, uint32_t element_end = 0,
              const std::string& read_seams = "", const std::string& assembly_seams = "")
        : read_id(read_id), type(type), contig_id(contig_id), read_strand(read_strand),
          out_clip(out_clip), in_clip(in_clip), read_clip_out(read_clip_out), read_clip_in(read_clip_in),
          span_start(span_start), span_end(span_end), read_span_start(read_span_start), read_span_end(read_span_end),
          element_contig(element_contig), element_strand(element_strand),
          element_start(element_start), element_end(element_end),
          read_seams(read_seams), assembly_seams(assembly_seams),
          has_alignment_X(false) {}
};

// worker class for detecting events in a single library
class RearrangeReadEvent {
private:
    const AlignmentStore& store;
    
    // parameters
    int max_margin;
    int min_element_length;
    int min_anchor_length;
    double max_anchor_mutations_percent;
    double max_element_mutation_percent;
    RearrangeVerify* verifier;
    ResolveSeams resolve_seams;
    std::string current_lib_id; // current library being processed
    
    // sequence references
    const AssemblySequences* assembly_sequences;
    const ReadSequences* read_sequences;
    
    // interval filtering data
    std::vector<Interval> intervals;
    std::set<uint32_t> relevant_read_indices;
    std::set<const Alignment*> relevant_alignments;

public:
    RearrangeReadEvent(const AlignmentStore& store,
                      RearrangeVerify* verifier = nullptr,
                      ResolveSeams resolve_seams = ResolveSeams::NO,
                      const AssemblySequences* assembly_seqs = nullptr,
                      const ReadSequences* read_seqs = nullptr,
                      int max_margin = 10,
                      int min_element_length = 50,
                      int min_anchor_length = 200,
                      double max_anchor_mutations_percent = 0.01,
                      double max_element_mutation_percent = 0.01);
    
    // main detection function - pass output containers by reference
    void detect_events(const std::string& lib_id, std::vector<ReadEvent>& output_events, 
                      std::map<std::string, size_t>& output_rejection_counts, const AlignmentStore& store,
                      const std::string& output_prefix = "");
    
    // interval filtering
    void set_intervals(const std::vector<Interval>& intervals);
    
    // set current library ID
    void set_current_lib(const std::string& lib_id) { current_lib_id = lib_id; }
    
    // utility methods
    size_t get_total_rejected_count(const std::map<std::string, size_t>& rejection_counts) const;
    
    // get read sequences for events (collects all read IDs and populates sequences map)
    void get_read_sequences(const std::string& lib_id, const std::vector<ReadEvent>& events, 
                           std::unordered_map<std::string, std::string>& sequences, const std::string& output_prefix = "");
    
    // resolve read seams for events (moved from ReadSequences)
    void resolve_read_seams(const std::string& lib_id, std::vector<ReadEvent>& events, 
                           const std::unordered_map<std::string, std::string>& sequences);
    
    
    // debug methods
    static void print_event(const ReadEvent& event, const AlignmentStore& store);
    void print_read_alignments(const std::vector<size_t>& alignment_indices,
                              const std::vector<Alignment>& alignments,
                              const std::vector<bool>& is_valid_anchor_vec) const;

private:
    // alignment validation
    bool is_valid_anchor(const Alignment& aln, bool verbose = false, size_t index = 0) const;
    bool is_valid_element(const Alignment& aln) const;
    
    // event testing functions
    EventTestResult test_trio_event(const Alignment& A, const Alignment& B, const Alignment& X, 
                                   bool is_read_reversed, ReadEvent& event,
                                   std::map<std::string, size_t>& rejection_counts,
                                   bool verbose = false);
    EventTestResult test_pair_event(const Alignment& A, const Alignment& B, 
                                   bool is_read_reversed, ReadEvent& event,
                                   std::map<std::string, size_t>& rejection_counts,
                                   bool verbose = false);
    
    // rejection tracking
    void add_rejection(EventTestResult result, std::map<std::string, size_t>& rejection_counts, bool verbose = false);
    
    // process single read functions
    bool process_read(const std::vector<size_t>& alignment_indices, 
                     const std::vector<Alignment>& alignments,
                     std::vector<ReadEvent>& output_events,
                     std::map<std::string, size_t>& rejection_counts,
                     bool verbose = false);
    bool process_read_trios(const std::vector<size_t>& alignment_indices,
                           const std::vector<Alignment>& alignments,
                           const std::vector<bool>& is_valid_anchor_vec,
                           std::vector<ReadEvent>& output_events,
                           std::map<std::string, size_t>& rejection_counts,
                           bool verbose = false);
    bool process_read_pairs(const std::vector<size_t>& alignment_indices,
                           const std::vector<Alignment>& alignments,
                           const std::vector<bool>& is_valid_anchor_vec,
                           std::vector<ReadEvent>& output_events,
                           std::map<std::string, size_t>& rejection_counts,
                           bool verbose = false);
    
    // helper functions
    std::string alignment_to_strand_string(bool is_reverse) const;
    bool alignments_overlap(const Alignment& A, const Alignment& B, bool verbose = false) const;
    ReadEvent create_event(const Alignment& A, const Alignment& B, const Alignment* X, 
                          RearrangementType type, bool is_read_reversed, bool verbose = false) const;
    
    // create seam intervals for seams between alignments
    std::vector<SeamInterval> get_seam_intervals(const Alignment& A, const Alignment& B, const Alignment* X,
                                               RearrangementType type, bool is_read_reversed, 
                                               const std::string& lib_id) const;
    
    // get assembly seams for events
    std::string get_assembly_seams(const Alignment& A, const Alignment& B, const Alignment* X,
                                 RearrangementType type) const;
    
    // helper functions for interval creation
    SeamInterval create_seam_interval(const Alignment& align1, const Alignment& align2,
                                    const std::string& lib_id, bool is_read_reversed) const;
    std::string extract_assembly_seam(const Alignment& align1, const Alignment& align2) const;
    
    // helper for assembly seam extraction
    std::string extract_seam(const std::string& sequence, uint32_t end_left, uint32_t start_right, bool use_placeholders = false) const;
    
    // helpers for read seam resolution (moved from ReadSequences)
    std::string extract_sequence_interval(const std::string& sequence, const SeamInterval& interval) const;
    
    // rejection tracking
    void add_rejection(const std::string& reason, std::map<std::string, size_t>& rejection_counts);
};

#endif // REARRANGE_READ_EVENT_H
