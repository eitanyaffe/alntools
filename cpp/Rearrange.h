#ifndef REARRANGE_H
#define REARRANGE_H

#include "alignment_store.h"
#include "utils.h"
#include "RearrangeErrors.h"
#include <map>
#include <memory>
#include <set>
#include <string>
#include <unordered_map>
#include <vector>

// Forward declaration
class VerifyRearrange;

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
    std::string strand;           // "+" or "-" for anchors A,B
    uint32_t out_clip;           // contig coordinates
    uint32_t in_clip;            // contig coordinates  
    uint32_t read_clip_out;      // read coordinates
    uint32_t read_clip_in;       // read coordinates
    std::string element_contig;   // contig of element X (empty for deletions)
    std::string element_strand;  // "+" or "-" for element X (empty for deletions)
    uint32_t element_start;      // contig coordinates (0 for deletions)
    uint32_t element_end;        // contig coordinates (0 for deletions)
    std::string left_shim;       // read sequence between A and X (empty for deletions)
    std::string right_shim;      // read sequence between X and B (empty for deletions)
    std::string middle_shim;     // read sequence between A and B (used for deletions)
    
    ReadEvent(const std::string& read_id = "", RearrangementType type = RearrangementType::LARGE_DELETE,
              const std::string& contig_id = "", const std::string& strand = "+",
              uint32_t out_clip = 0, uint32_t in_clip = 0,
              uint32_t read_clip_out = 0, uint32_t read_clip_in = 0,
              const std::string& element_contig = "", const std::string& element_strand = "",
              uint32_t element_start = 0, uint32_t element_end = 0,
              const std::string& left_shim = "", const std::string& right_shim = "", 
              const std::string& middle_shim = "")
        : read_id(read_id), type(type), contig_id(contig_id), strand(strand),
          out_clip(out_clip), in_clip(in_clip), read_clip_out(read_clip_out), read_clip_in(read_clip_in),
          element_contig(element_contig), element_strand(element_strand),
          element_start(element_start), element_end(element_end),
          left_shim(left_shim), right_shim(right_shim), middle_shim(middle_shim) {}
};

struct AggregatedEvent {
    std::string event_id;
    RearrangementType type;
    std::string contig_id;
    std::string strand;
    uint32_t out_clip;
    uint32_t in_clip;
    std::string element_contig;
    std::string element_strand;
    uint32_t element_start;
    uint32_t element_end;
    uint32_t read_count;
    uint32_t read_coverage;
    std::string left_shim;
    std::string right_shim;
    std::string middle_shim;
    
    AggregatedEvent(const std::string& event_id = "", RearrangementType type = RearrangementType::LARGE_DELETE,
                   const std::string& contig_id = "", const std::string& strand = "+",
                   uint32_t out_clip = 0, uint32_t in_clip = 0,
                   const std::string& element_contig = "", const std::string& element_strand = "",
                   uint32_t element_start = 0, uint32_t element_end = 0,
                   uint32_t read_count = 0, uint32_t read_coverage = 0,
                   const std::string& left_shim = "", const std::string& right_shim = "", 
                   const std::string& middle_shim = "")
        : event_id(event_id), type(type), contig_id(contig_id), strand(strand),
          out_clip(out_clip), in_clip(in_clip),
          element_contig(element_contig), element_strand(element_strand),
          element_start(element_start), element_end(element_end),
          read_count(read_count), read_coverage(read_coverage),
          left_shim(left_shim), right_shim(right_shim), middle_shim(middle_shim) {}
    
    // create aggregation key for grouping identical events
    std::string create_key() const;
};

enum class SupportType {
    A,    // single alignment contains both clip_out and clip_in
    AB,   // two alignments A-B with gap containing clip_out and clip_in
    AXB   // three alignments A-X-B with gap between A and B containing clip_out and clip_in
};

// add operator<< for SupportType
inline std::ostream& operator<<(std::ostream& os, const SupportType& type)
{
    switch (type) {
    case SupportType::A:
        os << "A";
        break;
    case SupportType::AB:
        os << "AB";
        break;
    case SupportType::AXB:
        os << "AXB";
        break;
    default:
        os << "unknown";
        break;
    }
    return os;
}

struct ReadSupport {
    std::string event_id;
    std::string read_id;
    std::string contig;
    uint32_t clip_out;
    uint32_t clip_in;
    SupportType support_type;
    uint32_t distance;  // abs(B.start - A.end), 0 for single alignment
    
    // alignment A (always present)
    uint32_t alignment_A_start;  // contig coordinates
    uint32_t alignment_A_end;
    std::string alignment_A_strand;
    uint32_t read_alignment_A_start;  // read coordinates
    uint32_t read_alignment_A_end;
    
    // alignment B (present for AB and AXB types, "none"/0 for A type)
    uint32_t alignment_B_start;  // contig coordinates  
    uint32_t alignment_B_end;
    std::string alignment_B_strand;
    uint32_t read_alignment_B_start;  // read coordinates
    uint32_t read_alignment_B_end;
    bool contains_event;  // true if this read originally detected the event
    
    ReadSupport(const std::string& event_id = "", const std::string& read_id = "",
               const std::string& contig = "", uint32_t clip_out = 0, uint32_t clip_in = 0,
               SupportType support_type = SupportType::A, uint32_t distance = 0,
               uint32_t alignment_A_start = 0, uint32_t alignment_A_end = 0, 
               const std::string& alignment_A_strand = "+",
               uint32_t read_alignment_A_start = 0, uint32_t read_alignment_A_end = 0,
               uint32_t alignment_B_start = 0, uint32_t alignment_B_end = 0,
               const std::string& alignment_B_strand = "none",
               uint32_t read_alignment_B_start = 0, uint32_t read_alignment_B_end = 0,
               bool contains_event = false)
        : event_id(event_id), read_id(read_id), contig(contig), 
          clip_out(clip_out), clip_in(clip_in), support_type(support_type), distance(distance),
          alignment_A_start(alignment_A_start), alignment_A_end(alignment_A_end), 
          alignment_A_strand(alignment_A_strand),
          read_alignment_A_start(read_alignment_A_start), read_alignment_A_end(read_alignment_A_end),
          alignment_B_start(alignment_B_start), alignment_B_end(alignment_B_end),
          alignment_B_strand(alignment_B_strand),
          read_alignment_B_start(read_alignment_B_start), read_alignment_B_end(read_alignment_B_end),
          contains_event(contains_event) {}
};

class Rearrange {
private:
    const AlignmentStore& store;
    
    // parameters
    int max_gap;
    int min_element_length;
    int min_anchor_length;
    double max_anchor_mutations_percent;
    double max_element_mutation_percent;
    VerifyRearrange* verifier;
    
    // results
    std::vector<ReadEvent> read_events;
    std::vector<ReadSupport> read_supports;
    
    // interval filtering data
    std::set<uint32_t> relevant_read_indices;
    std::set<const Alignment*> relevant_alignments;
    
    // contig-to-events mapping for efficient coverage lookup
    std::map<std::string, std::vector<std::pair<std::string, std::pair<uint32_t, uint32_t>>>> contig_to_events;
    
    // error tracking
    RearrangeErrors error_tracker;
    
    // read sequences for shim extraction
    std::unordered_map<std::string, std::string> read_sequences;   

public:
    Rearrange(const AlignmentStore& store,
       VerifyRearrange* verifier = nullptr,
       int max_gap = 10,
       int min_element_length = 50,
       int min_anchor_length = 200,
       double max_anchor_mutations_percent = 0.01,
       double max_element_mutation_percent = 0.01);

    // public access for RearrangeManager
    std::vector<AggregatedEvent> aggregated_events;
    
    // multi-library support methods (public for RearrangeManager)
    void get_events();
    void get_coverage();
    void aggregate_events(const std::map<std::string, std::string>& key_to_event_id);
    std::set<std::string> get_event_keys() const;
    
    // interval filtering methods
    void set_reads_and_alignments(const std::vector<Interval>& intervals);
    
    // accessor for error tracking (for R interface)
    const RearrangeErrors& get_error_tracker() const;
    
    // read sequence loading methods
    void load_read_sequences_from_file(const std::string& filename);
    void load_read_sequences_from_map(const std::unordered_map<std::string, std::string>& sequences);
    void clear_read_sequences();
        
    void write_to_csv(const std::string& ofn_prefix);
    
private:
    // alignment validation
    bool is_valid_anchor(const Alignment& aln) const;
    bool is_valid_element(const Alignment& aln) const;
    
    // new event testing functions
    EventTestResult test_trio_event(const Alignment& A, const Alignment& B, const Alignment& X, 
                                   bool is_read_reversed, ReadEvent& event);
    EventTestResult test_pair_event(const Alignment& A, const Alignment& B, 
                                   bool is_read_reversed, ReadEvent& event);
    
    // process single read functions
    bool process_read(const std::vector<size_t>& alignment_indices, 
                     const std::vector<Alignment>& alignments,
                     size_t& tested_reads, size_t& found_events);
    bool process_read_trios(const std::vector<size_t>& alignment_indices,
                           const std::vector<Alignment>& alignments,
                           const std::vector<bool>& is_valid_anchor_vec,
                           size_t& tested_reads, size_t& found_events, bool is_verbose = false);
    bool process_read_pairs(const std::vector<size_t>& alignment_indices,
                           const std::vector<Alignment>& alignments,
                           const std::vector<bool>& is_valid_anchor_vec,
                           size_t& tested_reads, size_t& found_events, bool is_verbose = false);
    
    
    // helper functions for event testing
    std::string extract_read_shim(const std::string& read_id, uint32_t start, uint32_t end, bool is_read_reversed) const;
    std::string alignment_to_strand_string(bool is_reverse) const;
    bool alignments_overlap(const Alignment& A, const Alignment& B) const;
    ReadEvent create_event(const Alignment& A, const Alignment& B, const Alignment* X, 
                          RearrangementType type, bool is_read_reversed) const;
    
    // debugging functions
    void print_event_debug(const std::string& context, const Alignment& S1, const Alignment& S2, 
                          const Alignment& S3, const Alignment* A, const Alignment* B, 
                          const Alignment* X, bool is_read_reversed) const;
    void print_pair_debug(const std::string& context, const Alignment& S1, const Alignment& S2,
                         const Alignment* A, const Alignment* B, bool is_read_reversed) const;
    
};

#endif // REARRANGE_H
