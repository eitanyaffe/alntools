#ifndef REARRANGE_READ_EVENT_H
#define REARRANGE_READ_EVENT_H

#include "alignment_store.h"
#include "utils.h"
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
    REJECT_READ_GAP_TOO_LARGE
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
    std::string element_contig;   // contig of element X (empty for deletions)
    std::string element_strand;  // "+" or "-" for element X (empty for deletions)
    uint32_t element_start;      // contig coordinates (0 for deletions)
    uint32_t element_end;        // contig coordinates (0 for deletions)
    std::string left_shim;       // read sequence between A and X (empty for deletions)
    std::string right_shim;      // read sequence between X and B (empty for deletions)
    std::string middle_shim;     // read sequence between A and B (used for deletions)
    
    ReadEvent(const std::string& read_id = "", RearrangementType type = RearrangementType::LARGE_DELETE,
              const std::string& contig_id = "", const std::string& read_strand = "+",
              uint32_t out_clip = 0, uint32_t in_clip = 0,
              uint32_t read_clip_out = 0, uint32_t read_clip_in = 0,
              const std::string& element_contig = "", const std::string& element_strand = "",
              uint32_t element_start = 0, uint32_t element_end = 0,
              const std::string& left_shim = "", const std::string& right_shim = "", 
              const std::string& middle_shim = "")
        : read_id(read_id), type(type), contig_id(contig_id), read_strand(read_strand),
          out_clip(out_clip), in_clip(in_clip), read_clip_out(read_clip_out), read_clip_in(read_clip_in),
          element_contig(element_contig), element_strand(element_strand),
          element_start(element_start), element_end(element_end),
          left_shim(left_shim), right_shim(right_shim), middle_shim(middle_shim) {}
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
    
    // read sequences for shim extraction
    std::unordered_map<std::string, std::string> read_sequences;
    
    // interval filtering data
    std::vector<Interval> intervals;
    std::set<uint32_t> relevant_read_indices;
    std::set<const Alignment*> relevant_alignments;

public:
    RearrangeReadEvent(const AlignmentStore& store,
                      RearrangeVerify* verifier = nullptr,
                      int max_margin = 10,
                      int min_element_length = 50,
                      int min_anchor_length = 200,
                      double max_anchor_mutations_percent = 0.01,
                      double max_element_mutation_percent = 0.01);
    
    // main detection function - pass output containers by reference
    void detect_events(std::vector<ReadEvent>& output_events, 
                      std::map<std::string, size_t>& output_rejection_counts);
    
    // interval filtering
    void set_intervals(const std::vector<Interval>& intervals);
    
    // read sequence management
    void load_read_sequences_from_file(const std::string& filename);
    void load_read_sequences_from_map(const std::unordered_map<std::string, std::string>& sequences);
    void clear_read_sequences();
    
    // utility methods
    size_t get_total_rejected_count(const std::map<std::string, size_t>& rejection_counts) const;
    
    // debug methods
    void print_event(const ReadEvent& event) const;

private:
    // alignment validation
    bool is_valid_anchor(const Alignment& aln) const;
    bool is_valid_element(const Alignment& aln) const;
    
    // event testing functions
    EventTestResult test_trio_event(const Alignment& A, const Alignment& B, const Alignment& X, 
                                   bool is_read_reversed, ReadEvent& event,
                                   std::map<std::string, size_t>& rejection_counts);
    EventTestResult test_pair_event(const Alignment& A, const Alignment& B, 
                                   bool is_read_reversed, ReadEvent& event,
                                   std::map<std::string, size_t>& rejection_counts);
    
    // rejection tracking
    void add_rejection(EventTestResult result, std::map<std::string, size_t>& rejection_counts);
    
    // process single read functions
    bool process_read(const std::vector<size_t>& alignment_indices, 
                     const std::vector<Alignment>& alignments,
                     std::vector<ReadEvent>& output_events,
                     std::map<std::string, size_t>& rejection_counts);
    bool process_read_trios(const std::vector<size_t>& alignment_indices,
                           const std::vector<Alignment>& alignments,
                           const std::vector<bool>& is_valid_anchor_vec,
                           std::vector<ReadEvent>& output_events,
                           std::map<std::string, size_t>& rejection_counts);
    bool process_read_pairs(const std::vector<size_t>& alignment_indices,
                           const std::vector<Alignment>& alignments,
                           const std::vector<bool>& is_valid_anchor_vec,
                           std::vector<ReadEvent>& output_events,
                           std::map<std::string, size_t>& rejection_counts);
    
    // helper functions
    std::string extract_read_shim(const std::string& read_id, uint32_t start, uint32_t end, bool is_read_reversed) const;
    std::string alignment_to_strand_string(bool is_reverse) const;
    bool alignments_overlap(const Alignment& A, const Alignment& B) const;
    ReadEvent create_event(const Alignment& A, const Alignment& B, const Alignment* X, 
                          RearrangementType type, bool is_read_reversed) const;
    
    // rejection tracking
    void add_rejection(const std::string& reason, std::map<std::string, size_t>& rejection_counts);
};

#endif // REARRANGE_READ_EVENT_H
