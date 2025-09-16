#ifndef REARRANGE_MANAGER_H
#define REARRANGE_MANAGER_H

#include "Rearrange.h"
#include "alignment_store.h"
#include "aln_types.h"
#include <map>
#include <memory>
#include <set>
#include <string>
#include <vector>

// Forward declaration
class VerifyRearrange;

// Cross-library data structures
struct SetRearrangementData {
    RearrangementType type;
    std::string contig_id;
    std::string strand;
    uint32_t out_clip;
    uint32_t in_clip;
    std::string element_contig;
    std::string element_strand;
    uint32_t element_start;
    uint32_t element_end;
    
    int total_support = 0;    // total reads across all libraries
    int total_coverage = 0;   // total coverage across all libraries  
    int library_count = 0;    // number of libraries with this event
    std::map<std::string, int> lib_support;   // lib_id -> read count
    std::map<std::string, int> lib_coverage;  // lib_id -> coverage count
    
    std::string create_key() const;
    double get_frequency() const {
        return (total_coverage > 0) ? (static_cast<double>(total_support) / total_coverage) : 0.0;
    }
};

struct RearrangementOutputRow {
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
    int library_count;
    int total_support;
    int total_coverage;
    double frequency;
    
    std::string create_key() const;
};

class RearrangeManager {
private:
    std::map<std::string, AlignmentStore> stores;
    std::map<std::string, std::unique_ptr<Rearrange>> rearrangers;
    std::vector<std::string> library_ids;
    std::vector<Interval> intervals;
    
    // parameters
    int max_gap;
    int min_element_length;
    int min_anchor_length;
    double max_mutations_percent;
    VerifyRearrange* verifier;
    
    // cross-library results
    std::map<std::string, SetRearrangementData> rearrangement_data;
    std::vector<RearrangementOutputRow> rearrangement_rows;
    
    // event key management
    std::set<std::string> all_event_keys;
    std::map<std::string, std::string> key_to_event_id;
    
public:
    RearrangeManager(const std::map<std::string, AlignmentStore>& stores,
                     const std::vector<Interval>& intervals,
                     VerifyRearrange* verifier = nullptr,
                     int max_gap = 10,
                     int min_element_length = 50,
                     int min_anchor_length = 200,
                     double max_mutations_percent = 0.01);
    
    void execute();
    void write_to_csv(const std::string& ofn_prefix);
    
    // getters for R interface
    const std::vector<RearrangementOutputRow>& get_rearrangement_rows() const { return rearrangement_rows; }
    const std::map<std::string, SetRearrangementData>& get_rearrangement_data() const { return rearrangement_data; }
    std::vector<std::string> get_library_ids() const { return library_ids; }
    
private:
    // two-pass execution flow
    void collect_all_event_keys();
    void sort_and_assign_event_ids();
    void process_libraries_with_event_ids();
    
    // cross-library aggregation
    void aggregate_across_libraries();
    void calculate_cross_library_coverage();
    
    // matrix output methods
    void write_events_file(const std::string& ofn_prefix);
    void write_support_matrix(const std::string& ofn_prefix);
    void write_coverage_matrix(const std::string& ofn_prefix);
};

#endif // REARRANGE_MANAGER_H
