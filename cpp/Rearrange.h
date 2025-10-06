#ifndef REARRANGE_H
#define REARRANGE_H

#include "alignment_store.h"
#include "RearrangeReadEvent.h"
#include "RearrangeEventGroup.h"
#include "RearrangeCoverage.h"
#include "RearrangeVerify.h"
#include "Sequences.h"
#include <map>
#include <string>
#include <vector>

// output row for read events table
struct ReadEventRow {
    std::string lib_id;
    std::string read_id;
    std::string event_id;
    RearrangementType type;
    std::string contig_id;
    std::string read_strand;
    uint32_t out_clip;
    uint32_t in_clip;
    uint32_t read_clip_out;
    uint32_t read_clip_in;
    uint32_t span_start;
    uint32_t span_end;
    uint32_t read_span_start;
    uint32_t read_span_end;
    std::string element_contig;
    std::string element_strand;
    uint32_t element_start;
    uint32_t element_end;
    std::string read_seams;      // seams in read coordinates and sequences
    std::string assembly_seams;  // seams in assembly coordinates and sequences
};

// main rearrangement detection class (similar interface to old RearrangeManager)
class Rearrange {
private:
    const std::map<std::string, AlignmentStore>& stores;
    std::vector<std::string> library_ids;
    const std::vector<Interval>& intervals;
    
    // parameters
    int max_margin;
    int min_element_length;
    int min_anchor_length;
    double max_anchor_mutations_percent;
    double max_element_mutation_percent;
    RearrangeVerify* verifier;
    bool write_per_library_files;
    ResolveSeams resolve_seams;
    
    // sequence objects
    const AssemblySequences* assembly_sequences;
    const ReadSequences* read_sequences;
    
    // output prefix for caching
    std::string output_prefix;
    
    // event grouper (reused across workflow)
    RearrangeEventGroup event_grouper;
    
    // map from lib_id to events
    std::map<std::string, std::vector<ReadEvent>> lib_events;

    // representative events (only container needed)
    std::vector<Event> rep_events;

    // groups of indices into all_events
    std::vector<std::vector<int>> event_groups;
    
    // results
    std::vector<ReadEventRow> read_event_rows;
    std::map<std::string, std::map<std::string, int>> support_matrix; // lib_id -> event_id -> count
    std::map<std::string, std::map<std::string, int>> coverage_matrix; // lib_id -> event_id -> count

public:
    Rearrange(const std::map<std::string, AlignmentStore>& stores,
              const std::vector<Interval>& intervals,
              RearrangeVerify* verifier = nullptr,
              ResolveSeams resolve_seams = ResolveSeams::NO,
              const AssemblySequences* assembly_seqs = nullptr,
              const ReadSequences* read_seqs = nullptr,
              int max_margin = 10,
              int min_element_length = 50,
              int min_anchor_length = 200,
              double max_anchor_mutations_percent = 0.01,
              double max_element_mutation_percent = 0.01,
              bool write_per_library_files = true,
              const std::string& output_prefix = "");
    
    // main execution function
    void execute();
        
    // output functions
    void write_to_csv(const std::string& ofn_prefix);
    
    // getters for R interface
    const std::vector<ReadEventRow>& get_read_event_rows() const { return read_event_rows; }
    const std::vector<Event>& get_events() const { return rep_events; }
    std::vector<std::string> get_library_ids() const { return library_ids; }
    const std::map<std::string, std::map<std::string, int>>& get_support_matrix() const { return support_matrix; }
    const std::map<std::string, std::map<std::string, int>>& get_coverage_matrix() const { return coverage_matrix; }
    
    // rejection statistics
    std::map<std::string, size_t> get_rejection_counts() const;

private:
    // workflow steps
    void detect_events_per_library();
    void group_and_select_representatives();
    void calculate_coverage_matrices();
    void prepare_output_tables();
    void calculate_event_statistics();
    
    // output file writers
    void write_read_events_file(const std::string& ofn_prefix);
    void write_events_file(const std::string& ofn_prefix);
    void write_support_matrix_file(const std::string& ofn_prefix);
    void write_coverage_matrix_file(const std::string& ofn_prefix);
};

#endif // REARRANGE_H