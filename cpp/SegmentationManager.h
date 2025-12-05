#ifndef SEGMENTATION_MANAGER_H
#define SEGMENTATION_MANAGER_H

#include "alignment_store.h"
#include "Segmentation.h"
#include <map>
#include <string>
#include <vector>

// structure to represent an aggregated breakpoint
struct AggregateBreakpoint {
    std::string breakpoint_id;
    std::string contig_id;
    uint32_t coord;
    std::string type;
    int read_support;
    double frequency;
    bool selected;
    bool is_segment_break;
    std::map<std::string, int> lib_support;
    std::map<std::string, int> lib_coverage;
    
    AggregateBreakpoint(const std::string& id = "", const std::string& contig = "",
                       uint32_t coord = 0, const std::string& type = "")
        : breakpoint_id(id), contig_id(contig), coord(coord), type(type),
          read_support(0), frequency(0.0), selected(false), is_segment_break(false) {}
};

// structure to represent a segment
struct Segment {
    std::string segment_id;
    std::string contig_id;
    uint32_t start;
    uint32_t end;
    uint32_t length;
    std::string start_breakpoint_id;
    std::string end_breakpoint_id;
    
    Segment(const std::string& id = "", const std::string& contig = "",
           uint32_t start = 0, uint32_t end = 0,
           const std::string& start_bp_id = "", const std::string& end_bp_id = "")
        : segment_id(id), contig_id(contig), start(start), end(end),
          length(end - start), start_breakpoint_id(start_bp_id), end_breakpoint_id(end_bp_id) {}
};

// main segmentation manager class
class SegmentationManager {
private:
    const std::map<std::string, AlignmentStore>& stores;
    std::vector<std::string> library_ids;
    
    // contig table
    std::string ifn_contig_table;
    std::map<std::string, uint32_t> contig_lengths;
    
    // parameters
    int max_margin;
    int min_anchor_length;
    int min_dangle_length;
    double max_anchor_mutations_percent;
    int min_alignment_distance;
    int min_breakpoint_read_support;
    double min_breakpoint_frequency;
    int min_segment_length;
    int max_segment_length;
    
    // results
    std::map<std::string, std::vector<ReadBreakpoint>> lib_breakpoints;
    std::vector<AggregateBreakpoint> aggregate_breakpoints;
    std::vector<Segment> segments;

public:
    SegmentationManager(const std::map<std::string, AlignmentStore>& stores,
                       const std::string& ifn_contig_table,
                       int max_margin = 20,
                       int min_anchor_length = 1000,
                       int min_dangle_length = 1000,
                       double max_anchor_mutations_percent = 0.001,
                       int min_alignment_distance = 200,
                       int min_breakpoint_read_support = 2,
                       double min_breakpoint_frequency = 0.2,
                       int min_segment_length = 200,
                       int max_segment_length = 500000);
    
    // main execution function
    void execute();
    
    // output functions
    void write_to_csv(const std::string& odir);
    
    // getters for results
    const std::map<std::string, std::vector<ReadBreakpoint>>& get_lib_breakpoints() const { 
        return lib_breakpoints; 
    }
    const std::vector<AggregateBreakpoint>& get_aggregate_breakpoints() const { 
        return aggregate_breakpoints; 
    }
    const std::vector<Segment>& get_segments() const { 
        return segments; 
    }

private:
    // workflow steps
    void load_contig_table();
    void detect_breakpoints_per_library();
    void aggregate_breakpoints_step();
    void calculate_coverage();
    void filter_breakpoints();
    void generate_segments();
    
    // helper functions
    void cluster_breakpoints(std::vector<ReadBreakpoint>& all_breakpoints);
    int calculate_breakpoint_coverage(const AggregateBreakpoint& bp, const std::string& lib_id) const;
    
    // segment generation helpers
    void add_artificial_coords(std::vector<uint32_t>& coords, uint32_t contig_length);
    void create_segments_from_coords(const std::vector<uint32_t>& coords, 
                                     const std::map<uint32_t, std::string>& coord_to_bp_id,
                                     const std::string& contig_id,
                                     uint32_t contig_length,
                                     int& segment_counter);
    
    // output file writers
    void write_read_breakpoints_file(const std::string& odir);
    void write_aggregate_breakpoints_file(const std::string& odir);
    void write_support_matrix_file(const std::string& odir);
    void write_coverage_matrix_file(const std::string& odir);
    void write_segments_file(const std::string& odir);
};

#endif // SEGMENTATION_MANAGER_H

