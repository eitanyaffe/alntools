#ifndef SEGMENTATION_H
#define SEGMENTATION_H

#include "alignment_store.h"
#include "utils.h"
#include <string>
#include <vector>
#include <map>

// structure to represent a read breakpoint
struct ReadBreakpoint {
    std::string lib_id;
    std::string read_id;
    std::string contig_id;
    uint32_t coord;
    std::string type;  // "dangle_left" or "dangle_right"
    uint32_t anchor_length;
    uint32_t anchor_mutations;
    uint32_t dangle_length;
    
    ReadBreakpoint(const std::string& lib_id = "", const std::string& read_id = "", 
                   const std::string& contig_id = "",
                   uint32_t coord = 0, const std::string& type = "",
                   uint32_t anchor_length = 0, uint32_t anchor_mutations = 0,
                   uint32_t dangle_length = 0)
        : lib_id(lib_id), read_id(read_id), contig_id(contig_id), coord(coord), type(type),
          anchor_length(anchor_length), anchor_mutations(anchor_mutations),
          dangle_length(dangle_length) {}
};

// worker class for detecting breakpoints in a single library
class Segmentation {
private:
    const AlignmentStore& store;
    
    // parameters
    int min_anchor_length;
    int min_dangle_length;
    double max_anchor_mutations_percent;
    int min_alignment_distance;
    
    std::string current_lib_id;

public:
    Segmentation(const AlignmentStore& store,
                int max_margin = 20,
                int min_anchor_length = 1000,
                int min_dangle_length = 1000,
                double max_anchor_mutations_percent = 0.001,
                int min_alignment_distance = 200);
    
    // main detection function
    void detect_breakpoints(const std::string& lib_id, 
                          std::vector<ReadBreakpoint>& output_breakpoints);

private:
    // check if alignment is valid anchor
    bool is_valid_anchor(const Alignment& aln) const;
    
    // check if alignment is valid dangle
    bool is_valid_dangle(const Alignment& dangle, const Alignment& anchor) const;
    
    // determine breakpoint type based on strand and end position
    std::string get_breakpoint_type(bool is_reverse, bool is_end) const;
    
    // get breakpoint coordinate from alignment
    uint32_t get_breakpoint_coord(const Alignment& aln, bool is_end) const;
    
    // process a single read
    void process_read(uint32_t read_idx, 
                     const std::vector<size_t>& alignment_indices,
                     const std::vector<Alignment>& alignments,
                     std::vector<ReadBreakpoint>& output_breakpoints);
    
    // check if two alignments are on different contigs or far apart
    bool is_valid_dangle_position(const Alignment& anchor, const Alignment& dangle) const;
};

#endif // SEGMENTATION_H

