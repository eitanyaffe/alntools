#pragma once

#include "alignment_store.h"
#include <string>
#include <vector>
#include <map>
#include <tuple>
#include <unordered_set>
#include <unordered_map>
#include <cstdint>

struct SegMatrixSegment {
    std::string id;
    std::string contig;
    uint32_t start;  // 1-based inclusive
    uint32_t end;    // 1-based inclusive
    uint32_t length;
    size_t index;
    uint32_t left_side_start;  // 0-based inclusive
    uint32_t left_side_end;    // 0-based exclusive
    uint32_t right_side_start; // 0-based inclusive
    uint32_t right_side_end;   // 0-based exclusive

    SegMatrixSegment()
        : start(0)
        , end(0)
        , length(0)
        , index(0)
        , left_side_start(0)
        , left_side_end(0)
        , right_side_start(0)
        , right_side_end(0)
    {
    }
};

struct ReadInterval {
    std::string segment_id;
    std::string contig;
    uint32_t contig_start;  // 0-based half-open
    uint32_t contig_end;    // 0-based half-open
    // estimated read-space coordinates derived from alignment
    uint32_t read_start_est;
    uint32_t read_end_est;
    bool is_reverse;
    bool covers_left_side;
    bool covers_right_side;
    bool is_short_segment;
    double left_mutation_percent;
    double right_mutation_percent;
    size_t segment_index;   // index in segments vector

    ReadInterval()
        : contig_start(0), contig_end(0),
          read_start_est(0), read_end_est(0),
          is_reverse(false), covers_left_side(false), covers_right_side(false),
          is_short_segment(false), left_mutation_percent(0.0), right_mutation_percent(0.0),
          segment_index(0) {}
};

class SegMatrix {
private:
    std::vector<SegMatrixSegment> segments;
    AlignmentStore& store;
    std::unordered_map<std::string, size_t> segment_index_map;

    // parameters
    double max_mutation_percent;
    uint32_t max_adjacency_distance;
    uint32_t max_margin;
    int min_indel_length;
    uint32_t side_length;
    uint32_t side_margin;
    uint32_t min_segment_length;

    std::map<std::tuple<std::string, std::string, std::string, std::string>, uint32_t> adjacency_matrix;
    std::map<std::tuple<std::string, std::string, std::string, std::string>, uint32_t> reach_matrix;

    struct SegmentStats {
        uint64_t total_read_count_left;
        uint64_t total_read_count_right;
        uint64_t associated_read_count_left;
        uint64_t associated_read_count_right;
        std::unordered_set<std::string> associated_segments_left;
        std::unordered_set<std::string> associated_segments_right;

        SegmentStats()
            : total_read_count_left(0)
            , total_read_count_right(0)
            , associated_read_count_left(0)
            , associated_read_count_right(0)
        {
        }
    };
    std::map<std::string, SegmentStats> segment_stats;

    struct SegmentReadContribution {
        bool total_left;
        bool total_right;
        bool associated_left;
        bool associated_right;

        SegmentReadContribution()
            : total_left(false)
            , total_right(false)
            , associated_left(false)
            , associated_right(false)
        {
        }
    };

    std::map<std::string, std::vector<size_t>> contig_to_segments;

    uint32_t debug_read_count;

    void load_segments(const std::string& ifn_segments);
    void load_library();
    void build_contig_to_segments_map();

    void process_read(uint32_t read_index);
    std::vector<const Alignment*> get_sorted_alignments_for_read(uint32_t read_index);
    std::vector<const Alignment*> filter_overlapping_alignments(const std::vector<const Alignment*>& alignments, bool debug_mode);

    std::vector<ReadInterval> intersect_alignment_with_segments(const Alignment& aln);
    void compute_side_coverage_and_mutations(ReadInterval& interval, const Alignment& aln);
    double compute_mutation_percent(const Alignment& aln, uint32_t contig_start, uint32_t contig_end);

    void process_interval_sequence(const std::vector<ReadInterval>& intervals, bool debug_mode);
    void create_association(const ReadInterval& exit_interval,
                            const ReadInterval& entry_interval,
                            bool debug_mode,
                            std::unordered_map<std::string, SegmentReadContribution>& read_contrib);
    bool process_interval_side(const ReadInterval& interval,
                               const std::string& direction,
                               bool debug_mode,
                               std::unordered_map<std::string, SegmentReadContribution>& read_contrib,
                               std::string& side_out,
                               double& mutation_percent_out);
    bool should_skip_mutation_threshold(double mutation_percent, const ReadInterval& interval,
                                        const std::string& direction, bool debug_mode) const;
    bool is_pointing_out(const ReadInterval& interval) const;
    bool is_pointing_in(const ReadInterval& interval) const;
    bool is_within_adjacency_distance(const ReadInterval& interval1, const ReadInterval& interval2) const;
    bool read_covers_seam(const std::vector<const Alignment*>& alignments) const;

    void print_read_debug_info(uint32_t read_index, const std::string& read_id,
                               const std::vector<const Alignment*>& alignments_before,
                               const std::vector<const Alignment*>& alignments_after,
                               const std::vector<ReadInterval>& intervals) const;
    void print_read_header(uint32_t read_index, const std::string& read_id) const;
    void print_interval(const ReadInterval& interval, size_t idx) const;
    void print_interval_sequence(const std::vector<ReadInterval>& intervals) const;
    void print_exit_detected(const ReadInterval& interval, double mutation_percent) const;
    void print_entry_detected(const ReadInterval& interval, double mutation_percent) const;
    void print_association(const ReadInterval& exit_interval, const ReadInterval& entry_interval, bool is_adjacent) const;
    void print_short_segment_skipped(const ReadInterval& interval) const;

    void write_matrices(const std::string& odir);
    void write_segment_summary(const std::string& odir);
    void write_matrix(const std::map<std::tuple<std::string, std::string, std::string, std::string>, uint32_t>& matrix,
                      const std::string& filename,
                      const std::string& matrix_name) const;

public:
    SegMatrix(AlignmentStore& store_ref);

    void compute(const std::string& ifn_segments,
                 const std::string& odir,
                 double max_mutation_percent,
                 uint32_t max_adjacency_distance,
                 uint32_t max_margin,
                 int min_indel_length,
                 uint32_t side_length,
                 uint32_t side_margin);
};
