#ifndef QUERYVARIANTS_H
#define QUERYVARIANTS_H

#include "alignment_store.h" // Includes aln_types.h indirectly
#include "utils.h"
#include <map>
#include <string>
#include <vector>

// Structure to represent a variant across multiple libraries
struct SetVariantData {
  std::string contig;
  uint32_t position; // 1-based coordinate
  std::string type; // "sub", "ins", "del", "left_clip", "right_clip"
  std::string sequence; // for mutations: nts, for clips: empty
  int total_support = 0; // total reads supporting this variant across all libs
  int total_coverage = 0; // total coverage at this position across all libs
  int library_count = 0; // number of libraries with this variant
  std::map<std::string, int> lib_support; // library_id -> read support
  std::map<std::string, int> lib_coverage; // library_id -> coverage
  
  // create unique key for variant
  std::string create_key() const;
  
  // calculate frequency
  double get_frequency() const {
    return (total_coverage > 0) ? (static_cast<double>(total_support) / total_coverage) : 0.0;
  }
};

// Structure representing a single row in the variants output table
struct VariantOutputRow {
  std::string variant_id;
  std::string contig;
  uint32_t coord; // 1-based
  std::string type;
  std::string sequence;
  std::string desc; // human-readable description (e.g., "A:T", "+ATG", "-CC")
  int library_count;
  int total_support;
  int total_coverage;
  double frequency;
};

class QueryVariants {
  private:
  const std::vector<Interval>& intervals;
  std::map<std::string, AlignmentStore> stores; // library_id -> store
  
  // QueryVariants-specific parameters
  int min_variants_variant_support;
  int min_variants_library_support;
  int min_variants_coverage_support;
  
  // Alignment filtering parameters (same as other queries)
  ClipMode clip_mode;
  int clip_margin;
  double min_mutations_percent;
  double max_mutations_percent;
  int min_alignment_length;
  int max_alignment_length;
  
  // Results storage
  std::map<std::string, SetVariantData> variant_data; // variant_key -> data
  std::vector<VariantOutputRow> variant_rows;
  std::vector<std::string> library_ids; // ordered list of library IDs
  
  // efficiency optimization - map from contig_index to intervals for efficient lookup
  std::map<uint32_t, std::vector<const Interval*>> contig_to_intervals_map;
  
  // process a single alignment from a library
  void process_alignment(const Alignment& aln, const AlignmentStore& store, 
                        const std::string& lib_id);
  
  // process mutations from an alignment
  void process_mutations(const Alignment& aln, const AlignmentStore& store, 
                        const std::string& lib_id);
  
  // process clips from an alignment
  void process_clips(const Alignment& aln, const AlignmentStore& store, 
                    const std::string& lib_id);
  
  // helper function to process a single clip
  void process_single_clip(uint32_t clip_pos, const std::string& clip_type, 
                          const AlignmentStore& store, const Alignment& aln, 
                          const std::string& lib_id);
  
  // calculate coverage for all variant positions
  void calculate_coverage();
  
  // build contig to intervals mapping for efficient lookup
  void build_contig_to_intervals_map(const AlignmentStore& store);
  
  // apply QueryVariants-specific filters and generate output rows
  void apply_filters_and_generate_rows();
  
  // write individual output files
  void write_variants_file(const std::string& ofn_prefix);
  void write_support_file(const std::string& ofn_prefix);
  void write_coverage_file(const std::string& ofn_prefix);

  public:
  QueryVariants(const std::vector<Interval>& intervals,
           const std::map<std::string, AlignmentStore>& stores,
           int min_variants_variant_support = 3,
           int min_variants_library_support = 1,
           int min_variants_coverage_support = 10,
           ClipMode clip_mode = ClipMode::ALL,
           int clip_margin = 10,
           double min_mutations_percent = 0.0,
           double max_mutations_percent = 10.0,
           int min_alignment_length = 0,
           int max_alignment_length = 0);

  // execute the query
  void execute();

  // write all output files
  void write_to_csv(const std::string& ofn_prefix);

  // get access to output rows (for R interface)
  const std::vector<VariantOutputRow>& get_variant_rows() const {
    return variant_rows;
  }
  
  const std::vector<std::string>& get_library_ids() const {
    return library_ids;
  }
  
  const std::map<std::string, SetVariantData>& get_variant_data() const {
    return variant_data;
  }
};

#endif // QUERYVARIANTS_H
