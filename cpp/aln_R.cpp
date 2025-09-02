// [[Rcpp::plugins(cpp17)]]

// OpenMP support - configured via environment variables in R startup
#ifdef _OPENMP
#include <omp.h>
#endif

#include "QueryBin.h"
#include "QueryConsensus.h"
#include "QueryFull.h"
#include "QueryPileup.h"
#include "alignment_store.h"
#include "paf_reader.h"
#include <Rcpp.h>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_set>
#include <vector>

using namespace std;
using namespace Rcpp;

////////////////////////////////////////////////////////////////////////////////
// Helper functions
////////////////////////////////////////////////////////////////////////////////

void vet_intervals(DataFrame df)
{
  // Check for each required column individually to provide specific error message
  if (!df.containsElementNamed("contig")) {
    stop("Intervals dataframe is missing required column: contig");
  }
  if (!df.containsElementNamed("start")) {
    stop("Intervals dataframe is missing required column: start");
  }
  if (!df.containsElementNamed("end")) {
    stop("Intervals dataframe is missing required column: end");
  }
  if (df.nrows() == 0) {
    stop("Intervals dataframe must contain at least one row");
  }
}

// Helper function to convert R DataFrame to C++ std::vector<Interval>
std::vector<Interval> Rcpp_DataFrame_to_Intervals(DataFrame df)
{
  // verify the dataframe has the required columns
  // Check if the dataframe has all required columns
  vet_intervals(df);

  CharacterVector contig = df["contig"];
  IntegerVector start = df["start"]; // Expecting 1-based
  IntegerVector end = df["end"]; // Expecting 1-based closed

  int n = df.nrows();
  std::vector<Interval> intervals;
  intervals.reserve(n);

  for (int i = 0; i < n; ++i) {
    // Convert R 1-based closed [start, end] to C++ 0-based half-open [start-1, end)
    if (start[i] <= 0) {
      stop("Interval start coordinates must be positive (1-based). Found %d at row %d", start[i], i + 1);
    }
    // Allow start == end for zero-length interval representation?
    if (end[i] < start[i]) {
      stop("Interval end coordinate must be >= start coordinate. Found start=%d, end=%d at row %d", start[i], end[i], i + 1);
    }
    intervals.emplace_back(as<std::string>(contig[i]), start[i] - 1, end[i]);
  }
  return intervals;
}

////////////////////////////////////////////////////////////////////////////////
// QueryByReadIds function
////////////////////////////////////////////////////////////////////////////////

// [[Rcpp::export]]
DataFrame aln_alignments_from_read_ids(
    XPtr<AlignmentStore> store_ptr,
    CharacterVector read_ids)
{
  // Validate the external pointer
  if (!store_ptr) {
    stop("Invalid AlignmentStore pointer provided.");
  }

  // Get reference to the AlignmentStore object (non-const since get_read_index isn't const)
  AlignmentStore& store = *store_ptr;

  // Get alignments
  const std::vector<Alignment>& all_alignments = store.get_alignments();

  // Collect read indices
  std::unordered_set<uint32_t> target_read_indices;
  for (int i = 0; i < read_ids.length(); ++i) {
    std::string read_id = as<std::string>(read_ids[i]);
    if (!store.has_read_id(read_id)) {
      cout << "Read ID '" << read_id << "' not found in alignment store" << endl;
      continue;
    }
    size_t read_index = store.get_read_index(read_id);
    target_read_indices.insert(read_index);
  }

  // Output vectors for DataFrame
  NumericVector out_aln_idx;
  CharacterVector out_aln_read_id;
  IntegerVector out_aln_read_length;
  CharacterVector out_aln_contig_id;
  IntegerVector out_aln_read_start;
  IntegerVector out_aln_read_end;
  IntegerVector out_aln_contig_start;
  IntegerVector out_aln_contig_end;
  LogicalVector out_aln_is_reverse;
  IntegerVector out_aln_num_mutations;

  // Collect alignments for target reads
  for (size_t i = 0; i < all_alignments.size(); ++i) {
    const Alignment& aln = all_alignments[i];

    if (target_read_indices.find(aln.read_index) != target_read_indices.end()) {
      // This alignment belongs to one of our target reads
      out_aln_idx.push_back(i);
      out_aln_read_id.push_back(store.get_read_id(aln.read_index));
      out_aln_read_length.push_back(store.get_reads()[aln.read_index].length);
      out_aln_contig_id.push_back(store.get_contig_id(aln.contig_index));
      out_aln_read_start.push_back(aln.read_start);
      out_aln_read_end.push_back(aln.read_end);
      out_aln_contig_start.push_back(aln.contig_start);
      out_aln_contig_end.push_back(aln.contig_end);
      out_aln_is_reverse.push_back(aln.is_reverse);
      out_aln_num_mutations.push_back(aln.get_mutation_count());
    }
  }

  return DataFrame::create(
      Named("alignment_index") = out_aln_idx,
      Named("read_id") = out_aln_read_id,
      Named("read_length") = out_aln_read_length,
      Named("contig_id") = out_aln_contig_id,
      Named("read_start") = out_aln_read_start,
      Named("read_end") = out_aln_read_end,
      Named("contig_start") = out_aln_contig_start,
      Named("contig_end") = out_aln_contig_end,
      Named("is_reverse") = out_aln_is_reverse,
      Named("num_mutations") = out_aln_num_mutations,
      Named("stringsAsFactors") = false);
}

////////////////////////////////////////////////////////////////////////////////
// Load AlignmentStore from file
////////////////////////////////////////////////////////////////////////////////

// [[Rcpp::export]]
XPtr<AlignmentStore> aln_load(std::string filepath)
{
  // Create a new AlignmentStore instance on the heap
  AlignmentStore* store = new AlignmentStore();

  try {
    // Attempt to load the data
    Rcout << "Loading AlignmentStore from: " << filepath << std::endl;
    store->load(filepath);

    // Create an external pointer managed by R's garbage collector
    XPtr<AlignmentStore> ptr(store, true);
    return ptr;
  } catch (const std::runtime_error& e) {
    // If load fails (e.g., file not found, assertion fails),
    // clean up the allocated store and signal an error to R.
    delete store;
    stop("Failed to load AlignmentStore: %s", e.what());
  } catch (...) {
    // Catch any other potential C++ exceptions
    delete store;
    stop("An unknown C++ error occurred during loading.");
  }
}

////////////////////////////////////////////////////////////////////////////////
// QueryBin functions
////////////////////////////////////////////////////////////////////////////////

// [[Rcpp::export]]
DataFrame aln_query_bin(
    XPtr<AlignmentStore> store_ptr,
    DataFrame intervals_df,
    int binsize,
    double seg_threshold = 0.2,
    double non_ref_threshold = 0.9,
    int num_threads = 0,
    std::string clip_mode_str = "all",
    int clip_margin = 10,
    double min_mutations_percent = 0.0,
    double max_mutations_percent = 10.0,
    int min_alignment_length = 0,
    int max_alignment_length = 0,
    int min_indel_length = 3)
{
  // Validate the external pointer
  if (!store_ptr) {
    stop("Invalid AlignmentStore pointer provided.");
  }

  // Get reference to the AlignmentStore object
  AlignmentStore& store = *store_ptr;

  // count short indels before query
  store.count_short_indels(min_indel_length);

  // Convert intervals
  std::vector<Interval> intervals = Rcpp_DataFrame_to_Intervals(intervals_df);

  // Convert clip_mode_str to ClipMode enum
  ClipMode clip_mode = string_to_clip_mode(clip_mode_str);
  
  std::vector<BinOutputRow> results;
  
  QueryBin queryBin(intervals, store, binsize, seg_threshold, non_ref_threshold, num_threads, clip_mode, clip_margin, min_mutations_percent, max_mutations_percent, min_alignment_length, max_alignment_length);
  queryBin.execute();
  results = queryBin.get_output_rows();

  // Convert results to R DataFrame
  CharacterVector out_contig;
  IntegerVector out_bin_start;
  IntegerVector out_bin_end;
  IntegerVector out_bin_length;
  IntegerVector out_sequenced_bp;
  IntegerVector out_read_count;
  IntegerVector out_mutation_count;
  NumericVector out_median_mutation_density;
  NumericVector out_seg_sites_density;
  NumericVector out_non_ref_sites_density;
  NumericVector out_seg_clip_density;
  NumericVector out_non_ref_clip_density;
  IntegerVector out_dist_none;
  IntegerVector out_dist_5;
  IntegerVector out_dist_4;
  IntegerVector out_dist_3;
  IntegerVector out_dist_2;
  IntegerVector out_dist_1_plus;

  for (const auto& row : results) {
    out_contig.push_back(row.contig);
    out_bin_start.push_back(row.bin_start);
    out_bin_end.push_back(row.bin_end);
    out_bin_length.push_back(row.bin_length);
    out_sequenced_bp.push_back(row.sequenced_basepairs);
    out_read_count.push_back(row.read_count);
    out_mutation_count.push_back(row.mutation_count);
    out_median_mutation_density.push_back(row.median_mutation_density);
    out_seg_sites_density.push_back(row.seg_sites_density);
    out_non_ref_sites_density.push_back(row.non_ref_sites_density);
    out_seg_clip_density.push_back(row.seg_clip_density);
    out_non_ref_clip_density.push_back(row.non_ref_clip_density);
    out_dist_none.push_back(row.dist_none);
    out_dist_5.push_back(row.dist_5);
    out_dist_4.push_back(row.dist_4);
    out_dist_3.push_back(row.dist_3);
    out_dist_2.push_back(row.dist_2);
    out_dist_1_plus.push_back(row.dist_1_plus);
  }

  return DataFrame::create(
      Named("contig") = out_contig,
      Named("start") = out_bin_start,
      Named("end") = out_bin_end,
      Named("length") = out_bin_length,
      Named("sequenced_bp") = out_sequenced_bp,
      Named("read_count") = out_read_count,
      Named("mutation_count") = out_mutation_count,
      Named("median_mutation_density") = out_median_mutation_density,
      Named("seg_sites_density") = out_seg_sites_density,
      Named("non_ref_sites_density") = out_non_ref_sites_density,
      Named("seg_clip_density") = out_seg_clip_density,
      Named("non_ref_clip_density") = out_non_ref_clip_density,
      Named("dist_none") = out_dist_none,
      Named("dist_5") = out_dist_5,
      Named("dist_4") = out_dist_4,
      Named("dist_3") = out_dist_3,
      Named("dist_2") = out_dist_2,
      Named("dist_1_plus") = out_dist_1_plus,
      Named("stringsAsFactors") = false // Good practice
  );
}

////////////////////////////////////////////////////////////////////////////////
// QueryConsensus functions
////////////////////////////////////////////////////////////////////////////////

// [[Rcpp::export]]
DataFrame aln_query_consensus(
    XPtr<AlignmentStore> store_ptr,
    DataFrame intervals_df,
    double consensus_threshold = 0.9,
    int num_threads = 0,
    std::string clip_mode_str = "all",
    int clip_margin = 10,
    double min_mutations_percent = 0.0,
    double max_mutations_percent = 10.0,
    int min_alignment_length = 0,
    int max_alignment_length = 0,
    int min_indel_length = 3)
{
  // Validate the external pointer
  if (!store_ptr) {
    stop("Invalid AlignmentStore pointer provided.");
  }

  // Get reference to the AlignmentStore object
  AlignmentStore& store = *store_ptr;

  // count short indels before query
  store.count_short_indels(min_indel_length);

  // Convert intervals
  std::vector<Interval> intervals = Rcpp_DataFrame_to_Intervals(intervals_df);

  // Convert clip_mode_str to ClipMode enum
  ClipMode clip_mode = string_to_clip_mode(clip_mode_str);
  
  std::vector<ConsensusOutputRow> results;
  
  QueryConsensus queryConsensus(intervals, store, consensus_threshold, num_threads, clip_mode, clip_margin, min_mutations_percent, max_mutations_percent, min_alignment_length, max_alignment_length);
  queryConsensus.execute();
  results = queryConsensus.get_output_rows();

  // Convert results to R DataFrame
  CharacterVector out_contig;
  IntegerVector out_position;
  CharacterVector out_variant_type;
  CharacterVector out_variant_desc;
  IntegerVector out_count;
  IntegerVector out_coverage;
  NumericVector out_frequency;

  for (const auto& row : results) {
    out_contig.push_back(row.contig);
    out_position.push_back(row.position);
    out_variant_type.push_back(row.variant_type);
    out_variant_desc.push_back(row.variant_desc);
    out_count.push_back(row.count);
    out_coverage.push_back(row.coverage);
    out_frequency.push_back(row.frequency);
  }

  return DataFrame::create(
      Named("contig") = out_contig,
      Named("position") = out_position,
      Named("variant_type") = out_variant_type,
      Named("variant_desc") = out_variant_desc,
      Named("count") = out_count,
      Named("coverage") = out_coverage,
      Named("frequency") = out_frequency,
      Named("stringsAsFactors") = false // Good practice
  );
}

////////////////////////////////////////////////////////////////////////////////
// QueryPileup functions
////////////////////////////////////////////////////////////////////////////////

// [[Rcpp::export]]
DataFrame aln_query_pileup(
    XPtr<AlignmentStore> store_ptr,
    DataFrame intervals_df,
    std::string report_mode_str,
    std::string clip_mode_str = "all",
    int clip_margin = 10,
    double min_mutations_percent = 0.0,
    double max_mutations_percent = 10.0,
    int min_alignment_length = 0,
    int max_alignment_length = 0,
    int min_indel_length = 3)
{
  // Validate the external pointer
  if (!store_ptr) {
    stop("Invalid AlignmentStore pointer provided.");
  }

  // Get reference to the AlignmentStore object
  AlignmentStore& store = *store_ptr;

  // count short indels before query
  store.count_short_indels(min_indel_length);

  // Convert intervals
  std::vector<Interval> intervals = Rcpp_DataFrame_to_Intervals(intervals_df);

  PileupReportMode report_mode = string_to_pileup_report_mode(report_mode_str);
  
  // Convert clip_mode_str to ClipMode enum
  ClipMode clip_mode = string_to_clip_mode(clip_mode_str);
  
  QueryPileup queryPileup(intervals, store, report_mode, clip_mode, clip_margin, min_mutations_percent, max_mutations_percent, min_alignment_length, max_alignment_length);

  // Run the steps
  queryPileup.execute();

  // Get the results
  const std::vector<PileupOutputRow>& results = queryPileup.get_output_rows();

  // Convert results to R DataFrame
  CharacterVector out_contig;
  IntegerVector out_position;
  CharacterVector out_variant;
  IntegerVector out_count;
  IntegerVector out_coverage;
  IntegerVector out_cumsum;

  for (const auto& row : results) {
    out_contig.push_back(row.contig);
    out_position.push_back(row.position);
    out_variant.push_back(row.variant);
    out_count.push_back(row.count);
    out_coverage.push_back(row.coverage);
    out_cumsum.push_back(row.cumsum);
  }

  return DataFrame::create(
      Named("contig") = out_contig,
      Named("position") = out_position,
      Named("variant") = out_variant,
      Named("count") = out_count,
      Named("coverage") = out_coverage,
      Named("cumsum") = out_cumsum,
      Named("stringsAsFactors") = false // Good practice
  );
}

////////////////////////////////////////////////////////////////////////////////
// QueryFull functions
////////////////////////////////////////////////////////////////////////////////

// [[Rcpp::export]]
List aln_query_full(
    XPtr<AlignmentStore> store_ptr,
    DataFrame intervals_df,
    std::string height_style_str = "by_coord_left",
    int max_alignments = 0,
    std::string clip_mode_str = "all",
    int clip_margin = 10,
    double min_mutations_percent = 0.0,
    double max_mutations_percent = 10.0,
    int min_alignment_length = 0,
    int max_alignment_length = 0,
    int min_indel_length = 3)
{
  // Validate the external pointer
  if (!store_ptr) {
    stop("Invalid AlignmentStore pointer provided.");
  }

  // Convert height_style_str to HeightStyle enum
  HeightStyle height_style;
  if (height_style_str == "by_mutations") {
    height_style = HeightStyle::BY_MUTATIONS;
  } else if (height_style_str == "by_coord_left") {
    height_style = HeightStyle::BY_COORD_LEFT;
  } else if (height_style_str == "by_coord_right") {
    height_style = HeightStyle::BY_COORD_RIGHT;
  } else {
    stop("Invalid height_style parameter. Must be 'by_coord_left', 'by_coord_right', or 'by_mutations'.");
  }

  // Get reference to the AlignmentStore object
  AlignmentStore& store = *store_ptr;

  // count short indels before query
  store.count_short_indels(min_indel_length);

  // Convert intervals
  std::vector<Interval> intervals = Rcpp_DataFrame_to_Intervals(intervals_df);

  // Convert clip_mode_str to ClipMode enum
  ClipMode clip_mode = string_to_clip_mode(clip_mode_str);
  
  QueryFull queryFull(intervals, store, height_style, max_alignments, clip_mode, clip_margin, min_mutations_percent, max_mutations_percent, min_alignment_length, max_alignment_length);

  // Run the steps
  queryFull.execute();

  // --- Create Alignments DataFrame ---
  const std::vector<FullOutputAlignments>& alignments = queryFull.get_output_alignments();
  NumericVector out_aln_idx;
  CharacterVector out_aln_read_id;
  IntegerVector out_aln_read_length;
  CharacterVector out_aln_contig_id;
  IntegerVector out_aln_read_start;
  IntegerVector out_aln_read_end;
  IntegerVector out_aln_contig_start;
  IntegerVector out_aln_contig_end;
  LogicalVector out_aln_is_reverse;
  CharacterVector out_aln_cs_tag;
  IntegerVector out_aln_height;
  IntegerVector out_aln_num_mutations;

  for (const auto& aln : alignments) {
    out_aln_idx.push_back(static_cast<double>(aln.alignment_index + 1)); // R numeric can hold uint64_t
    out_aln_read_id.push_back(aln.read_id);
    out_aln_read_length.push_back(aln.read_length);
    out_aln_contig_id.push_back(aln.contig_id);
    out_aln_read_start.push_back(aln.read_start);
    out_aln_read_end.push_back(aln.read_end);
    out_aln_contig_start.push_back(aln.contig_start);
    out_aln_contig_end.push_back(aln.contig_end);
    out_aln_is_reverse.push_back(aln.is_reverse);
    out_aln_cs_tag.push_back(aln.cs_tag);
    out_aln_height.push_back(aln.height);
    out_aln_num_mutations.push_back(aln.num_mutations);
  }

  DataFrame alignments_df = DataFrame::create(
      Named("alignment_index") = out_aln_idx,
      Named("read_id") = out_aln_read_id,
      Named("read_length") = out_aln_read_length,
      Named("contig_id") = out_aln_contig_id,
      Named("read_start") = out_aln_read_start,
      Named("read_end") = out_aln_read_end,
      Named("contig_start") = out_aln_contig_start,
      Named("contig_end") = out_aln_contig_end,
      Named("is_reverse") = out_aln_is_reverse,
      Named("cs_tag") = out_aln_cs_tag,
      Named("mutation_count") = out_aln_num_mutations,
      Named("height") = out_aln_height,
      Named("stringsAsFactors") = false);

  // --- Create Mutations DataFrame ---
  const std::vector<FullOutputMutations>& mutations = queryFull.get_output_mutations();
  NumericVector out_mut_aln_idx;
  CharacterVector out_mut_read_id;
  CharacterVector out_mut_contig_id;
  CharacterVector out_mut_type; // Convert enum to string
  IntegerVector out_mut_position;
  CharacterVector out_mut_desc;
  IntegerVector out_mut_height;

  // lookup table for mutation types to avoid stringstream overhead
  static const std::vector<std::string> mutation_type_lookup = {"SUB", "INS", "DEL"};

  for (const auto& mut : mutations) {
    out_mut_aln_idx.push_back(static_cast<double>(mut.alignment_index + 1)); // R numeric can hold uint64_t
    out_mut_read_id.push_back(mut.read_id);
    out_mut_contig_id.push_back(mut.contig_id);
    // convert MutationType enum to string using lookup table
    out_mut_type.push_back(mutation_type_lookup[static_cast<int>(mut.type)]);
    out_mut_position.push_back(mut.position);
    out_mut_desc.push_back(mut.desc);
    out_mut_height.push_back(mut.height);
  }

  DataFrame mutations_df = DataFrame::create(
      Named("alignment_index") = out_mut_aln_idx,
      Named("read_id") = out_mut_read_id,
      Named("contig_id") = out_mut_contig_id,
      Named("type") = out_mut_type,
      Named("position") = out_mut_position,
      Named("desc") = out_mut_desc,
      Named("height") = out_mut_height,
      Named("stringsAsFactors") = false);

  // --- Create Reads DataFrame ---
  const std::vector<FullOutputReads>& reads = queryFull.get_output_reads();
  CharacterVector out_read_id;
  CharacterVector out_read_contig_id;
  IntegerVector out_read_length;
  IntegerVector out_span_start;
  IntegerVector out_span_end;
  IntegerVector out_total_aligned_length;
  IntegerVector out_num_alignments;
  IntegerVector out_read_num_mutations;
  IntegerVector out_read_height;
  LogicalVector out_read_reversed;

  for (const auto& read : reads) {
    out_read_id.push_back(read.read_id);
    out_read_contig_id.push_back(read.contig_id);
    out_read_length.push_back(read.read_length);
    out_span_start.push_back(read.span_start);
    out_span_end.push_back(read.span_end);
    out_total_aligned_length.push_back(read.total_aligned_length);
    out_num_alignments.push_back(read.num_alignments);
    out_read_num_mutations.push_back(read.num_mutations);
    out_read_height.push_back(read.height);
    out_read_reversed.push_back(read.read_reversed);
  }

  DataFrame reads_df = DataFrame::create(
      Named("read_id") = out_read_id,
      Named("contig_id") = out_read_contig_id,
      Named("read_length") = out_read_length,
      Named("span_start") = out_span_start,
      Named("span_end") = out_span_end,
      Named("total_aligned_length") = out_total_aligned_length,
      Named("num_alignments") = out_num_alignments,
      Named("num_mutations") = out_read_num_mutations,
      Named("height") = out_read_height,
      Named("is_reverse") = out_read_reversed,
      Named("stringsAsFactors") = false);

  return List::create(
      Named("alignments") = alignments_df,
      Named("mutations") = mutations_df,
      Named("reads") = reads_df);
}

////////////////////////////////////////////////////////////////////////////////
// Construct AlignmentStore from PAF file
////////////////////////////////////////////////////////////////////////////////

// [[Rcpp::export]]
XPtr<AlignmentStore> aln_construct(
    std::string paf_file,
    int max_reads = 0)
{
  // Create a new AlignmentStore instance on the heap
  AlignmentStore* store = new AlignmentStore();

  try {
    // Create a PafReader to handle the file
    PafReader reader;

    // Read the PAF file without verification
    Rcout << "reading PAF file: " << paf_file << "\n";
    reader.read_paf(paf_file, *store, max_reads, false, true);

    // Organize alignments after loading
    store->organize_alignments();

    // Create an external pointer managed by R's garbage collector
    XPtr<AlignmentStore> ptr(store, true);
    return ptr;

  } catch (const std::runtime_error& e) {
    // Clean up allocated store and signal an error to R
    delete store;
    stop("failed to construct AlignmentStore: %s", e.what());
  } catch (...) {
    delete store;
    stop("an unknown C++ error occurred during construction");
  }
}

////////////////////////////////////////////////////////////////////////////////
// Get all contigs
////////////////////////////////////////////////////////////////////////////////

// [[Rcpp::export]]
DataFrame aln_get_contigs(XPtr<AlignmentStore> store_ptr)
{
  // validate the external pointer
  if (!store_ptr) {
    stop("invalid AlignmentStore pointer provided");
  }

  // get reference to the AlignmentStore object
  const AlignmentStore& store = *store_ptr;

  // get all contigs
  const std::vector<Contig>& contigs = store.get_contigs();

  // prepare output vectors
  CharacterVector out_contig_id;
  IntegerVector out_length;

  // fill output vectors
  for (const auto& contig : contigs) {
    out_contig_id.push_back(contig.id);
    out_length.push_back(contig.length);
  }

  return DataFrame::create(
      Named("contig_id") = out_contig_id,
      Named("length") = out_length,
      Named("stringsAsFactors") = false);
}

// [[Rcpp::export]]
DataFrame aln_get_reads(XPtr<AlignmentStore> store_ptr)
{
  // validate the external pointer
  if (!store_ptr) {
    stop("invalid AlignmentStore pointer provided");
  }

  // get reference to the AlignmentStore object
  const AlignmentStore& store = *store_ptr;

  // get all reads
  const std::vector<Read>& reads = store.get_reads();

  // prepare output vectors
  CharacterVector out_read_id;
  IntegerVector out_length;

  // fill output vectors
  for (const auto& read : reads) {
    out_read_id.push_back(read.id);
    out_length.push_back(read.length);
  }

  return DataFrame::create(
      Named("read_id") = out_read_id,
      Named("length") = out_length,
      Named("stringsAsFactors") = false);
}

// [[Rcpp::export]]
void aln_save(XPtr<AlignmentStore> store_ptr, std::string filepath)
{
  // Validate the external pointer
  if (!store_ptr) {
    stop("invalid AlignmentStore pointer provided");
  }

  try {
    // Save the alignment store to a file
    Rcout << "saving AlignmentStore to: " << filepath << "\n";
    store_ptr->save(filepath);

  } catch (const std::runtime_error& e) {
    stop("failed to save AlignmentStore: %s", e.what());
  } catch (...) {
    stop("an unknown C++ error occurred during saving");
  }
}

////////////////////////////////////////////////////////////////////////////////
// Break position detection
////////////////////////////////////////////////////////////////////////////////

// [[Rcpp::export]]
DataFrame aln_find_breaks(
    XPtr<AlignmentStore> store_ptr,
    int window_size,
    double p_threshold,
    int min_reads = 1)
{
  // validate the external pointer
  if (!store_ptr) {
    stop("invalid AlignmentStore pointer provided");
  }

  // validate parameters
  if (window_size <= 0) {
    stop("window size must be positive");
  }
  if (p_threshold <= 0.0 || p_threshold > 1.0) {
    stop("p threshold must be between 0 and 1");
  }
  if (min_reads <= 0) {
    stop("min reads must be positive");
  }

  // get reference to the AlignmentStore object
  const AlignmentStore& store = *store_ptr;

  try {
    // find break positions
    Rcout << "finding break positions with window size " << window_size 
          << ", p threshold " << p_threshold 
          << ", and min reads " << min_reads << "\n";
    
    auto results = store.find_break_positions(window_size, p_threshold, min_reads);

    // prepare output vectors
    CharacterVector out_contig_id;
    IntegerVector out_position;
    CharacterVector out_orientation;
    IntegerVector out_t;
    NumericVector out_e;
    NumericVector out_enrichment;
    NumericVector out_pval;
    NumericVector out_qval;

    // fill output vectors
    for (const auto& result : results) {
      out_contig_id.push_back(result.contig_id);
      out_position.push_back(result.position);
      out_orientation.push_back(result.orientation);
      out_t.push_back(result.t);
      out_e.push_back(result.e);
      out_enrichment.push_back(result.enrichment);
      out_pval.push_back(result.pval);
      out_qval.push_back(result.qval);
    }

    Rcout << "found " << results.size() << " significant break positions\n";

    return DataFrame::create(
        Named("contig") = out_contig_id,
        Named("position") = out_position,
        Named("orientation") = out_orientation,
        Named("t") = out_t,
        Named("e") = out_e,
        Named("enrichment") = out_enrichment,
        Named("pval") = out_pval,
        Named("qval") = out_qval,
        Named("stringsAsFactors") = false);

  } catch (const std::runtime_error& e) {
    stop("failed to find break positions: %s", e.what());
  } catch (...) {
    stop("an unknown C++ error occurred during break position detection");
  }
}
