#include "QueryBin.h"
#include "QueryFull.h"
#include "QueryPileup.h"
#include "alignment_store.h"
#include <Rcpp.h>
#include <sstream>
#include <stdexcept>
#include <string>
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
    int binsize)
{
  // Validate the external pointer
  if (!store_ptr) {
    stop("Invalid AlignmentStore pointer provided.");
  }

  // Get reference to the AlignmentStore object
  const AlignmentStore& store = *store_ptr;

  // Convert intervals
  std::vector<Interval> intervals = Rcpp_DataFrame_to_Intervals(intervals_df);

  QueryBin queryBin(intervals, store, binsize);

  // Run the steps
  queryBin.execute();

  // Get the results
  const std::vector<BinOutputRow>& results = queryBin.get_output_rows();

  // Convert results to R DataFrame
  CharacterVector out_contig;
  IntegerVector out_bin_start;
  IntegerVector out_bin_end;
  IntegerVector out_bin_length;
  IntegerVector out_sequenced_bp; // Using IntegerVector for R 'integer' type
  IntegerVector out_mutation_count;

  for (const auto& row : results) {
    out_contig.push_back(row.contig);
    out_bin_start.push_back(row.bin_start);
    out_bin_end.push_back(row.bin_end);
    out_bin_length.push_back(row.bin_length);
    out_sequenced_bp.push_back(row.sequenced_basepairs);
    out_mutation_count.push_back(row.mutation_count);
  }

  return DataFrame::create(
      Named("contig") = out_contig,
      Named("start") = out_bin_start,
      Named("end") = out_bin_end,
      Named("length") = out_bin_length,
      Named("read_count") = out_sequenced_bp,
      Named("mutation_count") = out_mutation_count,
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
    std::string report_mode_str)
{
  // Validate the external pointer
  if (!store_ptr) {
    stop("Invalid AlignmentStore pointer provided.");
  }

  // Get reference to the AlignmentStore object
  const AlignmentStore& store = *store_ptr;

  // Convert intervals
  std::vector<Interval> intervals = Rcpp_DataFrame_to_Intervals(intervals_df);

  PileupReportMode report_mode = string_to_pileup_report_mode(report_mode_str);
  QueryPileup queryPileup(intervals, store, report_mode);

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
    DataFrame intervals_df)
{
  // Validate the external pointer
  if (!store_ptr) {
    stop("Invalid AlignmentStore pointer provided.");
  }

  // Get reference to the AlignmentStore object
  const AlignmentStore& store = *store_ptr;

  // Convert intervals
  std::vector<Interval> intervals = Rcpp_DataFrame_to_Intervals(intervals_df);

  QueryFull queryFull(intervals, store);

  // Run the steps
  queryFull.execute();

  // --- Create Alignments DataFrame ---
  const std::vector<FullOutputAlignments>& alignments = queryFull.get_output_alignments();
  NumericVector out_aln_idx;
  CharacterVector out_aln_read_id;
  CharacterVector out_aln_contig_id;
  IntegerVector out_aln_read_start;
  IntegerVector out_aln_read_end;
  IntegerVector out_aln_contig_start;
  IntegerVector out_aln_contig_end;
  LogicalVector out_aln_is_reverse;
  CharacterVector out_aln_cs_tag;

  for (const auto& aln : alignments) {
    out_aln_idx.push_back(static_cast<double>(aln.alignment_index + 1)); // R numeric can hold uint64_t
    out_aln_read_id.push_back(aln.read_id);
    out_aln_contig_id.push_back(aln.contig_id);
    out_aln_read_start.push_back(aln.read_start);
    out_aln_read_end.push_back(aln.read_end);
    out_aln_contig_start.push_back(aln.contig_start);
    out_aln_contig_end.push_back(aln.contig_end);
    out_aln_is_reverse.push_back(aln.is_reverse);
    out_aln_cs_tag.push_back(aln.cs_tag);
  }

  DataFrame alignments_df = DataFrame::create(
      Named("alignment_index") = out_aln_idx,
      Named("read_id") = out_aln_read_id,
      Named("contig_id") = out_aln_contig_id,
      Named("read_start") = out_aln_read_start,
      Named("read_end") = out_aln_read_end,
      Named("contig_start") = out_aln_contig_start,
      Named("contig_end") = out_aln_contig_end,
      Named("is_reverse") = out_aln_is_reverse,
      Named("cs_tag") = out_aln_cs_tag,
      Named("stringsAsFactors") = false);

  // --- Create Mutations DataFrame ---
  const std::vector<FullOutputMutations>& mutations = queryFull.get_output_mutations();
  NumericVector out_mut_aln_idx;
  CharacterVector out_mut_read_id;
  CharacterVector out_mut_contig_id;
  CharacterVector out_mut_type; // Convert enum to string
  IntegerVector out_mut_position;
  CharacterVector out_mut_desc;

  for (const auto& mut : mutations) {
    out_mut_aln_idx.push_back(static_cast<double>(mut.alignment_index + 1)); // R numeric can hold uint64_t
    out_mut_read_id.push_back(mut.read_id);
    out_mut_contig_id.push_back(mut.contig_id);
    // Convert MutationType enum to string for R
    std::stringstream ss;
    ss << mut.type; // Use the overloaded operator<< from aln_types.h
    out_mut_type.push_back(ss.str());
    out_mut_position.push_back(mut.position);
    out_mut_desc.push_back(mut.desc);
  }

  DataFrame mutations_df = DataFrame::create(
      Named("alignment_index") = out_mut_aln_idx,
      Named("read_id") = out_mut_read_id,
      Named("contig_id") = out_mut_contig_id,
      Named("type") = out_mut_type,
      Named("position") = out_mut_position,
      Named("desc") = out_mut_desc,
      Named("stringsAsFactors") = false);

  return List::create(
      Named("alignments") = alignments_df,
      Named("mutations") = mutations_df);
}