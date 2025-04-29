#include "QueryBin.h" // Include QueryBin header
#include "alignment_store.h"
#include <Rcpp.h>
#include <stdexcept> // To catch potential exceptions from AlignmentStore::load
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
