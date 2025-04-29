#include "alignment_store.h"
#include <Rcpp.h>
#include <stdexcept> // To catch potential exceptions from AlignmentStore::load
#include <string>

// [[Rcpp::export]]
Rcpp::XPtr<AlignmentStore> aln_load(std::string filepath)
{
  // Create a new AlignmentStore instance on the heap
  AlignmentStore* store = new AlignmentStore();

  try {
    // Attempt to load the data
    Rcpp::Rcout << "Loading AlignmentStore from: " << filepath << std::endl;
    store->load(filepath);
    Rcpp::Rcout << "Load successful." << std::endl;

    // Create an external pointer managed by R's garbage collector
    Rcpp::XPtr<AlignmentStore> ptr(store, true);
    return ptr;
  } catch (const std::runtime_error& e) {
    // If load fails (e.g., file not found, assertion fails),
    // clean up the allocated store and signal an error to R.
    delete store;
    Rcpp::stop("Failed to load AlignmentStore: %s", e.what());
  } catch (...) {
    // Catch any other potential C++ exceptions
    delete store;
    Rcpp::stop("An unknown C++ error occurred during loading.");
  }
}