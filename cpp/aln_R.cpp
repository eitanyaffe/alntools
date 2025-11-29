// [[Rcpp::plugins(cpp17)]]

// OpenMP support - configured via environment variables in R startup
#ifdef _OPENMP
#include <omp.h>
#endif

#include "QueryBin.h"
#include "QueryConsensus.h"
#include "QueryFull.h"
#include "QueryPileup.h"
#include "QueryVariants.h"
#include "Rearrange.h"
#include "RearrangeVerify.h"
#include "Sequences.h"
#include "Homologs.h"
#include "alignment_store.h"
#include "paf_reader.h"
#include "utils.h"
#include <Rcpp.h>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_set>
#include <vector>

using namespace std;
using namespace Rcpp;

// Helper functions to extract R data to C++ containers
bool extract_genes_from_dataframe(const DataFrame& gene_df, std::vector<Gene>& gene_vector) {
    // check required columns
    if (!gene_df.containsElementNamed("gene") || !gene_df.containsElementNamed("contig") ||
        !gene_df.containsElementNamed("start") || !gene_df.containsElementNamed("end") ||
        !gene_df.containsElementNamed("strand")) {
        Rcpp::Rcerr << "error: gene DataFrame must contain columns: gene, contig, start, end, strand" << std::endl;
        return false;
    }
    
    // extract columns
    CharacterVector gene_ids = gene_df["gene"];
    CharacterVector contigs = gene_df["contig"];
    IntegerVector starts = gene_df["start"];
    IntegerVector ends = gene_df["end"];
    CharacterVector strands = gene_df["strand"];
    
    // optional description column
    CharacterVector descriptions;
    bool has_desc = gene_df.containsElementNamed("desc");
    if (has_desc) {
        descriptions = gene_df["desc"];
    }
    
    int nrows = gene_df.nrows();
    gene_vector.clear();
    gene_vector.reserve(nrows);
    
    for (int i = 0; i < nrows; i++) {
        std::string gene_id = as<std::string>(gene_ids[i]);
        std::string contig = as<std::string>(contigs[i]);
        uint32_t start = static_cast<uint32_t>(starts[i]);
        uint32_t end = static_cast<uint32_t>(ends[i]);
        std::string strand_str = as<std::string>(strands[i]);
        
        char strand = '+';
        if (!strand_str.empty()) {
            strand = strand_str[0];
        }
        
        std::string desc = "";
        if (has_desc && i < descriptions.size()) {
            desc = as<std::string>(descriptions[i]);
        }
        
        // create gene object and add to vector
        gene_vector.emplace_back(gene_id, contig, start, end, strand, desc);
    }
    
    return true;
}

bool extract_sequences_from_list(const List& seq_list, std::unordered_map<std::string, std::string>& sequences_map) {
    // check if list has names
    CharacterVector names = seq_list.names();
    if (names.size() == 0) {
        Rcpp::Rcerr << "error: reference sequences List must be named (names = contig IDs)" << std::endl;
        return false;
    }
    
    int n_sequences = seq_list.size();
    sequences_map.clear();
    
    for (int i = 0; i < n_sequences; i++) {
        std::string contig_id = as<std::string>(names[i]);
        
        // sequences can be either character vectors or lists (from seqinr)
        SEXP seq_element = seq_list[i];
        std::string sequence;
        
        if (TYPEOF(seq_element) == STRSXP) {
            // character vector - join all elements
            CharacterVector seq_vec = as<CharacterVector>(seq_element);
            for (int j = 0; j < seq_vec.size(); j++) {
                sequence += as<std::string>(seq_vec[j]);
            }
        } else if (TYPEOF(seq_element) == VECSXP) {
            // list (seqinr format) - extract sequence
            List seq_obj = as<List>(seq_element);
            if (seq_obj.containsElementNamed("seq")) {
                CharacterVector seq_vec = seq_obj["seq"];
                for (int j = 0; j < seq_vec.size(); j++) {
                    sequence += as<std::string>(seq_vec[j]);
                }
            } else {
                // assume the list itself contains the sequence
                for (int j = 0; j < seq_obj.size(); j++) {
                    sequence += as<std::string>(seq_obj[j]);
                }
            }
        } else {
            Rcpp::Rcerr << "warning: unsupported sequence format for contig " << contig_id << ", skipping" << std::endl;
            continue;
        }
        
        // convert to uppercase
        std::transform(sequence.begin(), sequence.end(), sequence.begin(), ::toupper);
        
        sequences_map[contig_id] = sequence;
    }
    
    return true;
}

////////////////////////////////////////////////////////////////////////////////
// Gene management functions
////////////////////////////////////////////////////////////////////////////////

// [[Rcpp::export]]
SEXP aln_genes(DataFrame gene_table, List reference_sequences, std::string codon_table_path) {
    // extract gene table from DataFrame to C++ vector
    std::vector<Gene> gene_vector;
    if (!extract_genes_from_dataframe(gene_table, gene_vector)) {
        stop("failed to extract gene table from DataFrame");
    }
    
    // extract reference sequences from List to C++ map
    std::unordered_map<std::string, std::string> sequences_map;
    if (!extract_sequences_from_list(reference_sequences, sequences_map)) {
        stop("failed to extract reference sequences from List");
    }
    
    // create genes object
    Genes* genes = new Genes();
    
    // load data into genes object
    if (!genes->load_gene_table_from_vector(gene_vector)) {
        delete genes;
        stop("failed to load gene table into Genes object");
    }
    
    if (!genes->load_codon_table(codon_table_path)) {
        delete genes;
        stop("failed to load codon table from %s", codon_table_path.c_str());
    }
    
    if (!genes->load_reference_sequences_from_map(sequences_map)) {
        delete genes;
        stop("failed to load reference sequences into Genes object");
    }
    
    // return as external pointer
    XPtr<Genes> genes_ptr(genes);
    return genes_ptr;
}

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
// QueryByReadId function - returns alignments and mutations for a single read
////////////////////////////////////////////////////////////////////////////////

// [[Rcpp::export]]
List aln_alignments_from_read_id(
    XPtr<AlignmentStore> store_ptr,
    std::string read_id,
    int min_indel_length = 3)
{
  if (!store_ptr) {
    stop("invalid AlignmentStore pointer provided");
  }

  AlignmentStore& store = *store_ptr;

  // count short indels for filtering (same as QueryFull)
  store.count_short_indels(min_indel_length);

  if (!store.has_read_id(read_id)) {
    stop("read ID '%s' not found in alignment store", read_id.c_str());
  }

  size_t read_index = store.get_read_index(read_id);
  uint32_t read_length = store.get_reads()[read_index].length;

  const std::vector<Alignment>& all_alignments = store.get_alignments();

  // collect alignments for this read
  struct AlignmentData {
    size_t index;
    const Alignment* aln;
  };
  std::vector<AlignmentData> read_alignments;

  for (size_t i = 0; i < all_alignments.size(); ++i) {
    const Alignment& aln = all_alignments[i];
    if (aln.read_index == read_index) {
      read_alignments.push_back({i, &aln});
    }
  }

  // build alignment output vectors
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

  // build mutation output vectors
  NumericVector out_mut_aln_idx;
  CharacterVector out_mut_type;
  IntegerVector out_mut_contig_coord;
  IntegerVector out_mut_read_coord;
  CharacterVector out_mut_desc;

  for (size_t i = 0; i < read_alignments.size(); ++i) {
    size_t aln_index = read_alignments[i].index;
    const Alignment* aln = read_alignments[i].aln;

    out_aln_idx.push_back(static_cast<double>(aln_index));
    out_aln_read_id.push_back(read_id);
    out_aln_read_length.push_back(static_cast<int>(read_length));
    out_aln_contig_id.push_back(store.get_contig_id(aln->contig_index));
    // convert to 1-based coordinates for R output
    out_aln_read_start.push_back(static_cast<int>(aln->read_start + 1));
    out_aln_read_end.push_back(static_cast<int>(aln->read_end));
    out_aln_contig_start.push_back(static_cast<int>(aln->contig_start + 1));
    out_aln_contig_end.push_back(static_cast<int>(aln->contig_end));
    out_aln_is_reverse.push_back(aln->is_reverse);
    out_aln_num_mutations.push_back(aln->get_mutation_count());

    // process mutations for this alignment
    int contig_start = static_cast<int>(aln->contig_start);
    int contig_end = static_cast<int>(aln->contig_end);
    int r_start = static_cast<int>(aln->read_start);
    int r_end = static_cast<int>(aln->read_end);

    for (uint32_t mut_idx : aln->mutations) {
      // skip short indels (same logic as QueryFull)
      if (should_skip_short_indel(store, aln->contig_index, mut_idx)) {
        continue;
      }

      const Mutation& mut = store.get_mutation(aln->contig_index, mut_idx);
      int contig_coord = static_cast<int>(mut.position);

      // linear extrapolation to read space
      int read_coord;
      if (aln->is_reverse) {
        // reverse strand: read_coord decreases as contig_coord increases
        double frac = static_cast<double>(contig_coord - contig_start) / 
                      static_cast<double>(contig_end - contig_start);
        read_coord = r_end - static_cast<int>(frac * (r_end - r_start));
      } else {
        // forward strand: read_coord increases with contig_coord
        double frac = static_cast<double>(contig_coord - contig_start) / 
                      static_cast<double>(contig_end - contig_start);
        read_coord = r_start + static_cast<int>(frac * (r_end - r_start));
      }

      // mutation type as string
      std::string type_str;
      switch (mut.type) {
        case MutationType::SUBSTITUTION: type_str = "SUB"; break;
        case MutationType::INSERTION: type_str = "INS"; break;
        case MutationType::DELETION: type_str = "DEL"; break;
        default: type_str = "UNK"; break;
      }

      out_mut_aln_idx.push_back(static_cast<double>(aln_index));
      out_mut_type.push_back(type_str);
      out_mut_contig_coord.push_back(contig_coord + 1); // 1-based for R
      out_mut_read_coord.push_back(read_coord);
      out_mut_desc.push_back(mut.to_string());
    }
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
      Named("num_mutations") = out_aln_num_mutations,
      Named("stringsAsFactors") = false);

  DataFrame mutations_df = DataFrame::create(
      Named("alignment_index") = out_mut_aln_idx,
      Named("type") = out_mut_type,
      Named("contig_coord") = out_mut_contig_coord,
      Named("read_coord") = out_mut_read_coord,
      Named("desc") = out_mut_desc,
      Named("stringsAsFactors") = false);

  return List::create(
      Named("alignments") = alignments_df,
      Named("mutations") = mutations_df,
      Named("read_length") = static_cast<int>(read_length));
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
    int min_consensus_coverage = 5,
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
  
  QueryConsensus queryConsensus(intervals, store, consensus_threshold, min_consensus_coverage, num_threads, clip_mode, clip_margin, min_mutations_percent, max_mutations_percent, min_alignment_length, max_alignment_length);
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
    int max_margin = 10,
    std::string chunk_type_str = "break_on_overlap",
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

  // Convert chunk_type_str to ChunkType enum
  ChunkType chunk_type;
  if (chunk_type_str == "read") {
    chunk_type = ChunkType::READ;
  } else if (chunk_type_str == "alignment") {
    chunk_type = ChunkType::ALIGNMENT;
  } else if (chunk_type_str == "break_on_overlap") {
    chunk_type = ChunkType::BREAK_ON_OVERLAP;
  } else if (chunk_type_str == "break_on_gap") {
    chunk_type = ChunkType::BREAK_ON_GAP;
  } else {
    stop("Invalid chunk_type parameter. Must be 'read', 'alignment', 'break_on_overlap', or 'break_on_gap'.");
  }

  // Get reference to the AlignmentStore object
  AlignmentStore& store = *store_ptr;

  // count short indels before query
  store.count_short_indels(min_indel_length);

  // Convert intervals
  std::vector<Interval> intervals = Rcpp_DataFrame_to_Intervals(intervals_df);

  // Convert clip_mode_str to ClipMode enum
  ClipMode clip_mode = string_to_clip_mode(clip_mode_str);
  
  QueryFull queryFull(intervals, store, height_style, max_alignments, clip_mode, clip_margin, min_mutations_percent, max_mutations_percent, min_alignment_length, max_alignment_length, max_margin, chunk_type);

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
  CharacterVector out_aln_chunk_id;

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
    out_aln_chunk_id.push_back(aln.chunk_id);
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
      Named("chunk_id") = out_aln_chunk_id,
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

  // --- Create Chunks DataFrame ---
  const std::vector<FullOutputChunk>& chunks = queryFull.get_output_chunks();
  CharacterVector out_chunk_id;
  CharacterVector out_chunk_read_id;
  CharacterVector out_chunk_contig_id;
  IntegerVector out_chunk_read_length;
  IntegerVector out_chunk_span_start;
  IntegerVector out_chunk_span_end;
  IntegerVector out_chunk_total_aligned_length;
  IntegerVector out_chunk_num_alignments;
  IntegerVector out_chunk_num_mutations;
  IntegerVector out_chunk_height;
  LogicalVector out_chunk_read_reversed;

  for (const auto& chunk : chunks) {
    out_chunk_id.push_back(chunk.chunk_id);
    out_chunk_read_id.push_back(chunk.read_id);
    out_chunk_contig_id.push_back(chunk.contig_id);
    out_chunk_read_length.push_back(chunk.read_length);
    out_chunk_span_start.push_back(chunk.span_start);
    out_chunk_span_end.push_back(chunk.span_end);
    out_chunk_total_aligned_length.push_back(chunk.total_aligned_length);
    out_chunk_num_alignments.push_back(chunk.num_alignments);
    out_chunk_num_mutations.push_back(chunk.num_mutations);
    out_chunk_height.push_back(chunk.height);
    out_chunk_read_reversed.push_back(chunk.read_reversed);
  }

  DataFrame chunks_df = DataFrame::create(
      Named("chunk_id") = out_chunk_id,
      Named("read_id") = out_chunk_read_id,
      Named("contig_id") = out_chunk_contig_id,
      Named("read_length") = out_chunk_read_length,
      Named("span_start") = out_chunk_span_start,
      Named("span_end") = out_chunk_span_end,
      Named("total_aligned_length") = out_chunk_total_aligned_length,
      Named("num_alignments") = out_chunk_num_alignments,
      Named("num_mutations") = out_chunk_num_mutations,
      Named("height") = out_chunk_height,
      Named("read_reversed") = out_chunk_read_reversed,
      Named("stringsAsFactors") = false);

  return List::create(
      Named("alignments") = alignments_df,
      Named("mutations") = mutations_df,
      Named("reads") = reads_df,
      Named("chunks") = chunks_df);
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

////////////////////////////////////////////////////////////////////////////////
// QueryVariants functions
////////////////////////////////////////////////////////////////////////////////

// [[Rcpp::export]]
List aln_query_variants(
    List store_list,
    DataFrame intervals_df,
    int min_variants_variant_support = 3,
    int min_variants_library_support = 1,
    int min_variants_coverage_support = 10,
    std::string clip_mode_str = "all",
    int clip_margin = 10,
    double min_mutations_percent = 0.0,
    double max_mutations_percent = 10.0,
    int min_alignment_length = 0,
    int max_alignment_length = 0,
    int min_indel_length = 3,
    SEXP genes = R_NilValue)
{
    // extract genes pointer if provided
    Genes* genes_ptr = nullptr;
    if (!Rf_isNull(genes)) {
        XPtr<Genes> genes_xptr(genes);
        genes_ptr = genes_xptr.get();
    }
    
    // validate store_list
    if (store_list.size() == 0) {
      stop("store_list must contain at least one AlignmentStore");
    }
    
    // extract library names (preserves order from R)
    CharacterVector lib_names = store_list.names();
    if (lib_names.size() != store_list.size()) {
      stop("all elements in store_list must be named");
    }
    
    // build ordered library_ids vector
    std::vector<std::string> library_ids;
    for (int i = 0; i < lib_names.size(); ++i) {
      library_ids.push_back(as<std::string>(lib_names[i]));
    }
    
    // convert to C++ stores map
    std::map<std::string, AlignmentStore> stores;
    for (int i = 0; i < store_list.size(); ++i) {
      std::string lib_id = as<std::string>(lib_names[i]);
      XPtr<AlignmentStore> store_ptr = store_list[i];
      
      if (!store_ptr) {
        stop("invalid AlignmentStore pointer for library: %s", lib_id.c_str());
      }
      
      // copy the store (since QueryVariants expects ownership)
      stores[lib_id] = *store_ptr;
      
      // count short indels for this store
      stores[lib_id].count_short_indels(min_indel_length);
    }
    
    // convert intervals
    std::vector<Interval> intervals = Rcpp_DataFrame_to_Intervals(intervals_df);
    
    // convert clip_mode_str to ClipMode enum
    ClipMode clip_mode = string_to_clip_mode(clip_mode_str);
    
    // execute QueryVariants
    QueryVariants queryVariants(intervals, stores, library_ids, min_variants_variant_support, min_variants_library_support, 
                     min_variants_coverage_support, clip_mode, clip_margin, min_mutations_percent, 
                     max_mutations_percent, min_alignment_length, max_alignment_length, genes_ptr);
    queryVariants.execute();
    
    // get results
    const std::vector<VariantOutputRow>& variant_rows = queryVariants.get_variant_rows();
    const std::vector<std::string>& library_ids_result = queryVariants.get_library_ids();
    const std::map<std::string, SetVariantData>& variant_data = queryVariants.get_variant_data();
    
    // convert variants table to R DataFrame
    CharacterVector out_variant_id;
    CharacterVector out_contig;
    IntegerVector out_coord;
    CharacterVector out_type;
    CharacterVector out_sequence;
    CharacterVector out_desc;
    IntegerVector out_library_count;
    IntegerVector out_total_support;
    IntegerVector out_total_coverage;
    NumericVector out_frequency;
    
    for (const auto& row : variant_rows) {
      out_variant_id.push_back(row.variant_id);
      out_contig.push_back(row.contig);
      out_coord.push_back(row.coord);
      out_type.push_back(row.type);
      out_sequence.push_back(row.sequence);
      out_desc.push_back(row.desc);
      out_library_count.push_back(row.library_count);
      out_total_support.push_back(row.total_support);
      out_total_coverage.push_back(row.total_coverage);
      out_frequency.push_back(row.frequency);
    }
    
    DataFrame variants_df = DataFrame::create(
        Named("variant_id") = out_variant_id,
        Named("contig") = out_contig,
        Named("coord") = out_coord,
        Named("type") = out_type,
        Named("sequence") = out_sequence,
        Named("desc") = out_desc,
        Named("library_count") = out_library_count,
        Named("total_support") = out_total_support,
        Named("total_coverage") = out_total_coverage,
        Named("frequency") = out_frequency,
        Named("stringsAsFactors") = false);
    
    // create support matrix
    IntegerMatrix support_matrix(variant_rows.size(), library_ids.size());
    CharacterVector support_rownames;
    CharacterVector support_colnames;
    
    for (size_t i = 0; i < variant_rows.size(); ++i) {
      const VariantOutputRow& row = variant_rows[i];
      support_rownames.push_back(row.variant_id);
      
      // find the original variant data
      std::string variant_key = row.contig + ":" + std::to_string(row.coord) + ":" + row.type + ":" + row.sequence;
      auto variant_it = variant_data.find(variant_key);
      
      for (size_t j = 0; j < library_ids_result.size(); ++j) {
        if (i == 0) {
          support_colnames.push_back(library_ids_result[j]);
        }
        
        int support = 0;
        if (variant_it != variant_data.end()) {
          const SetVariantData& variant = variant_it->second;
          auto support_it = variant.lib_support.find(library_ids_result[j]);
          if (support_it != variant.lib_support.end()) {
            support = support_it->second;
          }
        }
        support_matrix(i, j) = support;
      }
    }
    
    rownames(support_matrix) = support_rownames;
    colnames(support_matrix) = support_colnames;
    
    // create coverage matrix
    IntegerMatrix coverage_matrix(variant_rows.size(), library_ids_result.size());
    
    for (size_t i = 0; i < variant_rows.size(); ++i) {
      const VariantOutputRow& row = variant_rows[i];
      
      // find the original variant data
      std::string variant_key = row.contig + ":" + std::to_string(row.coord) + ":" + row.type + ":" + row.sequence;
      auto variant_it = variant_data.find(variant_key);
      
      for (size_t j = 0; j < library_ids_result.size(); ++j) {
        int coverage = 0;
        if (variant_it != variant_data.end()) {
          const SetVariantData& variant = variant_it->second;
          auto coverage_it = variant.lib_coverage.find(library_ids_result[j]);
          if (coverage_it != variant.lib_coverage.end()) {
            coverage = coverage_it->second;
          }
        }
        coverage_matrix(i, j) = coverage;
      }
    }
    
    rownames(coverage_matrix) = support_rownames;
    colnames(coverage_matrix) = support_colnames;
    
    // create gene annotation dataframes if genes were used
    DataFrame genic_df = DataFrame::create();
    DataFrame intergenic_df = DataFrame::create();
    
    if (genes_ptr != nullptr) {
      // get gene annotation data
      const std::vector<GenicRow>& genic_rows = queryVariants.get_genic_rows();
      const std::vector<IntergenicRow>& intergenic_rows = queryVariants.get_intergenic_rows();
      
      // create genic dataframe
      if (!genic_rows.empty()) {
        CharacterVector genic_row_id;
        CharacterVector genic_gene_id;
        CharacterVector genic_gene_desc;
        IntegerVector genic_aa_coord;
        CharacterVector genic_variant_codon;
        CharacterVector genic_ref_codon;
        CharacterVector genic_variant_type;
        CharacterVector genic_mutation_desc;
        
        for (const auto& row : genic_rows) {
          genic_row_id.push_back(row.row_id);
          genic_gene_id.push_back(row.gene_id);
          genic_gene_desc.push_back(row.gene_desc);
          genic_aa_coord.push_back(row.aa_coord);
          genic_variant_codon.push_back(row.variant_codon);
          genic_ref_codon.push_back(row.ref_codon);
          genic_variant_type.push_back(row.variant_type);
          genic_mutation_desc.push_back(row.mutation_desc);
        }
        
        genic_df = DataFrame::create(
            Named("row_id") = genic_row_id,
            Named("gene_id") = genic_gene_id,
            Named("gene_desc") = genic_gene_desc,
            Named("aa_coord") = genic_aa_coord,
            Named("variant_codon") = genic_variant_codon,
            Named("ref_codon") = genic_ref_codon,
            Named("variant_type") = genic_variant_type,
            Named("mutation_desc") = genic_mutation_desc,
            Named("stringsAsFactors") = false);
      }
      
      // create intergenic dataframe
      if (!intergenic_rows.empty()) {
        CharacterVector intergenic_row_id;
        CharacterVector intergenic_gene_left;
        CharacterVector intergenic_gene_right;
        CharacterVector intergenic_orientation_left;
        CharacterVector intergenic_orientation_right;
        IntegerVector intergenic_distance_left;
        IntegerVector intergenic_distance_right;
        
        for (const auto& row : intergenic_rows) {
          intergenic_row_id.push_back(row.row_id);
          intergenic_gene_left.push_back(row.gene_left);
          intergenic_gene_right.push_back(row.gene_right);
          intergenic_orientation_left.push_back(row.orientation_left);
          intergenic_orientation_right.push_back(row.orientation_right);
          intergenic_distance_left.push_back(row.distance_left);
          intergenic_distance_right.push_back(row.distance_right);
        }
        
        intergenic_df = DataFrame::create(
            Named("row_id") = intergenic_row_id,
            Named("gene_left") = intergenic_gene_left,
            Named("gene_right") = intergenic_gene_right,
            Named("orientation_left") = intergenic_orientation_left,
            Named("orientation_right") = intergenic_orientation_right,
            Named("distance_left") = intergenic_distance_left,
            Named("distance_right") = intergenic_distance_right,
            Named("stringsAsFactors") = false);
      }
    }
    
    // return as named list
    return List::create(
        Named("variants") = variants_df,
        Named("support") = support_matrix,
        Named("coverage") = coverage_matrix,
        Named("library_ids") = CharacterVector(library_ids.begin(), library_ids.end()),
        Named("genic") = genic_df,
        Named("intergenic") = intergenic_df);
}

////////////////////////////////////////////////////////////////////////////////
// Rearrangement Analysis
////////////////////////////////////////////////////////////////////////////////

// Helper function to convert RearrangementType to string
std::string rearrangement_type_to_string(RearrangementType type) {
    switch (type) {
        case RearrangementType::LARGE_INSERT: return "large_insert";
        case RearrangementType::LARGE_INVERT: return "large_invert";
        case RearrangementType::LARGE_DELETE: return "large_delete";
        default: return "unknown";
    }
}

// [[Rcpp::export]]
List aln_rearrange(
    List store_list,
    SEXP intervals_df = R_NilValue,
    int max_margin = 10,
    int min_element_length = 50,
    int min_anchor_length = 200,
    double max_anchor_mutations_percent = 10.0,
    double max_element_mutation_percent = 100.0,
    int min_indel_length = 3,
    bool should_verify = false,
    SEXP reference_contigs = R_NilValue,
    std::string resolve_seams = "no")
{
    // validate store_list
    if (store_list.size() == 0) {
        stop("store_list must contain at least one AlignmentStore");
    }
    
    // extract library names (preserves order from R)
    CharacterVector lib_names = store_list.names();
    if (lib_names.size() != store_list.size()) {
        stop("all elements in store_list must be named");
    }
    
    // build ordered library_ids vector
    std::vector<std::string> library_ids;
    for (int i = 0; i < lib_names.size(); ++i) {
        library_ids.push_back(as<std::string>(lib_names[i]));
    }
    
    // prepare stores map - modify original stores in-place to avoid copying
    std::map<std::string, AlignmentStore> stores_map;
    for (int i = 0; i < store_list.size(); ++i) {
        std::string lib_id = as<std::string>(lib_names[i]);
        XPtr<AlignmentStore> store_ptr = store_list[i];
        
        if (!store_ptr) {
            stop("invalid AlignmentStore pointer for library: %s", lib_id.c_str());
        }
        
        // prepare the store by ensuring required indices are built
        store_ptr->count_short_indels(min_indel_length);
        
        if (!store_ptr->is_read_alignment_index_built()) {
            store_ptr->init_read_alignment_index();
        }
        
        // add to map (copying is required to build the map from individual pointers)
        stores_map[lib_id] = *store_ptr;
    }
    
    // convert intervals
    std::vector<Interval> intervals;
    if (!Rf_isNull(intervals_df)) {
        DataFrame df = as<DataFrame>(intervals_df);
        intervals = Rcpp_DataFrame_to_Intervals(df);
    }
    
    // parse resolve_seams parameter
    ResolveSeams resolve_seams_mode = string_to_resolve_seams(resolve_seams);
    
    // create sequence objects
    std::unique_ptr<AssemblySequences> assembly_sequences;
    std::unique_ptr<ReadSequences> read_sequences_obj;
    
    // load assembly sequences if needed
    if (should_verify || resolve_seams_mode == ResolveSeams::REFERENCE_ONLY || resolve_seams_mode == ResolveSeams::COMPLETE) {
        if (Rf_isNull(reference_contigs)) {
            stop("reference_contigs required for verification or reference seam resolution");
        }
        
        List contigs_list = as<List>(reference_contigs);
        if (contigs_list.size() == 0) {
            stop("reference_contigs must not be empty");
        }
        
        std::unordered_map<std::string, std::string> contig_sequences;
        if (!extract_sequences_from_list(contigs_list, contig_sequences)) {
            stop("failed to extract contig sequences from reference_contigs");
        }
        
        assembly_sequences = std::make_unique<AssemblySequences>();
        assembly_sequences->load_from_map(contig_sequences);
    }
    
    // read sequences not supported through R interface
    if (resolve_seams_mode == ResolveSeams::READS_ONLY || resolve_seams_mode == ResolveSeams::COMPLETE) {
        stop("read sequence resolution not supported through R interface - use resolve_seams='no' or 'reference_only'");
    }
    
    // create verifier if needed
    RearrangeVerify* verifier = nullptr;
    std::unique_ptr<RearrangeVerify> verify_obj;
    if (should_verify) {
        verify_obj = std::make_unique<RearrangeVerify>(assembly_sequences.get());
        verifier = verify_obj.get();
    }
    
    // create Rearrange (disable file writing for R interface, no caching)
    Rearrange manager(stores_map, library_ids, intervals, verifier, resolve_seams_mode, assembly_sequences.get(), read_sequences_obj.get(),
                     max_margin, min_element_length, min_anchor_length, max_anchor_mutations_percent, max_element_mutation_percent, false, "");
    
    // execute Rearrange
    manager.execute();
    
    // get results
    const std::vector<Event>& events = manager.get_events();
    const std::vector<ReadEventRow>& read_event_rows = manager.get_read_event_rows();
    const std::vector<std::string>& library_ids_result = manager.get_library_ids();
    const std::map<std::string, std::map<std::string, int>>& support_matrix = manager.get_support_matrix();
    const std::map<std::string, std::map<std::string, int>>& coverage_matrix = manager.get_coverage_matrix();
    
    // convert events table to R DataFrame
    CharacterVector out_event_id;
    CharacterVector out_type;
    CharacterVector out_contig;
    IntegerVector out_out_clip;
    IntegerVector out_in_clip;
    CharacterVector out_element_contig;
    CharacterVector out_element_strand;
    IntegerVector out_element_start;
    IntegerVector out_element_end;
    IntegerVector out_library_count;
    IntegerVector out_total_support;
    IntegerVector out_total_coverage;
    NumericVector out_frequency;
    CharacterVector out_read_seams;
    CharacterVector out_assembly_seams;
    
    for (const auto& event : events) {
        out_event_id.push_back(event.event_id);
        out_type.push_back(rearrangement_type_to_string(event.type));
        out_contig.push_back(event.contig_id);
        out_out_clip.push_back(event.out_clip);
        out_in_clip.push_back(event.in_clip);
        out_element_contig.push_back(event.element_contig);
        out_element_strand.push_back(event.element_strand);
        out_element_start.push_back(event.element_start);
        out_element_end.push_back(event.element_end);

        // use pre-calculated statistics from event
        out_library_count.push_back(event.library_count);
        out_total_support.push_back(event.total_support);
        out_total_coverage.push_back(event.total_coverage);
        out_frequency.push_back(event.frequency);
        out_read_seams.push_back(event.read_seams);
        out_assembly_seams.push_back(event.assembly_seams);
    }
    
    DataFrame events_df = DataFrame::create(
        Named("event_id") = out_event_id,
        Named("type") = out_type,
        Named("contig") = out_contig,
        Named("out_clip") = out_out_clip,
        Named("in_clip") = out_in_clip,
        Named("element_contig") = out_element_contig,
        Named("element_strand") = out_element_strand,
        Named("element_start") = out_element_start,
        Named("element_end") = out_element_end,
        Named("library_count") = out_library_count,
        Named("total_support") = out_total_support,
        Named("total_coverage") = out_total_coverage,
        Named("frequency") = out_frequency,
        Named("read_seams") = out_read_seams,
        Named("assembly_seams") = out_assembly_seams,
        Named("stringsAsFactors") = false);
    
    // convert read events table to R DataFrame
    CharacterVector read_out_lib_id;
    CharacterVector read_out_read_id;
    CharacterVector read_out_event_id;
    CharacterVector read_out_type;
    CharacterVector read_out_contig_id;
    CharacterVector read_out_read_strand;
    IntegerVector read_out_out_clip;
    IntegerVector read_out_in_clip;
    IntegerVector read_out_read_clip_out;
    IntegerVector read_out_read_clip_in;
    IntegerVector read_out_span_start;
    IntegerVector read_out_span_end;
    IntegerVector read_out_read_span_start;
    IntegerVector read_out_read_span_end;
    CharacterVector read_out_element_contig;
    CharacterVector read_out_element_strand;
    IntegerVector read_out_element_start;
    IntegerVector read_out_element_end;
    CharacterVector read_out_read_seams;
    CharacterVector read_out_assembly_seams;
    
    for (const ReadEventRow& row : read_event_rows) {
        read_out_lib_id.push_back(row.lib_id);
        read_out_read_id.push_back(row.read_id);
        read_out_event_id.push_back(row.event_id);
        read_out_type.push_back(rearrangement_type_to_string(row.type));
        read_out_contig_id.push_back(row.contig_id);
        read_out_read_strand.push_back(row.read_strand);
        read_out_out_clip.push_back(row.out_clip);
        read_out_in_clip.push_back(row.in_clip);
        read_out_read_clip_out.push_back(row.read_clip_out);
        read_out_read_clip_in.push_back(row.read_clip_in);
        read_out_span_start.push_back(row.span_start);
        read_out_span_end.push_back(row.span_end);
        read_out_read_span_start.push_back(row.read_span_start);
        read_out_read_span_end.push_back(row.read_span_end);
        read_out_element_contig.push_back(row.element_contig);
        read_out_element_strand.push_back(row.element_strand);
        read_out_element_start.push_back(row.element_start);
        read_out_element_end.push_back(row.element_end);
        read_out_read_seams.push_back(row.read_seams);
        read_out_assembly_seams.push_back(row.assembly_seams);
    }
    
    DataFrame read_events_df = DataFrame::create(
        Named("lib_id") = read_out_lib_id,
        Named("read_id") = read_out_read_id,
        Named("event_id") = read_out_event_id,
        Named("type") = read_out_type,
        Named("contig_id") = read_out_contig_id,
        Named("read_strand") = read_out_read_strand,
        Named("out_clip") = read_out_out_clip,
        Named("in_clip") = read_out_in_clip,
        Named("read_clip_out") = read_out_read_clip_out,
        Named("read_clip_in") = read_out_read_clip_in,
        Named("span_start") = read_out_span_start,
        Named("span_end") = read_out_span_end,
        Named("read_span_start") = read_out_read_span_start,
        Named("read_span_end") = read_out_read_span_end,
        Named("element_contig") = read_out_element_contig,
        Named("element_strand") = read_out_element_strand,
        Named("element_start") = read_out_element_start,
        Named("element_end") = read_out_element_end,
        Named("read_seams") = read_out_read_seams,
        Named("assembly_seams") = read_out_assembly_seams,
        Named("stringsAsFactors") = false);
    
    // create support matrix
    IntegerMatrix support_matrix_r(events.size(), library_ids.size());
    CharacterVector support_rownames;
    CharacterVector support_colnames;
    
    for (size_t i = 0; i < events.size(); ++i) {
        const Event& event = events[i];
        support_rownames.push_back(event.event_id);
        
        for (size_t j = 0; j < library_ids_result.size(); ++j) {
            if (i == 0) {
                support_colnames.push_back(library_ids_result[j]);
            }
            
            int support = 0;
            auto lib_it = support_matrix.find(library_ids_result[j]);
            if (lib_it != support_matrix.end()) {
                auto event_it = lib_it->second.find(event.event_id);
                if (event_it != lib_it->second.end()) {
                    support = event_it->second;
                }
            }
            support_matrix_r(i, j) = support;
        }
    }
    
    rownames(support_matrix_r) = support_rownames;
    colnames(support_matrix_r) = support_colnames;
    
    // create coverage matrix
    IntegerMatrix coverage_matrix_r(events.size(), library_ids_result.size());
    
    for (size_t i = 0; i < events.size(); ++i) {
        const Event& event = events[i];
        
        for (size_t j = 0; j < library_ids_result.size(); ++j) {
            int coverage = 0;
            auto lib_it = coverage_matrix.find(library_ids_result[j]);
            if (lib_it != coverage_matrix.end()) {
                auto event_it = lib_it->second.find(event.event_id);
                if (event_it != lib_it->second.end()) {
                    coverage = event_it->second;
                }
            }
            coverage_matrix_r(i, j) = coverage;
        }
    }
    
    rownames(coverage_matrix_r) = support_rownames;
    colnames(coverage_matrix_r) = support_colnames;
    
    // create rejection summary - simplified for now
    DataFrame rejections_df = DataFrame::create(
        Named("library_id") = CharacterVector(),
        Named("rejection_reason") = CharacterVector(),
        Named("count") = IntegerVector(),
        Named("stringsAsFactors") = false);
    
    // return as named list
    return List::create(
        Named("events") = events_df,
        Named("read_events") = read_events_df,
        Named("support") = support_matrix_r,
        Named("coverage") = coverage_matrix_r,
        Named("rejections") = rejections_df,
        Named("library_ids") = CharacterVector(library_ids.begin(), library_ids.end()));
}

// [[Rcpp::export]]
DataFrame homologs_search(List assembly_sequences,
                         const std::string& query_contig,
                         int query_start,
                         int query_end,
                         int k = 21,
                         int num_kmers = 10,
                         double threshold = 80.0,
                         int num_threads = 0) {
    
    try {
        // extract assembly sequences from List to C++ map
        std::unordered_map<std::string, std::string> sequences_map;
        if (!extract_sequences_from_list(assembly_sequences, sequences_map)) {
            stop("failed to extract assembly sequences from List");
        }
        
        // load sequences into AssemblySequences object
        AssemblySequences assembly;
        assembly.load_from_map(sequences_map);
        
        Homologs homologs(num_threads);
        std::vector<ContigRegion> regions = homologs.search_homologs(
            assembly, query_contig, 
            static_cast<uint32_t>(query_start), 
            static_cast<uint32_t>(query_end),
            static_cast<uint32_t>(k), 
            static_cast<uint32_t>(num_kmers), 
            threshold);
        
        // convert results to R vectors
        CharacterVector assembly_vec;
        CharacterVector contig_vec;
        IntegerVector start_vec;
        IntegerVector end_vec;
        CharacterVector desc_vec;
        CharacterVector id_vec;
        IntegerVector length_vec;
        IntegerVector kmer_count_vec;
        NumericVector coverage_vec;
        
        for (size_t i = 0; i < regions.size(); i++) {
            const ContigRegion& region = regions[i];
            
            assembly_vec.push_back("assembly");
            contig_vec.push_back(region.contig);
            start_vec.push_back(static_cast<int>(region.start));
            end_vec.push_back(static_cast<int>(region.end));
            length_vec.push_back(static_cast<int>(region.length));
            kmer_count_vec.push_back(static_cast<int>(region.kmer_count));
            coverage_vec.push_back(region.coverage_percent);
            
            // create description with length and coverage info
            std::string desc = "kmer_match_length_" + std::to_string(region.length) + 
                              "_coverage_" + std::to_string(static_cast<int>(region.coverage_percent)) + "pct";
            desc_vec.push_back(desc);
            
            std::string region_id = "homolog_" + std::to_string(i + 1);
            id_vec.push_back(region_id);
        }
        
        // create and return DataFrame
        return DataFrame::create(
            Named("assembly") = assembly_vec,
            Named("contig") = contig_vec,
            Named("start") = start_vec,
            Named("end") = end_vec,
            Named("desc") = desc_vec,
            Named("id") = id_vec,
            Named("length") = length_vec,
            Named("kmer_count") = kmer_count_vec,
            Named("coverage_percent") = coverage_vec,
            Named("stringsAsFactors") = false);
            
    } catch (const std::exception& e) {
        stop("homologs search failed: " + std::string(e.what()));
    }
}
