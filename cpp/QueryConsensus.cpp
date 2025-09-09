#include "QueryConsensus.h"
#include <algorithm> // For std::min/max
#include <cassert>
#include <chrono>
#include <cmath> // For std::isfinite
#include <fstream>
#include <iomanip> // For std::setprecision
#include <iostream>
#include <map>
#include <set>
#include <string>
#include <vector>
#include <thread>
#include <tuple>
#ifdef _OPENMP
#include <omp.h>
#endif

using namespace std;

QueryConsensus::QueryConsensus(
    const std::vector<Interval>& intervals,
    const AlignmentStore& store,
    double consensus_threshold,
    int min_consensus_coverage,
    int num_threads,
    ClipMode clip_mode,
    int clip_margin,
    double min_mutations_percent,
    double max_mutations_percent,
    int min_alignment_length,
    int max_alignment_length)
    : QueryBase(intervals, clip_mode, clip_margin, min_mutations_percent, max_mutations_percent, min_alignment_length, max_alignment_length)
    , store(store)
    , consensus_threshold(consensus_threshold)
    , min_consensus_coverage(min_consensus_coverage)
    , num_threads(num_threads)
{
  // validate and correct consensus_threshold
  if (consensus_threshold <= 0.0 || consensus_threshold >= 1.0 || !std::isfinite(consensus_threshold)) {
    cerr << "warning: invalid consensus_threshold (" << consensus_threshold << "), setting to default 0.9" << endl;
    this->consensus_threshold = 0.9;
  }
  
  // validate and correct min_consensus_coverage
  if (min_consensus_coverage < 1) {
    cerr << "warning: invalid min_consensus_coverage (" << min_consensus_coverage << "), setting to default 5" << endl;
    this->min_consensus_coverage = 5;
  }
  
  // set number of threads
  if (this->num_threads <= 0) {
    this->num_threads = std::thread::hardware_concurrency();
    if (this->num_threads <= 0) this->num_threads = 1;
  }
  
#ifdef _OPENMP
  omp_set_num_threads(this->num_threads);
#endif
}

void QueryConsensus::build_contig_to_intervals_map()
{
  contig_to_intervals_map.clear();
  
  for (const auto& interval : intervals) {
    // skip intervals for contigs that don't exist in the store
    if (!store.has_contig_id(interval.contig)) {
      continue;
    }
    
    uint32_t contig_index = store.get_contig_index(interval.contig);
    contig_to_intervals_map[contig_index].push_back(&interval);
  }
  
  cout << "built contig to intervals map for " << contig_to_intervals_map.size() << " contigs" << endl;
}

void QueryConsensus::merge_variant_data(const std::map<std::tuple<uint32_t, uint32_t, std::string>, VariantData>& local_data)
{
  for (const auto& entry : local_data) {
    const auto& key = entry.first;
    const auto& local_variant = entry.second;
    
    auto it = variant_results.find(key);
    if (it != variant_results.end()) {
      // merge into existing variant
      VariantData& global_variant = it->second;
      global_variant.count += local_variant.count;
      // coverage will be calculated separately
    } else {
      // create new variant entry
      variant_results[key] = local_variant;
    }
  }
}

void QueryConsensus::calc_variant_coverage()
{
  // collect all positions that need coverage calculation
  std::map<std::pair<uint32_t, uint32_t>, int> position_coverage;
  
  for (const auto& entry : variant_results) {
    uint32_t contig_index = std::get<0>(entry.first);
    uint32_t position = std::get<1>(entry.first);
    position_coverage[{contig_index, position}] = 0;
  }
  
  cout << "calculating coverage for " << position_coverage.size() << " variant positions" << endl;
  
  // calculate coverage for each position
  for (auto& pos_entry : position_coverage) {
    uint32_t contig_index = pos_entry.first.first;
    uint32_t position = pos_entry.first.second;
    
    // create a small interval around this position
    std::string contig_id = store.get_contig_id(contig_index);
    Interval pos_interval(contig_id, position, position + 1);
    auto alignment_refs = store.get_alignments_in_interval(pos_interval);
    
    int coverage = 0;
    for (const auto& alignment_ref : alignment_refs) {
      const Alignment& aln = alignment_ref.get();
      
      // apply alignment filtering
      if (!passes_filter(aln, store)) {
        continue;
      }
      
      // check if alignment covers this position
      if (aln.contig_index == contig_index && 
          aln.contig_start <= position && 
          aln.contig_end > position) {
        coverage++;
      }
    }
    
    pos_entry.second = coverage;
  }
  
  // update variant results with coverage
  for (auto& entry : variant_results) {
    uint32_t contig_index = std::get<0>(entry.first);
    uint32_t position = std::get<1>(entry.first);
    entry.second.coverage = position_coverage[{contig_index, position}];
  }
}

void QueryConsensus::process_single_alignment(const Alignment& aln, std::map<std::tuple<uint32_t, uint32_t, std::string>, VariantData>& target_variant_results)
{
  // apply alignment filtering
  if (!passes_filter(aln, store)) {
    return;
  }
  
  // get intervals for this alignment's contig only
  auto contig_intervals_it = contig_to_intervals_map.find(aln.contig_index);
  if (contig_intervals_it == contig_to_intervals_map.end()) {
    // no intervals for this contig, nothing to process
    return;
  }
  
  const std::vector<const Interval*>& relevant_intervals = contig_intervals_it->second;

  // process mutations for this alignment
  for (uint32_t mutation_index : aln.mutations) {
    // skip short indels
    if (store.is_short_indel(aln.contig_index, mutation_index)) {
      continue;
    }

    // fetch the mutation object
    const Mutation& mutation = store.get_mutation(aln.contig_index, mutation_index);
    uint32_t mutation_contig_pos = mutation.position;

    // check if this mutation falls within any of our relevant intervals
    for (const Interval* interval : relevant_intervals) {
      // check if mutation is within this interval
      if (mutation_contig_pos >= interval->start && mutation_contig_pos < interval->end) {
        // create variant key
        std::string variant_key = std::to_string(static_cast<int>(mutation.type)) + "_" + mutation.nts;
        
        auto key = std::make_tuple(aln.contig_index, mutation_contig_pos, variant_key);
        auto it = target_variant_results.find(key);
        if (it != target_variant_results.end()) {
          it->second.count++;
        } else {
          target_variant_results[key] = VariantData(1, 0); // coverage calculated later
        }
        
        break; // only count the mutation once even if it's in multiple overlapping intervals
      }
    }
  }
}

void QueryConsensus::aggregate_data()
{
  variant_results.clear();
  
  // build read-to-alignments index for LOCAL_ALIGN filtering
  if (clip_mode == ClipMode::LOCAL_ALIGN) {
    const_cast<AlignmentStore&>(store).init_read_alignment_index();
  }
  
  // build the contig to intervals mapping first
  build_contig_to_intervals_map();

  // count intervals with and without contigs
  int intervals_with_contigs = 0;
  int intervals_without_contigs = 0;

  for (const auto& interval : intervals) {
    if (!store.has_contig_id(interval.contig)) {
      intervals_without_contigs++;
    } else {
      intervals_with_contigs++;
    }
  }

  // report interval summary
  cout << "intervals summary: " << intervals_with_contigs << " with alignments, " 
       << intervals_without_contigs << " without alignments" << endl;

  // collect all unique alignments across all intervals
  std::set<const Alignment*> processed_alignments_set;
  
  for (const auto& interval : intervals) {
    // skip if contig doesn't exist
    if (!store.has_contig_id(interval.contig)) {
      continue;
    }
    
    // get alignments overlapping this interval
    std::vector<std::reference_wrapper<const Alignment>> alignments = store.get_alignments_in_interval(interval);
    
    for (const auto& alignment_ref : alignments) {
      const auto& aln = alignment_ref.get();
      processed_alignments_set.insert(&aln);
    }
  }

  // convert to vector for indexed access in OpenMP
  std::vector<const Alignment*> processed_alignments(processed_alignments_set.begin(), processed_alignments_set.end());

  // start timing
  auto start_time = std::chrono::high_resolution_clock::now();

  // calculate progress reporting thresholds (every 10%)
  size_t total_alignments = processed_alignments.size();
  size_t progress_step = std::max(static_cast<size_t>(1), total_alignments / 10 / num_threads);
  
  // process each alignment only once (parallelized)
#ifdef _OPENMP
  #pragma omp parallel
  {
    #pragma omp single
    {
      cout << "using " << omp_get_num_threads() << " threads for " << processed_alignments.size() << " alignments" << endl;
    }
    // thread-local variant results
    std::map<std::tuple<uint32_t, uint32_t, std::string>, VariantData> local_variant_results;
    
    size_t thread0_counter = 0;
    #pragma omp for
    for (size_t i = 0; i < processed_alignments.size(); ++i) {
      const Alignment* aln_ptr = processed_alignments[i];
      process_single_alignment(*aln_ptr, local_variant_results);
      
      // estimated progress reporting (only from thread 0)
      if (omp_get_thread_num() == 0) {
        ++thread0_counter;
        if (thread0_counter % progress_step == 0) {
          double estimated_progress = static_cast<double>(thread0_counter * omp_get_num_threads()) / total_alignments;
          cout << "progress: " << std::fixed << std::setprecision(1) << estimated_progress << " done" << endl;
        }
      }
    }
    
    // merge local results into global variant_results
    #pragma omp critical
    {
      merge_variant_data(local_variant_results);
    }
  }
  
  // calculate coverage once after all threads are merged
  cout << "calculating variant coverage" << endl;
  calc_variant_coverage();
#else
  cout << "processing " << processed_alignments.size() << " alignments sequentially" << endl;
  for (size_t i = 0; i < processed_alignments.size(); ++i) {
    const Alignment* aln_ptr = processed_alignments[i];
    process_single_alignment(*aln_ptr, variant_results);
    
    // progress reporting every 10%
    if ((i + 1) % progress_step == 0) {
      double progress = 100.0 * (i + 1) / total_alignments;
      cout << "progress: " << std::fixed << std::setprecision(0) << progress << "% (" << (i + 1) << "/" << total_alignments << " alignments)" << endl;
    }
  }
  
  // calculate coverage for sequential processing
  cout << "calculating variant coverage" << endl;
  calc_variant_coverage();
#endif

  // end timing and report
  auto end_time = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
  cout << "consensus analysis completed in " << duration.count() << " ms" << endl;
}

void QueryConsensus::generate_output_rows()
{
  output_rows.clear();

  cout << "filtering variants with consensus threshold " << consensus_threshold 
       << " and min coverage " << min_consensus_coverage << endl;
  
  int total_variants = 0;
  int consensus_variants = 0;
  int low_coverage_variants = 0;
  
  for (const auto& entry : variant_results) {
    const auto& key = entry.first;
    const auto& data = entry.second;

    uint32_t contig_index = std::get<0>(key);
    uint32_t position = std::get<1>(key);
    std::string variant_key = std::get<2>(key);

    total_variants++;
    
    // calculate frequency
    double frequency = (data.coverage > 0) ? static_cast<double>(data.count) / data.coverage : 0.0;
    
    // filter by coverage first
    if (data.coverage < min_consensus_coverage) {
      low_coverage_variants++;
      continue;
    }
    
    // filter by consensus threshold
    if (frequency >= consensus_threshold) {
      consensus_variants++;
      
      string contig_id = store.get_contig_id(contig_index);
      
      // parse variant key to get type and description
      size_t underscore_pos = variant_key.find('_');
      if (underscore_pos == string::npos) continue;
      
      string type_str = variant_key.substr(0, underscore_pos);
      string nts = variant_key.substr(underscore_pos + 1);
      
      // convert type number to string
      string variant_type;
      string variant_desc;
      int type_int = std::stoi(type_str);
      
      if (type_int == static_cast<int>(MutationType::SUBSTITUTION)) {
        variant_type = "SUB";
        if (nts.length() >= 2) {
          variant_desc = string(1, nts[0]) + ":" + string(1, nts[1]);
        } else {
          variant_desc = nts;
        }
      } else if (type_int == static_cast<int>(MutationType::INSERTION)) {
        variant_type = "INS";
        variant_desc = "+" + nts;
      } else if (type_int == static_cast<int>(MutationType::DELETION)) {
        variant_type = "DEL";
        variant_desc = "-" + nts;
      } else {
        variant_type = "UNK";
        variant_desc = nts;
      }

      output_rows.push_back({ contig_id, position + 1, variant_type, variant_desc,
          data.count, data.coverage, frequency });
    }
  }
  
  cout << "found " << consensus_variants << " consensus variants out of " << total_variants 
       << " total variants (" << low_coverage_variants << " filtered by low coverage)" << endl;
}

// function to write the generated rows to a file
void QueryConsensus::write_to_csv(const std::string& ofn_prefix)
{
  string filename = ofn_prefix + "_consensus.tsv";
  cout << "writing consensus data rows to " << filename << endl;
  ofstream ofs(filename);

  if (!ofs.is_open()) {
    cerr << "error: could not open file " << filename << endl;
    exit(1);
  }

  // write header with columns
  ofs << "contig\tposition\tvariant_type\tvariant_desc\tcount\tcoverage\tfrequency\n";

  for (const auto& row : output_rows) {
    ofs << row.contig << "\t"
        << row.position << "\t"
        << row.variant_type << "\t"
        << row.variant_desc << "\t"
        << row.count << "\t"
        << row.coverage << "\t"
        << std::fixed << std::setprecision(4) << row.frequency << "\n";
  }

  ofs.close();
}

void QueryConsensus::execute()
{
  aggregate_data();
  generate_output_rows();
}
