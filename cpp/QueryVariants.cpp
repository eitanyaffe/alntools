#include "QueryVariants.h"
#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>

using namespace std;

// Helper function to create human-readable mutation description
std::string create_desc(const std::string& type, const std::string& sequence) {
  if (type == "sub") {
    // substitution: sequence like "AG" means A→G, format as "A:G"
    if (sequence.length() == 2) {
      return sequence.substr(0, 1) + ":" + sequence.substr(1, 1);
    } else {
      return "SUB"; // fallback
    }
  } else if (type == "ins") {
    // insertion: sequence contains inserted nucleotides, format as "+XXX"
    return "+" + sequence;
  } else if (type == "del") {
    // deletion: sequence contains deleted nucleotides, format as "-XXX"
    return "-" + sequence;
  } else if (type == "left_clip" || type == "right_clip") {
    // clips have no sequence, return clip identifier
    return "CLIP";
  } else {
    return "UNK";
  }
}

std::string SetVariantData::create_key() const
{
  return contig + ":" + std::to_string(position) + ":" + type + ":" + sequence;
}

QueryVariants::QueryVariants(
    const std::vector<Interval>& intervals,
    const std::map<std::string, AlignmentStore>& stores,
    int min_variants_variant_support,
    int min_variants_library_support,
    int min_variants_coverage_support,
    ClipMode clip_mode,
    int clip_margin,
    double min_mutations_percent,
    double max_mutations_percent,
    int min_alignment_length,
    int max_alignment_length,
    Genes* genes)
    : QueryBase(intervals, clip_mode, clip_margin, min_mutations_percent, max_mutations_percent, 
                min_alignment_length, max_alignment_length)
    , stores(stores)
    , min_variants_variant_support(min_variants_variant_support)
    , min_variants_library_support(min_variants_library_support)
    , min_variants_coverage_support(min_variants_coverage_support)
    , genes(genes)
{
  // extract ordered library IDs
  for (const auto& entry : stores) {
    library_ids.push_back(entry.first);
  }
  std::sort(library_ids.begin(), library_ids.end());
}

void QueryVariants::build_contig_to_intervals_map(const AlignmentStore& store)
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

void QueryVariants::execute()
{
  cout << "executing QueryVariants with " << stores.size() << " libraries and " 
       << intervals.size() << " intervals" << endl;
  
  // process each library
  for (const auto& store_entry : stores) {
    const string& lib_id = store_entry.first;
    const AlignmentStore& store = store_entry.second;
    
    cout << "processing library: " << lib_id << endl;
    
    // build read-to-alignments index for LOCAL_ALIGN filtering if needed
    init_local_align_if_needed(const_cast<AlignmentStore&>(store), clip_mode);
    
    // build the contig to intervals mapping for this store
    build_contig_to_intervals_map(store);
    
    // process each interval
    for (const Interval& interval : intervals) {
      // skip intervals for contigs that don't exist in this store
      if (!store.has_contig_id(interval.contig)) {
        continue;
      }
      
      // get alignments overlapping this interval
      auto alignment_refs = store.get_alignments_intersecting_interval(interval);
      
      for (const auto& alignment_ref : alignment_refs) {
        const Alignment& aln = alignment_ref.get();
        process_alignment(aln, store, lib_id);
      }
    }
  }
  
  cout << "found " << variant_data.size() << " unique variants before filtering" << endl;
  
  // calculate coverage for all variant positions
  calculate_coverage();
  
  // apply filters and generate output rows
  apply_filters_and_generate_rows();
  
  cout << "after filtering: " << variant_rows.size() << " variants remain" << endl;
}

void QueryVariants::process_alignment(const Alignment& aln, const AlignmentStore& store, 
                                const std::string& lib_id)
{
  // apply alignment filtering first
  if (!passes_filter(aln, store)) {
    return;
  }
  
  // process mutations from this alignment
  process_mutations(aln, store, lib_id);
  
  // process clips from this alignment
  process_clips(aln, store, lib_id);
}

void QueryVariants::process_mutations(const Alignment& aln, const AlignmentStore& store, 
                                const std::string& lib_id)
{
  // process each mutation in this alignment
  for (uint32_t mutation_index : aln.mutations) {
    // skip short indels
    if (should_skip_short_indel(store, aln.contig_index, mutation_index)) {
      continue;
    }
    
    // fetch the mutation object
    const Mutation& mutation = store.get_mutation(aln.contig_index, mutation_index);
    uint32_t mutation_pos = mutation.position;
    
    // get intervals for this alignment's contig only
    auto contig_intervals_it = contig_to_intervals_map.find(aln.contig_index);
    if (contig_intervals_it == contig_to_intervals_map.end()) {
      // no intervals for this contig, skip this mutation
      continue;
    }
    
    const std::vector<const Interval*>& relevant_intervals = contig_intervals_it->second;
    
    // check if mutation falls within any of the relevant intervals
    bool in_interval = false;
    for (const Interval* interval : relevant_intervals) {
      if (mutation_pos >= interval->start && mutation_pos < interval->end) {
        in_interval = true;
        break;
      }
    }
    
    if (!in_interval) continue;
    
    // create variant data
    SetVariantData variant;
    variant.contig = store.get_contig_id(aln.contig_index);
    variant.position = mutation_pos + 1; // convert to 1-based
    variant.sequence = mutation.nts;
    
    // set type based on mutation type
    switch (mutation.type) {
      case MutationType::SUBSTITUTION:
        variant.type = "sub";
        break;
      case MutationType::INSERTION:
        variant.type = "ins";
        break;
      case MutationType::DELETION:
        variant.type = "del";
        break;
      default:
        continue; // skip unknown types
    }
    
    // add to variant data
    string variant_key = variant.create_key();
    variant_data[variant_key].contig = variant.contig;
    variant_data[variant_key].position = variant.position;
    variant_data[variant_key].type = variant.type;
    variant_data[variant_key].sequence = variant.sequence;
    variant_data[variant_key].lib_support[lib_id]++;
  }
}

void QueryVariants::process_single_clip(uint32_t clip_pos, const std::string& clip_type, 
                                  const AlignmentStore& store, const Alignment& aln, 
                                  const std::string& lib_id)
{
  // get intervals for this alignment's contig only
  auto contig_intervals_it = contig_to_intervals_map.find(aln.contig_index);
  if (contig_intervals_it == contig_to_intervals_map.end()) {
    // no intervals for this contig, nothing to process
    return;
  }
  
  const std::vector<const Interval*>& relevant_intervals = contig_intervals_it->second;
  
  // check if clip position is within any of the relevant intervals
  bool in_interval = false;
  for (const Interval* interval : relevant_intervals) {
    if (clip_pos >= interval->start && clip_pos < interval->end) {
      in_interval = true;
      break;
    }
  }
  
  if (in_interval) {
    SetVariantData variant;
    variant.contig = store.get_contig_id(aln.contig_index);
    variant.position = clip_pos + 1; // convert to 1-based
    variant.type = clip_type;
    variant.sequence = "";
    
    string variant_key = variant.create_key();
    variant_data[variant_key].contig = variant.contig;
    variant_data[variant_key].position = variant.position;
    variant_data[variant_key].type = variant.type;
    variant_data[variant_key].sequence = variant.sequence;
    variant_data[variant_key].lib_support[lib_id]++;
  }
}

void QueryVariants::process_clips(const Alignment& aln, const AlignmentStore& store, 
                            const std::string& lib_id)
{
  uint32_t read_length = store.get_reads()[aln.read_index].length;
  bool is_left_clipped = (aln.read_start > static_cast<uint32_t>(clip_margin));
  bool is_right_clipped = (aln.read_end < (read_length - static_cast<uint32_t>(clip_margin)));
  
  // process left clip
  if (is_left_clipped) {
    process_single_clip(aln.contig_start, "left_clip", store, aln, lib_id);
  }
  
  // process right clip
  if (is_right_clipped) {
    process_single_clip(aln.contig_end, "right_clip", store, aln, lib_id);
  }
}

void QueryVariants::calculate_coverage()
{
  cout << "calculating coverage for " << variant_data.size() << " variant positions" << endl;
  
  // for each variant position, calculate coverage across all libraries
  for (auto& variant_entry : variant_data) {
    SetVariantData& variant = variant_entry.second;
    
    // process each library
    for (const auto& store_entry : stores) {
      const string& lib_id = store_entry.first;
      const AlignmentStore& store = store_entry.second;
      
      // skip if contig doesn't exist in this store
      if (!store.has_contig_id(variant.contig)) {
        variant.lib_coverage[lib_id] = 0;
        continue;
      }
      
      uint32_t pos_0based = variant.position - 1; // convert to 0-based
      
      // create interval for this position
      Interval pos_interval(variant.contig, pos_0based, pos_0based + 1);
      auto alignment_refs = store.get_alignments_intersecting_interval(pos_interval);
      
      int coverage = 0;
      for (const auto& alignment_ref : alignment_refs) {
        const Alignment& aln = alignment_ref.get();
        
        // apply alignment filtering
        if (!passes_filter(aln, store)) {
          continue;
        }
        
        // check if alignment covers this position
        if (aln.contig_start <= pos_0based && aln.contig_end > pos_0based) {
          coverage++;
        }
      }
      
      variant.lib_coverage[lib_id] = coverage;
    }
  }
}

void QueryVariants::apply_filters_and_generate_rows()
{
  cout << "applying filters: min_variant_support=" << min_variants_variant_support 
       << ", min_library_support=" << min_variants_library_support 
       << ", min_coverage_support=" << min_variants_coverage_support << endl;
  
  // collect variants that pass filters
  vector<pair<string, SetVariantData*>> passing_variants;
  
  for (auto& variant_entry : variant_data) {
    SetVariantData& variant = variant_entry.second;
    
    // calculate totals
    variant.total_support = 0;
    variant.total_coverage = 0;
    variant.library_count = 0;
    
    for (const string& lib_id : library_ids) {
      int support = variant.lib_support[lib_id];
      int coverage = variant.lib_coverage[lib_id];
      
      variant.total_support += support;
      variant.total_coverage += coverage;
      
      if (support > 0) {
        variant.library_count++;
      }
    }
    
    // apply filters
    if (variant.total_support >= min_variants_variant_support &&
        variant.library_count >= min_variants_library_support &&
        variant.total_coverage >= min_variants_coverage_support) {
      passing_variants.push_back(make_pair(variant_entry.first, &variant));
    }
  }
  
  cout << passing_variants.size() << " variants passed filters" << endl;
  
  // sort by contig, position, then coverage (descending)
  sort(passing_variants.begin(), passing_variants.end(), 
       [](const pair<string, SetVariantData*>& a, const pair<string, SetVariantData*>& b) {
         if (a.second->contig != b.second->contig) {
           return a.second->contig < b.second->contig;
         }
         if (a.second->position != b.second->position) {
           return a.second->position < b.second->position;
         }
         return a.second->total_coverage > b.second->total_coverage;
       });
  
  // generate output rows with sequential variant IDs
  variant_rows.clear();
  variant_rows.reserve(passing_variants.size());
  
  for (size_t i = 0; i < passing_variants.size(); ++i) {
    const SetVariantData* variant = passing_variants[i].second;
    
    VariantOutputRow row;
    row.variant_id = "v" + std::to_string(i + 1);
    row.contig = variant->contig;
    row.coord = variant->position;
    row.type = variant->type;
    row.sequence = variant->sequence;
    row.desc = create_desc(variant->type, variant->sequence);
    row.library_count = variant->library_count;
    row.total_support = variant->total_support;
    row.total_coverage = variant->total_coverage;
    row.frequency = variant->get_frequency();
    row.is_genic = false;

    // add gene annotation if genes object is available
    if (genes != nullptr) {
      // annotate_position now returns whether the variant is genic
      row.is_genic = annotate_position(row.variant_id, variant->contig, variant->position, 
                                      variant->type, variant->sequence);
    }
    
    variant_rows.push_back(row);
  }
}

void QueryVariants::write_to_csv(const std::string& ofn_prefix)
{
  write_variants_file(ofn_prefix);
  write_support_file(ofn_prefix);
  write_coverage_file(ofn_prefix);
  
  // write gene annotation files if genes were used
  if (genes != nullptr) {
    write_genic_file(ofn_prefix);
    write_intergenic_file(ofn_prefix);
  }
}

void QueryVariants::write_variants_file(const std::string& ofn_prefix)
{
  string filename = ofn_prefix + "_variants.tsv";
  cout << "writing variants data to " << filename << endl;
  
  ofstream file;
  safe_open_file_for_writing(filename, file);
  
  // write header
  file << "variant_id\tcontig\tcoord\ttype\tsequence\tdesc\tlibrary_count\ttotal_support\ttotal_coverage\tfrequency\tis_genic" << endl;
  
  // write data rows
  for (const auto& row : variant_rows) {
    file << row.variant_id << "\t"
         << row.contig << "\t"
         << row.coord << "\t"
         << row.type << "\t"
         << row.sequence << "\t"
         << row.desc << "\t"
         << row.library_count << "\t"
         << row.total_support << "\t"
         << row.total_coverage << "\t"
         << fixed << setprecision(6) << row.frequency << "\t"
         << (row.is_genic ? "true" : "false") << endl;
  }
  
  file.close();
}

void QueryVariants::write_support_file(const std::string& ofn_prefix)
{
  string filename = ofn_prefix + "_support.tsv";
  cout << "writing support data to " << filename << endl;
  
  ofstream file;
  safe_open_file_for_writing(filename, file);
  
  // write header
  file << "variant_id";
  for (const string& lib_id : library_ids) {
    file << "\t" << lib_id;
  }
  file << endl;
  
  // write data rows
  for (size_t i = 0; i < variant_rows.size(); ++i) {
    const VariantOutputRow& row = variant_rows[i];
    file << row.variant_id;
    
    // find the original variant data
    string variant_key = row.contig + ":" + std::to_string(row.coord) + ":" + row.type + ":" + row.sequence;
    auto variant_it = variant_data.find(variant_key);
    
    if (variant_it != variant_data.end()) {
      const SetVariantData& variant = variant_it->second;
      for (const string& lib_name : library_ids) {
        auto support_it = variant.lib_support.find(lib_name);
        int support = (support_it != variant.lib_support.end()) ? support_it->second : 0;
        file << "\t" << support;
      }
    } else {
      // fallback: write zeros
      for (size_t j = 0; j < library_ids.size(); ++j) {
        file << "\t0";
      }
    }
    file << endl;
  }
  
  file.close();
}

void QueryVariants::write_coverage_file(const std::string& ofn_prefix)
{
  string filename = ofn_prefix + "_coverage.tsv";
  cout << "writing coverage data to " << filename << endl;
  
  ofstream file;
  safe_open_file_for_writing(filename, file);
  
  // write header
  file << "variant_id";
  for (const string& lib_id : library_ids) {
    file << "\t" << lib_id;
  }
  file << endl;
  
  // write data rows
  for (size_t i = 0; i < variant_rows.size(); ++i) {
    const VariantOutputRow& row = variant_rows[i];
    file << row.variant_id;
    
    // find the original variant data
    string variant_key = row.contig + ":" + std::to_string(row.coord) + ":" + row.type + ":" + row.sequence;
    auto variant_it = variant_data.find(variant_key);
    
    if (variant_it != variant_data.end()) {
      const SetVariantData& variant = variant_it->second;
      for (const string& lib_id : library_ids) {
        auto coverage_it = variant.lib_coverage.find(lib_id);
        int coverage = (coverage_it != variant.lib_coverage.end()) ? coverage_it->second : 0;
        file << "\t" << coverage;
      }
    } else {
      // fallback: write zeros
      for (size_t j = 0; j < library_ids.size(); ++j) {
        file << "\t0";
      }
    }
    file << endl;
  }
  
  file.close();
}


bool QueryVariants::annotate_position(const std::string& variant_id, const std::string& contig, uint32_t position,
                                     const std::string& mutation_type, const std::string& mutation_sequence)
{
  if (genes == nullptr) {
    return false; // no annotation needed, not genic
  }
  
  GeneAnnotation annotation = genes->annotate_variant(contig, position, mutation_type, mutation_sequence);
  
  if (annotation.loc == "genic") {
    GenicRow genic_row;
    genic_row.row_id = variant_id;
    genic_row.gene_id = annotation.gene_id;
    genic_row.gene_desc = annotation.gene_desc;
    genic_row.aa_coord = annotation.aa_coord;
    genic_row.variant_codon = annotation.variant_codon;
    genic_row.ref_codon = annotation.ref_codon;
    genic_row.variant_type = annotation.variant_type;
    genic_row.mutation_desc = annotation.mutation_desc;
    genic_rows.push_back(genic_row);
    return true; // variant is genic
  } else {
    IntergenicRow intergenic_row;
    intergenic_row.row_id = variant_id;
    intergenic_row.gene_left = annotation.gene_left;
    intergenic_row.gene_right = annotation.gene_right;
    intergenic_row.orientation_left = annotation.orientation_left;
    intergenic_row.orientation_right = annotation.orientation_right;
    intergenic_row.distance_left = annotation.distance_left;
    intergenic_row.distance_right = annotation.distance_right;
    intergenic_rows.push_back(intergenic_row);
    return false; // variant is intergenic
  }
}

void QueryVariants::write_genic_file(const std::string& ofn_prefix, const std::string& suffix)
{
  if (genic_rows.empty()) {
    cout << "no genic entries to write" << endl;
    return;
  }
  
  string filename = ofn_prefix + "_genic" + suffix + ".tsv";
  cout << "writing genic data to " << filename << endl;
  
  ofstream file;
  safe_open_file_for_writing(filename, file);
  
  // write header
  file << "row_id\tgene_id\tgene_desc\taa_coord\tvariant_codon\tref_codon\tvariant_type\tmutation_desc" << endl;
  
  // write data rows
  for (const auto& row : genic_rows) {
    file << row.row_id << "\t"
         << row.gene_id << "\t"
         << row.gene_desc << "\t"
         << row.aa_coord << "\t"
         << row.variant_codon << "\t"
         << row.ref_codon << "\t"
         << row.variant_type << "\t"
         << row.mutation_desc << endl;
  }
  
  file.close();
}

void QueryVariants::write_intergenic_file(const std::string& ofn_prefix, const std::string& suffix)
{
  if (intergenic_rows.empty()) {
    cout << "no intergenic entries to write" << endl;
    return;
  }
  
  string filename = ofn_prefix + "_intergenic" + suffix + ".tsv";
  cout << "writing intergenic data to " << filename << endl;
  
  ofstream file;
  safe_open_file_for_writing(filename, file);
  
  // write header
  file << "row_id\tgene_left\tgene_right\torientation_left\torientation_right\tdistance_left\tdistance_right" << endl;
  
  // write data rows
  for (const auto& row : intergenic_rows) {
    file << row.row_id << "\t"
         << row.gene_left << "\t"
         << row.gene_right << "\t"
         << row.orientation_left << "\t"
         << row.orientation_right << "\t"
         << row.distance_left << "\t"
         << row.distance_right << endl;
  }
  
  file.close();
}
