#include "Params.h"
#include "QueryBin.h"
#include "QueryConsensus.h"
#include "QueryFull.h"
#include "QueryPileup.h"
#include "QueryVariants.h"
#include "alignment_store.h"
#include "utils.h"
#include <fstream>
#include <filesystem>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

using namespace std;

// convert string to HeightStyle enum
HeightStyle string_to_height_style(const std::string& style_str)
{
  if (style_str == "by_mutations") {
    return HeightStyle::BY_MUTATIONS;
  } else if (style_str == "by_coord_left") {
    return HeightStyle::BY_COORD_LEFT;
  } else if (style_str == "by_coord_right") {
    return HeightStyle::BY_COORD_RIGHT;
  } else {
    // default to BY_COORD_LEFT
    return HeightStyle::BY_COORD_LEFT;
  }
}

// convert string to ChunkType enum
ChunkType string_to_chunk_type(const std::string& chunk_str)
{
  if (chunk_str == "read") {
    return ChunkType::READ;
  } else if (chunk_str == "alignment") {
    return ChunkType::ALIGNMENT;
  } else if (chunk_str == "break_on_overlap") {
    return ChunkType::BREAK_ON_OVERLAP;
  } else if (chunk_str == "break_on_gap") {
    return ChunkType::BREAK_ON_GAP;
  } else {
    // default to BREAK_ON_OVERLAP
    return ChunkType::BREAK_ON_OVERLAP;
  }
}

void query_params(const char* name, int argc, char** argv, Parameters& params)
{
  params.add_parser("ifn_aln", new ParserFilename("input ALN file"), false);
  params.add_parser("ifn_intervals", new ParserFilename("input table with query contig intervals"), true);
  params.add_parser("odir", new ParserFilename("output directory"), true);
  params.add_parser("mode", new ParserString("query mode (full, pileup, bin, consensus, variants)", "full"), true);
  params.add_parser("pileup_mode", new ParserString("pileup report mode (all, covered, mutated)", "covered"), false);
  params.add_parser("binsize", new ParserInteger("bin size for 'bin' mode", 100), false);
  params.add_parser("seg_threshold", new ParserDouble("segregating sites threshold for 'bin' mode", 0.2), false);
  params.add_parser("non_ref_threshold", new ParserDouble("non-reference sites threshold for 'bin' mode", 0.9), false);
  params.add_parser("height_style", new ParserString("read height style for 'full' mode (by_coord_left, by_coord_right, by_mutations)", "by_coord_left"), false);
  params.add_parser("chunk_type", new ParserString("chunk definition style for 'full' mode (read, alignment, break_on_overlap, break_on_gap)", "break_on_overlap"), false);
  params.add_parser("max_alignments", new ParserInteger("maximum number of alignments to return per interval (0 for no limit)", 0), false);
  params.add_parser("clip_mode", new ParserString("clipping mode (all, complete, allow_one_side_clip, only_one_side_clipped, only_two_side_clipped, only_clipped)", "all"), false);
  params.add_parser("clip_margin", new ParserInteger("clipping margin in bases (default 10)", 10), false);
  params.add_parser("min_mutations_percent", new ParserDouble("minimum mutations percentage (default 0.0)", 0.0), false);
  params.add_parser("max_mutations_percent", new ParserDouble("maximum mutations percentage (default 10.0)", 10.0), false);
  params.add_parser("min_alignment_length", new ParserInteger("minimum alignment length in read coordinates (default 0)", 0), false);
  params.add_parser("max_alignment_length", new ParserInteger("maximum alignment length in read coordinates (default 0, no limit)", 0), false);
  params.add_parser("max_margin", new ParserInteger("maximum margin between alignments in same chunk for 'full' mode (default 10)", 10), false);
  params.add_parser("min_indel_length", new ParserInteger("minimum indel length to include in mutation density calculations (default 3)", 3), false);
  params.add_parser("num_threads", new ParserInteger("number of threads for 'bin' and 'consensus' modes (0 for auto)", 0), false);
  params.add_parser("consensus_threshold", new ParserDouble("consensus threshold for 'consensus' mode (default 0.9)", 0.9), false);
  params.add_parser("min_consensus_coverage", new ParserInteger("minimum coverage for 'consensus' mode (default 5)", 5), false);
  params.add_parser("ifn_libraries", new ParserFilename("input table with library definitions for 'variants' mode"), false);
  params.add_parser("min_variants_variant_support", new ParserInteger("minimum variant support across all libraries for 'variants' mode (default 3)", 3), false);
  params.add_parser("min_variants_library_support", new ParserInteger("minimum number of libraries with variant for 'variants' mode (default 1)", 1), false);
  params.add_parser("min_variants_coverage_support", new ParserInteger("minimum total coverage for 'variants' mode (default 10)", 10), false);
  params.add_parser("ifn_gene_table", new ParserFilename("input gene table file (tab-delimited with gene, contig, start, end, strand columns)"), false);
  params.add_parser("ifn_codon_table", new ParserFilename("input codon table file or name (e.g., table11)"), false);
  params.add_parser("ifn_reference_fasta", new ParserFilename("input reference FASTA file for codon analysis (required when use_genes=true)"), false);
  params.add_parser("use_genes", new ParserBoolean("use gene annotation (requires gene table, codon table, and reference FASTA)", false), false);

  if (argc == 1) {
    params.usage(name);
    exit(1);
  }

  // read command line params
  params.read(argc, argv);
  params.parse();
  params.verify_mandatory();

  // Validate mode
  string mode = params.get_string("mode");
  if (mode != "full" && mode != "pileup" && mode != "bin" && mode != "consensus" && mode != "variants") {
    cerr << "error: invalid mode specified: " << mode << ". Must be 'full', 'pileup', 'bin', 'consensus', or 'variants'." << endl;
    exit(1);
  }
  
  // Validate mode-specific parameters
  if (mode == "variants") {
    string ifn_libraries = params.get_string("ifn_libraries");
    if (ifn_libraries.empty()) {
      cerr << "error: ifn_libraries parameter is required for 'variants' mode." << endl;
      exit(1);
    }
  } else {
    // For non-variants modes, ifn_aln is required
    string ifn_aln = params.get_string("ifn_aln");
    if (ifn_aln.empty()) {
      cerr << "error: ifn_aln parameter is required for modes other than 'variants'." << endl;
      exit(1);
    }
  }

  // If mode is 'bin', binsize must be positive
  if (mode == "bin") {
    int binsize = params.get_int("binsize");
    if (binsize <= 0) {
      cerr << "error: binsize must be a positive integer for mode 'bin'." << endl;
      exit(1);
    }
  }

  // Validate height_style and chunk_type if mode is 'full'
  if (mode == "full") {
    string height_style = params.get_string("height_style");
    if (height_style != "by_coord_left" && height_style != "by_coord_right" && 
        height_style != "by_mutations") {
      cerr << "error: invalid height_style specified: " << height_style << ". Must be 'by_coord_left', 'by_coord_right', or 'by_mutations'." << endl;
      exit(1);
    }
    
    string chunk_type = params.get_string("chunk_type");
    if (chunk_type != "read" && chunk_type != "alignment" && 
        chunk_type != "break_on_overlap" && chunk_type != "break_on_gap") {
      cerr << "error: invalid chunk_type specified: " << chunk_type << ". Must be 'read', 'alignment', 'break_on_overlap', or 'break_on_gap'." << endl;
      exit(1);
    }
    
    // validate clip_mode parameter
    string clip_mode_str = params.get_string("clip_mode");
    if (clip_mode_str != "all" && clip_mode_str != "complete" && clip_mode_str != "allow_one_side_clip" && 
        clip_mode_str != "only_one_side_clipped" && clip_mode_str != "only_two_side_clipped" && 
        clip_mode_str != "only_clipped" && clip_mode_str != "local_align") {
      cerr << "error: invalid clip_mode specified: " << clip_mode_str << ". Must be 'all', 'complete', 'allow_one_side_clip', 'only_one_side_clipped', 'only_two_side_clipped', 'only_clipped', or 'local_align'." << endl;
      exit(1);
    }
  }
  
  // validate gene annotation parameters
  bool use_genes = params.get_bool("use_genes");
  if (use_genes && mode == "variants") {
    string ifn_gene_table = params.get_string("ifn_gene_table");
    string ifn_codon_table = params.get_string("ifn_codon_table");
    string ifn_reference_fasta = params.get_string("ifn_reference_fasta");
    
    if (ifn_gene_table.empty() || ifn_codon_table.empty() || ifn_reference_fasta.empty()) {
      cerr << "error: use_genes=true requires all three parameters: ifn_gene_table, ifn_codon_table, and ifn_reference_fasta" << endl;
      exit(1);
    }
  } else if (use_genes && mode != "variants") {
    cerr << "warning: gene annotation is only supported for variants mode, ignoring use_genes=true" << endl;
  }

  params.print(cout);
}

int query_main(const char* name, int argc, char** argv)
{
  Parameters params;
  query_params(name, argc, argv, params);

  string ifn_aln = params.get_string("ifn_aln");
  string ifn_intervals = params.get_string("ifn_intervals");
  string odir = params.get_string("odir");
  string mode = params.get_string("mode");
  int binsize = params.get_int("binsize"); // Will be 0 if not specified or mode is not 'bin'
  double seg_threshold = params.get_double("seg_threshold");
  double non_ref_threshold = params.get_double("non_ref_threshold");
  double consensus_threshold = params.get_double("consensus_threshold");
  int min_consensus_coverage = params.get_int("min_consensus_coverage");

  // Get pileup mode string and convert to enum
  PileupReportMode pileup_mode = string_to_pileup_report_mode(params.get_string("pileup_mode"));

  // Get height style and convert to enum
  HeightStyle height_style = string_to_height_style(params.get_string("height_style"));
  
  // Get chunk type and convert to enum
  ChunkType chunk_type = string_to_chunk_type(params.get_string("chunk_type"));
  
  // get max alignments parameter
  int max_alignments = params.get_int("max_alignments");
  
  // get alignment filtering parameters
  ClipMode clip_mode = string_to_clip_mode(params.get_string("clip_mode"));
  int clip_margin = params.get_int("clip_margin");
  double min_mutations_percent = params.get_double("min_mutations_percent");
  double max_mutations_percent = params.get_double("max_mutations_percent");
  int min_alignment_length = params.get_int("min_alignment_length");
  int max_alignment_length = params.get_int("max_alignment_length");
  int max_margin = params.get_int("max_margin");
  int min_indel_length = params.get_int("min_indel_length");
  
  // get gene annotation parameters
  string ifn_gene_table = params.get_string("ifn_gene_table");
  string ifn_codon_table = params.get_string("ifn_codon_table");
  string ifn_reference_fasta = params.get_string("ifn_reference_fasta");
  bool use_genes = params.get_bool("use_genes");

  cout << "query command called:" << endl;
  cout << "  ifn_aln: " << ifn_aln << endl;
  cout << "  ifn_intervals: " << ifn_intervals << endl;
  cout << "  odir: " << odir << endl;
  cout << "  mode: " << mode << endl;
  if (mode == "bin") {
    cout << "  binsize: " << binsize << endl;
    cout << "  seg_threshold: " << seg_threshold << endl;
    cout << "  non_ref_threshold: " << non_ref_threshold << endl;
  }
  if (mode == "consensus") {
    cout << "  consensus_threshold: " << consensus_threshold << endl;
    cout << "  min_consensus_coverage: " << min_consensus_coverage << endl;
  }
  if (mode == "pileup") {
    cout << "  pileup_mode: " << params.get_string("pileup_mode") << endl;
  }
  if (mode == "full") {
    cout << "  height_style: " << params.get_string("height_style") << endl;
    cout << "  max_alignments: " << max_alignments << endl;
  }
  // print filtering parameters for all modes
  cout << "  clip_mode: " << params.get_string("clip_mode") << endl;
  cout << "  clip_margin: " << clip_margin << endl;
  cout << "  min_mutations_percent: " << min_mutations_percent << endl;
  cout << "  max_mutations_percent: " << max_mutations_percent << endl;
  cout << "  min_alignment_length: " << min_alignment_length << endl;
  cout << "  max_alignment_length: " << max_alignment_length << endl;
  cout << "  min_indel_length: " << min_indel_length << endl;
  cout << "  use_genes: " << (use_genes ? "true" : "false") << endl;
  if (use_genes) {
    cout << "  ifn_gene_table: " << ifn_gene_table << endl;
    cout << "  ifn_codon_table: " << ifn_codon_table << endl;
    if (!ifn_reference_fasta.empty()) {
      cout << "  ifn_reference_fasta: " << ifn_reference_fasta << endl;
    }
  }

  vector<Interval> intervals;
  read_intervals(ifn_intervals, intervals);
  cout << "found " << intervals.size() << " intervals in " << ifn_intervals << endl;

  cout << "executing mode: " << mode << endl;
  
  if (mode == "variants") {
    if (odir.empty()) {
      cerr << "error: odir must be provided for variants mode" << endl;
      exit(1);
    }
    {
      std::error_code ec;
      std::filesystem::create_directories(odir, ec);
    }
    // for variants mode, load multiple ALN files based on libraries table
    string ifn_libraries = params.get_string("ifn_libraries");
    int min_variants_variant_support = params.get_int("min_variants_variant_support");
    int min_variants_library_support = params.get_int("min_variants_library_support");
    int min_variants_coverage_support = params.get_int("min_variants_coverage_support");
    
    // read libraries table
    map<string, string> library_files;
    vector<string> library_ids;
    read_library_table(ifn_libraries, library_files, library_ids);
    
    // load all ALN stores in library_ids order
    map<string, AlignmentStore> stores;
    for (const string& lib_id : library_ids) {
      auto it = library_files.find(lib_id);
      if (it != library_files.end()) {
        const string& aln_file = it->second;
        
        cout << "loading library " << lib_id << " from " << aln_file << endl;
        stores[lib_id].load(aln_file);
        stores[lib_id].count_short_indels(min_indel_length);
      }
    }
    
    // create genes object if gene annotation is requested
    Genes genes;
    Genes* genes_ptr = nullptr;
    if (use_genes) {
      cout << "loading gene annotation with reference sequences..." << endl;
      
      if (!genes.load_gene_table(ifn_gene_table)) {
        cerr << "error: failed to load gene table from " << ifn_gene_table << endl;
        exit(1);
      }
      
      if (!genes.load_codon_table(ifn_codon_table)) {
        cerr << "error: failed to load codon table from " << ifn_codon_table << endl;
        exit(1);
      }
      
      if (!genes.load_reference_sequences(ifn_reference_fasta)) {
        cerr << "error: failed to load reference sequences from " << ifn_reference_fasta << endl;
        exit(1);
      }
      
      cout << "successfully loaded genes (" << genes.get_gene_count() << " total) with reference sequences" << endl;
      genes_ptr = &genes;
    }
    
    QueryVariants queryVariants(intervals, stores, library_ids, min_variants_variant_support, min_variants_library_support, 
                     min_variants_coverage_support, clip_mode, clip_margin, min_mutations_percent, 
                     max_mutations_percent, min_alignment_length, max_alignment_length, genes_ptr);
    queryVariants.execute();
    queryVariants.write_to_csv(odir);
    
  } else {
    // for other modes, load single ALN file
    AlignmentStore store;
    store.load(ifn_aln);

    // count short indels for filtering
    store.count_short_indels(min_indel_length);
    
    if (odir.empty()) { cerr << "error: odir must be provided" << endl; exit(1); }
    {
      std::error_code ec;
      std::filesystem::create_directories(odir, ec);
    }
    if (mode == "full") {
      QueryFull queryFull(intervals, store, height_style, max_alignments, clip_mode, clip_margin, min_mutations_percent, max_mutations_percent, min_alignment_length, max_alignment_length, max_margin, chunk_type);
      queryFull.execute();
      queryFull.write_to_csv(odir);
    } else if (mode == "pileup") {
      QueryPileup queryPileup(intervals, store, pileup_mode, clip_mode, clip_margin, min_mutations_percent, max_mutations_percent, min_alignment_length, max_alignment_length);
      queryPileup.execute();
      queryPileup.write_to_csv(odir);
    } else if (mode == "bin") {
      int num_threads = params.get_int("num_threads");
      QueryBin queryBin(intervals, store, binsize, seg_threshold, non_ref_threshold, num_threads, clip_mode, clip_margin, min_mutations_percent, max_mutations_percent, min_alignment_length, max_alignment_length);
      queryBin.execute();
      queryBin.write_to_csv(odir);
    } else if (mode == "consensus") {
      int num_threads = params.get_int("num_threads");
      QueryConsensus queryConsensus(intervals, store, consensus_threshold, min_consensus_coverage, num_threads, clip_mode, clip_margin, min_mutations_percent, max_mutations_percent, min_alignment_length, max_alignment_length);
      queryConsensus.execute();
      queryConsensus.write_to_csv(odir);
    }
  }

  return 0;
}