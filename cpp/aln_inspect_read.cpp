#include <cstddef>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

#include "Params.h"
#include "alignment_store.h"
#include "aln_types.h"
#include "utils.h"

using namespace std;

void inspect_read_command(const string& aln_file, const string& read_id)
{
  double size_mb = get_file_size_mb(aln_file);
  cout << "loading alignment file " << aln_file << " (" << size_mb << " MB)" << endl;

  AlignmentStore store;
  store.load(aln_file);

  if (!store.has_read_id(read_id)) {
    cerr << "error: read ID '" << read_id << "' not found in alignment file" << endl;
    exit(1);
  }

  size_t read_index = store.get_read_index(read_id);
  const Read& read = store.get_reads()[read_index];

  store.init_read_alignment_index();
  const auto& read_to_alignment_indices = store.get_read_to_alignment_indices();
  auto it = read_to_alignment_indices.find(static_cast<uint32_t>(read_index));
  
  vector<size_t> alignment_indices;
  if (it != read_to_alignment_indices.end()) {
    alignment_indices = it->second;
  }

  size_t total_mutations = 0;
  for (size_t aln_idx : alignment_indices) {
    const Alignment& aln = store.get_alignments()[aln_idx];
    total_mutations += aln.mutations.size();
  }

  cout << "\n=== Read Summary ===" << endl;
  cout << "Read ID: " << read.id << endl;
  cout << "Length: " << read.length << " bp" << endl;
  cout << "Number of alignments: " << alignment_indices.size() << endl;
  cout << "Total mutations: " << total_mutations << endl;

  if (alignment_indices.empty()) {
    cout << "\nno alignments found for this read" << endl;
    return;
  }

  cout << "\n=== Alignments ===" << endl;
  cout << "alignment_index\tcontig\tread_start\tread_end\tcontig_start\tcontig_end\tis_reverse\tmutations" << endl;
  
  for (size_t aln_idx : alignment_indices) {
    const Alignment& aln = store.get_alignments()[aln_idx];
    const string& contig_id = store.get_contig_id(aln.contig_index);
    
    cout << aln_idx << "\t"
         << contig_id << "\t"
         << aln.read_start << "\t"
         << aln.read_end << "\t"
         << aln.contig_start << "\t"
         << aln.contig_end << "\t"
         << (aln.is_reverse ? "true" : "false") << "\t"
         << aln.mutations.size() << endl;
  }

  cout << "\n=== Mutations ===" << endl;
  cout << "alignment_index\tcontig\tmutation_index\ttype\tposition\tdescription" << endl;
  
  for (size_t aln_idx : alignment_indices) {
    const Alignment& aln = store.get_alignments()[aln_idx];
    const string& contig_id = store.get_contig_id(aln.contig_index);
    
    for (size_t mut_idx = 0; mut_idx < aln.mutations.size(); ++mut_idx) {
      uint32_t mutation_index = aln.mutations[mut_idx];
      const Mutation& mutation = store.get_mutation(aln.contig_index, mutation_index);
      
      string type_str;
      switch (mutation.type) {
      case MutationType::SUBSTITUTION:
        type_str = "SUB";
        break;
      case MutationType::INSERTION:
        type_str = "INS";
        break;
      case MutationType::DELETION:
        type_str = "DEL";
        break;
      }
      
      cout << aln_idx << "\t"
           << contig_id << "\t"
           << mutation_index << "\t"
           << type_str << "\t"
           << mutation.position << "\t"
           << mutation.to_string() << endl;
    }
  }
}

void inspect_read_params(const char* name, int argc, char** argv, Parameters& params)
{
  params.add_parser("ifn_aln", new ParserFilename("input ALN file"), true);
  params.add_parser("read_id", new ParserString("read ID to inspect"), true);

  if (argc == 1) {
    params.usage(name);
    exit(1);
  }

  params.read(argc, argv);
  params.parse();
  params.verify_mandatory();
  params.print(cout);
}

int inspect_read_main(const char* name, int argc, char** argv)
{
  Parameters params;
  inspect_read_params(name, argc, argv, params);

  string ifn_aln = params.get_string("ifn_aln");
  string read_id = params.get_string("read_id");

  inspect_read_command(ifn_aln, read_id);

  return 0;
}

