#include <cstddef>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "Params.h"
#include "alignment_store.h"
#include "aln_types.h"
#include "utils.h"

using namespace std;

void info_command(const string& aln_file)
{
  double size_mb = get_file_size_mb(aln_file);
  cout << "loading alignment file " << aln_file << " (" << size_mb << " MB)" << endl;

  AlignmentStore store;
  store.load(aln_file);

  cout << "Total alignments: " << store.get_alignment_count() << "\n";
  cout << "Total reads: " << store.get_read_count() << "\n";

  size_t total_length = 0;
  for (const auto& aln : store.get_alignments()) {
    total_length += (aln.read_end - aln.read_start);
  }
  double avg_length = store.get_alignment_count() > 0 ? static_cast<double>(total_length) / store.get_alignment_count() : 0;

  cout << "Average alignment length: " << avg_length << " bp\n";

  size_t total_mutations = 0;
  for (const auto& aln : store.get_alignments()) {
    total_mutations += aln.mutations.size();
  }

  double avg_mutations = store.get_alignment_count() > 0 ? static_cast<double>(total_mutations) / store.get_alignment_count() : 0;

  cout << "Total mutations: " << total_mutations << "\n";
  cout << "Average mutations per alignment: " << avg_mutations << "\n";
}

void print_read_info(const string& read_id, AlignmentStore& store)
{
  cout << string(60, '=') << "\n";

  if (!store.has_read_id(read_id)) {
    cout << "read not found: " << read_id << "\n";
    return;
  }

  size_t read_index = store.get_read_index(read_id);
  const Read& read = store.get_reads()[read_index];

  cout << "read: " << read.id << "\n";
  cout << "length: " << read.length << " bp\n";

  const auto& read_to_alns = store.get_read_to_alignment_indices();
  auto it = read_to_alns.find(static_cast<uint32_t>(read_index));
  if (it == read_to_alns.end()) {
    cout << "alignment count: 0\n";
    return;
  }

  const vector<size_t>& aln_indices = it->second;
  cout << "alignment count: " << aln_indices.size() << "\n";

  for (size_t i = 0; i < aln_indices.size(); ++i) {
    const Alignment& aln = store.get_alignments()[aln_indices[i]];
    const string& contig_id = store.get_contig_id(aln.contig_index);

    cout << "\nalignment " << (i + 1) << ":\n";
    cout << "  contig: " << contig_id << "\n";
    cout << "  contig coords: " << aln.contig_start << "-" << aln.contig_end << "\n";
    cout << "  read coords: " << aln.read_start << "-" << aln.read_end << "\n";
    cout << "  strand: " << (aln.is_reverse ? "-" : "+") << "\n";
    cout << "  mutation count: " << aln.mutations.size() << "\n";

    for (size_t j = 0; j < aln.mutations.size(); ++j) {
      const Mutation& mut = store.get_mutation(aln.contig_index, aln.mutations[j]);
      cout << "  mutation " << (j + 1) << ": " << mut.type
           << " pos=" << mut.position
           << " " << mut.to_string() << "\n";
    }
  }
}

void reads_info_command(const string& aln_file, const string& ids_str)
{
  double size_mb = get_file_size_mb(aln_file);
  cout << "loading alignment file " << aln_file << " (" << size_mb << " MB)" << endl;

  AlignmentStore store;
  store.load(aln_file);
  store.init_read_alignment_index();

  vector<string> read_ids;
  istringstream ss(ids_str);
  string token;
  while (getline(ss, token, ',')) {
    if (!token.empty())
      read_ids.push_back(token);
  }

  for (const string& read_id : read_ids)
    print_read_info(read_id, store);

  cout << string(60, '=') << "\n";
}

void info_params(const char* name, int argc, char** argv, Parameters& params)
{
  params.add_parser("ifn",      new ParserFilename("input ALN file"), true);
  params.add_parser("read_ids", new ParserString("comma-separated list of read IDs to inspect (optional)", ""), false);

  if (argc == 1) {
    params.usage(name);
    exit(1);
  }

  params.read(argc, argv);
  params.parse();
  params.verify_mandatory();
  params.print(cout);
}

int info_main(const char* name, int argc, char** argv)
{
  Parameters params;
  info_params(name, argc, argv, params);

  string ifn      = params.get_string("ifn");
  string read_ids = params.get_string("read_ids");

  if (!read_ids.empty())
    reads_info_command(ifn, read_ids);
  else
    info_command(ifn);

  return 0;
}
