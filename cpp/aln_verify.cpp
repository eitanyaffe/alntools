#include <iostream>
#include <fstream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include "alignment_store.h"
#include "Params.h"

using namespace std;

// Function to read FASTA file and create a map of contig IDs to sequences
unordered_map<string, string> read_fasta(const string& filename,
					 const unordered_set<string>& contig_ids)
{
    cout << "Reading FASTA file: " << filename << endl;
    unordered_map<string, string> contigs;
    ifstream file(filename);
    string line, id, sequence;
    while (getline(file, line)) {
        if (line[0] == '>') {
            if (!id.empty() && contig_ids.find(id) != contig_ids.end()) {
                contigs[id] = sequence;
            }
            id = line.substr(1); // Remove '>'
            sequence.clear();
        } else {
            sequence += line;
        }
    }
    if (!id.empty() && contig_ids.find(id) != contig_ids.end()) {
        contigs[id] = sequence;
    }
    return contigs;
}

// Function to read FASTQ file and create a map of read IDs to sequences
unordered_map<string, string> read_fastq(const string& filename,
					 const unordered_set<string>& read_ids)
{
    cout << "Reading FASTQ file: " << filename << endl;
    unordered_map<string, string> reads;
    ifstream file(filename);
    string line, id, sequence;
    while (getline(file, line)) {
        if (line[0] == '@') {
            id = line;
            getline(file, sequence); // Read the sequence line
            if (read_ids.find(id) != read_ids.end()) {
                reads[id] = sequence;
            }
            getline(file, line); // Skip the '+' line
            getline(file, line); // Skip the quality line
        }
    }
    return reads;
}

// Function to apply mutations to a contig fragment
string apply_mutations(const string& contig_fragment,
		       const vector<Mutation>& mutations,
		       uint32_t contig_start)
{
    string mutated_fragment = contig_fragment;
    int offset = 0; // Track changes in length due to insertions/deletions

    for (const auto& mutation : mutations) {
        int adjusted_position = mutation.position - contig_start + offset;
        switch (mutation.type) {
            case MutationType::SUBSTITUTION:
                mutated_fragment[adjusted_position] = mutation.query_bases[0];
                break;
            case MutationType::INSERTION:
                mutated_fragment.insert(adjusted_position, mutation.query_bases);
                offset += mutation.query_bases.size();
                break;
            case MutationType::DELETION:
                mutated_fragment.erase(adjusted_position, mutation.target_bases.size());
                offset -= mutation.target_bases.size();
                break;
        }
    }
    return mutated_fragment;
}

// Main verification function
void verify_command(const string& ifn_aln,
		    const string& ifn_reads,
		    const string& ifn_contigs)
{
    AlignmentStore store;
    cout << "Reading alignment file: " << ifn_aln << "\n";
    store.load(ifn_aln);

    // Collect contig and read IDs from the store
    vector<string> contig_ids, read_ids;
    for (const auto& alignment : store.get_alignments()) {
        contig_ids.push_back(store.get_contig_id(alignment.contig_index));
        read_ids.push_back(store.get_read_id(alignment.read_index));
    }

    unordered_set<string> contig_set(contig_ids.begin(), contig_ids.end());
    unordered_set<string> read_set(read_ids.begin(), read_ids.end());

    auto contigs = read_fasta(ifn_contigs, contig_set);
    auto reads = read_fastq(ifn_reads, read_set);

    int bad_alignment_count = 0;
    for (const auto& alignment : store.get_alignments()) {
        const string& contig_id = store.get_contig_id(alignment.contig_index);
        const string& read_id = store.get_read_id(alignment.read_index);

        if (contigs.find(contig_id) == contigs.end() || reads.find(read_id) == reads.end()) {
            cerr << "Error: Contig or read not found for alignment.\n";
            abort();
        }

        string contig_fragment = contigs[contig_id].substr(alignment.contig_start,
							   alignment.contig_end - alignment.contig_start);
        string mutated_contig = apply_mutations(contig_fragment, alignment.mutations, alignment.contig_start);
        string read_segment = reads[read_id].substr(alignment.read_start,
						    alignment.read_end - alignment.read_start);

        bool mismatch_found = false;
        for (size_t i = 0; i < read_segment.size(); ++i) {
            if (mutated_contig[i] != read_segment[i]) {
                cout << "Mismatch in alignment:\n";
                cout << "Expected " << read_segment[i] << " at position " <<
		  i << " in read but found " << mutated_contig[i] << "\n";
                cout << "Mutations: ";
                for (const auto& mutation : alignment.mutations) {
                    cout << mutation.type << " at " << mutation.position << "; ";
                }
                cout << "\n";
                mismatch_found = true;
                break;
            }
        }

        if (mismatch_found) {
            bad_alignment_count++;
            if (bad_alignment_count >= 10) {
                cerr << "Too many bad alignments. Exiting.\n";
                exit(-1);
            }
        }
    }

    cout << "Verification complete. Total alignments processed: " << store.get_alignment_count() << "\n";
} 

void verify_params(const char* name, int argc, char **argv, Parameters& params)
{
  params.add_parser("ifn_aln", new ParserFilename("input ALN file"), true);
  params.add_parser("ifn_reads", new ParserFilename("input reads, FASTQ"), true);
  params.add_parser("ifn_contigs", new ParserFilename("input contigs, FASTA"), true);
  
  if (argc == 1) {
    params.usage(name);
    exit(1);
  }

  // read command line params
  params.read(argc, argv);
  params.parse();
  params.verify_mandatory();
  params.print(cout);
}

int verify_main(const char* name, int argc, char **argv)
{
  Parameters params;
  verify_params(name, argc, argv, params);
  
  string ifn_aln = params.get_string("ifn_aln");
  string ifn_reads = params.get_string("ifn_reads");
  string ifn_contigs = params.get_string("ifn_contigs");

  verify_command(ifn_aln, ifn_reads, ifn_contigs);
  
  return 0;
}
