#include <iostream>
#include <string>
#include <cstring>
#include <fstream>
#include "aln_types.h"
#include "paf_reader.h"
#include "alignment_store.h"
#include "utils.h"

using std::cout;
using std::cerr;
using std::endl;
using std::ofstream;
using std::stod;

void print_usage() 
{
    cout << "Usage: alntools <command> [options]\n\n"
              << "Commands:\n"
              << "  construct <paf_file> <aln_file>    Construct alignment file from PAF\n"
              << "  info <aln_file>                    Show basic stats of alignment file\n"
              << "  save <aln_file> <output_prefix>    Save alignment file to tab-delim files\n\n"
              << "Example:\n"
              << "  alntools construct input.paf output.aln\n"
              << "  alntools info output.aln\n"
              << "  alntools save output.aln output_prefix\n";
}

void construct_alignment(const string& paf_file, const string& aln_file) 
{
    cout << "Constructing alignment file from PAF...\n";
    cout << "Reading PAF file: " << paf_file << "\n";
    
    PafReader reader;
    AlignmentStore store;
    
    cout << "Processing alignments...\n";
    if (!reader.read_paf(paf_file, store)) {
        cerr << "Error reading PAF file: " << paf_file << endl;
        return;
    }
    
    cout << "Writing alignment file: " << aln_file << "\n";
    store.save(aln_file);
    
    cout << "Done! Processed " << store.get_alignment_count() << " alignments\n";
}

void show_info(const string& aln_file) 
{
    cout << "Reading alignment file: " << aln_file << "\n";
    
    AlignmentStore store;
    store.load(aln_file);
    
    cout << "\nAlignment File Statistics:\n";
    cout << "-------------------------\n";
    cout << "Total alignments: " << store.get_alignment_count() << "\n";
    cout << "Total reads: " << store.get_read_count() << "\n";
    
    // Calculate average alignment length
    size_t total_length = 0;
    for (const auto& aln : store.get_alignments()) {
        total_length += aln.alignment_length;
    }
    double avg_length = store.get_alignment_count() > 0 ? 
        static_cast<double>(total_length) / store.get_alignment_count() : 0;
    
    cout << "Average alignment length: " << avg_length << " bp\n";
    cout << "-------------------------\n";
}

void save_to_tab_delimited(const string& aln_file, const string& output_prefix) 
{
    cout << "Reading alignment file: " << aln_file << "\n";
    
    AlignmentStore store;
    store.load(aln_file);
    
    string reads_file = output_prefix + "_reads.txt";
    string alignments_file = output_prefix + "_alignments.txt";
    
    cout << "Saving reads to: " << reads_file << "\n";
    cout << "Saving alignments to: " << alignments_file << "\n";
    
    // Save reads
    ofstream reads_out(reads_file);
    reads_out << "ReadID\tLength\tStrand\n";
    for (const auto& read : store.get_reads()) {
        reads_out << read.id << "\t" 
                 << read.length << "\t"
                 << (read.is_reverse ? "-" : "+") << "\n";
    }
    
    // Save alignments
    ofstream alns_out(alignments_file);
    alns_out << "ReadID\tTargetID\tStart\tEnd\tLength\tStrand\tScore\n";
    for (const auto& aln : store.get_alignments()) {
        alns_out << aln.query_name << "\t"
                << aln.target_name << "\t"
                << aln.query_start << "\t"
                << aln.query_end << "\t"
                << aln.alignment_length << "\t"
                << (aln.is_reverse ? "-" : "+") << "\t"
                << aln.alignment_score << "\n";
    }
    
    cout << "Done! Saved " << store.get_read_count() << " reads and "
              << store.get_alignment_count() << " alignments\n";
}

int main(int argc, char* argv[]) 
{
    if (argc < 2) {
        print_usage();
        return 1;
    }
    
    string command = argv[1];
    
    try {
        if (command == "construct") {
            massert(argc == 4, "construct command requires PAF file and output ALN file");
            construct_alignment(argv[2], argv[3]);
        }
        else if (command == "info") {
            massert(argc == 3, "info command requires ALN file");
            show_info(argv[2]);
        }
        else if (command == "save") {
            massert(argc == 4, "save command requires ALN file and output prefix");
            save_to_tab_delimited(argv[2], argv[3]);
        }
        else {
            massert(false, "Unknown command '%s'", command.c_str());
        }
    }
    catch (const std::exception& e) {
        massert(false, "Error: %s", e.what());
    }
    
    return 0;
} 