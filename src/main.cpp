#include <iostream>
#include <string>
#include <cstring>
#include <fstream>
#include "aln_types.h"
#include "paf_reader.h"
#include "alignment_store.h"

void print_usage() 
{
    std::cout << "Usage: alntools <command> [options]\n\n"
              << "Commands:\n"
              << "  construct <paf_file> <aln_file>    Construct alignment file from PAF\n"
              << "  info <aln_file>                    Show basic stats of alignment file\n"
              << "  save <aln_file> <output_prefix>    Save alignment file to tab-delim files\n\n"
              << "Example:\n"
              << "  alntools construct input.paf output.aln\n"
              << "  alntools info output.aln\n"
              << "  alntools save output.aln output_prefix\n";
}

void construct_alignment(const std::string& paf_file, const std::string& aln_file) 
{
    std::cout << "Constructing alignment file from PAF...\n";
    std::cout << "Reading PAF file: " << paf_file << "\n";
    
    PafReader reader(paf_file);
    AlignmentStore store;
    
    std::cout << "Processing alignments...\n";
    while (auto aln = reader.next()) {
        store.add_alignment(*aln);
    }
    
    std::cout << "Writing alignment file: " << aln_file << "\n";
    store.save(aln_file);
    
    std::cout << "Done! Processed " << store.get_alignment_count() << " alignments\n";
}

void show_info(const std::string& aln_file) 
{
    std::cout << "Reading alignment file: " << aln_file << "\n";
    
    AlignmentStore store;
    store.load(aln_file);
    
    std::cout << "\nAlignment File Statistics:\n";
    std::cout << "-------------------------\n";
    std::cout << "Total alignments: " << store.get_alignment_count() << "\n";
    std::cout << "Total reads: " << store.get_read_count() << "\n";
    
    // Calculate average alignment length
    size_t total_length = 0;
    for (const auto& aln : store.get_alignments()) {
        total_length += aln.alignment_length;
    }
    double avg_length = store.get_alignment_count() > 0 ? 
        static_cast<double>(total_length) / store.get_alignment_count() : 0;
    
    std::cout << "Average alignment length: " << avg_length << " bp\n";
    std::cout << "-------------------------\n";
}

void save_to_tab_delimited(const std::string& aln_file, const std::string& output_prefix) 
{
    std::cout << "Reading alignment file: " << aln_file << "\n";
    
    AlignmentStore store;
    store.load(aln_file);
    
    std::string reads_file = output_prefix + "_reads.txt";
    std::string alignments_file = output_prefix + "_alignments.txt";
    
    std::cout << "Saving reads to: " << reads_file << "\n";
    std::cout << "Saving alignments to: " << alignments_file << "\n";
    
    // Save reads
    std::ofstream reads_out(reads_file);
    reads_out << "ReadID\tLength\tStrand\n";
    for (const auto& read : store.get_reads()) {
        reads_out << read.id << "\t" 
                 << read.length << "\t"
                 << (read.is_reverse ? "-" : "+") << "\n";
    }
    
    // Save alignments
    std::ofstream alns_out(alignments_file);
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
    
    std::cout << "Done! Saved " << store.get_read_count() << " reads and "
              << store.get_alignment_count() << " alignments\n";
}

int main(int argc, char* argv[]) 
{
    if (argc < 2) {
        print_usage();
        return 1;
    }
    
    std::string command = argv[1];
    
    try {
        if (command == "construct") {
            if (argc != 4) {
                std::cerr << "Error: construct command requires PAF file and output ALN file\n";
                print_usage();
                return 1;
            }
            construct_alignment(argv[2], argv[3]);
        }
        else if (command == "info") {
            if (argc != 3) {
                std::cerr << "Error: info command requires ALN file\n";
                print_usage();
                return 1;
            }
            show_info(argv[2]);
        }
        else if (command == "save") {
            if (argc != 4) {
                std::cerr << "Error: save command requires ALN file and output prefix\n";
                print_usage();
                return 1;
            }
            save_to_tab_delimited(argv[2], argv[3]);
        }
        else {
            std::cerr << "Error: Unknown command '" << command << "'\n";
            print_usage();
            return 1;
        }
    }
    catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }
    
    return 0;
} 