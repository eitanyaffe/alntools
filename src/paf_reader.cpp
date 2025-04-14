#include "../include/paf_reader.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <unordered_map>

PafReader::PafReader(const std::string& filename) : filename_(filename) 
{
    // Empty implementation for now
}

std::unique_ptr<Alignment> PafReader::next() 
{
    // Empty implementation for now
    return nullptr;
}

bool PafReader::read_paf(const std::string& filename, AlignmentStore& store) 
{
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Failed to open file: " << filename << std::endl;
        return false;
    }

    // Maps for tracking existing contigs and reads
    std::unordered_map<std::string, uint32_t> contig_indices;
    std::unordered_map<std::string, uint32_t> read_indices;

    std::string line;
    while (std::getline(file, line)) {
        std::vector<std::string> fields;
        split_line(line, '\t', fields);

        if (fields.size() < 12) {
            std::cerr << "Warning: Skipping malformed line with fewer than 12 fields\n";
            continue;
        }

        // Get read information
        std::string read_id = fields[0];
        uint64_t read_length = std::stoull(fields[1]);
        
        // Get contig information
        std::string contig_name = fields[5];
        uint64_t contig_length = std::stoull(fields[6]);
        
        // Get alignment information
        uint64_t read_start = std::stoull(fields[2]);
        uint64_t read_end = std::stoull(fields[3]);
        bool is_reverse = (fields[4] == "-");
        uint64_t contig_start = std::stoull(fields[7]);
        uint64_t contig_end = std::stoull(fields[8]);
        uint32_t alignment_score = fields.size() > 11 ? std::stoul(fields[11]) : 0;

        // Check if contig already exists, if not add it
        uint32_t contig_index;
        if (contig_indices.find(contig_name) == contig_indices.end()) {
            store.add_contig(Contig(contig_name, contig_length));
            contig_index = store.get_contigs().size() - 1;
            contig_indices[contig_name] = contig_index;
        } else {
            contig_index = contig_indices[contig_name];
        }

        // Check if read already exists, if not add it
        uint32_t read_index;
        if (read_indices.find(read_id) == read_indices.end()) {
            store.add_read(Read(read_id, read_length));
            read_index = store.get_reads().size() - 1;
            read_indices[read_id] = read_index;
        } else {
            read_index = read_indices[read_id];
        }

        // Create alignment
        Alignment alignment(read_id, contig_name,
                          read_length, contig_length,
                          read_start, read_end,
                          contig_start, contig_end,
                          read_end - read_start, alignment_score,
                          is_reverse);

        // Look for cg:Z and cs:Z tags
        std::string cigar_string;
        std::string cs_string;
        
        for (size_t i = 12; i < fields.size(); ++i) {
            if (fields[i].substr(0, 5) == "cg:Z:") {
                cigar_string = fields[i].substr(5);
            } else if (fields[i].substr(0, 5) == "cs:Z:") {
                cs_string = fields[i].substr(5);
            }
        }

        // Parse mutations from cs string
        if (!cs_string.empty()) {
            parse_cs_string(cs_string, alignment);
        }

        // Compare CIGAR and cs strings if both are available
        if (!cigar_string.empty() && !cs_string.empty()) {
            compare_cigar_and_cs(cigar_string, cs_string, read_id);
        }

        // Add alignment to store
        store.add_alignment(alignment);
        
        // Update read's alignment indices
        Read& read = store.get_reads()[read_index];
        read.alignment_indices.push_back(store.get_alignments().size() - 1);
    }

    file.close();
    return true;
}

void PafReader::split_line(const std::string& line, char delimiter, std::vector<std::string>& fields) 
{
    // Empty implementation for now
}

void PafReader::parse_cs_string(const std::string& cs_string, Alignment& alignment) 
{
    // Empty implementation for now
}

void PafReader::compare_cigar_and_cs(const std::string& cigar, const std::string& cs, const std::string& read_id) 
{
    // Empty implementation for now
} 