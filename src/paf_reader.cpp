#include "../include/paf_reader.h"
#include "../include/utils.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <unordered_map>
#include <regex>
#include <stdexcept>
#include <string_view>

using std::ifstream;
using std::ofstream;
using std::cerr;
using std::endl;
using std::unordered_map;
using std::stoull;
using std::stoul;
using std::regex;
using std::smatch;
using std::regex_match;
using std::stringstream;
using std::invalid_argument;
using std::string_view;

bool PafReader::read_paf(const string& filename, AlignmentStore& store) 
{
    ifstream file(filename);
    massert(file.is_open(), "Failed to open file: %s", filename.c_str());

    // Maps for tracking existing contigs and reads
    unordered_map<string, uint32_t> contig_indices;
    unordered_map<string, uint32_t> read_indices;

    string line;
    size_t line_number = 0;
    
    while (std::getline(file, line)) {
        line_number++;
        if (line_number % 10 == 0) {
            std::cout << "Processed " << line_number << " lines..." << std::endl;
        }
        vector<string> fields;
        split_line(line, '\t', fields);

        massert(fields.size() >= 12, "Malformed line %zu with fewer than 12 fields: %s", 
                line_number, line.c_str());

        // Get read information
        string read_id = fields[0];
        uint64_t read_length = stoull(fields[1]);
        
        // Get contig information
        string contig_name = fields[5];
        uint64_t contig_length = stoull(fields[6]);
        
        // Get alignment information
        uint64_t read_start = stoull(fields[2]);
        uint64_t read_end = stoull(fields[3]);
        bool is_reverse = (fields[4] == "-");
        uint64_t contig_start = stoull(fields[7]);
        uint64_t contig_end = stoull(fields[8]);
        uint32_t alignment_score = fields.size() > 11 ? stoul(fields[11]) : 0;

        // Validate numeric values
        massert(read_end > read_start, "Invalid read coordinates on line %zu: end (%zu) <= start (%zu)", 
                line_number, read_end, read_start);
        
        massert(contig_end > contig_start, "Invalid contig coordinates on line %zu: end (%zu) <= start (%zu)", 
                line_number, contig_end, contig_start);

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

        // Look for cs:Z tags
        string cs_string;
        
        if (fields.size() > 12) {
            for (size_t i = 12; i < fields.size(); ++i) {
                if (fields[i].substr(0, 5) == "cs:Z:") {
                    cs_string = fields[i].substr(5);
                }
            }
        }

        // Parse mutations from cs string
        if (!cs_string.empty()) {
            parse_cs_string(cs_string, alignment);
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

void PafReader::split_line(const string& line, char delimiter, vector<string>& fields) 
{
    fields.clear();
    
    // Reserve space for fields to avoid reallocations
    // Most PAF lines have around 12-15 fields
    fields.reserve(15);
    
    // Use string_view for efficient substring operations
    string_view line_view(line);
    size_t start = 0;
    size_t end = 0;
    
    while ((end = line_view.find(delimiter, start)) != string_view::npos) {
        fields.emplace_back(line_view.substr(start, end - start));
        start = end + 1;
    }
    
    // Add the last field
    if (start < line_view.size()) {
        fields.emplace_back(line_view.substr(start));
    }
}

void PafReader::parse_cs_string(const string& cs_string, Alignment& alignment) 
{
    // Clear any existing mutations
    alignment.clear_mutations();
    
    // Current position in the reference sequence
    uint32_t ref_pos = alignment.target_start;
    
    // Split the cs string by ':'
    vector<string> segments;
    split_line(cs_string, ':', segments);
    
    for (const auto& segment : segments) 
    {
        if (segment.empty()) 
        {
            continue; // Skip empty segments
        }
        
        char action = segment[0];
        
        switch (action) 
        {
            case '*': // Substitution
            {
                massert(segment.length() == 3, "Invalid substitution segment length: %zu", segment.length());
                char ref_base = segment[1];
                char query_base = segment[2];
                alignment.add_mutation(Mutation(MutationType::SUBSTITUTION, ref_pos, string(1, query_base), string(1, ref_base)));
                ref_pos++;
                break;
            }
            case '+': // Insertion
            {
                string inserted_bases = segment.substr(1);
                alignment.add_mutation(Mutation(MutationType::INSERTION, ref_pos, inserted_bases));
                // ref_pos doesn't change for insertions
                break;
            }
            case '-': // Deletion
            {
                string deleted_bases = segment.substr(1);
                alignment.add_mutation(Mutation(MutationType::DELETION, ref_pos, "", deleted_bases));
                ref_pos += deleted_bases.length();
                break;
            }
            default:
                // Skip other actions
                break;
        }
    }
}