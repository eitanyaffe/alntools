#include "paf_reader.h"
#include "utils.h"
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

void PafReader::read_paf(const string& filename, AlignmentStore& store) 
{
    ifstream file(filename);
    massert(file.is_open(), "Failed to open file: %s", filename.c_str());

    string line;
    size_t line_number = 0;
    size_t mutation_count = 0;

    while (std::getline(file, line)) {
        line_number++;
        if (line_number % 10000 == 0) {
            std::cout << "Processed " << line_number << " alignments..." << std::endl;
        }
        vector<string> fields;
        split_line(line, '\t', fields);

        massert(fields.size() >= 12, "Malformed line %zu with fewer than 12 fields: %s", 
                line_number, line.c_str());

        // Get read information
        string read_id = fields[0];
        uint64_t read_length = stoull(fields[1]);
        bool is_reverse = fields[4] == "-"; // Assuming the reverse flag is in the 5th field

        // Use AlignmentStore to get or add read index
        size_t read_index = store.add_or_get_read_index(read_id, read_length);

        // Get contig information
        string contig_id = fields[5]; // Assuming the contig ID is in the 6th field
        uint64_t contig_length = stoull(fields[6]); // Assuming the contig length is in the 7th field

        // Use AlignmentStore to get or add contig index
        size_t contig_index = store.add_or_get_contig_index(contig_id, contig_length);

        // Get alignment information
        uint64_t read_start = stoull(fields[2]);
        uint64_t read_end = stoull(fields[3]);
        uint64_t contig_start = stoull(fields[7]);
        uint64_t contig_end = stoull(fields[8]);

        // Validate numeric values
        massert(read_end > read_start, "Invalid read coordinates on line %zu: end (%zu) <= start (%zu)", 
                line_number, read_end, read_start);
        
        massert(contig_end > contig_start, "Invalid contig coordinates on line %zu: end (%zu) <= start (%zu)", 
                line_number, contig_end, contig_start);

        // Create alignment
        Alignment alignment(read_index, contig_index,
                          contig_start, contig_end,
                          read_start, read_end,
                          is_reverse);

        // Look for cs:Z tags      
        if (fields.size() > 12) {
            for (size_t i = 12; i < fields.size(); ++i) {
                if (fields[i].substr(0, 5) == "cs:Z:") {
                    string cs_string = fields[i].substr(5);
                    mutation_count += add_mutations(cs_string, alignment);
                }
            }
        }      

        // Add alignment to store
        store.add_alignment(alignment);
    }
    std::cout << "Total mutations found: " << mutation_count << "\n";

    file.close();
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

void PafReader::parse_cs_string(const std::string& cs_string,
                                std::vector<char>& actions,
                                std::vector<std::string>& values) 
{
    std::string current;
    char action = '\0';

    for (size_t i = 0; i < cs_string.size(); ++i) {
        char c = cs_string[i];
        if (c == ':' || c == '=' || c == '*' || c == '+' || c == '-' || c == '~') {
            if (!current.empty()) {
                massert(action != '\0', "Invalid action: %c", action);
                actions.push_back(action);
                values.push_back(current);
            }
            action = c;
            current.clear();
        } else {
            current += c;
        }
    }

    // Push the last segment
    if (!current.empty()) {
        massert(action != '\0', "Invalid action: %c", action);
        actions.push_back(action);
        values.push_back(current);
    }
}

size_t PafReader::add_mutations(const string& cs_string, Alignment& alignment) 
{
    // Clear any existing mutations
    alignment.clear_mutations();
    
    // Current position in the reference sequence
    uint32_t ref_pos = alignment.contig_start;
    size_t mutation_count = 0;

    std::vector<char> actions;
    std::vector<std::string> values;
    parse_cs_string(cs_string, actions, values);

    for (size_t i = 0; i < actions.size(); ++i) {
        char action = actions[i];
        string segment = values[i];
         switch (action) {
            case '*': // Substitution
            {
                massert(segment.length() == 2, "Invalid substitution segment length: %zu", segment.length());
                char ref_base = segment[0];
                char query_base = segment[1];
                alignment.add_mutation(Mutation(MutationType::SUBSTITUTION, ref_pos,
						string(1, query_base), string(1, ref_base)));
                ref_pos++;
                mutation_count++;
                break;
            }
            case '+': // Insertion
            {
                alignment.add_mutation(Mutation(MutationType::INSERTION, ref_pos, segment));
                ref_pos += segment.length();
                mutation_count++;
                break;
            }
            case '-': // Deletion
            {
                alignment.add_mutation(Mutation(MutationType::DELETION, ref_pos, "", segment));
                ref_pos += segment.length();
                mutation_count++;
                break;
            }
            default:
                // Skip other actions
                break;
        }
    }
    
    return mutation_count;
}
