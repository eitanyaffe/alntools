#include "paf_reader.h"
#include "utils.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <unordered_map>
#include <regex>
#include <stdexcept>
#include <string_view>

using namespace std;

void PafReader::load_reads_contigs(const string& ifn_reads,
				   const string& ifn_contigs)
{
  unordered_set<string> empty_set;
  
  read_fasta(ifn_contigs, empty_set, m_contigs);
  read_fastq(ifn_reads, empty_set, m_reads);

  cout << "Loaded " << m_contigs.size() << " contigs and " << m_reads.size() << " reads" << endl;
}

bool PafReader::verify_alignment(Alignment& alignment,
				 const string& read_id,
				 const string& contig_id,
				 bool quit_on_error)
{
   cout << "==================\n"
         << "Read: " << read_id << " [" << alignment.read_start << "," << alignment.read_end << "] "
         << "  Contig: " << contig_id << " [" << alignment.contig_start << "," << alignment.contig_end << "] "
         << "  Is reverse: " << (alignment.is_reverse ? "yes" : "no") << "\n";

    // Count mutations by type
    auto count_mutations_by_type = [](const vector<Mutation>& mutations) {
        size_t subs = 0, ins = 0, dels = 0;
        for (const auto& mut : mutations) {
            switch (mut.type) {
                case MutationType::SUBSTITUTION: subs++; break;
                case MutationType::INSERTION: ins++; break;
                case MutationType::DELETION: dels++; break;
            }
        }
        return make_tuple(subs, ins, dels);
    };

    auto [num_subs, num_ins, num_dels] = count_mutations_by_type(alignment.mutations);
    cout << "Mutations - Substitutions: " << num_subs 
         << ", Insertions: " << num_ins
         << ", Deletions: " << num_dels << "\n";

    massert(m_contigs.find(contig_id) != m_contigs.end(),
            "Error: Contig '%s' not found in FASTA file.", contig_id.c_str());
    massert(m_reads.find(read_id) != m_reads.end(),
            "Error: Read '%s' not found in FASTQ file.", read_id.c_str());

    cout << "mutating contig with " << alignment.mutations.size() << " mutations" << endl;

    massert(m_contigs.find(contig_id) != m_contigs.end(), "contig %s not found in map", contig_id.c_str());
    massert(m_reads.find(read_id) != m_reads.end(), "read %s not found in map", read_id.c_str());

    string contig_fragment = m_contigs[contig_id].substr(alignment.contig_start,
							 alignment.contig_end - alignment.contig_start);

    string mutated_contig = apply_mutations(contig_fragment, alignment.mutations, alignment.contig_start, quit_on_error);
    string read_segment = m_reads[read_id].substr(alignment.read_start,
						  alignment.read_end - alignment.read_start);

    if (alignment.is_reverse)
      mutated_contig = reverse_complement(mutated_contig);

    if(read_segment.size() != mutated_contig.size()) {
      printf("read segment length (%zu) does not match mutated contig length (%zu)", 
	     read_segment.size(), mutated_contig.size());
      cout << "read_segment='" << read_segment << "'" << endl;
      cout << "mutated_contig='" << mutated_contig << "'" << endl;
      return false;
    }

    bool mismatch_found = false;
    for (size_t i = 0; i < read_segment.size(); ++i)
    {
      if (mutated_contig[i] != read_segment[i])
      {
 
        // Calculate start and end indices for the segment
        size_t start = (i >= 5) ? i - 8 : 0;
        size_t end = (i + 5 < read_segment.size()) ? i + 8 : read_segment.size() - 1;

        cout << "Mismatch found, fragment coordinate=" << i << endl; 
        cout << "read        : " << read_segment.substr(start, end - start + 1) << "\n";
        cout << "contig_mut  : " << mutated_contig.substr(start, end - start + 1) << "\n";
        cout << "contig_orig : " << contig_fragment.substr(start, end - start + 1) << "\n";
 
        mismatch_found = true;
	break;
      }
    }

    return !mismatch_found;
}

void PafReader::read_paf(const string& filename, AlignmentStore& store, int max_reads, bool should_verify, bool quit_on_error)
{
    ifstream file(filename);
    massert(file.is_open(), "Failed to open file: %s", filename.c_str());

    string line;
    size_t line_number = 0;
    size_t mutation_count = 0;
    size_t bad_alignment_count = 0;

    while (std::getline(file, line)) {
        line_number++;
        if (line_number % 10000 == 0)
            std::cout << "Processed " << line_number << " alignments..." << std::endl;
	if (line_number > (size_t)max_reads && max_reads != 0)
	  break;
	
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

	if (should_verify) {
	  bool is_valid = verify_alignment(alignment, read_id, contig_id, quit_on_error);
          if (!is_valid) {
            bad_alignment_count++;
            if (bad_alignment_count >= 10) {
              cout << "reached maximum number of bad alignments (10), stopping" << endl;
              break;
            }
          }
	}
        // Add alignment to store
        store.add_alignment(alignment);
    }
    std::cout << "Total mutations found: " << mutation_count << "\n";
    if (should_verify) {
      std::cout << "found " << bad_alignment_count << " bad alignments" << std::endl;
    }

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
                char ref_base = toupper(segment[0]);
                char read_base = toupper(segment[1]);
                alignment.add_mutation(Mutation(MutationType::SUBSTITUTION, ref_pos,
						string(1, read_base), string(1, ref_base)));
                ref_pos++;
                mutation_count++;
                break;
            }
            case '+': // Insertion
            {
	            string insertion_bases = to_upper(segment);
                alignment.add_mutation(Mutation(MutationType::INSERTION, ref_pos, insertion_bases, ""));
                mutation_count++;
                break;
            }
            case '-': // Deletion
            {
                string deleted_bases = to_upper(segment);
                alignment.add_mutation(Mutation(MutationType::DELETION, ref_pos, "", deleted_bases));
                ref_pos += deleted_bases.length();
                mutation_count++;
                break;
            }
            case ':': // Identity
            {
                char* end;
                long val = strtol(segment.c_str(), &end, 10);
                if (end == segment.c_str() || *end != '\0' || val < 0) {
                    mexit("Failed to convert segment '%s' to valid positive integer", segment.c_str());
                }
                ref_pos += val;
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
