#include "utils.h"
#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>

#include "aln_types.h"

using namespace std;

void massert(bool cond, const char* fmt, ...)
{
  if (cond) {
    return;
  }

  fprintf(stderr, "Error: ");

  va_list argp;
  va_start(argp, fmt);
  vfprintf(stderr, fmt, argp);
  va_end(argp);

  fprintf(stderr, "\n");
  exit(-1);
}

void mexit(const char *fmt, ...)
{
  fprintf(stderr, "Error: ");

  va_list argp;
  va_start(argp, fmt);
  vfprintf(stderr, fmt, argp);
  va_end(argp);

  fprintf(stderr, "\n");
  exit(-1);
}

// Function to convert a string to uppercase
string to_upper(string str)
{
  std::transform(str.begin(), str.end(), str.begin(), ::toupper);
  return str;
}

std::string reverse_complement(std::string seq) {
  std::string result = seq;
  int N = seq.length();
  for (int i=0; i<N; i++) {
    char c = seq[N-i-1], r;
    switch( c )
      {
      case 'A': r = 'T'; break;
      case 'G': r = 'C'; break;
      case 'C': r = 'G'; break;
      case 'T': r = 'A'; break;
      default: r = c;
      }
    result[i] = r;
  }
  return(result);
}

// Function to read FASTA file and create a map of contig IDs to sequences
void read_fasta(const string &filename,
		const unordered_set<string> &contig_ids,
		unordered_map<string, string>& contigs)
{
  cout << "Reading FASTA file: " << filename << endl;
  ifstream file(filename);
  string line, id, sequence;
  while (getline(file, line))
    {
      if (line[0] == '>')
	{
	  if (!id.empty() && (contig_ids.empty() || contig_ids.find(id) != contig_ids.end()))
	    {
	      contigs[id] = sequence;
	    }
	  id = line.substr(1, line.find_first_of(" ") - 1); // Extract ID before first space/tab
	  sequence.clear();
	}
      else
	{
	  sequence += line;
	}
    }
  if (!id.empty() && (contig_ids.empty() || contig_ids.find(id) != contig_ids.end()))
    {
      contigs[id] = sequence;
    }
}

// Function to read FASTQ file and create a map of read IDs to sequences
void read_fastq(const string &filename,
		const unordered_set<string> &read_ids,
		unordered_map<string, string>& reads)
{
  cout << "Reading FASTQ file: " << filename << endl;
  ifstream file(filename);
  string line, id, sequence;
  while (getline(file, line))
    {
      if (line[0] == '@')
	{
	  id = line.substr(1);
	  getline(file, sequence); // Read the sequence line
	  if (read_ids.empty() || read_ids.find(id) != read_ids.end())
	    reads[id] = sequence;
	  getline(file, line); // Skip the '+' line
	  getline(file, line); // Skip the quality line
	}
    }
}

// Function to apply mutations to a contig fragment
string apply_mutations(const string &contig_fragment,
                       const vector<Mutation> &mutations,
                       uint32_t contig_start, bool quit_on_error)
{
  string result;
  size_t current_pos = 0;
  bool error_found = false;
  
  // Create a copy of mutations and sort by position
  vector<Mutation> sorted_mutations = mutations;
  std::sort(sorted_mutations.begin(), sorted_mutations.end(), 
            [](const Mutation& a, const Mutation& b) {
              return a.position < b.position;
            });
    
  // Process each mutation in order
  for (const auto &mutation : sorted_mutations)
    {
      // Convert absolute contig position to position within fragment
      size_t fragment_pos = mutation.position - contig_start;
      
      // Verify position is within bounds
      massert(mutation.position >= contig_start && fragment_pos < contig_fragment.size(),
	      "mutation position %u is outside fragment bounds", mutation.position);
        
      // Copy unchanged sequence up to this mutation
      if (fragment_pos > current_pos) {
	result.append(contig_fragment.substr(current_pos, fragment_pos - current_pos));
      }
        
      string ref_nts = to_upper(mutation.ref_nts);
      string read_nts = to_upper(mutation.read_nts);

      // Apply mutation based on type
      switch (mutation.type)
        {
	case MutationType::SUBSTITUTION:
   
	  // Verify reference bases match expected
	  if(contig_fragment.substr(fragment_pos, ref_nts.size()) != ref_nts) {
	    printf("reference bases at position %zu do not match expected. expected: %s, found: %s\n",
		   fragment_pos, ref_nts.c_str(), 
		   contig_fragment.substr(fragment_pos, ref_nts.size()).c_str());
	    error_found = true;
	  }
	  // Verify substitution lengths match
	  if (read_nts.size() != ref_nts.size()) {
	    printf("substitution lengths do not match at position %zu. read: %zu, ref: %zu\n",
		   fragment_pos, read_nts.size(), ref_nts.size());
	    error_found = true;
	  }
	  // Add the substituted bases
	  result.append(read_nts);
                
	  // Move position past the substituted bases
	  current_pos = fragment_pos + ref_nts.size();
	  break;
                
	case MutationType::INSERTION:
	  // Add inserted bases
	  result.append(read_nts);
                
	  // Position doesn't change since nothing is consumed from reference
	  current_pos = fragment_pos;
	  break;
                
	case MutationType::DELETION:
	  // Verify reference bases match expected
	  if(contig_fragment.substr(fragment_pos, ref_nts.size()) != ref_nts) {
	    printf("reference bases at position %zu do not match expected for deletion. expected: %s, found: %s\n",
		   fragment_pos, ref_nts.c_str(), 
		   contig_fragment.substr(fragment_pos, ref_nts.size()).c_str());
	    error_found = true;
	  }
                
	  // Skip over the deleted bases (don't add them to result)
	  current_pos = fragment_pos + ref_nts.size();
	  break;
        }
    }
    
  // Append any remaining reference sequence
  if (current_pos < contig_fragment.size()) {
    result.append(contig_fragment.substr(current_pos));
  }

  if (quit_on_error && error_found) {
    cout << "error found, quiting" << endl;
  }
    
  return result;
}
