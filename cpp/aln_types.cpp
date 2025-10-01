#include "aln_types.h"
#include "alignment_store.h"
#include "utils.h"

#include <sstream>
#include <string>

// Function to print Alignment with read and contig names
void print_alignment(const Alignment& align, const AlignmentStore& store, const std::string& label)
{
  const std::string& read_id = store.get_read_id(align.read_index);
  const std::string& contig_id = store.get_contig_id(align.contig_index);
  
  if (!label.empty()) {
    std::cerr << label << ": ";
  }
  
  std::cerr << "Alignment(read=\"" << read_id << "\" contig=\"" << contig_id << "\" "
            << "contig_start=" << align.contig_start << " contig_end=" << align.contig_end
            << " read_start=" << align.read_start << " read_end=" << align.read_end 
            << " is_reverse=" << align.is_reverse << ")" << std::endl;
}

// Create a unique string key for this mutation on a given contig
std::string Mutation::create_key(uint32_t contig_index) const
{
  std::ostringstream oss;
  oss << contig_index << "_" << position << "_";

  switch (type) {
  case MutationType::SUBSTITUTION:
    massert(nts.length() == 2, "SUB mutation nts length is not 2 for cs tag generation");
    oss << "SUB_" << nts;
    break;
  case MutationType::INSERTION:
    oss << "INS_" << nts;
    break;
  case MutationType::DELETION:
    oss << "DEL_" << nts;
    break;
  }
  return oss.str();
}