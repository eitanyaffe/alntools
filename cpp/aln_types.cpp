#include "aln_types.h"
#include "utils.h"

#include <sstream>
#include <string>

// Implementation of any non-inline functions from aln_types.h would go here
// Currently all functions are inline, so this file is empty but kept for future use

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