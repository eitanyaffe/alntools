#include <assert.h>
#include <fstream>
#include <iostream>
#include <map>
#include <math.h>
#include <set>
#include <sstream>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <vector>
#ifdef _OPENMP
#include <omp.h>
#endif

int construct_main(const char* name, int argc, char** argv);
int info_main(const char* name, int argc, char** argv);
int extract_main(const char* name, int argc, char** argv);
int verify_main(const char* name, int argc, char** argv);
int query_main(const char* name, int argc, char** argv);
int coverage_main(const char* name, int argc, char** argv);
int breaks_main(const char* name, int argc, char** argv);
int rearrange_main(const char* name, int argc, char** argv);
int homologs_main(const char* name, int argc, char** argv);
int segments_main(const char* name, int argc, char** argv);
int cov_matrix_main(const char* name, int argc, char** argv);
int csegment_coverage_main(const char* name, int argc, char** argv);
int seg_matrix_main(const char* name, int argc, char** argv);
int inspect_read_main(const char* name, int argc, char** argv);
int get_read_ids_main(const char* name, int argc, char** argv);
int cov_intervals_main(const char* name, int argc, char** argv);
int get_local_deletions_main(const char* name, int argc, char** argv);

using namespace std;

void usage(const char* name)
{
  fprintf(stderr, "alntools: \n");
  fprintf(stderr, "usage: %s <command> [options]\n", name);
  fprintf(stderr, "commands:\n");
  fprintf(stderr, "  construct: Construct ALN file from PAF file\n");
  fprintf(stderr, "  info: Show basic info and stats for ALN file\n");
  fprintf(stderr, "  extract: Save ALN file to tab-delimited tables\n");
  fprintf(stderr, "  verify: verify ALN file using reads and contigs\n");
  fprintf(stderr, "  query: query ALN file\n");
  fprintf(stderr, "  coverage: generate read and contig alignment coverage statistics\n");
  fprintf(stderr, "  breaks: find positions with excessive read start/end events\n");
  fprintf(stderr, "  rearrange: detect genome rearrangements from read alignments\n");
  fprintf(stderr, "  segments: detect assembly breakpoints and generate segments\n");
  fprintf(stderr, "  cov_matrix: generate coverage matrix for segments across multiple libraries\n");
  fprintf(stderr, "  csegment_coverage: generate unique read coverage for csegments across multiple libraries\n");
  fprintf(stderr, "  seg_matrix: generate segment adjacency and reach matrices from read alignments\n");
  fprintf(stderr, "  homologs: find homologous regions using kmer-based search\n");
  fprintf(stderr, "  inspect_read: print detailed information about a specific read\n");
  fprintf(stderr, "  get_read_ids: assign each read to the bin with the longest segment intersection\n");
  fprintf(stderr, "  cov_intervals: compute intervalcoverage (fraction covered) per segment interval\n");
  fprintf(stderr, "  get_local_deletions: find reads with a deletion matching a query interval\n");
  fprintf(stderr, "\n");
#ifdef _OPENMP
  fprintf(stderr, "thread support: enabled (OpenMP available, max threads: %d)\n", omp_get_max_threads());
#else
  fprintf(stderr, "thread support: disabled (compiled without OpenMP)\n");
#endif
}

int main(int argc, char** argv)
{
  if (argc == 1) {
    usage(argv[0]);
    exit(1);
  }
  string command(argv[1]);
  string name = string(argv[0]) + " " + command;

  int rc = 0;
  if (command == "construct") {
    rc = construct_main(name.c_str(), argc - 1, argv + 1);
  } else if (command == "info") {
    rc = info_main(name.c_str(), argc - 1, argv + 1);
  } else if (command == "extract") {
    rc = extract_main(name.c_str(), argc - 1, argv + 1);
  } else if (command == "verify") {
    rc = verify_main(name.c_str(), argc - 1, argv + 1);
  } else if (command == "query") {
    rc = query_main(name.c_str(), argc - 1, argv + 1);
  } else if (command == "coverage") {
    rc = coverage_main(name.c_str(), argc - 1, argv + 1);
  } else if (command == "breaks") {
    rc = breaks_main(name.c_str(), argc - 1, argv + 1);
  } else if (command == "rearrange") {
    rc = rearrange_main(name.c_str(), argc - 1, argv + 1);
  } else if (command == "segments") {
    rc = segments_main(name.c_str(), argc - 1, argv + 1);
  } else if (command == "cov_matrix") {
    rc = cov_matrix_main(name.c_str(), argc - 1, argv + 1);
  } else if (command == "csegment_coverage") {
    rc = csegment_coverage_main(name.c_str(), argc - 1, argv + 1);
  } else if (command == "seg_matrix") {
    rc = seg_matrix_main(name.c_str(), argc - 1, argv + 1);
  } else if (command == "homologs") {
    rc = homologs_main(name.c_str(), argc - 1, argv + 1);
  } else if (command == "inspect_read") {
    rc = inspect_read_main(name.c_str(), argc - 1, argv + 1);
  } else if (command == "get_read_ids") {
    rc = get_read_ids_main(name.c_str(), argc - 1, argv + 1);
  } else if (command == "cov_intervals") {
    rc = cov_intervals_main(name.c_str(), argc - 1, argv + 1);
  } else if (command == "get_local_deletions") {
    rc = get_local_deletions_main(name.c_str(), argc - 1, argv + 1);
  } else {
    printf("unknown command: %s\n", command.c_str());
    usage(argv[0]);
    exit(1);
  }

  return rc;
}
