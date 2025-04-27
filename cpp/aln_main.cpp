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

int construct_main(const char* name, int argc, char** argv);
int info_main(const char* name, int argc, char** argv);
int extract_main(const char* name, int argc, char** argv);
int verify_main(const char* name, int argc, char** argv);
int query_main(const char* name, int argc, char** argv);

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
  } else {
    printf("unknown command: %s\n", command.c_str());
    usage(argv[0]);
    exit(1);
  }

  return rc;
}
