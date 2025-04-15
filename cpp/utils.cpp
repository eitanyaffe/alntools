#include "utils.h"

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
