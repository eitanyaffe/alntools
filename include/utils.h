#pragma once

#include <cstdarg>
#include <cstdio>
#include <cstdlib>

void massert(bool cond, const char* fmt, ...);
void mexit(const char *fmt, ...);
