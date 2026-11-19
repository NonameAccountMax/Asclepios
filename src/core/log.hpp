#pragma once
#include <cstdarg>
#include <cstdio>
#include <string>

namespace abq {

struct Log {
  static FILE*& dat() { static FILE* f = nullptr; return f; }
  static int& warnings() { static int n = 0; return n; }
  static bool& quiet() { static bool q = false; return q; }
};

inline void vlogf(const char* prefix, const char* fmt, va_list ap) {
  char buf[4096];
  vsnprintf(buf, sizeof(buf), fmt, ap);
  if (!Log::quiet()) std::fprintf(stdout, "%s%s\n", prefix, buf);
  if (Log::dat()) { std::fprintf(Log::dat(), "%s%s\n", prefix, buf); std::fflush(Log::dat()); }
}
inline void info(const char* fmt, ...) { va_list ap; va_start(ap, fmt); vlogf("", fmt, ap); va_end(ap); }
inline void warn(const char* fmt, ...) {
  ++Log::warnings();
  va_list ap; va_start(ap, fmt); vlogf("***WARNING: ", fmt, ap); va_end(ap);
}

}  // namespace abq
