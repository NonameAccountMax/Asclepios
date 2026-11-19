// Minimal binary archive for restart files (same binary, same machine: raw bytes of POD data).
#pragma once
#include <cstdint>
#include <cstdio>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>

namespace abq {

class BinOut {
 public:
  explicit BinOut(const std::string& path) : f_(std::fopen(path.c_str(), "wb")) {
    if (!f_) throw std::runtime_error("cannot write " + path);
  }
  ~BinOut() { if (f_) std::fclose(f_); }
  template <class T> void pod(const T& x) {
    static_assert(std::is_trivially_copyable<T>::value, "POD only");
    std::fwrite(&x, sizeof(T), 1, f_);
  }
  template <class T> void vec(const std::vector<T>& v) {
    static_assert(std::is_trivially_copyable<T>::value, "POD only");
    uint64_t n = v.size();
    pod(n);
    if (n) std::fwrite(v.data(), sizeof(T), n, f_);
  }
  void tag(const char* s) { uint64_t n = std::char_traits<char>::length(s); pod(n); std::fwrite(s, 1, n, f_); }

 private:
  std::FILE* f_;
};

class BinIn {
 public:
  explicit BinIn(const std::string& path) : f_(std::fopen(path.c_str(), "rb")) {
    if (!f_) throw std::runtime_error("cannot read " + path);
  }
  ~BinIn() { if (f_) std::fclose(f_); }
  template <class T> void pod(T& x) {
    if (std::fread(&x, sizeof(T), 1, f_) != 1) throw std::runtime_error("restart file truncated");
  }
  template <class T> void vec(std::vector<T>& v) {
    uint64_t n;
    pod(n);
    v.resize(n);
    if (n && std::fread(v.data(), sizeof(T), n, f_) != n) throw std::runtime_error("restart file truncated");
  }
  std::string read_tag() {
    uint64_t n;
    pod(n);
    std::string got(n, ' ');
    if (n && std::fread(&got[0], 1, n, f_) != n) throw std::runtime_error("restart file truncated");
    return got;
  }
  // checks a section tag written by BinOut::tag
  void tag(const char* s) {
    uint64_t n;
    pod(n);
    std::string got(n, ' ');
    if (n && std::fread(&got[0], 1, n, f_) != n) throw std::runtime_error("restart file truncated");
    if (got != s) throw std::runtime_error("restart file: expected section '" + std::string(s) + "', found '" + got + "'");
  }

 private:
  std::FILE* f_;
};

}  // namespace abq
