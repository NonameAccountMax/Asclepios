// Results: VTK XML unstructured grids (one .vtu per frame, collected in a .pvd), energy history
// (CSV) and an Abaqus-like status file (.sta).
#pragma once
#include <cstdio>
#include <string>
#include <vector>

namespace abq {

struct Analysis;

class OutputWriter {
 public:
  OutputWriter(Analysis& A, const std::string& job);
  ~OutputWriter();
  void frame(int step, const std::string& step_name, double step_time, double total_time);
  void history(int step, int inc, double step_time, double total_time, double dt, double ke);
  void sta_header(int step, const std::string& name, bool explicit_step);
  void rf_history(double total_time, const double r[3]);
  void mon_history(double total_time, int label, const double* u6);
  void sta_line(const char* fmt, ...);
  void close();
  int frames() const { return nframes_; }

 private:
  Analysis& A_;
  std::string job_, dir_;
  FILE* pvd_ = nullptr;
  FILE* hist_ = nullptr;
  FILE* sta_ = nullptr;
  FILE* rfh_ = nullptr;
  FILE* monh_ = nullptr;
  int nframes_ = 0;
  std::vector<std::string> pvd_entries_;
  void write_pvd();
};

}  // namespace abq
