// asclepios: Abaqus-like lite finite element solver (explicit dynamics + static implicit).
#include <omp.h>

#include <chrono>
#include <cstdio>
#include <cstring>
#include <stdexcept>
#include <string>

#include "core/binio.hpp"
#include "core/log.hpp"
#include "model/model.hpp"
#include "output/output.hpp"
#include "solver/solvers.hpp"

using namespace abq;

static void usage() {
  std::fprintf(stderr,
               "usage: asclepios job=<name> [input=<file.inp>] [cpus=<n>] [save=<dt>] [restart=<file.rst>]\n"
               "       asclepios <file.inp> [-j <cpus>] [-o <job>]\n"
               "  save=<dt>        write <job>.rst every dt of total time (explicit steps)\n"
               "  restart=<file>   resume from a restart file written by the same model and binary\n");
}

int main(int argc, char** argv) {
  std::string input, job;
  int cpus = 0;
  double save_dt = 0;
  std::string restart;
  for (int i = 1; i < argc; ++i) {
    std::string a = argv[i];
    if (a.rfind("job=", 0) == 0) job = a.substr(4);
    else if (a.rfind("save=", 0) == 0) save_dt = std::stod(a.substr(5));
    else if (a.rfind("restart=", 0) == 0) restart = a.substr(8);
    else if (a.rfind("input=", 0) == 0) input = a.substr(6);
    else if (a.rfind("cpus=", 0) == 0) cpus = std::stoi(a.substr(5));
    else if (a == "-j" && i + 1 < argc) cpus = std::stoi(argv[++i]);
    else if (a == "-o" && i + 1 < argc) job = argv[++i];
    else if (a == "-q") Log::quiet() = true;
    else if (a == "-h" || a == "--help") { usage(); return 0; }
    else if (input.empty()) input = a;
    else { usage(); return 1; }
  }
  if (input.empty() && !job.empty()) input = job + ".inp";
  if (input.empty()) { usage(); return 1; }
  if (job.empty()) {
    job = input;
    size_t sl = job.find_last_of('/');
    if (sl != std::string::npos) job = job.substr(sl + 1);
    if (job.size() > 4 && job.substr(job.size() - 4) == ".inp") job = job.substr(0, job.size() - 4);
  }
  if (cpus <= 0) cpus = std::max(1, omp_get_num_procs() / 2);
  Log::dat() = std::fopen((job + ".dat").c_str(), "w");
  auto t0 = std::chrono::steady_clock::now();
  try {
    info("asclepios — job %s, input %s, %d threads", job.c_str(), input.c_str(), cpus);
    Model m = build_model(input);
    info("model: %d nodes, %d elements, %zu steps", m.nn(), m.ne(), m.steps.size());
    Analysis A(m);
    A.job = job;
    A.nthreads = cpus;
    A.setup();
    A.writer = std::make_unique<OutputWriter>(A, job);
    A.save_interval = save_dt;
    A.next_save = save_dt;
    size_t s0 = 0;
    if (!restart.empty()) {
      BinIn in(restart);
      in.tag("asclepios-restart");
      int si;
      in.pod(si);
      if (si < 0 || si >= (int)m.steps.size() || m.steps[si].proc != Step::EXPLICIT)
        throw std::runtime_error("restart file: invalid step");
      A.restart_in = restart;
      A.resume_step = si;
      s0 = si;
    }
    for (size_t s = s0; s < m.steps.size(); ++s) {
      const Step& st = m.steps[s];
      if (st.proc == Step::EXPLICIT) run_explicit_step(A, st, (int)s);
      else if (st.riks) run_riks_step(A, st, (int)s);
      else if (st.proc == Step::HEAT) run_heat_step(A, st, (int)s);
      else if (st.proc != Step::STATIC || st.perturbation) run_linear_step(A, st, (int)s);
      else run_static_step(A, st, (int)s);
    }
    double el = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    info("\nTHE ANALYSIS HAS COMPLETED SUCCESSFULLY  (%d increments, %.1f s, %d warnings, %d frames in %s.pvd)",
         A.total_increments, el, Log::warnings(), A.writer->frames(), job.c_str());
    A.writer->close();
  } catch (std::exception& e) {
    std::fprintf(stderr, "***ERROR: %s\n", e.what());
    if (Log::dat()) std::fprintf(Log::dat(), "***ERROR: %s\n", e.what());
    return 2;
  }
  if (Log::dat()) std::fclose(Log::dat());
  return 0;
}
