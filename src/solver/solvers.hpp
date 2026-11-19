#pragma once
#include "solver/analysis.hpp"

namespace abq {
void run_explicit_step(Analysis& A, const Step& st, int si);
void run_static_step(Analysis& A, const Step& st, int si);
}  // namespace abq
