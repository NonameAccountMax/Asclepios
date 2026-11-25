#pragma once
#include "solver/analysis.hpp"

namespace abq {
void run_explicit_step(Analysis& A, const Step& st, int si);
void run_static_step(Analysis& A, const Step& st, int si);
// connector stops (stops.cpp): presence of supported stops; kinematic correction of the accelerations a so that the
// predicted configuration X + u + dpred respects them (c: displacement per unit acceleration correction), forces in fk
bool has_connector_stops(const Model& m);
void connector_stops(const Analysis& A, double c, std::vector<double>& dpred, std::vector<double>& a, std::vector<double>& fk);
}  // namespace abq
