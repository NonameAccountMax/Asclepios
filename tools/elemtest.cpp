// Element check: finite-difference stiffness of the first element of a deck; symmetry and eigenvalues.
#include <Eigen/Dense>
#include <cstdio>
#include "element/element.hpp"
#include "model/model.hpp"
using namespace abq;
int main(int argc, char** argv) {
  Model m = build_model(argv[1]);
  auto groups = make_element_groups(m);
  auto& G = *groups[0];
  int nen = G.nen, nd = 6 * nen;
  std::vector<double> X(3 * nen), du(nd, 0.0), f0(nd), f1(nd);
  for (int I = 0; I < nen; ++I) for (int k = 0; k < 3; ++k) X[3 * I + k] = m.X[G.conn[I]][k];
  ElemIn in; in.X = X.data(); in.xo = X.data(); in.du = du.data(); in.dt = 1e-7; in.dynamic = false; in.nlgeom = argc > 2;
  ElemOut out; out.f = f0.data();
  G.compute(0, in, out, false);
  Eigen::MatrixXd K(nd, nd);
  double h = 1e-7;
  for (int j = 0; j < nd; ++j) {
    du[j] = h; out.f = f1.data(); G.compute(0, in, out, false);
    du[j] = -h; std::vector<double> f2(nd); out.f = f2.data(); G.compute(0, in, out, false);
    du[j] = 0;
    for (int i = 0; i < nd; ++i) K(i, j) = (f1[i] - f2[i]) / (2 * h);
  }
  double asym = (K - K.transpose()).norm() / K.norm();
  printf("type %s nen %d  |K-K^T|/|K| = %.3e  dt_stable %.4e\n", eltype_name(G.type), nen, asym, out.dt_stable);
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(0.5 * (K + K.transpose()));
  printf("eigenvalues:"); for (int i = 0; i < nd; ++i) printf(" %.3e", es.eigenvalues()[i]); printf("\n");
  std::vector<double> mt(nen), mr(nen);
  G.mass(0, X.data(), mt.data(), mr.data());
  Eigen::VectorXd minv(nd);
  for (int I = 0; I < nen; ++I) for (int k = 0; k < 6; ++k) { double mm = k < 3 ? mt[I] : mr[I]; minv(6*I+k) = mm > 0 ? 1/std::sqrt(mm) : 0; }
  Eigen::MatrixXd A = minv.asDiagonal() * (0.5*(K+K.transpose())) * minv.asDiagonal();
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es2(A);
  double lmax = es2.eigenvalues()[nd-1];
  printf("mass t %.3e rot %.3e ; omega_max %.4e -> dt_crit %.4e (element reports %.4e)\n", mt[0], mr[0], std::sqrt(lmax), 2/std::sqrt(lmax), out.dt_stable);
  Eigen::VectorXd v = es2.eigenvectors().col(nd-1);
  printf("critical mode:"); for (int i = 0; i < nd; ++i) printf(" %.2f", v(i)); printf("\n");
  if (asym > 1e-6) {
    for (int i = 0; i < nd; ++i) { for (int j = 0; j < nd; ++j) printf("%9.2e ", K(i,j)-K(j,i)); printf("\n"); }
  }
}
