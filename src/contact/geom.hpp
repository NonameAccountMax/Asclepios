// Geometric helpers shared by the contact algorithms.
#pragma once
#include "core/math.hpp"

namespace abq {

// barycentric coordinates of the projection of p on the plane of the triangle abc (outside [0, 1] beyond its edges)
inline void plane_coords(const Vec3& p, const Vec3& a, const Vec3& b, const Vec3& c, double w[3]) {
  const Vec3 e1 = b - a, e2 = c - a, d = p - a;
  const double d11 = dot(e1, e1), d12 = dot(e1, e2), d22 = dot(e2, e2), det = d11 * d22 - d12 * d12;
  if (!(det > 0)) { w[0] = 1; w[1] = w[2] = 0; return; }
  const double p1 = dot(d, e1), p2 = dot(d, e2);
  w[1] = (d22 * p1 - d12 * p2) / det;
  w[2] = (d11 * p2 - d12 * p1) / det;
  w[0] = 1 - w[1] - w[2];
}


// closest point on a triangle (Ericson), barycentric weights; interior = projection inside
inline Vec3 closest_on_tri(const Vec3& p, const Vec3& a, const Vec3& b, const Vec3& c, double w[3], bool& interior) {
  Vec3 ab = b - a, ac = c - a, ap = p - a;
  double d1 = dot(ab, ap), d2 = dot(ac, ap);
  interior = false;
  if (d1 <= 0 && d2 <= 0) { w[0] = 1; w[1] = w[2] = 0; return a; }
  Vec3 bp = p - b;
  double d3 = dot(ab, bp), d4 = dot(ac, bp);
  if (d3 >= 0 && d4 <= d3) { w[1] = 1; w[0] = w[2] = 0; return b; }
  double vc = d1 * d4 - d3 * d2;
  if (vc <= 0 && d1 >= 0 && d3 <= 0) { double v = d1 / (d1 - d3); w[0] = 1 - v; w[1] = v; w[2] = 0; return a + ab * v; }
  Vec3 cp = p - c;
  double d5 = dot(ab, cp), d6 = dot(ac, cp);
  if (d6 >= 0 && d5 <= d6) { w[2] = 1; w[0] = w[1] = 0; return c; }
  double vb = d5 * d2 - d1 * d6;
  if (vb <= 0 && d2 >= 0 && d6 <= 0) { double v = d2 / (d2 - d6); w[0] = 1 - v; w[2] = v; w[1] = 0; return a + ac * v; }
  double va = d3 * d6 - d5 * d4;
  if (va <= 0 && (d4 - d3) >= 0 && (d5 - d6) >= 0) {
    double v = (d4 - d3) / ((d4 - d3) + (d5 - d6));
    w[1] = 1 - v; w[2] = v; w[0] = 0;
    return b + (c - b) * v;
  }
  double den = 1.0 / (va + vb + vc);
  double v = vb * den, ww = vc * den;
  w[0] = 1 - v - ww; w[1] = v; w[2] = ww;
  interior = true;
  return a + ab * v + ac * ww;
}


}  // namespace abq
