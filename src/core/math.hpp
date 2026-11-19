// Small fixed-size linear algebra used everywhere (double precision).
#pragma once
#include <array>
#include <cmath>
#include <algorithm>

namespace abq {

struct Vec3 {
  double x = 0, y = 0, z = 0;
  Vec3() = default;
  constexpr Vec3(double a, double b, double c) : x(a), y(b), z(c) {}
  double& operator[](int i) { return (&x)[i]; }
  double operator[](int i) const { return (&x)[i]; }
  Vec3 operator+(const Vec3& o) const { return {x + o.x, y + o.y, z + o.z}; }
  Vec3 operator-(const Vec3& o) const { return {x - o.x, y - o.y, z - o.z}; }
  Vec3 operator-() const { return {-x, -y, -z}; }
  Vec3 operator*(double s) const { return {x * s, y * s, z * s}; }
  Vec3 operator/(double s) const { return {x / s, y / s, z / s}; }
  Vec3& operator+=(const Vec3& o) { x += o.x; y += o.y; z += o.z; return *this; }
  Vec3& operator-=(const Vec3& o) { x -= o.x; y -= o.y; z -= o.z; return *this; }
  Vec3& operator*=(double s) { x *= s; y *= s; z *= s; return *this; }
};
inline Vec3 operator*(double s, const Vec3& v) { return v * s; }
inline double dot(const Vec3& a, const Vec3& b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
inline Vec3 cross(const Vec3& a, const Vec3& b) {
  return {a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x};
}
inline double norm(const Vec3& a) { return std::sqrt(dot(a, a)); }
inline double norm2(const Vec3& a) { return dot(a, a); }
inline Vec3 normalized(const Vec3& a) {
  double n = norm(a);
  return n > 0 ? a / n : Vec3{0, 0, 0};
}
inline Vec3 load3(const double* p) { return {p[0], p[1], p[2]}; }
inline void store3(double* p, const Vec3& v) { p[0] = v.x; p[1] = v.y; p[2] = v.z; }
inline void add3(double* p, const Vec3& v) { p[0] += v.x; p[1] += v.y; p[2] += v.z; }

// 3x3 matrix, row-major. Used for rotations / local frames: rows = local axes.
struct Mat3 {
  double m[3][3] = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
  static Mat3 identity() { Mat3 r; r.m[0][0] = r.m[1][1] = r.m[2][2] = 1; return r; }
  static Mat3 rows(const Vec3& a, const Vec3& b, const Vec3& c) {
    Mat3 r;
    for (int j = 0; j < 3; ++j) { r.m[0][j] = a[j]; r.m[1][j] = b[j]; r.m[2][j] = c[j]; }
    return r;
  }
  Vec3 row(int i) const { return {m[i][0], m[i][1], m[i][2]}; }
  Vec3 col(int j) const { return {m[0][j], m[1][j], m[2][j]}; }
  Vec3 operator*(const Vec3& v) const {
    return {m[0][0] * v.x + m[0][1] * v.y + m[0][2] * v.z, m[1][0] * v.x + m[1][1] * v.y + m[1][2] * v.z,
            m[2][0] * v.x + m[2][1] * v.y + m[2][2] * v.z};
  }
  Mat3 operator*(const Mat3& o) const {
    Mat3 r;
    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j) r.m[i][j] = m[i][0] * o.m[0][j] + m[i][1] * o.m[1][j] + m[i][2] * o.m[2][j];
    return r;
  }
  Mat3 operator+(const Mat3& o) const { Mat3 r; for (int i = 0; i < 9; ++i) (&r.m[0][0])[i] = (&m[0][0])[i] + (&o.m[0][0])[i]; return r; }
  Mat3 operator-(const Mat3& o) const { Mat3 r; for (int i = 0; i < 9; ++i) (&r.m[0][0])[i] = (&m[0][0])[i] - (&o.m[0][0])[i]; return r; }
  Mat3 operator*(double s) const { Mat3 r; for (int i = 0; i < 9; ++i) (&r.m[0][0])[i] = (&m[0][0])[i] * s; return r; }
  Mat3 T() const { Mat3 r; for (int i = 0; i < 3; ++i) for (int j = 0; j < 3; ++j) r.m[i][j] = m[j][i]; return r; }
  Vec3 Tmul(const Vec3& v) const {  // transpose(M) * v
    return {m[0][0] * v.x + m[1][0] * v.y + m[2][0] * v.z, m[0][1] * v.x + m[1][1] * v.y + m[2][1] * v.z,
            m[0][2] * v.x + m[1][2] * v.y + m[2][2] * v.z};
  }
  double det() const {
    return m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1]) - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0]) +
           m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);
  }
  Mat3 inverse() const {
    Mat3 r;
    double d = det();
    double id = 1.0 / d;
    r.m[0][0] = (m[1][1] * m[2][2] - m[1][2] * m[2][1]) * id;
    r.m[0][1] = (m[0][2] * m[2][1] - m[0][1] * m[2][2]) * id;
    r.m[0][2] = (m[0][1] * m[1][2] - m[0][2] * m[1][1]) * id;
    r.m[1][0] = (m[1][2] * m[2][0] - m[1][0] * m[2][2]) * id;
    r.m[1][1] = (m[0][0] * m[2][2] - m[0][2] * m[2][0]) * id;
    r.m[1][2] = (m[0][2] * m[1][0] - m[0][0] * m[1][2]) * id;
    r.m[2][0] = (m[1][0] * m[2][1] - m[1][1] * m[2][0]) * id;
    r.m[2][1] = (m[0][1] * m[2][0] - m[0][0] * m[2][1]) * id;
    r.m[2][2] = (m[0][0] * m[1][1] - m[0][1] * m[1][0]) * id;
    return r;
  }
};

// Rotation matrix from a rotation vector (Rodrigues).
inline Mat3 rotation_from_vector(const Vec3& th) {
  double a = norm(th);
  Mat3 R = Mat3::identity();
  if (a < 1e-14) {
    // first order
    R.m[0][1] = -th.z; R.m[0][2] = th.y; R.m[1][0] = th.z; R.m[1][2] = -th.x; R.m[2][0] = -th.y; R.m[2][1] = th.x;
    return R;
  }
  Vec3 k = th / a;
  double c = std::cos(a), s = std::sin(a), v = 1 - c;
  R.m[0][0] = c + k.x * k.x * v;       R.m[0][1] = k.x * k.y * v - k.z * s; R.m[0][2] = k.x * k.z * v + k.y * s;
  R.m[1][0] = k.y * k.x * v + k.z * s; R.m[1][1] = c + k.y * k.y * v;       R.m[1][2] = k.y * k.z * v - k.x * s;
  R.m[2][0] = k.z * k.x * v - k.y * s; R.m[2][1] = k.z * k.y * v + k.x * s; R.m[2][2] = c + k.z * k.z * v;
  return R;
}

// Hughes-Winget incremental rotation for a spin increment W = skew(w) (w = rotation vector increment).
// Q = (I - W/2)^-1 (I + W/2), exactly orthogonal.
inline Mat3 hughes_winget(const Vec3& w) {
  double w2 = norm2(w);
  double c = 2.0 / (4.0 + w2);
  Mat3 W;
  W.m[0][1] = -w.z; W.m[0][2] = w.y; W.m[1][0] = w.z; W.m[1][2] = -w.x; W.m[2][0] = -w.y; W.m[2][1] = w.x;
  // Q = I + c*(2W + W^2)  (closed form, Hughes & Winget 1980)
  Mat3 W2 = W * W;
  Mat3 Q = Mat3::identity();
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j) Q.m[i][j] += c * (2.0 * W.m[i][j] + W2.m[i][j]);
  return Q;
}

// Symmetric tensor in Abaqus Voigt order: 11, 22, 33, 12, 13, 23 (tensor shear components, not engineering).
struct Sym6 {
  double v[6] = {0, 0, 0, 0, 0, 0};
  double& operator[](int i) { return v[i]; }
  double operator[](int i) const { return v[i]; }
  double trace() const { return v[0] + v[1] + v[2]; }
  Mat3 full() const {
    Mat3 r;
    r.m[0][0] = v[0]; r.m[1][1] = v[1]; r.m[2][2] = v[2];
    r.m[0][1] = r.m[1][0] = v[3]; r.m[0][2] = r.m[2][0] = v[4]; r.m[1][2] = r.m[2][1] = v[5];
    return r;
  }
  static Sym6 from(const Mat3& a) {
    Sym6 s;
    s.v[0] = a.m[0][0]; s.v[1] = a.m[1][1]; s.v[2] = a.m[2][2];
    s.v[3] = 0.5 * (a.m[0][1] + a.m[1][0]); s.v[4] = 0.5 * (a.m[0][2] + a.m[2][0]); s.v[5] = 0.5 * (a.m[1][2] + a.m[2][1]);
    return s;
  }
};
inline double mises(const Sym6& s) {
  double a = s[0] - s[1], b = s[1] - s[2], c = s[2] - s[0];
  return std::sqrt(0.5 * (a * a + b * b + c * c) + 3.0 * (s[3] * s[3] + s[4] * s[4] + s[5] * s[5]));
}
// Rotate symmetric tensor: R * S * R^T
inline Sym6 rotate(const Mat3& R, const Sym6& s) { return Sym6::from(R * s.full() * R.T()); }

// Eigen-decomposition of a symmetric 2x2 (a, b; b, c): returns eigenvalues l1>=l2 and angle of first eigenvector.
inline void eig2sym(double a, double b, double c, double& l1, double& l2, double& cs, double& sn) {
  double m = 0.5 * (a + c), d = std::sqrt(0.25 * (a - c) * (a - c) + b * b);
  l1 = m + d; l2 = m - d;
  double th = 0.5 * std::atan2(2 * b, a - c);
  cs = std::cos(th); sn = std::sin(th);
}

// Symmetric 3x3 eigenvalues (Jacobi), used only for output (LE).
inline void eig3sym(const Mat3& A, double ev[3], Mat3* vecs = nullptr) {
  double a[3][3];
  double v[3][3] = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
  for (int i = 0; i < 3; ++i) for (int j = 0; j < 3; ++j) a[i][j] = A.m[i][j];
  for (int sweep = 0; sweep < 30; ++sweep) {
    double off = std::fabs(a[0][1]) + std::fabs(a[0][2]) + std::fabs(a[1][2]);
    if (off < 1e-15 * (std::fabs(a[0][0]) + std::fabs(a[1][1]) + std::fabs(a[2][2]) + 1e-300)) break;
    for (int p = 0; p < 2; ++p)
      for (int q = p + 1; q < 3; ++q) {
        if (std::fabs(a[p][q]) < 1e-300) continue;
        double th = 0.5 * (a[q][q] - a[p][p]) / a[p][q];
        double t = (th >= 0 ? 1.0 : -1.0) / (std::fabs(th) + std::sqrt(th * th + 1));
        double c = 1 / std::sqrt(t * t + 1), s = t * c;
        for (int k = 0; k < 3; ++k) {
          double akp = a[k][p], akq = a[k][q];
          a[k][p] = c * akp - s * akq; a[k][q] = s * akp + c * akq;
        }
        for (int k = 0; k < 3; ++k) {
          double apk = a[p][k], aqk = a[q][k];
          a[p][k] = c * apk - s * aqk; a[q][k] = s * apk + c * aqk;
        }
        for (int k = 0; k < 3; ++k) {
          double vkp = v[k][p], vkq = v[k][q];
          v[k][p] = c * vkp - s * vkq; v[k][q] = s * vkp + c * vkq;
        }
      }
  }
  for (int i = 0; i < 3; ++i) ev[i] = a[i][i];
  if (vecs) for (int i = 0; i < 3; ++i) for (int j = 0; j < 3; ++j) vecs->m[i][j] = v[i][j];
}

}  // namespace abq
