// General contact (explicit): penalty node-to-facet and edge-to-edge with shell/beam thickness offsets, softened
// pressure-overclosure laws, exclusions and Coulomb friction.
#pragma once
#include <vector>

#include "core/binio.hpp"
#include "model/model.hpp"

namespace abq {

struct Analysis;
struct Energies;

class GeneralContact {
 public:
  explicit GeneralContact(Analysis& A);
  ~GeneralContact();
  // (re)define from a step's *CONTACT INCLUSIONS; returns false if nothing to do
  bool define(const Step& st);
  bool active() const { return active_; }
  // penalty forces for the end-of-increment configuration (X + u + du); dt = current increment
  void compute(double dt, std::vector<double>& f, Energies& en, int increment);
  int ncontacts() const { return ncontacts_; }
  // implicit procedures: contact points of the configuration X + u + du, without state: for each slave node the
  // deepest penetrating facet (shells two-sided, on the side of the start of the increment; solid faces and rigid
  // facets one-sided), with the nodes carrying the force and their weights; beam segments and edges not included
  struct Point { int slave; int nm; int mn[4]; double w[4]; Vec3 n; double pen; int master; double area; };
  // near > 0: also the points within that distance of contact (pen > -near), for the tangent of the implicit solver
  void implicit_points(std::vector<Point>& pts, double near = 0);
  double friction() const;
  const SurfaceInteraction* softened() const;  // softened pressure-overclosure law, null for hard contact
  double clearance() const;
  double mean_facet_length() const;
  // stable increment bound of the current contacts (1/W of the stiffest pair; huge without contact)
  double dt_limit() const;
  // restart (after define() of the same step)
  void save(BinOut& o) const;
  void load(BinIn& in);

 private:
  struct Impl;
  Impl* p_ = nullptr;
  void prepare();
  Analysis& A_;
  bool active_ = false;
  int ncontacts_ = 0;
};

}  // namespace abq
