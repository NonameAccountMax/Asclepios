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
  // stable increment bound of the current contacts (1/W of the stiffest pair; huge without contact)
  double dt_limit() const;
  // restart (after define() of the same step)
  void save(BinOut& o) const;
  void load(BinIn& in);

 private:
  struct Impl;
  Impl* p_ = nullptr;
  Analysis& A_;
  bool active_ = false;
  int ncontacts_ = 0;
};

}  // namespace abq
