// Contact pairs (explicit, *CONTACT PAIR): kinematic hard contact by prediction and correction of the increment.
#pragma once
#include <vector>

#include "model/model.hpp"

namespace abq {

struct Analysis;

class ContactPairs {
 public:
  explicit ContactPairs(Analysis& A);
  ~ContactPairs();
  // (re)define from a step's *CONTACT PAIR options; returns false if there is none
  bool define(const Step& st);
  bool active() const { return active_; }
  // corrects the accelerations a so that the slave nodes do not penetrate the master surfaces in the predicted
  // configuration X + u + du (du: displacement from the stored u to the predicted position, moved by c x a correction
  // x); fk receives the contact forces of the corrections (mass x correction, on prescribed degrees of freedom the
  // force taken by the reaction)
  void correct(double c, std::vector<double>& du, std::vector<double>& a, std::vector<double>& fk);
  int ncontacts() const { return ncontacts_; }

 private:
  struct Impl;
  Impl* p_ = nullptr;
  Analysis& A_;
  bool active_ = false;
  int ncontacts_ = 0;
};

}  // namespace abq
