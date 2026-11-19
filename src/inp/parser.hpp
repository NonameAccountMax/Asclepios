// Abaqus-style keyword input reader: comments, *INCLUDE, *PARAMETER substitution.
#pragma once
#include <map>
#include <string>
#include <vector>

namespace abq {

struct Card {
  std::string keyword;                         // upper case, single spaces, without '*'
  std::map<std::string, std::string> params;   // UPPER key -> value as written (trimmed); "" if no '='
  std::vector<std::vector<std::string>> data;  // one vector of trimmed tokens per data line
  std::string file;
  int line = 0;

  bool has(const std::string& k) const { return params.count(k) != 0; }
  std::string get(const std::string& k, const std::string& def = "") const {
    auto it = params.find(k);
    return it == params.end() ? def : it->second;
  }
  std::string getu(const std::string& k, const std::string& def = "") const;  // upper-cased value
  double getd(const std::string& k, double def) const;
  int geti(const std::string& k, int def) const;
  std::string where() const { return file + ":" + std::to_string(line); }
};

std::vector<Card> parse_inp(const std::string& path);

// helpers shared by the model builder
std::string upper(std::string s);
std::string trim(const std::string& s);
double to_double(const std::string& s, const Card* c = nullptr);
int to_int(const std::string& s, const Card* c = nullptr);
bool is_number(const std::string& s);

}  // namespace abq
