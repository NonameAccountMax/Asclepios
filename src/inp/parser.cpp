#include "inp/parser.hpp"

#include <cctype>
#include <cmath>
#include <fstream>
#include <sstream>
#include <stdexcept>

#include "core/log.hpp"

namespace abq {

std::string upper(std::string s) {
  for (auto& c : s) c = (char)std::toupper((unsigned char)c);
  return s;
}
std::string trim(const std::string& s) {
  size_t a = 0, b = s.size();
  while (a < b && std::isspace((unsigned char)s[a])) ++a;
  while (b > a && std::isspace((unsigned char)s[b - 1])) --b;
  return s.substr(a, b - a);
}
bool is_number(const std::string& s) {
  if (s.empty()) return false;
  char* end = nullptr;
  std::strtod(s.c_str(), &end);
  return end && *end == '\0';
}
double to_double(const std::string& s, const Card* c) {
  std::string t = trim(s);
  if (t.empty()) return 0.0;
  // Fortran style exponents 1.0D-3
  for (auto& ch : t) if (ch == 'd' || ch == 'D') ch = 'e';
  char* end = nullptr;
  double v = std::strtod(t.c_str(), &end);
  if (!end || *end != '\0') throw std::runtime_error("bad number '" + s + "'" + (c ? " at " + c->where() : ""));
  return v;
}
int to_int(const std::string& s, const Card* c) {
  double v = to_double(s, c);
  return (int)std::lround(v);
}

std::string Card::getu(const std::string& k, const std::string& def) const { return upper(get(k, def)); }
double Card::getd(const std::string& k, double def) const { return has(k) ? to_double(get(k), this) : def; }
int Card::geti(const std::string& k, int def) const { return has(k) ? to_int(get(k), this) : def; }

namespace {

// --- tiny expression evaluator for *PARAMETER (numbers, + - * / **, parentheses, names) ---
struct Expr {
  const std::string& s;
  size_t p = 0;
  const std::map<std::string, double>& vars;
  Expr(const std::string& str, const std::map<std::string, double>& v) : s(str), vars(v) {}
  void ws() { while (p < s.size() && std::isspace((unsigned char)s[p])) ++p; }
  double primary() {
    ws();
    if (p < s.size() && s[p] == '(') { ++p; double v = sum(); ws(); if (p < s.size() && s[p] == ')') ++p; return v; }
    if (p < s.size() && (s[p] == '-' || s[p] == '+')) { char c = s[p++]; double v = power(); return c == '-' ? -v : v; }
    size_t b = p;
    if (p < s.size() && (std::isalpha((unsigned char)s[p]) || s[p] == '_')) {
      while (p < s.size() && (std::isalnum((unsigned char)s[p]) || s[p] == '_')) ++p;
      std::string name = upper(s.substr(b, p - b));
      auto it = vars.find(name);
      if (it == vars.end()) throw std::runtime_error("unknown parameter '" + name + "'");
      return it->second;
    }
    char* end = nullptr;
    double v = std::strtod(s.c_str() + p, &end);
    if (end == s.c_str() + p) throw std::runtime_error("bad parameter expression '" + s + "'");
    p = end - s.c_str();
    return v;
  }
  double power() {
    double v = primary();
    ws();
    if (p + 1 < s.size() && s[p] == '*' && s[p + 1] == '*') { p += 2; return std::pow(v, power()); }
    return v;
  }
  double prod() {
    double v = power();
    for (;;) {
      ws();
      if (p < s.size() && s[p] == '*' && !(p + 1 < s.size() && s[p + 1] == '*')) { ++p; v *= power(); }
      else if (p < s.size() && s[p] == '/') { ++p; v /= power(); }
      else return v;
    }
  }
  double sum() {
    double v = prod();
    for (;;) {
      ws();
      if (p < s.size() && s[p] == '+') { ++p; v += prod(); }
      else if (p < s.size() && s[p] == '-') { ++p; v -= prod(); }
      else return v;
    }
  }
};

std::string fmt_param(double v) {
  char buf[64];
  std::snprintf(buf, sizeof(buf), "%.17g", v);
  return buf;
}

struct Reader {
  std::vector<Card> cards;
  std::map<std::string, double> params;
  bool in_param_card = false;
  int depth = 0;

  std::string subst(const std::string& line) {
    if (line.find('<') == std::string::npos) return line;
    std::string out;
    for (size_t i = 0; i < line.size();) {
      if (line[i] == '<') {
        size_t j = line.find('>', i);
        if (j != std::string::npos) {
          std::string name = upper(trim(line.substr(i + 1, j - i - 1)));
          auto it = params.find(name);
          if (it != params.end()) { out += fmt_param(it->second); i = j + 1; continue; }
        }
      }
      out += line[i++];
    }
    return out;
  }

  static std::vector<std::string> split_commas(const std::string& s) {
    std::vector<std::string> t;
    std::string cur;
    for (char c : s) {
      if (c == ',') { t.push_back(trim(cur)); cur.clear(); }
      else cur += c;
    }
    t.push_back(trim(cur));
    return t;
  }

  void read(const std::string& path) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("cannot open input file '" + path + "'");
    if (++depth > 10) throw std::runtime_error("*INCLUDE nested too deeply");
    std::string dir;
    size_t sl = path.find_last_of('/');
    if (sl != std::string::npos) dir = path.substr(0, sl + 1);
    std::string raw;
    int lineno = 0;
    std::string pending_kw;  // keyword line continued by trailing comma
    int pending_line = 0;
    while (std::getline(in, raw)) {
      ++lineno;
      if (!raw.empty() && raw.back() == '\r') raw.pop_back();
      std::string line = trim(raw);
      if (line.empty()) continue;
      if (line.size() >= 2 && line[0] == '*' && line[1] == '*') continue;
      if (!pending_kw.empty()) {
        pending_kw += line;
        if (line.back() == ',') continue;
        start_card(pending_kw, path, pending_line, dir);
        pending_kw.clear();
        continue;
      }
      if (line[0] == '*') {
        if (line.back() == ',' && line.find(',') != line.size() - 1) {
          // "*Keyword, a=1," continuation (only if parameters exist already)
          pending_kw = line; pending_line = lineno; continue;
        }
        start_card(line, path, lineno, dir);
        continue;
      }
      if (cards.empty()) continue;  // data before any keyword: ignore
      if (in_param_card) {
        std::string l = line;
        size_t eq = l.find('=');
        if (eq == std::string::npos) continue;
        std::string name = upper(trim(l.substr(0, eq)));
        std::string expr = trim(l.substr(eq + 1));
        Expr e(expr, params);
        params[name] = e.sum();
        continue;
      }
      cards.back().data.push_back(split_commas(subst(line)));
    }
    --depth;
  }

  void start_card(const std::string& kwline_raw, const std::string& path, int lineno, const std::string& dir) {
    std::string kwline = subst(kwline_raw);
    auto toks = split_commas(kwline.substr(1));
    Card c;
    // normalize keyword: upper, collapse spaces
    std::string kw;
    bool sp = false;
    for (char ch : toks[0]) {
      if (std::isspace((unsigned char)ch)) { sp = true; continue; }
      if (sp && !kw.empty()) kw += ' ';
      sp = false;
      kw += (char)std::toupper((unsigned char)ch);
    }
    c.keyword = kw;
    for (size_t i = 1; i < toks.size(); ++i) {
      if (toks[i].empty()) continue;
      size_t eq = toks[i].find('=');
      if (eq == std::string::npos) {
        std::string k = upper(trim(toks[i]));
        // collapse internal whitespace in flags ("ALL EXTERIOR")
        c.params[k] = "";
      } else {
        std::string k = upper(trim(toks[i].substr(0, eq)));
        std::string k2;
        bool s2 = false;
        for (char ch : k) { if (ch == ' ') { s2 = true; continue; } if (s2 && !k2.empty()) k2 += ' '; s2 = false; k2 += ch; }
        c.params[k2] = trim(toks[i].substr(eq + 1));
      }
    }
    c.file = path;
    c.line = lineno;
    in_param_card = false;
    if (c.keyword == "INCLUDE") {
      std::string inc = c.get("INPUT");
      if (inc.size() >= 2 && (inc.front() == '"' || inc.front() == '\'')) inc = inc.substr(1, inc.size() - 2);
      std::string full = (inc.size() && inc[0] == '/') ? inc : dir + inc;
      std::ifstream test(full);
      if (!test) full = inc;  // relative to the working directory
      read(full);
      return;
    }
    if (c.keyword == "PARAMETER") { in_param_card = true; return; }
    cards.push_back(std::move(c));
  }
};

}  // namespace

std::vector<Card> parse_inp(const std::string& path) {
  Reader r;
  r.read(path);
  return std::move(r.cards);
}

}  // namespace abq
