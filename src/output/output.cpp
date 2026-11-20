#include "output/output.hpp"

#include <sys/stat.h>

#include <cstdarg>
#include <cstdint>
#include <cstring>

#include "core/log.hpp"
#include "solver/analysis.hpp"

namespace abq {

namespace {

std::string base64(const unsigned char* data, size_t n) {
  static const char* T = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";
  std::string out;
  out.reserve((n + 2) / 3 * 4);
  size_t i = 0;
  for (; i + 2 < n; i += 3) {
    uint32_t v = (data[i] << 16) | (data[i + 1] << 8) | data[i + 2];
    out += T[(v >> 18) & 63]; out += T[(v >> 12) & 63]; out += T[(v >> 6) & 63]; out += T[v & 63];
  }
  if (i < n) {
    uint32_t v = data[i] << 16;
    if (i + 1 < n) v |= data[i + 1] << 8;
    out += T[(v >> 18) & 63]; out += T[(v >> 12) & 63];
    out += (i + 1 < n) ? T[(v >> 6) & 63] : '=';
    out += '=';
  }
  return out;
}

template <class T>
void write_array(FILE* f, const char* type, const char* name, int ncomp, const std::vector<T>& v) {
  std::fprintf(f, "<DataArray type=\"%s\" Name=\"%s\" NumberOfComponents=\"%d\" format=\"binary\">\n", type, name, ncomp);
  std::vector<unsigned char> buf(4 + v.size() * sizeof(T));
  uint32_t nb = (uint32_t)(v.size() * sizeof(T));
  std::memcpy(buf.data(), &nb, 4);
  if (!v.empty()) std::memcpy(buf.data() + 4, v.data(), v.size() * sizeof(T));
  std::string s = base64(buf.data(), buf.size());
  std::fwrite(s.data(), 1, s.size(), f);
  std::fprintf(f, "\n</DataArray>\n");
}

int vtk_type(ElType t, int nen) {
  switch (t) {
    case ElType::C3D8R: case ElType::C3D8: case ElType::C3D8I: return 12;
    case ElType::C3D4: return 10;
    case ElType::S4R: case ElType::S4: return 9;
    case ElType::S3R: case ElType::S3RS: return 5;
    default: return 3;  // line
  }
}

}  // namespace

OutputWriter::OutputWriter(Analysis& A, const std::string& job) : A_(A), job_(job) {
  dir_ = job + "_vtk";
  mkdir(dir_.c_str(), 0755);
  hist_ = std::fopen((job + "_history.csv").c_str(), "w");
  if (hist_) std::fprintf(hist_, "step,increment,step_time,total_time,dt,ALLKE,ALLIE,ALLWK,ALLVD,ALLPD,ALLAE,ALLDMD,ALLFD,ALLCW,ETOTAL\n");
  sta_ = std::fopen((job + ".sta").c_str(), "w");
}

OutputWriter::~OutputWriter() { close(); }

void OutputWriter::close() {
  if (hist_) { std::fclose(hist_); hist_ = nullptr; }
  if (sta_) { std::fclose(sta_); sta_ = nullptr; }
  if (rfh_) { std::fclose(rfh_); rfh_ = nullptr; }
  if (monh_) { std::fclose(monh_); monh_ = nullptr; }
}

// total reaction of the node set of the history request, one line per requested increment (<job>_rf.csv)
void OutputWriter::rf_history(double total_time, const double r[3]) {
  if (!rfh_) {
    rfh_ = std::fopen((job_ + "_rf.csv").c_str(), "w");
    if (!rfh_) return;
    std::fprintf(rfh_, "total_time,RF1,RF2,RF3\n");
  }
  std::fprintf(rfh_, "%.12g,%.9g,%.9g,%.9g\n", total_time, r[0], r[1], r[2]);
}

// displacements and rotations of the nodes of the history request (<job>_mon.csv, one line per node and increment)
void OutputWriter::mon_history(double total_time, int label, const double* u6) {
  if (!monh_) {
    monh_ = std::fopen((job_ + "_mon.csv").c_str(), "w");
    if (!monh_) return;
    std::fprintf(monh_, "total_time,node,U1,U2,U3,UR1,UR2,UR3\n");
  }
  std::fprintf(monh_, "%.12g,%d,%.9g,%.9g,%.9g,%.9g,%.9g,%.9g\n", total_time, label, u6[0], u6[1], u6[2], u6[3], u6[4], u6[5]);
}

void OutputWriter::sta_header(int step, const std::string& name, bool explicit_step) {
  if (!sta_) return;
  std::fprintf(sta_, "\n STEP %d  (%s)\n", step + 1, name.c_str());
  if (explicit_step)
    std::fprintf(sta_, "  STEP     INCREMENT     TOTAL     WALLCLOCK   STABLE     CRITICAL    KINETIC      TOTAL      PERCENT\n"
                       "                  TIME       TIME       TIME     INCREMENT  ELEMENT      ENERGY       ENERGY     CHNG MASS\n");
  else
    std::fprintf(sta_, "  STEP  INC ATT SEVERE EQUIL TOTAL  TOTAL      STEP       INC OF\n"
                       "              DISCON ITERS ITERS  TIME/     TIME/LPF   TIME/LPF\n");
  std::fflush(sta_);
}

void OutputWriter::sta_line(const char* fmt, ...) {
  if (!sta_) return;
  va_list ap;
  va_start(ap, fmt);
  std::vfprintf(sta_, fmt, ap);
  va_end(ap);
  std::fputc('\n', sta_);
  std::fflush(sta_);
}

void OutputWriter::history(int step, int inc, double st, double tt, double dt, double ke) {
  if (!hist_) return;
  const Energies& E = A_.en;
  std::fprintf(hist_, "%d,%d,%.9g,%.9g,%.6g,%.9g,%.9g,%.9g,%.9g,%.9g,%.9g,%.9g,%.9g,%.9g,%.9g\n", step + 1, inc, st, tt, dt, ke,
               E.internal, E.external, E.viscous, E.plastic, E.artificial, E.damping, E.friction, E.mass_scaling_work,
               E.total(ke));
  std::fflush(hist_);
}

void OutputWriter::frame(int step, const std::string& step_name, double st, double tt) {
  const Model& m = A_.m;
  const int nn = A_.nn;
  char fname[512];
  const size_t js = job_.find_last_of('/');
  const std::string jbase = js == std::string::npos ? job_ : job_.substr(js + 1);
  std::snprintf(fname, sizeof(fname), "%s/%s_%05d.vtu", dir_.c_str(), jbase.c_str(), nframes_);
  FILE* f = std::fopen(fname, "w");
  if (!f) { warn("cannot write %s", fname); return; }
  // cells
  std::vector<int32_t> conn, offs;
  std::vector<uint8_t> types;
  std::vector<int32_t> labels, etype;
  std::vector<std::vector<float>> cdata(N_OUT);
  int off = 0;
  for (auto& gp : A_.groups) {
    auto& G = *gp;
    for (int e = 0; e < G.size(); ++e) {
      int nen = G.nen;
      int vt = vtk_type(G.type, nen);
      int nvis = vt == 3 ? 2 : nen;
      for (int I = 0; I < nvis; ++I) conn.push_back(G.conn[e * nen + I]);
      off += nvis;
      offs.push_back(off);
      types.push_back((uint8_t)vt);
      labels.push_back(m.elems[G.ids[e]].label);
      etype.push_back((int32_t)G.type);
      for (int c = 0; c < N_OUT; ++c) cdata[c].push_back(G.out[e * N_OUT + c]);
    }
  }
  int ncell = (int)types.size();
  std::fprintf(f, "<?xml version=\"1.0\"?>\n<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\" header_type=\"UInt32\">\n");
  std::fprintf(f, "<UnstructuredGrid>\n<Piece NumberOfPoints=\"%d\" NumberOfCells=\"%d\">\n", nn, ncell);
  std::vector<float> pts(3 * nn), U(3 * nn), V(3 * nn), RF(3 * nn), UR(3 * nn), RM(3 * nn);
  std::vector<float> T(nn);
  std::vector<int32_t> nl(nn);
  for (int i = 0; i < nn; ++i) {
    for (int k = 0; k < 3; ++k) {
      pts[3 * i + k] = (float)(A_.X[3 * i + k] + A_.u[6 * i + k]);
      U[3 * i + k] = (float)A_.u[6 * i + k];
      V[3 * i + k] = (float)A_.v[6 * i + k];
      RF[3 * i + k] = (float)A_.rf[6 * i + k];
      RM[3 * i + k] = (float)A_.rf[6 * i + 3 + k];
      UR[3 * i + k] = (float)A_.u[6 * i + 3 + k];
    }
    T[i] = (float)A_.temp[i];
    nl[i] = m.node_label[i];
  }
  std::fprintf(f, "<Points>\n");
  write_array(f, "Float32", "Points", 3, pts);
  std::fprintf(f, "</Points>\n<Cells>\n");
  write_array(f, "Int32", "connectivity", 1, conn);
  write_array(f, "Int32", "offsets", 1, offs);
  write_array(f, "UInt8", "types", 1, types);
  std::fprintf(f, "</Cells>\n<PointData Vectors=\"U\">\n");
  write_array(f, "Float32", "U", 3, U);
  write_array(f, "Float32", "V", 3, V);
  write_array(f, "Float32", "RF", 3, RF);
  write_array(f, "Float32", "UR", 3, UR);
  write_array(f, "Float32", "RM", 3, RM);
  if (A_.has_temp) write_array(f, "Float32", "NT11", 1, T);
  write_array(f, "Int32", "NodeLabel", 1, nl);
  std::fprintf(f, "</PointData>\n<CellData>\n");
  for (int c = 0; c < N_OUT; ++c) write_array(f, "Float32", out_channel_name(c), 1, cdata[c]);
  write_array(f, "Int32", "ElementLabel", 1, labels);
  write_array(f, "Int32", "ElementType", 1, etype);
  std::fprintf(f, "</CellData>\n</Piece>\n</UnstructuredGrid>\n</VTKFile>\n");
  std::fclose(f);
  char entry[600];
  const char* base = std::strrchr(fname, '/');
  const size_t ds = dir_.find_last_of('/');
  const std::string dbase = ds == std::string::npos ? dir_ : dir_.substr(ds + 1);
  std::snprintf(entry, sizeof(entry), "<DataSet timestep=\"%.9g\" group=\"\" part=\"0\" file=\"%s/%s\"/>", tt, dbase.c_str(), base ? base + 1 : fname);
  pvd_entries_.push_back(entry);
  ++nframes_;
  write_pvd();
}

void OutputWriter::write_pvd() {
  FILE* f = std::fopen((job_ + ".pvd").c_str(), "w");
  if (!f) return;
  std::fprintf(f, "<?xml version=\"1.0\"?>\n<VTKFile type=\"Collection\" version=\"0.1\">\n<Collection>\n");
  for (auto& e : pvd_entries_) std::fprintf(f, "%s\n", e.c_str());
  std::fprintf(f, "</Collection>\n</VTKFile>\n");
  std::fclose(f);
}

}  // namespace abq
