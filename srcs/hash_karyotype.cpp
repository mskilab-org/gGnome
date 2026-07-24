#include <Rcpp.h>
#include <vector>
#include <algorithm>
#include <utility>
#include <cstdint>
#include <cstdio>
using namespace Rcpp;

// Canonical-karyotype hash, C++ replacement for hash_snodelist (R/graphstats.R).
//   1. For each circular walk: take lex-min over the full rotation+RC orbit
//      by comparing Booth(W) with Booth(rc(W)) and keeping the smaller.
//   2. For each linear walk: replace with min(W, rc(W)) lex-wise.
//      rc(w)[i] = -w[n-1-i].
//   3. Build a multiset that doubles each entry with its rc (kept from the
//      original hash_snodelist; redundant under (1)+(2) but cheap).
//   4. Sort the multiset lex.
//   5. FNV-1a 64-bit hash over the bytes; return as a 16-char hex string.
//
// The earlier version of this code, and R's sort_snodes today, compute
// rc(Booth(W)) instead of Booth(rc(W)) in step (1). The two differ by a
// cyclic rotation in general, so circular karyotypes presented in different
// rotations of the same orbit ended up with different canonical forms. Fixed
// here and in R's sort_snodes; both paths now produce the genuine canonical
// representative of the rotation+RC orbit.

// ----------------------------------------------------------------------------
// Booth's algorithm: index of the lex-min cyclic rotation of s.
// Linear time, no auxiliary suffix array.
// ----------------------------------------------------------------------------
static int booth_rotation_start(const std::vector<int>& s) {
  const int n = (int)s.size();
  if (n <= 1) return 0;
  std::vector<int> s2(2 * n);
  for (int i = 0; i < n; i++) { s2[i] = s[i]; s2[i + n] = s[i]; }
  int i = 0, j = 1, k = 0;
  while (i < n && j < n && k < n) {
    int a = s2[i + k], b = s2[j + k];
    if (a == b) {
      k++;
    } else if (a > b) {
      i = i + k + 1;
      if (i <= j) i = j + 1;
      k = 0;
    } else {
      j = j + k + 1;
      if (j <= i) j = i + 1;
      k = 0;
    }
  }
  int start = std::min(i, j);
  if (start >= n) start -= n;
  return start;
}

static void apply_booth_rotation(std::vector<int>& s) {
  const int n = (int)s.size();
  if (n <= 1) return;
  int start = booth_rotation_start(s);
  if (start > 0 && start < n) {
    std::rotate(s.begin(), s.begin() + start, s.end());
  }
}

// rc(s)[i] = -s[n-1-i]
static std::vector<int> rc_walk(const std::vector<int>& s) {
  std::vector<int> r(s.size());
  for (size_t i = 0; i < s.size(); i++) r[i] = -s[s.size() - 1 - i];
  return r;
}

// ----------------------------------------------------------------------------
// FNV-1a 64-bit byte stream hash.
// ----------------------------------------------------------------------------
struct FNV1a64 {
  uint64_t h;
  FNV1a64() : h(0xcbf29ce484222325ULL) {}
  void update_byte(uint8_t b) { h = (h ^ (uint64_t)b) * 0x100000001b3ULL; }
  void update_int(int v) {
    const uint8_t* p = reinterpret_cast<const uint8_t*>(&v);
    for (size_t i = 0; i < sizeof(v); i++) update_byte(p[i]);
  }
};

// ----------------------------------------------------------------------------
// Public entry point.
// ----------------------------------------------------------------------------

// [[Rcpp::depends(Rcpp)]]
// [[Rcpp::export]]
std::string hash_karyotype_cpp(List snode_id, LogicalVector circular) {
  const int K = snode_id.size();
  if (K == 0) return std::string("0000000000000000");
  if (K != circular.size())
    Rcpp::stop("hash_karyotype_cpp: snode_id and circular have different lengths");

  // 1. Pull walks; canonicalize each.
  //    Circular: pick min(Booth(W), Booth(rc(W))) over the rotation+RC orbit.
  //    Linear:   pick min(W, rc(W)) lex.
  std::vector<std::vector<int>> walks(K);
  std::vector<bool> circ(K);
  for (int k = 0; k < K; k++) {
    IntegerVector iv = snode_id[k];
    std::vector<int> w(iv.begin(), iv.end());
    circ[k] = (bool)circular[k];

    if (circ[k]) {
      std::vector<int> w_rc = rc_walk(w);
      apply_booth_rotation(w);                 // w    = Booth(W)
      apply_booth_rotation(w_rc);              // w_rc = Booth(rc(W))
      walks[k] = (w_rc < w) ? std::move(w_rc) : std::move(w);
    } else {
      std::vector<int> w_rc = rc_walk(w);
      walks[k] = (w_rc < w) ? std::move(w_rc) : std::move(w);
    }
  }

  // 2. Build doubled multiset: (walk, circ) + (rc(walk), circ).
  std::vector<std::pair<std::vector<int>, bool>> entries;
  entries.reserve(2 * K);
  for (int k = 0; k < K; k++) {
    entries.emplace_back(walks[k], circ[k]);
    entries.emplace_back(rc_walk(walks[k]), circ[k]);
  }

  // 3. Sort lex.
  std::sort(entries.begin(), entries.end(),
    [](const std::pair<std::vector<int>, bool>& a,
       const std::pair<std::vector<int>, bool>& b) {
      if (a.first != b.first) return a.first < b.first;
      return a.second < b.second;
    });

  // 4. Hash.
  FNV1a64 hasher;
  for (const auto& e : entries) {
    hasher.update_byte(0xFE);                       // entry separator
    hasher.update_byte(e.second ? (uint8_t)'C' : (uint8_t)'L');
    for (int v : e.first) hasher.update_int(v);
  }

  // 5. Format as 16-char lowercase hex.
  char buf[17];
  std::snprintf(buf, sizeof(buf), "%016llx", (unsigned long long)hasher.h);
  return std::string(buf);
}
