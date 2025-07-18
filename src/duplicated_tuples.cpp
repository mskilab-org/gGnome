#include <Rcpp.h>
#include <unordered_set>
using namespace Rcpp;

// [[Rcpp::export]]
LogicalVector duplicated_tuples_source(List dots) {
  int n = as<NumericVector>(dots[0]).size();
  int m = dots.size();

  LogicalVector is_duplicate(n, false);
  std::unordered_set<double> seen;

  for (int i = 0; i < n; ++i) {
    bool any_seen = false;

    for (int j = 0; j < m; ++j) {
      NumericVector col = dots[j];
      double val = col[i];
      if (seen.count(val)) {
        any_seen = true;
        break;
      }
    }

    if (any_seen) {
      is_duplicate[i] = true;
    } else {
      for (int j = 0; j < m; ++j) {
        seen.insert(as<NumericVector>(dots[j])[i]);
      }
    }
  }

  return is_duplicate;
}