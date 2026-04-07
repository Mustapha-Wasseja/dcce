// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// Single-unit OLS: returns list(b, V, e, r2, df_resid, sigma2, rank_ok)
// X: T x K design matrix (already includes constant if desired)
// y: T x 1 response vector
// Column names, if any, are attached on the R side (see .run_unit_loop).
// [[Rcpp::export(.unit_ols_cpp)]]
List unit_ols_cpp(const arma::mat& X, const arma::vec& y) {
  int T_val = X.n_rows;
  int K     = X.n_cols;
  int df_resid = T_val - K;

  if (df_resid <= 0) {
    arma::vec b_na(K, fill::value(NA_REAL));
    arma::mat V_na(K, K, fill::value(NA_REAL));
    arma::vec e_na(T_val, fill::value(NA_REAL));
    return List::create(
      _["b"]        = wrap(b_na),
      _["V"]        = wrap(V_na),
      _["e"]        = wrap(e_na),
      _["r2"]       = NA_REAL,
      _["df_resid"] = df_resid,
      _["sigma2"]   = NA_REAL,
      _["rank_ok"]  = false
    );
  }

  // Normal equations
  arma::mat XtX = X.t() * X;
  arma::vec Xty = X.t() * y;
  arma::vec b(K);
  arma::mat XtX_inv;
  bool ok = true;

  // Try symmetric positive-definite solve first
  bool solved = arma::solve(b, XtX, Xty, arma::solve_opts::likely_sympd);
  if (!solved) {
    // Fall back to pseudoinverse for rank-deficient systems
    XtX_inv = arma::pinv(XtX);
    b = XtX_inv * Xty;
    ok = false;
  } else {
    // Invert via solve(XtX, I) for consistency with R implementation
    arma::mat I_k = arma::eye<arma::mat>(K, K);
    if (!arma::solve(XtX_inv, XtX, I_k, arma::solve_opts::likely_sympd)) {
      XtX_inv = arma::pinv(XtX);
      ok = false;
    }
  }

  arma::vec e      = y - X * b;
  double    rss    = arma::dot(e, e);
  double    ybar   = arma::mean(y);
  arma::vec y_dm   = y - ybar;
  double    tss    = arma::dot(y_dm, y_dm);
  double    r2     = (tss > 0) ? 1.0 - rss / tss : NA_REAL;
  double    sigma2 = rss / df_resid;
  arma::mat V      = sigma2 * XtX_inv;

  return List::create(
    _["b"]        = wrap(b),
    _["V"]        = wrap(V),
    _["e"]        = wrap(e),
    _["r2"]       = r2,
    _["df_resid"] = df_resid,
    _["sigma2"]   = sigma2,
    _["rank_ok"]  = ok
  );
}


// Batch OLS: runs unit_ols_cpp for each element of a named list of
// list(X, y) pairs, returning a list of the per-unit results keyed by
// the same names.
// [[Rcpp::export(.batch_ols_cpp)]]
List batch_ols_cpp(const List& panel_list) {
  int N = panel_list.size();
  List results(N);
  CharacterVector unit_names = panel_list.names();

  for (int i = 0; i < N; i++) {
    List unit_data = panel_list[i];
    arma::mat X = as<arma::mat>(unit_data["X"]);
    arma::vec y = as<arma::vec>(unit_data["y"]);
    results[i]  = unit_ols_cpp(X, y);
  }
  if (unit_names.size() == N) {
    results.names() = unit_names;
  }
  return results;
}
