#' Static CCE Estimator Internals
#'
#' Internal functions for the Pesaran (2006) Common Correlated Effects
#' estimator (both Mean Group and Pooled variants).
#'
#' @name cce_estimator
#' @keywords internal
NULL


# ──────────────────────────────────────────────────────────────────────────────
# .partial_out — remove CSA component from variables
# ──────────────────────────────────────────────────────────────────────────────

#' Partial out CSAs from regressors
#'
#' Projects out the CSA (cross-sectional average) component from a matrix
#' of variables using the Frisch-Waugh-Lovell projection:
#' `M_Z X = X - Z (Z'Z)^{-1} Z' X`
#'
#' @param X Numeric matrix of variables to be partialled.
#' @param Z Numeric matrix of CSA regressors (including intercept if present).
#' @return Numeric matrix of residuals after partialling out Z.
#' @keywords internal
.partial_out <- function(X, Z) {
  X <- as.matrix(X)
  Z <- as.matrix(Z)

  ZtZ <- crossprod(Z)
  rc <- rcond(ZtZ)
  if (is.na(rc) || rc < .Machine$double.eps^0.5) {
    ZtZ_inv <- .pinv(ZtZ)
  } else {
    ZtZ_inv <- solve(ZtZ)
  }

  # M_Z X = X - Z * (Z'Z)^{-1} * Z'X
  X - Z %*% (ZtZ_inv %*% crossprod(Z, X))
}


# ──────────────────────────────────────────────────────────────────────────────
# .pooled_vcov_pesaran — Pesaran (2006) non-parametric VCE for pooled coefs
# ──────────────────────────────────────────────────────────────────────────────

#' Pesaran (2006) non-parametric pooled variance
#'
#' For pooled (CCEP) coefficients, computes the non-parametric VCE:
#' `V_p = (1/N^2) sum_i (b_i - b_p)(b_i - b_p)'`
#' This is the same formula as MG variance but applied to the deviation
#' from the pooled estimate rather than the MG estimate.
#'
#' @param coef_list List of unit-level coefficient vectors.
#' @param b_pooled Numeric vector: the pooled coefficient.
#' @return Variance-covariance matrix.
#' @keywords internal
.pooled_vcov_pesaran <- function(coef_list, b_pooled) {
  N <- length(coef_list)
  B <- do.call(rbind, coef_list)
  B_centered <- sweep(B, 2, b_pooled)
  V <- crossprod(B_centered) / N^2
  colnames(V) <- rownames(V) <- names(b_pooled)
  V
}
