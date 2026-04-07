#' CS-ARDL Estimator Internals
#'
#' Internal functions for the Cross-Sectionally augmented ARDL (CS-ARDL)
#' estimator of Chudik et al. (2016). Estimates an ARDL(py, px) model
#' with CSAs and recovers long-run coefficients via delta method.
#'
#' @name csardl_estimator
#' @keywords internal
NULL


#' CS-ARDL long-run coefficient recovery
#'
#' Recovers long-run coefficients from ARDL(py, px) estimates:
#' `w2 = sum(b2_l, l=0..px) / (1 - sum(b1_l, l=1..py))`
#'
#' @param b Numeric vector: full coefficient vector from ARDL regression.
#' @param py Integer: number of lags of y.
#' @param px Integer: number of lags of x.
#' @param y_lag_names Character: names of lagged y coefficients in `b`.
#' @param x_lag_names Character: names of x and lagged x coefficients in `b`.
#' @param V Numeric matrix: variance-covariance of `b`.
#' @return A list with `lr_coef` and `lr_vcov`.
#' @keywords internal
.csardl_recover_lr <- function(b, py, px, y_lag_names, x_lag_names, V) {

  # g(b) = sum(b_x[0:px]) / (1 - sum(b_y[1:py]))
  g <- function(bvec) {
    b_y <- bvec[y_lag_names]
    b_x <- bvec[x_lag_names]
    sum(b_x) / (1 - sum(b_y))
  }

  # Extract relevant sub-vector and sub-matrix
  all_names <- c(y_lag_names, x_lag_names)
  b_sub <- b[all_names]
  V_sub <- V[all_names, all_names]

  dm <- .delta_method(g, b_sub, V_sub)

  list(lr_coef = dm$estimate, lr_vcov = dm$vcov)
}
