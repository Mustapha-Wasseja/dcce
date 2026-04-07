#' PMG / Error Correction Model Internals
#'
#' Internal functions for the Pooled Mean Group (PMG) estimator in
#' error-correction form (Shin, Pesaran & Smith, 1999).
#'
#' @name pmg_estimator
#' @keywords internal
NULL


#' Recover long-run coefficients from ECM
#'
#' Given the speed of adjustment `phi` and the coefficients on levels of
#' long-run variables `omega`, recovers structural long-run coefficients
#' `w = -omega / phi` with delta-method standard errors.
#'
#' @param b_phi Numeric: coefficient on the lagged level of y (speed of adjustment).
#' @param b_omega Numeric vector: coefficients on levels of long-run x variables.
#' @param vcov_joint Numeric matrix: joint variance-covariance of (phi, omega).
#' @param phi_idx Integer: index of phi in the coefficient vector.
#' @param omega_idx Integer vector: indices of omega in the coefficient vector.
#' @return A list with `lr_coef` (long-run coefficients) and `lr_vcov`.
#' @keywords internal
.pmg_recover_lr <- function(b_phi, b_omega, vcov_joint, phi_idx, omega_idx) {
  # Long-run: w_k = -omega_k / phi
  lr_coef <- -b_omega / b_phi
  n_lr <- length(b_omega)

  # Delta method: g(phi, omega) = -omega/phi
  # dg/dphi = omega/phi^2
  # dg/domega_k = -1/phi
  b_full <- numeric(length(phi_idx) + length(omega_idx))
  full_idx <- c(phi_idx, omega_idx)

  g <- function(b) {
    phi_val <- b[1]
    omega_val <- b[2:length(b)]
    -omega_val / phi_val
  }

  b_input <- c(b_phi, b_omega)
  V_input <- vcov_joint[c(phi_idx, omega_idx), c(phi_idx, omega_idx)]

  dm <- .delta_method(g, b_input, V_input)

  list(lr_coef = lr_coef, lr_vcov = dm$vcov)
}
