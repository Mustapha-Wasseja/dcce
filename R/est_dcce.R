#' Dynamic CCE Estimator Internals
#'
#' Internal functions specific to the Dynamic CCE (DCCE) estimator of
#' Chudik and Pesaran (2015), including jackknife bias correction and
#' collinearity checks on the augmented CSA matrix.
#'
#' @name dcce_estimator
#' @keywords internal
NULL


# ──────────────────────────────────────────────────────────────────────────────
# .check_csa_rank — collinearity check for augmented CSA matrix
# ──────────────────────────────────────────────────────────────────────────────

#' Check rank of augmented CSA matrix
#'
#' Verifies that the CSA matrix (with lags) is full column rank. Warns but
#' continues if rank deficient (MG estimates remain consistent per
#' Chudik & Pesaran 2015, p.398).
#'
#' @param Z Numeric matrix: the CSA matrix for a given unit.
#' @param unit_id Character: unit identifier for warning messages.
#' @return Logical: `TRUE` if full rank, `FALSE` otherwise.
#' @keywords internal
.check_csa_rank <- function(Z, unit_id = "") {
  if (ncol(Z) == 0L) return(TRUE)
  rk <- Matrix::rankMatrix(Z)
  if (rk < ncol(Z)) {
    cli::cli_warn(
      "CSA matrix is rank deficient for unit {.val {unit_id}} (rank {rk} < {ncol(Z)} columns). MG estimates remain consistent."
    )
    return(FALSE)
  }
  TRUE
}


# ──────────────────────────────────────────────────────────────────────────────
# .jackknife_bias_correction — half-panel jackknife
# ──────────────────────────────────────────────────────────────────────────────

#' Half-panel jackknife bias correction
#'
#' Implements the Chudik & Pesaran (2015) jackknife bias correction:
#' `b_jk = 2 * b_full - 0.5 * (b_first + b_second)`
#' where `b_first` and `b_second` are MG estimates from the first and
#' second half of the time dimension, respectively.
#'
#' @param panel A panel data.frame with all CSA columns already appended.
#' @param y_name Character: dependent variable name.
#' @param x_names Character: regressor names.
#' @param csa_colnames Character: CSA column names.
#' @param unit_var Character: unit variable name.
#' @param time_var Character: time variable name.
#' @param constant Logical: include intercept.
#' @param trend Logical: include trend.
#' @param b_full Numeric vector: full-sample MG coefficient.
#' @return A list with `b_jk` (bias-corrected coefficients).
#' @keywords internal
.jackknife_bias_correction <- function(panel, y_name, x_names, csa_colnames,
                                       unit_var, time_var,
                                       constant, trend, b_full) {
  times <- sort(unique(panel[[time_var]]))
  T_total <- length(times)
  T_half <- floor(T_total / 2)

  times_first  <- times[1:T_half]
  times_second <- times[(T_half + 1):T_total]

  # Estimate on first half
  b_first <- .estimate_half(panel, y_name, x_names, csa_colnames,
                            unit_var, time_var, constant, trend, times_first)

  # Estimate on second half
  b_second <- .estimate_half(panel, y_name, x_names, csa_colnames,
                             unit_var, time_var, constant, trend, times_second)

  # Jackknife correction
  b_jk <- 2 * b_full - 0.5 * (b_first + b_second)

  list(b_jk = b_jk, b_first = b_first, b_second = b_second)
}


#' Estimate MG on a subset of time periods
#'
#' @keywords internal
.estimate_half <- function(panel, y_name, x_names, csa_colnames,
                           unit_var, time_var, constant, trend, time_subset) {

  panel_sub <- panel[panel[[time_var]] %in% time_subset, , drop = FALSE]
  units <- unique(panel_sub[[unit_var]])

  coef_names_structural <- x_names
  if (constant) coef_names_structural <- c("(Intercept)", coef_names_structural)
  if (trend) coef_names_structural <- c(coef_names_structural, "(trend)")

  coef_list <- list()
  for (u in units) {
    idx <- which(panel_sub[[unit_var]] == u)
    yi <- panel_sub[[y_name]][idx]
    Xi <- .build_unit_X(panel_sub, idx, x_names, csa_colnames, constant, trend)

    complete <- stats::complete.cases(cbind(yi, Xi))
    yi <- yi[complete]
    Xi <- Xi[complete, , drop = FALSE]

    if (length(yi) <= ncol(Xi)) next

    fit_i <- .unit_ols(yi, Xi)
    coef_list[[as.character(u)]] <- fit_i$b[coef_names_structural]
  }

  if (length(coef_list) == 0L) return(rep(NA_real_, length(coef_names_structural)))
  .mg_aggregate(coef_list)
}
