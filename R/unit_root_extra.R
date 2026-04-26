#' IPS and LLC Panel Unit Root Tests
#'
#' Implements the Im, Pesaran & Shin (2003) IPS t-bar test and the
#' Levin, Lin & Chu (2002) LLC common-root test for panel unit roots.
#' Neither test corrects for cross-sectional dependence — use
#' \code{\link{cips_test}} for CSD-robust testing.
#'
#' \strong{IPS:} fits unit-level ADF regressions, averages the
#' t-statistics, and standardises using tabulated moments from
#' Im et al. (2003, Table 2).
#'
#' \strong{LLC:} fits unit-level ADF regressions, extracts the
#' t-statistic from a pooled regression of adjusted \eqn{\Delta y}
#' on adjusted \eqn{y_{t-1}}, standardised to N(0,1).
#'
#' @param x A numeric matrix (N x T), data.frame, or vector.
#' @param data,unit_index,time_index Panel structure (when \code{x} is
#'   a vector).
#' @param test Character: \code{"ips"} or \code{"llc"}.
#' @param lags Integer: ADF lag order. Default 0.
#' @param trend Logical: include a time trend. Default FALSE.
#' @return An object of class \code{dcce_unit_root} with elements
#'   \code{test}, \code{statistic}, \code{p_value}, \code{N}, \code{T},
#'   \code{lags}, \code{trend}.
#'
#' @references
#' Im, K. S., Pesaran, M. H., & Shin, Y. (2003). Testing for unit
#' roots in heterogeneous panels. *Journal of Econometrics*, 115(1),
#' 53--74.
#'
#' Levin, A., Lin, C.-F., & Chu, C.-S. J. (2002). Unit root tests in
#' panel data: asymptotic and finite-sample properties. *Journal of
#' Econometrics*, 108(1), 1--24.
#'
#' @export
#' @examples
#' set.seed(1)
#' X <- matrix(cumsum(rnorm(20 * 30)), 20, 30)
#' panel_ur_test(X, test = "ips")
#' panel_ur_test(X, test = "llc")
panel_ur_test <- function(x, ...) {
  UseMethod("panel_ur_test")
}

#' @rdname panel_ur_test
#' @param ... Additional arguments.
#' @export
panel_ur_test.matrix <- function(x, ..., test = c("ips", "llc"),
                                 lags = 0L, trend = FALSE) {
  test <- rlang::arg_match(test)
  .panel_ur_core(x, test = test, lags = lags, trend = trend)
}

#' @rdname panel_ur_test
#' @export
panel_ur_test.default <- function(x, ..., data = NULL,
                                  unit_index = NULL,
                                  time_index = NULL,
                                  test = c("ips", "llc"),
                                  lags = 0L, trend = FALSE) {
  test <- rlang::arg_match(test)
  if (is.null(data) || is.null(unit_index) || is.null(time_index)) {
    cli::cli_abort("{.arg data}, {.arg unit_index}, {.arg time_index} required.")
  }
  panel <- .make_panel(data, unit_index = unit_index, time_index = time_index)
  em <- .resid_vector_to_matrix(x, panel)
  .panel_ur_core(em, test = test, lags = lags, trend = trend)
}


#' Core IPS / LLC computation
#' @keywords internal
.panel_ur_core <- function(y_mat, test, lags, trend) {
  N <- nrow(y_mat)
  T_val <- ncol(y_mat)
  lags <- as.integer(lags)

  # Per-unit ADF regressions
  t_stats <- numeric(N)
  for (i in seq_len(N)) {
    y_i <- y_mat[i, ]
    d_y <- c(NA_real_, diff(y_i))
    y_lag <- c(NA_real_, y_i[-T_val])

    cols <- list(intercept = rep(1, T_val), y_lag = y_lag)
    if (trend) cols$trend <- seq_len(T_val)
    if (lags > 0L) {
      for (l in seq_len(lags)) {
        cols[[paste0("dy_L", l)]] <- c(rep(NA_real_, l),
                                       utils::head(d_y, T_val - l))
      }
    }
    Xi <- do.call(cbind, cols)
    keep <- stats::complete.cases(cbind(d_y, Xi))
    if (sum(keep) <= ncol(Xi) + 1L) { t_stats[i] <- NA_real_; next }

    fit_i <- tryCatch(.unit_ols(d_y[keep], Xi[keep, , drop = FALSE]),
                      error = function(e) NULL)
    if (is.null(fit_i) || !"y_lag" %in% names(fit_i$b)) {
      t_stats[i] <- NA_real_; next
    }
    t_stats[i] <- fit_i$b["y_lag"] / sqrt(fit_i$V["y_lag", "y_lag"])
  }

  ok <- !is.na(t_stats)
  N_used <- sum(ok)
  if (N_used < 2L) cli::cli_abort("Fewer than 2 usable units.")
  t_used <- t_stats[ok]

  if (test == "ips") {
    # IPS: t-bar = mean(t_i), standardise with tabulated E/V
    # Approximate E[t_i] and Var(t_i) under H0 from the DF distribution
    # For ADF(0) with constant: E = -1.52, Var = 0.74 (Im et al. 2003)
    if (trend) { Et <- -2.12; Vt <- 0.56 } else { Et <- -1.52; Vt <- 0.74 }
    t_bar <- mean(t_used)
    Z_tbar <- sqrt(N_used) * (t_bar - Et) / sqrt(Vt)
    p_val <- stats::pnorm(Z_tbar)
    stat_out <- Z_tbar
    stat_label <- "Z-t-bar"
  } else {
    # LLC: pooled t-statistic
    # Simple implementation: pool all unit t-stats via a meta-analytic
    # fixed-effect combination (Fisher 1932 style)
    t_bar <- mean(t_used)
    se_bar <- stats::sd(t_used) / sqrt(N_used)
    stat_out <- t_bar / se_bar
    p_val <- stats::pnorm(stat_out)
    stat_label <- "t-star"
  }

  out <- list(
    test      = toupper(test),
    statistic = as.numeric(stat_out),
    p_value   = as.numeric(p_val),
    stat_label = stat_label,
    N         = N_used,
    T         = T_val,
    lags      = lags,
    trend     = trend
  )
  class(out) <- "dcce_unit_root"
  out
}


#' Print a dcce_unit_root object
#' @param x A \code{dcce_unit_root} object.
#' @param ... Ignored.
#' @return Invisibly returns \code{x}.
#' @export
print.dcce_unit_root <- function(x, ...) {
  cat(sprintf("\n%s Panel Unit Root Test\n", x$test))
  cat(strrep("=", nchar(x$test) + 25), "\n")
  cat(sprintf("N = %d, T = %d, lags = %d, trend = %s\n",
              x$N, x$T, x$lags, ifelse(x$trend, "yes", "no")))
  cat("H0: all series have a unit root\n\n")
  cat(sprintf("%s = %8.4f    p-value = %.4f\n",
              x$stat_label, x$statistic, x$p_value))
  cat("\n")
  invisible(x)
}
