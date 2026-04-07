#' Pesaran CIPS Panel Unit Root Test
#'
#' Implements the cross-sectionally augmented IPS (CIPS) panel unit root
#' test of Pesaran (2007), which allows for cross-sectional dependence
#' through a single unobserved common factor. For each unit, a
#' cross-sectionally augmented Dickey-Fuller (CADF) regression is run:
#'
#' \deqn{\Delta y_{it} = a_i + b_i y_{i,t-1} + c_i \bar{y}_{t-1}
#'                     + d_i \Delta \bar{y}_t
#'                     + \sum_{j=1}^{p} \rho_{ij} \Delta y_{i,t-j}
#'                     + \sum_{j=0}^{p} \delta_{ij} \Delta \bar{y}_{t-j}
#'                     + u_{it}.}
#'
#' The CIPS statistic is the cross-sectional average of the unit-level
#' t-statistics for \eqn{b_i = 0} (the CADF statistic). Critical values
#' come from Pesaran (2007, Table II(b), constant case) and Pesaran
#' (2007, Table II(c), constant + trend case). The null hypothesis is
#' that all series contain a unit root.
#'
#' @param x A numeric vector, numeric matrix (N x T), data.frame, or
#'   `dcce_fit` object. If a vector, `data`, `unit_index`, and `time_index`
#'   must also be supplied.
#' @param data A data.frame containing the panel structure (when `x` is
#'   a vector).
#' @param unit_index Character: name of the unit variable in `data`.
#' @param time_index Character: name of the time variable in `data`.
#' @param lags Integer: number of lags of \eqn{\Delta y} to include
#'   in the CADF regression. Default 0 (pure CADF without augmentation).
#' @param trend Logical: include a linear time trend? Default `FALSE`.
#' @return An object of class `dcce_cips` with elements:
#'   \describe{
#'     \item{statistic}{The CIPS statistic (truncated cross-sectional average).}
#'     \item{p_value}{Approximate p-value from Pesaran (2007) critical values.}
#'     \item{unit_stats}{Per-unit truncated CADF t-statistics.}
#'     \item{N}{Number of units.}
#'     \item{T}{Time dimension.}
#'     \item{lags}{Number of augmentation lags.}
#'     \item{trend}{Whether a trend was included.}
#'   }
#' @references
#' Pesaran, M. H. (2007). A simple panel unit root test in the presence
#' of cross-section dependence. *Journal of Applied Econometrics*,
#' 22(2), 265-312.
#'
#' @export
#' @examples
#' set.seed(1)
#' N <- 20; T <- 30
#' f <- cumsum(rnorm(T))
#' X <- matrix(NA, N, T)
#' for (i in seq_len(N)) X[i, ] <- cumsum(rnorm(T)) + 0.5 * f
#' cips_test(X)
cips_test <- function(x, ...) {
  UseMethod("cips_test")
}

#' @rdname cips_test
#' @param ... Additional arguments passed to methods.
#' @export
cips_test.matrix <- function(x, ..., lags = 0L, trend = FALSE) {
  .cips_core(x, lags = lags, trend = trend)
}

#' @rdname cips_test
#' @export
cips_test.dcce_fit <- function(x, ..., lags = 0L, trend = FALSE) {
  em <- .resid_list_to_matrix(x$resid_list)
  .cips_core(em, lags = lags, trend = trend)
}

#' @rdname cips_test
#' @export
cips_test.default <- function(x, ..., data = NULL, unit_index = NULL,
                              time_index = NULL, lags = 0L, trend = FALSE) {
  if (!is.numeric(x)) {
    cli::cli_abort("{.arg x} must be numeric, a matrix, or a {.cls dcce_fit}.")
  }
  if (is.null(data) || is.null(unit_index) || is.null(time_index)) {
    cli::cli_abort(
      "{.arg data}, {.arg unit_index}, and {.arg time_index} are required."
    )
  }
  panel <- .make_panel(data, unit_index = unit_index, time_index = time_index)
  em <- .resid_vector_to_matrix(x, panel)
  .cips_core(em, lags = lags, trend = trend)
}


#' Core CIPS computation
#' @keywords internal
.cips_core <- function(y_mat, lags = 0L, trend = FALSE) {
  N <- nrow(y_mat)
  T_val <- ncol(y_mat)
  lags <- as.integer(lags)

  if (T_val < 3L + lags) {
    cli::cli_abort("Insufficient time periods for CIPS with lags = {lags}.")
  }

  # Cross-sectional means across units at each time period
  ybar <- colMeans(y_mat, na.rm = TRUE)
  d_ybar <- c(NA_real_, diff(ybar))

  # Per-unit CADF regressions
  t_stats <- numeric(N)
  for (i in seq_len(N)) {
    y_i <- y_mat[i, ]
    d_y <- c(NA_real_, diff(y_i))
    y_lag <- c(NA_real_, y_i[-T_val])
    ybar_lag <- c(NA_real_, ybar[-T_val])

    X_cols <- list(
      y_lag    = y_lag,
      ybar_lag = ybar_lag,
      d_ybar   = d_ybar
    )

    if (trend) {
      X_cols$trend <- seq_len(T_val)
    }

    # Add lagged Delta y and lagged Delta ybar
    if (lags > 0L) {
      for (j in seq_len(lags)) {
        lag_dy <- c(rep(NA_real_, j), utils::head(d_y, T_val - j))
        X_cols[[paste0("dy_L", j)]] <- lag_dy
        lag_dybar <- c(rep(NA_real_, j), utils::head(d_ybar, T_val - j))
        X_cols[[paste0("dybar_L", j)]] <- lag_dybar
      }
    }

    Xi <- do.call(cbind, X_cols)
    Xi <- cbind(intercept = 1, Xi)

    keep <- stats::complete.cases(cbind(d_y, Xi))
    if (sum(keep) <= ncol(Xi) + 1L) {
      t_stats[i] <- NA_real_
      next
    }

    fit_i <- tryCatch(
      .unit_ols(d_y[keep], Xi[keep, , drop = FALSE]),
      error = function(e) NULL
    )
    if (is.null(fit_i)) {
      t_stats[i] <- NA_real_
      next
    }

    # t-stat on y_lag (the CADF unit root coefficient)
    if ("y_lag" %in% names(fit_i$b)) {
      b_y <- fit_i$b["y_lag"]
      se_y <- sqrt(fit_i$V["y_lag", "y_lag"])
      t_stats[i] <- as.numeric(b_y / se_y)
    } else {
      t_stats[i] <- NA_real_
    }
  }

  # Truncate extreme values (Pesaran 2007 p.275)
  truncate_val <- function(tt, case_const = TRUE, case_trend = FALSE) {
    if (is.na(tt)) return(NA_real_)
    if (case_trend) {
      K1 <- -6.19; K2 <- 2.61
    } else {
      K1 <- -6.12; K2 <- 4.16
    }
    max(K1, min(K2, tt))
  }
  t_stats_trunc <- vapply(t_stats, truncate_val,
                          case_const = !trend, case_trend = trend,
                          FUN.VALUE = numeric(1))

  cips_stat <- mean(t_stats_trunc, na.rm = TRUE)

  # Approximate p-value from Pesaran (2007) critical values
  # Table II(b) [constant] and II(c) [constant + trend]
  crit_const <- c("1%" = -2.50, "5%" = -2.27, "10%" = -2.15)
  crit_trend <- c("1%" = -3.06, "5%" = -2.83, "10%" = -2.72)
  cv <- if (trend) crit_trend else crit_const

  # Linear interpolation on the inverse CDF for rough p-value
  if (cips_stat <= cv["1%"]) {
    p_val <- 0.01
  } else if (cips_stat <= cv["5%"]) {
    p_val <- 0.01 + (cv["1%"] - cips_stat) / (cv["1%"] - cv["5%"]) * (0.05 - 0.01)
    p_val <- min(max(p_val, 0.01), 0.05)
  } else if (cips_stat <= cv["10%"]) {
    p_val <- 0.05 + (cv["5%"] - cips_stat) / (cv["5%"] - cv["10%"]) * (0.10 - 0.05)
    p_val <- min(max(p_val, 0.05), 0.10)
  } else {
    p_val <- 0.5  # non-rejection; p > 0.10
  }

  out <- list(
    statistic  = as.numeric(cips_stat),
    p_value    = as.numeric(p_val),
    unit_stats = t_stats_trunc,
    critical   = cv,
    N          = N,
    T          = T_val,
    lags       = lags,
    trend      = trend
  )
  class(out) <- "dcce_cips"
  out
}


#' Print a dcce_cips object
#'
#' @param x A `dcce_cips` object.
#' @param ... Ignored.
#' @return Invisibly returns `x`.
#' @export
print.dcce_cips <- function(x, ...) {
  cat("\nPesaran (2007) CIPS Panel Unit Root Test\n")
  cat("========================================\n")
  cat(sprintf("N = %d, T = %d, lags = %d, trend = %s\n",
              x$N, x$T, x$lags, ifelse(x$trend, "yes", "no")))
  cat(sprintf("\nH0: all series have a unit root\n"))
  cat(sprintf("CIPS statistic = %8.4f    approx p-value = %.4f\n",
              x$statistic, x$p_value))
  cat("\nCritical values (Pesaran 2007):\n")
  cat(sprintf("  1%%: %6.2f    5%%: %6.2f    10%%: %6.2f\n",
              x$critical["1%"], x$critical["5%"], x$critical["10%"]))
  cat("\n")
  invisible(x)
}
