#' Dumitrescu-Hurlin Panel Granger Causality Test
#'
#' Tests whether \eqn{x} Granger-causes \eqn{y} in a heterogeneous panel,
#' following Dumitrescu & Hurlin (2012). For each unit \eqn{i} a bivariate
#' VAR-type regression is fitted:
#' \deqn{y_{it} = \alpha_i + \sum_{k=1}^{K} \gamma_{ik} y_{i,t-k}
#'                         + \sum_{k=1}^{K} \beta_{ik} x_{i,t-k}
#'                         + u_{it},}
#' and the individual Wald statistic for
#' \eqn{H_0^{(i)}: \beta_{i1} = \cdots = \beta_{iK} = 0} is computed.
#'
#' The panel statistics are:
#' \describe{
#'   \item{W-bar}{The cross-sectional average of the N unit-level Wald
#'     statistics.}
#'   \item{Z-bar}{Standardised version:
#'     \eqn{\tilde{Z} = \sqrt{N/(2K)} (W\text{-bar} - K)
#'     \xrightarrow{d} \mathcal{N}(0,1)}.}
#'   \item{Z-bar tilde}{Small-sample adjusted version using
#'     \eqn{E[W_i]} and \eqn{Var(W_i)} from the \eqn{F(K, T-3K-1)}
#'     distribution.}
#' }
#'
#' @param data A panel data.frame.
#' @param unit_index Character: unit identifier column.
#' @param time_index Character: time identifier column.
#' @param y Character: name of the dependent variable.
#' @param x Character: name of the potential Granger-causing variable.
#' @param lags Integer: lag order \eqn{K} for both variables. Default 1.
#'
#' @return An object of class \code{dcce_granger} with elements:
#'   \describe{
#'     \item{W_bar}{Mean Wald statistic across units.}
#'     \item{Z_bar}{Standardised statistic (large-T).}
#'     \item{Z_bar_tilde}{Small-sample adjusted statistic.}
#'     \item{p_value_Z}{p-value for Z-bar.}
#'     \item{p_value_Zt}{p-value for Z-bar tilde.}
#'     \item{unit_wald}{Named vector of per-unit Wald statistics.}
#'     \item{N}{Number of units used.}
#'     \item{T_bar}{Average time-series length.}
#'     \item{lags}{Lag order.}
#'   }
#'
#' @references
#' Dumitrescu, E.-I., & Hurlin, C. (2012). Testing for Granger
#' non-causality in heterogeneous panels. *Economic Modelling*,
#' 29(4), 1450--1460.
#'
#' @export
#' @examples
#' data(pwt8)
#' gc <- granger_test(
#'   data = pwt8, unit_index = "country", time_index = "year",
#'   y = "d_log_rgdpo", x = "log_ck", lags = 1
#' )
#' print(gc)
granger_test <- function(data,
                         unit_index,
                         time_index,
                         y,
                         x,
                         lags = 1L) {
  lags <- as.integer(lags)
  if (lags < 1L) cli::cli_abort("{.arg lags} must be >= 1.")

  panel <- .make_panel(data, unit_index = unit_index, time_index = time_index)
  unit_var <- attr(panel, "unit_var")
  time_var <- attr(panel, "time_var")

  if (!y %in% names(panel)) cli::cli_abort("Variable {.val {y}} not found.")
  if (!x %in% names(panel)) cli::cli_abort("Variable {.val {x}} not found.")

  units <- unique(panel[[unit_var]])
  N <- length(units)

  wald_vec <- rep(NA_real_, N)
  names(wald_vec) <- as.character(units)
  T_vec <- integer(N)

  for (i in seq_along(units)) {
    u <- units[i]
    idx <- which(panel[[unit_var]] == u)
    yi <- panel[[y]][idx]
    xi <- panel[[x]][idx]
    Ti <- length(yi)

    if (Ti <= 3L * lags + 1L) next
    T_vec[i] <- Ti

    # Build regression matrices
    # Unrestricted: y_t ~ const + y_{t-1}..y_{t-K} + x_{t-1}..x_{t-K}
    # Restricted:   y_t ~ const + y_{t-1}..y_{t-K}
    dep   <- yi[(lags + 1L):Ti]
    n_eff <- length(dep)
    X_r   <- matrix(1, nrow = n_eff, ncol = 1L + lags)
    X_ur  <- matrix(1, nrow = n_eff, ncol = 1L + 2L * lags)

    colnames_r  <- c("const", paste0("y_L", seq_len(lags)))
    colnames_ur <- c(colnames_r, paste0("x_L", seq_len(lags)))
    colnames(X_r)  <- colnames_r
    colnames(X_ur) <- colnames_ur

    for (k in seq_len(lags)) {
      y_lag <- yi[(lags + 1L - k):(Ti - k)]
      x_lag <- xi[(lags + 1L - k):(Ti - k)]
      X_r[,  1L + k]         <- y_lag
      X_ur[, 1L + k]         <- y_lag
      X_ur[, 1L + lags + k]  <- x_lag
    }

    # Drop rows with any NA (from first-differenced series, etc.)
    complete <- stats::complete.cases(cbind(dep, X_ur))
    dep_c  <- dep[complete]
    X_ur_c <- X_ur[complete, , drop = FALSE]
    X_r_c  <- X_r[complete, , drop = FALSE]
    n_eff  <- length(dep_c)

    if (n_eff <= ncol(X_ur_c) + 1L) next

    # OLS: unrestricted and restricted
    fit_ur <- .unit_ols(dep_c, X_ur_c)
    fit_r  <- .unit_ols(dep_c, X_r_c)

    rss_ur <- sum(fit_ur$e^2)
    rss_r  <- sum(fit_r$e^2)

    # Wald = ((RSS_r - RSS_ur) / K) / (RSS_ur / (T - 3K - 1))
    df1 <- lags
    df2 <- n_eff - 2L * lags - 1L
    if (df2 <= 0L) next
    wald_vec[i] <- ((rss_r - rss_ur) / df1) / (rss_ur / df2)
  }

  keep <- !is.na(wald_vec)
  if (sum(keep) < 2L) {
    cli::cli_abort(
      "Fewer than 2 units could be estimated. Reduce {.arg lags} or check data."
    )
  }
  wald_vec <- wald_vec[keep]
  T_vec    <- T_vec[keep]
  N_used   <- length(wald_vec)
  T_bar    <- mean(T_vec)
  K        <- lags

  # W-bar: cross-sectional average of unit Wald statistics
  W_bar <- mean(wald_vec)

  # Z-bar (large-T standardisation)
  # Under H0: E[W_i] = K, Var(W_i) = 2K (chi-sq approximation)
  Z_bar <- sqrt(N_used / (2 * K)) * (W_bar - K)
  p_Z   <- 2 * stats::pnorm(-abs(Z_bar))

  # Z-bar tilde (small-sample adjusted)
  # Uses the exact moments of the F(K, T-3K-1) distribution
  # E[F] = df2/(df2-2), Var[F] = 2*df2^2*(df1+df2-2)/(df1*(df2-2)^2*(df2-4))
  df2_vec <- T_vec - 3L * K - 1L
  EW <- rep(NA_real_, N_used)
  VW <- rep(NA_real_, N_used)
  for (j in seq_len(N_used)) {
    d2 <- df2_vec[j]
    if (d2 > 4L) {
      EW[j] <- K * d2 / (d2 - 2)
      VW[j] <- 2 * K * d2^2 * (K + d2 - 2) / ((d2 - 2)^2 * (d2 - 4))
    }
  }
  ok <- !is.na(EW) & !is.na(VW) & VW > 0
  if (sum(ok) >= 2L) {
    EW_bar <- mean(EW[ok])
    VW_bar <- mean(VW[ok]) / sum(ok)   # variance of the mean
    Z_bar_tilde <- (W_bar - EW_bar) / sqrt(VW_bar)
    p_Zt <- 2 * stats::pnorm(-abs(Z_bar_tilde))
  } else {
    Z_bar_tilde <- NA_real_
    p_Zt <- NA_real_
  }

  out <- list(
    W_bar       = W_bar,
    Z_bar       = Z_bar,
    Z_bar_tilde = Z_bar_tilde,
    p_value_Z   = p_Z,
    p_value_Zt  = p_Zt,
    unit_wald   = wald_vec,
    N           = N_used,
    T_bar       = T_bar,
    lags        = K,
    y           = y,
    x           = x
  )
  class(out) <- "dcce_granger"
  out
}


#' Print a dcce_granger object
#'
#' @param x A \code{dcce_granger} object.
#' @param ... Ignored.
#' @return Invisibly returns \code{x}.
#' @export
print.dcce_granger <- function(x, ...) {
  cat("\nDumitrescu-Hurlin Panel Granger Causality Test\n")
  cat("==============================================\n")
  cat(sprintf("H0: %s does not Granger-cause %s\n", x$x, x$y))
  cat(sprintf("N = %d, T_bar = %.1f, lags = %d\n\n", x$N, x$T_bar, x$lags))
  cat(sprintf("W-bar          = %8.4f\n", x$W_bar))
  cat(sprintf("Z-bar          = %8.4f    p-value = %.4f\n",
              x$Z_bar, x$p_value_Z))
  if (!is.na(x$Z_bar_tilde)) {
    cat(sprintf("Z-bar tilde    = %8.4f    p-value = %.4f\n",
                x$Z_bar_tilde, x$p_value_Zt))
  }
  cat("\n")
  invisible(x)
}
