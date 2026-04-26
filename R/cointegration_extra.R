#' Pedroni and Kao Panel Cointegration Tests
#'
#' Implements simplified versions of the Pedroni (1999, 2004) and
#' Kao (1999) residual-based panel cointegration tests.
#'
#' \strong{Pedroni:} Fits unit-level OLS regressions of \eqn{y} on
#' \eqn{x_1, \ldots, x_K}, extracts the residuals, runs an ADF
#' regression on each unit's residual series, and averages the
#' t-statistics. Reports the group-mean t-statistic (\code{group_t})
#' and the group-mean rho-statistic (\code{group_rho}).
#'
#' \strong{Kao:} Similar but pools the residuals rather than averaging
#' unit-level statistics. Reports a single ADF t-statistic on the
#' demeaned pooled residuals.
#'
#' Both tests have null hypothesis \strong{no cointegration}. A large
#' negative statistic rejects the null.
#'
#' @param data A panel data.frame.
#' @param unit_index Character: unit identifier column.
#' @param time_index Character: time identifier column.
#' @param formula Two-sided formula in levels: \code{y ~ x1 + x2}.
#' @param test Character: \code{"pedroni"} or \code{"kao"}.
#' @param lags Integer: ADF lag order for the residual regression.
#'   Default 1.
#'
#' @return An object of class \code{dcce_cointegration_extra} with
#'   elements \code{test}, \code{statistics} (a tibble), \code{N},
#'   \code{T_bar}, \code{lags}.
#'
#' @references
#' Pedroni, P. (1999). Critical values for cointegration tests in
#' heterogeneous panels with multiple regressors. *Oxford Bulletin of
#' Economics and Statistics*, 61(S1), 653--670.
#'
#' Pedroni, P. (2004). Panel cointegration: asymptotic and finite
#' sample properties of pooled time series tests with an application
#' to the PPP hypothesis. *Econometric Theory*, 20(3), 597--625.
#'
#' Kao, C. (1999). Spurious regression and residual-based tests for
#' cointegration in panel data. *Journal of Econometrics*, 90(1),
#' 1--44.
#'
#' @export
#' @examples
#' data(pwt8)
#' panel_coint_test(
#'   data = pwt8, unit_index = "country", time_index = "year",
#'   formula = log_rgdpo ~ log_hc + log_ck,
#'   test = "pedroni", lags = 1
#' )
panel_coint_test <- function(data,
                             unit_index,
                             time_index,
                             formula,
                             test = c("pedroni", "kao"),
                             lags = 1L) {
  test <- rlang::arg_match(test)
  lags <- as.integer(lags)
  call <- match.call()

  panel <- .make_panel(data, unit_index = unit_index, time_index = time_index)
  parsed <- .parse_dcce_formula(formula, panel, unit_index, time_index)
  y_name  <- parsed$y_name
  x_names <- parsed$x_names
  unit_var <- attr(panel, "unit_var")
  time_var <- attr(panel, "time_var")
  units <- unique(panel[[unit_var]])
  N <- length(units)

  # Step 1: unit-level OLS on levels, extract residuals
  resid_list <- list()
  for (u in units) {
    idx <- which(panel[[unit_var]] == u)
    yi <- panel[[y_name]][idx]
    Xi <- cbind(1, as.matrix(panel[idx, x_names, drop = FALSE]))
    colnames(Xi) <- c("const", x_names)
    keep <- stats::complete.cases(cbind(yi, Xi))
    if (sum(keep) < ncol(Xi) + 2L) next
    fit_i <- tryCatch(.unit_ols(yi[keep], Xi[keep, , drop = FALSE]),
                      error = function(e) NULL)
    if (!is.null(fit_i)) {
      resid_list[[as.character(u)]] <- fit_i$e
    }
  }

  if (length(resid_list) < 2L) {
    cli::cli_abort("Too few units with valid OLS residuals.")
  }

  if (test == "pedroni") {
    stats_out <- .pedroni_core(resid_list, lags = lags)
  } else {
    stats_out <- .kao_core(resid_list, lags = lags)
  }

  T_bar <- mean(vapply(resid_list, length, integer(1)))

  out <- list(
    test       = toupper(test),
    statistics = stats_out,
    N          = length(resid_list),
    T_bar      = T_bar,
    lags       = lags,
    call       = call
  )
  class(out) <- "dcce_cointegration_extra"
  out
}


#' Pedroni group-mean statistics
#' @keywords internal
.pedroni_core <- function(resid_list, lags) {
  N <- length(resid_list)
  t_stats <- numeric(N)
  rho_stats <- numeric(N)

  for (i in seq_len(N)) {
    e <- resid_list[[i]]
    Ti <- length(e)
    if (Ti < lags + 4L) { t_stats[i] <- NA_real_; rho_stats[i] <- NA_real_; next }

    d_e <- diff(e)
    e_lag <- e[-Ti]

    cols <- list(e_lag = e_lag[-1])
    if (lags > 0L) {
      for (l in seq_len(lags)) {
        if (l + 1L <= length(d_e)) {
          cols[[paste0("de_L", l)]] <- c(rep(NA_real_, l),
                                         utils::head(d_e, length(d_e) - l))[-1]
        }
      }
    }
    Xi <- do.call(cbind, cols)
    dep <- d_e[-1]
    keep <- seq_len(min(length(dep), nrow(Xi)))
    keep <- keep[stats::complete.cases(cbind(dep[keep], Xi[keep, , drop = FALSE]))]
    if (length(keep) < ncol(Xi) + 1L) {
      t_stats[i] <- NA_real_; rho_stats[i] <- NA_real_; next
    }

    fit_i <- tryCatch(.unit_ols(dep[keep], Xi[keep, , drop = FALSE]),
                      error = function(e) NULL)
    if (is.null(fit_i) || !"e_lag" %in% names(fit_i$b)) {
      t_stats[i] <- NA_real_; rho_stats[i] <- NA_real_; next
    }
    t_stats[i]   <- fit_i$b["e_lag"] / sqrt(fit_i$V["e_lag", "e_lag"])
    rho_stats[i] <- (Ti - 1L) * fit_i$b["e_lag"]
  }

  ok <- !is.na(t_stats)
  N_ok <- sum(ok)
  group_t   <- mean(t_stats[ok])
  group_rho <- mean(rho_stats[ok])

  # Standardise: under H0, E[t_i] approx -1.65 (Pedroni 2004)
  z_t   <- sqrt(N_ok) * (group_t - (-1.65)) / 1
  z_rho <- sqrt(N_ok) * group_rho   # already centred at 0 under H0

  tibble::tibble(
    test      = c("Group t", "Group rho"),
    statistic = c(z_t, z_rho),
    p_value   = stats::pnorm(c(z_t, z_rho))
  )
}


#' Kao ADF statistic
#' @keywords internal
.kao_core <- function(resid_list, lags) {
  # Pool all residuals, fit a single ADF
  e_all <- unlist(resid_list)
  T_per <- vapply(resid_list, length, integer(1))
  N <- length(resid_list)

  # Build lagged and differenced series respecting unit boundaries
  dep_all <- numeric(0)
  X_all   <- numeric(0)

  for (i in seq_len(N)) {
    e <- resid_list[[i]]
    Ti <- length(e)
    if (Ti < lags + 3L) next
    d_e <- diff(e)
    e_lag <- e[-Ti]
    dep <- d_e[-1]
    x_lag <- e_lag[-1]

    keep <- seq_len(min(length(dep), length(x_lag)))
    dep_all <- c(dep_all, dep[keep])
    X_all   <- c(X_all, x_lag[keep])
  }

  if (length(dep_all) < 3L) {
    return(tibble::tibble(test = "ADF", statistic = NA_real_,
                          p_value = NA_real_))
  }

  Xi <- cbind(e_lag = X_all)
  fit <- tryCatch(.unit_ols(dep_all, Xi), error = function(e) NULL)

  if (is.null(fit)) {
    return(tibble::tibble(test = "ADF", statistic = NA_real_,
                          p_value = NA_real_))
  }

  t_stat <- fit$b["e_lag"] / sqrt(fit$V["e_lag", "e_lag"])

  tibble::tibble(
    test      = "ADF",
    statistic = as.numeric(t_stat),
    p_value   = stats::pnorm(as.numeric(t_stat))
  )
}


#' Print a dcce_cointegration_extra object
#' @param x A \code{dcce_cointegration_extra} object.
#' @param ... Ignored.
#' @return Invisibly returns \code{x}.
#' @export
print.dcce_cointegration_extra <- function(x, ...) {
  cat(sprintf("\n%s Panel Cointegration Test\n", x$test))
  cat(strrep("=", nchar(x$test) + 30), "\n")
  cat(sprintf("N = %d, T_bar = %.1f, lags = %d\n", x$N, x$T_bar, x$lags))
  cat("H0: no cointegration\n\n")
  for (i in seq_len(nrow(x$statistics))) {
    cat(sprintf("  %-12s  statistic = %9.4f   p-value = %.4f\n",
                x$statistics$test[i],
                x$statistics$statistic[i],
                x$statistics$p_value[i]))
  }
  cat("\nLarge negative statistics reject the null.\n\n")
  invisible(x)
}
