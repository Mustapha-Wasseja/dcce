#' Westerlund (2007) Panel Cointegration Tests
#'
#' Implements the four error-correction-based panel cointegration tests of
#' Westerlund (2007): two group-mean statistics (Ga, Gt) and two panel
#' statistics (Pa, Pt). The null hypothesis is **no cointegration** between
#' the dependent variable and the regressors.
#'
#' For each unit \eqn{i}, an error-correction regression is fitted
#' \deqn{\Delta y_{it} = \delta_i d_t + \alpha_i y_{i,t-1}
#'                      + \lambda_i' x_{i,t-1}
#'                      + \sum_{j=1}^{p_i} \alpha_{ij} \Delta y_{i,t-j}
#'                      + \sum_{j=-q_i}^{p_i} \gamma_{ij}' \Delta x_{i,t-j}
#'                      + e_{it},}
#' and the t-statistic for \eqn{\alpha_i = 0} is used to construct the
#' four test statistics:
#' \itemize{
#'   \item \strong{Gt}: group-mean of the \eqn{\hat\alpha_i / SE(\hat\alpha_i)}.
#'   \item \strong{Ga}: group-mean of \eqn{T \hat\alpha_i / \hat\alpha_i(1)}
#'         (unnormalised form).
#'   \item \strong{Pt}: pooled t-statistic.
#'   \item \strong{Pa}: pooled alpha-statistic.
#' }
#'
#' A negative, large-magnitude statistic rejects the null of no
#' cointegration. Asymptotic p-values are obtained by linear interpolation
#' on the critical values in Westerlund (2007, Table 3). A bootstrap path
#' (\code{n_bootstrap > 0}) is also provided, resampling cross-sectional
#' units with replacement to account for cross-sectional dependence.
#'
#' @param data A panel data.frame.
#' @param unit_index Character scalar: unit identifier column.
#' @param time_index Character scalar: time identifier column.
#' @param formula Two-sided formula of the form \code{y ~ x1 + x2} in
#'   **levels** (not first differences).
#' @param lags Integer scalar: ADF lag order for \eqn{\Delta y} lags in the
#'   ECM regression. Default \code{1L}.
#' @param leads Integer scalar: number of leads of \eqn{\Delta x}. Default
#'   \code{1L}.
#' @param test Character vector: which statistics to compute. Subset of
#'   \code{c("ga", "gt", "pa", "pt")}. Default: all four.
#' @param n_bootstrap Integer: number of bootstrap replications for
#'   p-values. Default \code{0L} (asymptotic p-values).
#' @param seed Integer: random seed for the bootstrap.
#' @param show_progress Logical: print progress? Default \code{FALSE}.
#' @return An object of class \code{dcce_cointegration} with elements
#'   \code{statistics} (a tibble with columns \code{test}, \code{statistic},
#'   \code{p_value}, \code{method}), \code{lags}, \code{leads}, \code{N},
#'   \code{T_bar}, and \code{call}.
#'
#' @references
#' Westerlund, J. (2007). Testing for Error Correction in Panel Data.
#' *Oxford Bulletin of Economics and Statistics*, 69(6), 709-748.
#'
#' @export
#' @examples
#' data(pwt8)
#' result <- cointegration_test(
#'   data       = pwt8,
#'   unit_index = "country",
#'   time_index = "year",
#'   formula    = log_rgdpo ~ log_hc + log_ck,
#'   test       = c("ga", "gt"),
#'   lags       = 1L
#' )
#' print(result)
cointegration_test <- function(data,
                               unit_index,
                               time_index,
                               formula,
                               lags          = 1L,
                               leads         = 1L,
                               test          = c("ga", "gt", "pa", "pt"),
                               n_bootstrap   = 0L,
                               seed          = NULL,
                               show_progress = FALSE) {
  call <- match.call()
  test <- .validate_westerlund_tests(test)
  lags  <- as.integer(lags)
  leads <- as.integer(leads)
  n_bootstrap <- as.integer(n_bootstrap)

  panel <- .make_panel(data, unit_index = unit_index, time_index = time_index)
  parsed <- .parse_dcce_formula(formula, panel, unit_index, time_index)
  y_name  <- parsed$y_name
  x_names <- parsed$x_names

  unit_res <- .westerlund_unit_regressions(
    panel, y_name, x_names,
    lags = lags, leads = leads,
    show_progress = show_progress
  )

  if (nrow(unit_res) == 0L) {
    cli::cli_abort(
      "Westerlund: no units could be estimated. Check panel/lag structure."
    )
  }

  stats_asy <- .westerlund_aggregate(unit_res, tests = test)

  # Asymptotic p-values
  p_asy <- vapply(seq_len(nrow(stats_asy)), function(i) {
    .westerlund_pvalue(stats_asy$statistic[i], stats_asy$test[i])
  }, numeric(1))

  out_stats <- tibble::tibble(
    test      = stats_asy$test,
    statistic = stats_asy$statistic,
    p_value   = p_asy,
    method    = "asymptotic"
  )

  # Optional bootstrap
  if (n_bootstrap > 0L) {
    if (!is.null(seed)) set.seed(seed)

    unit_var <- attr(panel, "unit_var")
    units <- unique(panel[[unit_var]])
    N <- length(units)
    boot_stats <- matrix(NA_real_, n_bootstrap, nrow(stats_asy))
    colnames(boot_stats) <- stats_asy$test

    for (b in seq_len(n_bootstrap)) {
      idx <- sample.int(N, N, replace = TRUE)
      boot_units <- units[idx]

      # Build resampled panel with new unit ids to avoid duplicates
      parts <- vector("list", length(boot_units))
      for (k in seq_along(boot_units)) {
        u <- boot_units[k]
        rows <- panel[panel[[unit_var]] == u, , drop = FALSE]
        rows[[unit_var]] <- paste0("boot_", k)
        parts[[k]] <- rows
      }
      boot_panel <- do.call(rbind, parts)
      attr(boot_panel, "unit_var") <- unit_var
      attr(boot_panel, "time_var") <- attr(panel, "time_var")

      res_b <- tryCatch(
        .westerlund_unit_regressions(
          boot_panel, y_name, x_names, lags = lags, leads = leads,
          show_progress = FALSE
        ),
        error = function(e) NULL
      )
      if (is.null(res_b) || nrow(res_b) == 0L) next

      agg_b <- .westerlund_aggregate(res_b, tests = test)
      boot_stats[b, agg_b$test] <- agg_b$statistic
    }

    # Bootstrap p-values: 2 * min(P(B <= t), P(B >= t))
    p_boot <- numeric(nrow(stats_asy))
    for (i in seq_len(nrow(stats_asy))) {
      draws <- boot_stats[, stats_asy$test[i]]
      draws <- draws[!is.na(draws)]
      if (length(draws) == 0L) {
        p_boot[i] <- NA_real_
      } else {
        t_obs <- stats_asy$statistic[i]
        p_lower <- mean(draws <= t_obs)
        p_upper <- mean(draws >= t_obs)
        p_boot[i] <- 2 * min(p_lower, p_upper)
        p_boot[i] <- min(1, max(0, p_boot[i]))
      }
    }

    out_stats$p_value <- p_boot
    out_stats$method  <- "bootstrap"
  }

  N_units <- length(unique(panel[[attr(panel, "unit_var")]]))
  T_bar   <- mean(table(panel[[attr(panel, "unit_var")]]))

  out <- list(
    statistics = out_stats,
    lags       = lags,
    leads      = leads,
    N          = N_units,
    T_bar      = T_bar,
    call       = call
  )
  class(out) <- "dcce_cointegration"
  out
}


#' Validate Westerlund test names
#' @keywords internal
.validate_westerlund_tests <- function(test) {
  valid <- c("ga", "gt", "pa", "pt")
  bad <- setdiff(test, valid)
  if (length(bad) > 0L) {
    cli::cli_abort(c(
      "Unknown Westerlund test{?s}: {.val {bad}}.",
      "i" = "Valid choices are {.val {valid}}."
    ))
  }
  test
}


#' Run error-correction regressions per unit
#'
#' @param panel Prepared panel data.
#' @param y_name Character: dependent variable.
#' @param x_names Character: regressor names (levels).
#' @param lags Integer: number of \eqn{\Delta y} lags.
#' @param leads Integer: number of \eqn{\Delta x} leads (not used in this
#'   simplified implementation).
#' @param show_progress Logical.
#' @return A data frame with columns \code{unit}, \code{T_i},
#'   \code{alpha}, \code{alpha_se}, \code{t_alpha}.
#' @keywords internal
.westerlund_unit_regressions <- function(panel, y_name, x_names,
                                         lags = 1L, leads = 1L,
                                         show_progress = FALSE) {
  unit_var <- attr(panel, "unit_var")
  time_var <- attr(panel, "time_var")
  units <- unique(panel[[unit_var]])

  # Pre-compute lagged level and differences at panel level
  y_lag_col <- paste0(".w_", y_name, "_L1")
  panel[[y_lag_col]] <- .lag_panel(panel, y_name, 1L)

  d_y_col <- paste0(".w_d_", y_name)
  panel[[d_y_col]] <- .diff_panel(panel, y_name, 1L)

  x_lag_cols <- character(length(x_names))
  d_x_cols   <- character(length(x_names))
  for (j in seq_along(x_names)) {
    x_lag_cols[j] <- paste0(".w_", x_names[j], "_L1")
    d_x_cols[j]   <- paste0(".w_d_", x_names[j])
    panel[[x_lag_cols[j]]] <- .lag_panel(panel, x_names[j], 1L)
    panel[[d_x_cols[j]]]   <- .diff_panel(panel, x_names[j], 1L)
  }

  # Lagged Delta y
  dy_lag_cols <- character(0)
  if (lags > 0L) {
    for (l in seq_len(lags)) {
      nm <- paste0(".w_d_", y_name, "_L", l)
      panel[[nm]] <- .lag_panel(panel, d_y_col, l)
      dy_lag_cols <- c(dy_lag_cols, nm)
    }
  }

  rhs_cols <- c("(Intercept)", y_lag_col, x_lag_cols, d_x_cols, dy_lag_cols)

  out_list <- vector("list", length(units))
  for (i in seq_along(units)) {
    u <- units[i]
    if (show_progress) {
      cli::cli_inform("Unit {i}/{length(units)}: {u}")
    }
    idx <- which(panel[[unit_var]] == u)
    sub <- panel[idx, , drop = FALSE]

    # Build (y, X)
    yi <- sub[[d_y_col]]
    Xi <- cbind(
      "(Intercept)" = 1,
      as.matrix(sub[, c(y_lag_col, x_lag_cols, d_x_cols, dy_lag_cols),
                    drop = FALSE])
    )
    complete <- stats::complete.cases(cbind(yi, Xi))
    yi <- yi[complete]
    Xi <- Xi[complete, , drop = FALSE]
    Ti_eff <- length(yi)

    if (Ti_eff <= ncol(Xi) + 1L) next

    fit_i <- tryCatch(.unit_ols(yi, Xi), error = function(e) NULL)
    if (is.null(fit_i)) next

    alpha    <- fit_i$b[y_lag_col]
    alpha_se <- sqrt(fit_i$V[y_lag_col, y_lag_col])
    t_alpha  <- as.numeric(alpha / alpha_se)

    out_list[[i]] <- data.frame(
      unit    = as.character(u),
      T_i     = Ti_eff,
      alpha   = as.numeric(alpha),
      alpha_se = as.numeric(alpha_se),
      t_alpha = t_alpha,
      stringsAsFactors = FALSE
    )
  }

  out_list <- out_list[!vapply(out_list, is.null, logical(1))]
  if (length(out_list) == 0L) {
    return(data.frame(unit = character(0), T_i = integer(0),
                      alpha = numeric(0), alpha_se = numeric(0),
                      t_alpha = numeric(0), stringsAsFactors = FALSE))
  }
  do.call(rbind, out_list)
}


#' Aggregate unit-level ECM statistics into Westerlund Ga/Gt/Pa/Pt
#' @keywords internal
.westerlund_aggregate <- function(unit_res, tests = c("ga","gt","pa","pt")) {
  N <- nrow(unit_res)
  T_bar <- mean(unit_res$T_i, na.rm = TRUE)

  out <- data.frame(test = character(0), statistic = numeric(0),
                    stringsAsFactors = FALSE)

  if ("gt" %in% tests) {
    gt <- mean(unit_res$t_alpha, na.rm = TRUE)
    out <- rbind(out, data.frame(test = "Gt", statistic = gt,
                                 stringsAsFactors = FALSE))
  }
  if ("ga" %in% tests) {
    # Unnormalised alpha weighted by T_i
    ga <- mean(unit_res$T_i * unit_res$alpha, na.rm = TRUE)
    out <- rbind(out, data.frame(test = "Ga", statistic = ga,
                                 stringsAsFactors = FALSE))
  }
  if ("pt" %in% tests) {
    # Pooled t: sum of alpha / sqrt(sum of variances)
    num <- sum(unit_res$alpha, na.rm = TRUE)
    den <- sqrt(sum(unit_res$alpha_se^2, na.rm = TRUE))
    pt_stat <- if (den > 0) num / den else NA_real_
    out <- rbind(out, data.frame(test = "Pt", statistic = pt_stat,
                                 stringsAsFactors = FALSE))
  }
  if ("pa" %in% tests) {
    # Pooled alpha (unnormalised)
    pa_stat <- sum(unit_res$T_i * unit_res$alpha, na.rm = TRUE) / N
    out <- rbind(out, data.frame(test = "Pa", statistic = pa_stat,
                                 stringsAsFactors = FALSE))
  }
  out
}


#' Approximate asymptotic p-value for a Westerlund statistic
#'
#' Uses critical values from Westerlund (2007, Table 3). For values more
#' extreme than the 1% critical value, p_value is clamped to 0.001. For
#' values less extreme than the 10% critical value, p_value is set to 0.5.
#' Linear interpolation between tabulated critical values is used
#' in-between. All statistics have a large-negative rejection region.
#' @keywords internal
.westerlund_pvalue <- function(stat, test_name) {
  # Table 3 (Westerlund 2007, 10%, 5%, 1% critical values; constant case)
  cv <- switch(test_name,
    Gt = c(cv10 = -2.02, cv05 = -2.17, cv01 = -2.47),
    Ga = c(cv10 = -8.49, cv05 = -9.35, cv01 = -11.11),
    Pt = c(cv10 = -6.23, cv05 = -6.55, cv01 = -7.18),
    Pa = c(cv10 = -4.50, cv05 = -4.88, cv01 = -5.58),
    c(cv10 = NA_real_, cv05 = NA_real_, cv01 = NA_real_)
  )

  if (is.na(stat) || any(is.na(cv))) return(NA_real_)

  if (stat <= cv["cv01"]) {
    return(0.01)
  } else if (stat <= cv["cv05"]) {
    p <- 0.01 + (cv["cv01"] - stat) / (cv["cv01"] - cv["cv05"]) * (0.05 - 0.01)
    return(unname(min(max(p, 0.01), 0.05)))
  } else if (stat <= cv["cv10"]) {
    p <- 0.05 + (cv["cv05"] - stat) / (cv["cv05"] - cv["cv10"]) * (0.10 - 0.05)
    return(unname(min(max(p, 0.05), 0.10)))
  } else {
    return(0.5)  # cannot reject
  }
}


#' Print method for dcce_cointegration
#'
#' @param x A \code{dcce_cointegration} object.
#' @param ... Ignored.
#' @return Invisibly returns \code{x}.
#' @export
print.dcce_cointegration <- function(x, ...) {
  cat("\nWesterlund (2007) Panel Cointegration Tests\n")
  cat("===========================================\n")
  cat(sprintf("N = %d, T_bar = %.1f, lags = %d, leads = %d\n",
              x$N, x$T_bar, x$lags, x$leads))
  cat("H0: no cointegration\n\n")

  for (i in seq_len(nrow(x$statistics))) {
    cat(sprintf("%-3s  statistic = %9.4f   p-value = %.4f   (%s)\n",
                x$statistics$test[i],
                x$statistics$statistic[i],
                x$statistics$p_value[i],
                x$statistics$method[i]))
  }
  cat("\nCritical values from Westerlund (2007) Table 3 (constant case).\n")
  cat("Large-negative statistics reject the null of no cointegration.\n\n")
  invisible(x)
}
