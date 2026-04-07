#' Cross-Sectional Average (CSA) Utilities
#'
#' Internal functions for computing cross-sectional averages and their lags,
#' used to approximate unobserved common factors in CCE-type estimators.
#'
#' @name csa_utils
#' @keywords internal
NULL

# ──────────────────────────────────────────────────────────────────────────────
# .build_csa — cross-sectional averages with lags
# ──────────────────────────────────────────────────────────────────────────────

#' Build cross-sectional averages
#'
#' Computes cross-sectional means of `vars` at each time period and appends
#' contemporaneous plus lagged CSAs as new columns. For unbalanced panels,
#' means are computed using only observed units at each period.
#'
#' @param panel A panel data.frame from `.make_panel()`.
#' @param vars Character vector of variable names to average.
#' @param lags Integer: number of lags of CSAs to include (0 = contemporaneous
#'   only). Can also be a named integer vector for per-variable lags.
#' @return The input `panel` with additional CSA columns appended.
#' @keywords internal
.build_csa <- function(panel, vars, lags = 0L) {

  unit_var <- attr(panel, "unit_var")
  time_var <- attr(panel, "time_var")
  times <- sort(unique(panel[[time_var]]))

  # Handle per-variable lags
  if (length(lags) == 1L && is.null(names(lags))) {
    lag_vec <- stats::setNames(rep(as.integer(lags), length(vars)), vars)
  } else {
    lag_vec <- as.integer(lags[vars])
    names(lag_vec) <- vars
  }

  # Compute cross-sectional means at each time period
  csa_means <- list()
  for (v in vars) {
    means <- tapply(panel[[v]], panel[[time_var]], mean, na.rm = TRUE)
    csa_means[[v]] <- means
  }

  # Merge contemporaneous CSAs back to panel
  for (v in vars) {
    col_name <- paste0("csa_", v)
    panel[[col_name]] <- as.numeric(csa_means[[v]][as.character(panel[[time_var]])])
  }

  # Build lagged CSAs
  max_lag <- max(lag_vec)
  if (max_lag > 0L) {
    for (v in vars) {
      p_v <- lag_vec[v]
      if (p_v == 0L) next
      means_vec <- csa_means[[v]]
      for (l in seq_len(p_v)) {
        col_name <- paste0("csa_", v, "_L", l)
        lagged_means <- rep(NA_real_, length(times))
        names(lagged_means) <- names(means_vec)
        if (l < length(times)) {
          lagged_means[(l + 1L):length(times)] <- means_vec[1L:(length(times) - l)]
        }
        panel[[col_name]] <- as.numeric(lagged_means[as.character(panel[[time_var]])])
      }
    }
  }

  # Preserve panel attributes
  attr(panel, "unit_var") <- unit_var
  attr(panel, "time_var") <- time_var
  panel
}


# ──────────────────────────────────────────────────────────────────────────────
# .build_global_csa — CSAs from a wider sample
# ──────────────────────────────────────────────────────────────────────────────

#' Build cross-sectional averages from a wider sample
#'
#' Computes cross-sectional means using `full_panel` (a broader sample) and
#' merges them into `panel` (the estimation sample). Useful when the estimation
#' sample is a subset of the available data.
#'
#' @param panel A panel data.frame (estimation sample).
#' @param full_panel A panel data.frame (wider sample for computing averages).
#' @param vars Character vector of variable names.
#' @param lags Integer or named integer vector of CSA lags.
#' @return `panel` with CSA columns appended (computed from `full_panel`).
#' @keywords internal
.build_global_csa <- function(panel, full_panel, vars, lags = 0L) {
  unit_var <- attr(panel, "unit_var")
  time_var <- attr(panel, "time_var")

  # Compute means from full_panel
  times <- sort(unique(full_panel[[time_var]]))

  if (length(lags) == 1L && is.null(names(lags))) {
    lag_vec <- stats::setNames(rep(as.integer(lags), length(vars)), vars)
  } else {
    lag_vec <- as.integer(lags[vars])
    names(lag_vec) <- vars
  }

  csa_means <- list()
  for (v in vars) {
    means <- tapply(full_panel[[v]], full_panel[[time_var]], mean, na.rm = TRUE)
    csa_means[[v]] <- means
  }

  # Merge into estimation panel
  for (v in vars) {
    col_name <- paste0("csa_", v)
    panel[[col_name]] <- as.numeric(csa_means[[v]][as.character(panel[[time_var]])])
  }

  max_lag <- max(lag_vec)
  if (max_lag > 0L) {
    for (v in vars) {
      p_v <- lag_vec[v]
      if (p_v == 0L) next
      means_vec <- csa_means[[v]]
      for (l in seq_len(p_v)) {
        col_name <- paste0("csa_", v, "_L", l)
        lagged_means <- rep(NA_real_, length(times))
        names(lagged_means) <- names(means_vec)
        if (l < length(times)) {
          lagged_means[(l + 1L):length(times)] <- means_vec[1L:(length(times) - l)]
        }
        panel[[col_name]] <- as.numeric(lagged_means[as.character(panel[[time_var]])])
      }
    }
  }

  attr(panel, "unit_var") <- unit_var
  attr(panel, "time_var") <- time_var
  panel
}


# ──────────────────────────────────────────────────────────────────────────────
# .build_cluster_csa — CSAs within clusters
# ──────────────────────────────────────────────────────────────────────────────

#' Build cluster-level cross-sectional averages
#'
#' Computes cross-sectional means within groups defined by column `by`.
#' Useful for models with group-specific common factors.
#'
#' @param panel A panel data.frame from `.make_panel()`.
#' @param vars Character vector of variable names.
#' @param lags Integer or named integer vector of CSA lags.
#' @param by Character: name of the grouping variable in `panel`.
#' @return `panel` with CSA columns appended (group-specific averages).
#' @keywords internal
.build_cluster_csa <- function(panel, vars, lags = 0L, by) {
  unit_var <- attr(panel, "unit_var")
  time_var <- attr(panel, "time_var")
  groups <- unique(panel[[by]])

  if (length(lags) == 1L && is.null(names(lags))) {
    lag_vec <- stats::setNames(rep(as.integer(lags), length(vars)), vars)
  } else {
    lag_vec <- as.integer(lags[vars])
    names(lag_vec) <- vars
  }

  # Initialize CSA columns
  for (v in vars) {
    panel[[paste0("csa_", v)]] <- NA_real_
    p_v <- lag_vec[v]
    if (p_v > 0L) {
      for (l in seq_len(p_v)) {
        panel[[paste0("csa_", v, "_L", l)]] <- NA_real_
      }
    }
  }

  # Compute per-group CSAs
  for (g in groups) {
    g_idx <- which(panel[[by]] == g)
    g_panel <- panel[g_idx, , drop = FALSE]
    g_times <- sort(unique(g_panel[[time_var]]))

    for (v in vars) {
      means <- tapply(g_panel[[v]], g_panel[[time_var]], mean, na.rm = TRUE)
      panel[[paste0("csa_", v)]][g_idx] <- as.numeric(means[as.character(g_panel[[time_var]])])

      p_v <- lag_vec[v]
      if (p_v > 0L) {
        for (l in seq_len(p_v)) {
          lagged_means <- rep(NA_real_, length(g_times))
          names(lagged_means) <- names(means)
          if (l < length(g_times)) {
            lagged_means[(l + 1L):length(g_times)] <- means[1L:(length(g_times) - l)]
          }
          panel[[paste0("csa_", v, "_L", l)]][g_idx] <-
            as.numeric(lagged_means[as.character(g_panel[[time_var]])])
        }
      }
    }
  }

  attr(panel, "unit_var") <- unit_var
  attr(panel, "time_var") <- time_var
  panel
}


# ──────────────────────────────────────────────────────────────────────────────
# .get_csa_colnames — extract CSA column names from panel
# ──────────────────────────────────────────────────────────────────────────────

#' Get CSA column names
#'
#' Returns all column names matching the CSA naming convention.
#'
#' @param panel A panel data.frame with CSA columns.
#' @return Character vector of CSA column names.
#' @keywords internal
.get_csa_colnames <- function(panel) {
  grep("^csa_", names(panel), value = TRUE)
}
