#' Augmented Mean Group (AMG) Estimator Internals
#'
#' Internal helpers for the Augmented Mean Group (AMG) estimator of
#' Eberhardt & Teal (2010) and Bond & Eberhardt (2013). AMG accounts for
#' cross-sectional dependence via a two-step procedure:
#'
#' \enumerate{
#'   \item Fit a pooled first-difference regression of \eqn{\Delta y_{it}}
#'         on \eqn{\Delta x_{it}} augmented with \eqn{T-1} time dummies.
#'   \item Extract the time-dummy coefficients as the Common Dynamic
#'         Process (CDP), a non-parametric proxy for unobserved
#'         common factors.
#'   \item Cumulate the CDP within each unit (back to the level) and add
#'         it as an extra regressor in a unit-level OLS on levels.
#'   \item Average the unit-level slopes (excluding the CDP) to obtain
#'         the AMG Mean Group estimate.
#' }
#'
#' This implementation uses base R throughout and does not add any new
#' package dependencies.
#'
#' @name amg_estimator
#' @keywords internal
NULL


#' Compute the Common Dynamic Process from a pooled FD regression
#'
#' Runs a pooled first-difference regression with time dummies and returns
#' a two-column data.frame mapping each time period to its estimated
#' common dynamic process value. The reference time period (the first
#' one after differencing) has CDP 0.
#'
#' @param panel A panel data.frame from \code{.make_panel()}.
#' @param y_name Character: dependent variable name.
#' @param x_names Character: regressor names (level variables).
#' @return A two-column data frame: \code{time} and \code{cdp}.
#' @keywords internal
.amg_common_dynamic_process <- function(panel, y_name, x_names) {
  unit_var <- attr(panel, "unit_var")
  time_var <- attr(panel, "time_var")

  # First differences of y and each x within units
  dy_name <- paste0(".amg_d_", y_name)
  panel[[dy_name]] <- .diff_panel(panel, y_name, 1L)

  dx_names <- character(length(x_names))
  for (j in seq_along(x_names)) {
    dx_names[j] <- paste0(".amg_d_", x_names[j])
    panel[[dx_names[j]]] <- .diff_panel(panel, x_names[j], 1L)
  }

  # Keep rows where all differenced variables are available
  keep <- stats::complete.cases(panel[, c(dy_name, dx_names, time_var)])
  pool <- panel[keep, , drop = FALSE]
  if (nrow(pool) == 0L) {
    cli::cli_abort("AMG: no observations available after first-differencing.")
  }

  # Distinct time periods present after differencing (first period dropped)
  times_post <- sort(unique(pool[[time_var]]))
  if (length(times_post) < 2L) {
    cli::cli_abort("AMG: too few time periods after differencing.")
  }

  # Build the pooled design matrix: intercept-free with time dummies
  # Reference category = first post-differencing period
  dummy_times <- times_post[-1L]
  D_mat <- matrix(0, nrow = nrow(pool), ncol = length(dummy_times))
  colnames(D_mat) <- paste0(".amg_yr_", format(dummy_times))
  for (j in seq_along(dummy_times)) {
    D_mat[, j] <- as.integer(pool[[time_var]] == dummy_times[j])
  }

  X_pool <- cbind(
    Intercept = 1,
    as.matrix(pool[, dx_names, drop = FALSE]),
    D_mat
  )
  y_pool <- pool[[dy_name]]

  fit_pool <- .unit_ols(y_pool, X_pool)
  cdp_coefs <- fit_pool$b[colnames(D_mat)]

  # CDP data frame — include reference time (CDP = 0)
  cdp_df <- data.frame(
    time = c(times_post[1L], dummy_times),
    cdp  = c(0, as.numeric(cdp_coefs)),
    stringsAsFactors = FALSE
  )
  names(cdp_df)[1L] <- time_var
  cdp_df
}


#' Attach the cumulative CDP to a panel as a new column
#'
#' Merges the CDP into the panel by time period (units before the first
#' differenced period receive NA), then cumulates within each unit to
#' recover the level-form process \eqn{\sum_{s \le t} \hat\mu_s}.
#'
#' @param panel A panel data.frame.
#' @param cdp_df Two-column data.frame from
#'   \code{.amg_common_dynamic_process()}.
#' @return The panel augmented with a \code{cdp_level} column.
#' @keywords internal
.amg_attach_cdp_level <- function(panel, cdp_df) {
  unit_var <- attr(panel, "unit_var")
  time_var <- attr(panel, "time_var")

  # Map time -> cdp via match (base R)
  cdp_lookup <- stats::setNames(cdp_df$cdp, as.character(cdp_df[[time_var]]))
  panel$cdp_inc <- as.numeric(
    cdp_lookup[as.character(panel[[time_var]])]
  )
  panel$cdp_inc[is.na(panel$cdp_inc)] <- 0

  # Cumulative sum within unit (sorted by time thanks to .make_panel)
  panel$cdp_level <- NA_real_
  units <- unique(panel[[unit_var]])
  for (u in units) {
    idx <- which(panel[[unit_var]] == u)
    panel$cdp_level[idx] <- cumsum(panel$cdp_inc[idx])
  }
  panel$cdp_inc <- NULL

  # Preserve attributes
  attr(panel, "unit_var") <- unit_var
  attr(panel, "time_var") <- time_var
  panel
}
