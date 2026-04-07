#' CS-DL Estimator Internals
#'
#' Internal functions for the Cross-Sectionally augmented Distributed Lag
#' (CS-DL) estimator of Chudik, Mohaddes, Pesaran & Raissi (2016). The CS-DL
#' regression is
#' \deqn{\Delta y_{it} = \alpha_i + w'_i x_{it}
#'                      + \sum_{\ell=0}^{p_x} \phi'_{i\ell} \Delta x_{i,t-\ell}
#'                      + \delta'_i \bar{z}_t + u_{it},}
#' where the long-run coefficient \eqn{w_i} is identified directly as the
#' coefficient on the level of \eqn{x_{it}}. This file provides helpers to
#' augment a formula with the required \eqn{\Delta x_{t-\ell}} terms before
#' the main \code{dcce()} pipeline runs.
#'
#' @name csdl_estimator
#' @keywords internal
NULL


#' Augment a CS-DL panel with first-difference lags of the regressors
#'
#' Given a panel and the list of base regressor names, attach columns for
#' the first difference of each regressor and its lags 1, ..., \code{p_x},
#' and also the first difference of the dependent variable (used as the
#' CS-DL LHS). The returned list contains the augmented panel plus the
#' vectors of column names that should be used on the LHS and RHS of the
#' CS-DL regression.
#'
#' @param panel A panel data.frame with \code{unit_var}/\code{time_var}
#'   attributes.
#' @param y_name Character scalar: dependent variable name.
#' @param x_names Character vector: base regressor names (levels).
#' @param p_x Integer: number of \eqn{\Delta x} lags to include (default 3).
#' @return A list with
#'   \describe{
#'     \item{panel}{Augmented panel with new columns.}
#'     \item{y_diff_name}{Name of the first-difference dependent variable.}
#'     \item{rhs_terms}{Character vector: level x + lagged \eqn{\Delta x}
#'       terms for the CS-DL regression.}
#'     \item{lr_names}{Character vector: column names whose coefficients are
#'       the long-run effects (one per base x).}
#'   }
#' @keywords internal
.csdl_augment <- function(panel, y_name, x_names, p_x = 3L) {
  unit_var <- attr(panel, "unit_var")
  time_var <- attr(panel, "time_var")

  # First difference of y as the LHS
  y_diff_col <- paste0("d.", y_name)
  if (!y_diff_col %in% names(panel)) {
    panel[[y_diff_col]] <- .diff_panel(panel, y_name, 1L)
  }

  rhs_terms <- character(0)
  lr_names <- character(0)

  for (xv in x_names) {
    # Level of x is a long-run regressor
    rhs_terms <- c(rhs_terms, xv)
    lr_names  <- c(lr_names, xv)

    # Contemporaneous Delta x
    d_col <- paste0("d.", xv)
    if (!d_col %in% names(panel)) {
      panel[[d_col]] <- .diff_panel(panel, xv, 1L)
    }
    rhs_terms <- c(rhs_terms, d_col)

    # Lagged Delta x (l = 1, ..., p_x)
    if (p_x > 0L) {
      for (l in seq_len(p_x)) {
        lag_col <- paste0("d.", xv, "_L", l)
        if (!lag_col %in% names(panel)) {
          panel[[lag_col]] <- .lag_panel(panel, d_col, l)
        }
        rhs_terms <- c(rhs_terms, lag_col)
      }
    }
  }

  # Preserve panel attributes
  attr(panel, "unit_var") <- unit_var
  attr(panel, "time_var") <- time_var

  list(
    panel       = panel,
    y_diff_name = y_diff_col,
    rhs_terms   = rhs_terms,
    lr_names    = lr_names
  )
}


#' Post-process a CS-DL fit to label long-run coefficients
#'
#' The coefficients on the level x terms are the long-run effects. This
#' helper extracts them and stores them in the fit's \code{lr_coef} and
#' \code{lr_vcov} fields so that downstream S3 methods can print a long-run
#' block.
#'
#' @param fit A \code{dcce_fit} object.
#' @param lr_names Character vector of x level names whose coefficients are
#'   the long-run effects.
#' @return Augmented fit object.
#' @keywords internal
.csdl_postprocess <- function(fit, lr_names) {
  keep <- lr_names[lr_names %in% names(fit$coefficients)]
  if (length(keep) == 0L) return(fit)

  fit$lr_coef <- fit$coefficients[keep]
  fit$lr_vcov <- fit$vcov[keep, keep, drop = FALSE]
  fit$lr_se   <- sqrt(diag(fit$lr_vcov))
  fit$lr_names_used <- keep
  fit
}
