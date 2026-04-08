#' High-Dimensional Fixed Effect Absorption
#'
#' Internal helpers to project out one or more grouping factors from a
#' set of numeric columns before the main \code{dcce()} unit loop runs.
#' When a single factor is supplied the within transformation reduces to
#' a simple group-mean demeaning; for two or more factors the alternating
#' projections algorithm of Guimaraes & Portugal (2010) / Correia (2016)
#' is used.
#'
#' The implementation is pure base R and requires no additional package
#' dependencies.
#'
#' @name absorb_internals
#' @keywords internal
NULL


#' Resolve absorb argument to a list of grouping vectors
#'
#' @param absorb Either \code{NULL}, a character vector of column names,
#'   or a one-sided formula like \code{~ industry + year}.
#' @param panel A panel data.frame.
#' @return A named list of grouping factors (one per variable) or
#'   \code{NULL} if no absorb requested.
#' @keywords internal
.absorb_resolve <- function(absorb, panel) {
  if (is.null(absorb)) return(NULL)

  if (inherits(absorb, "formula")) {
    if (length(absorb) != 2L) {
      cli::cli_abort("{.arg absorb} must be a one-sided formula.")
    }
    vars <- all.vars(absorb)
  } else if (is.character(absorb)) {
    vars <- absorb
  } else {
    cli::cli_abort(
      "{.arg absorb} must be NULL, a character vector, or a one-sided formula."
    )
  }

  bad <- setdiff(vars, names(panel))
  if (length(bad) > 0L) {
    cli::cli_abort(
      "Absorb variable{?s} {.val {bad}} not found in the panel."
    )
  }

  groups <- lapply(vars, function(v) as.factor(panel[[v]]))
  names(groups) <- vars
  groups
}


#' Within transform: demean numeric columns by a single factor
#'
#' @param X Numeric vector or matrix.
#' @param g A factor of length matching the rows of X.
#' @return Demeaned X.
#' @keywords internal
.absorb_demean_one <- function(X, g) {
  if (is.null(dim(X))) {
    means <- stats::ave(X, g, FUN = function(v) mean(v, na.rm = TRUE))
    return(X - means)
  }
  out <- X
  for (j in seq_len(ncol(X))) {
    col_j <- X[, j]
    means <- stats::ave(col_j, g, FUN = function(v) mean(v, na.rm = TRUE))
    out[, j] <- col_j - means
  }
  out
}


#' Alternating projections for multi-factor demeaning
#'
#' Iteratively projects X off each grouping factor until convergence.
#' Correia (2016) shows this converges geometrically for any finite
#' number of factors; in practice 10-30 iterations suffice for the
#' tolerances used below.
#'
#' @param X Numeric vector or matrix.
#' @param groups A list of factors (from \code{.absorb_resolve()}).
#' @param tol Numeric: convergence tolerance on max absolute change.
#' @param max_iter Integer: hard iteration cap.
#' @return Within-transformed X.
#' @keywords internal
.absorb_demean_multi <- function(X, groups, tol = 1e-10, max_iter = 200L) {
  if (length(groups) == 0L) return(X)
  if (length(groups) == 1L) return(.absorb_demean_one(X, groups[[1L]]))

  out <- X
  for (it in seq_len(max_iter)) {
    prev <- out
    for (g in groups) {
      out <- .absorb_demean_one(out, g)
    }
    diff <- max(abs(out - prev), na.rm = TRUE)
    if (!is.finite(diff) || diff < tol) break
  }
  out
}


#' Apply absorb projection to a panel
#'
#' Demeans the dependent variable, the regressors, and any user-supplied
#' CSA source variables using the given grouping factors. Returns the
#' panel with the affected numeric columns overwritten by their within
#' transforms. Non-numeric columns and columns not listed in the target
#' set are left untouched.
#'
#' @param panel A panel data.frame from \code{.make_panel()}.
#' @param target_vars Character vector: column names to be demeaned.
#' @param groups A list of grouping factors from \code{.absorb_resolve()}.
#' @return The panel with the target columns demeaned.
#' @keywords internal
.absorb_apply <- function(panel, target_vars, groups) {
  if (is.null(groups) || length(groups) == 0L) return(panel)

  unit_var <- attr(panel, "unit_var")
  time_var <- attr(panel, "time_var")

  target_vars <- intersect(target_vars, names(panel))
  for (v in target_vars) {
    if (!is.numeric(panel[[v]])) next
    panel[[v]] <- .absorb_demean_multi(panel[[v]], groups)
  }

  attr(panel, "unit_var") <- unit_var
  attr(panel, "time_var") <- time_var
  panel
}
