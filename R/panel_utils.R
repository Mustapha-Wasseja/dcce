#' Panel Data Utilities
#'
#' Internal functions for constructing, validating, and manipulating panel data.
#' Exported formula helpers `L()`, `D()`, and `Lrange()` allow lag/difference
#' operations inside `dcce()` formulas.
#'
#' @name panel_utils
#' @keywords internal
NULL

# ──────────────────────────────────────────────────────────────────────────────
# .pinv — Moore-Penrose pseudoinverse via SVD
# ──────────────────────────────────────────────────────────────────────────────

#' Moore-Penrose pseudoinverse
#'
#' @param M A numeric matrix.
#' @param tol Tolerance for zero singular values.
#' @return The pseudoinverse of `M`.
#' @keywords internal
.pinv <- function(M, tol = .Machine$double.eps^0.5) {
  s <- svd(M)
  d <- s$d
  pos <- d > max(tol * d[1L], 0)
  if (!any(pos)) return(matrix(0, ncol(M), nrow(M)))
  s$v[, pos, drop = FALSE] %*% (1 / d[pos] * t(s$u[, pos, drop = FALSE]))
}


# ──────────────────────────────────────────────────────────────────────────────
# .make_panel — validate and prepare panel data
# ──────────────────────────────────────────────────────────────────────────────

#' Prepare a panel data frame
#'
#' Validates the data and index columns, sorts by unit then time, and attaches
#' `"unit_var"` and `"time_var"` attributes. Accepts either the new
#' `unit_index`/`time_index` arguments or the legacy `index` argument.
#'
#' @param data A data.frame.
#' @param unit_index Character scalar: name of the unit identifier column.
#' @param time_index Character scalar: name of the time identifier column.
#' @param index Character vector of length 2 (legacy): `c(unit_col,
#'   time_col)`. Ignored if `unit_index` and `time_index` are provided.
#' @return A sorted data.frame with attributes `unit_var` and `time_var`.
#' @keywords internal
.make_panel <- function(data,
                        unit_index = NULL,
                        time_index = NULL,
                        index = NULL) {

  if (!is.data.frame(data)) {
    cli::cli_abort("{.arg data} must be a data frame, not {.cls {class(data)}}.")
  }

  # Resolve unit_var / time_var from new or legacy arguments
  # Support positional call: .make_panel(data, c("id","t")) where the

  # length-2 vector lands in unit_index
  if (!is.null(unit_index) && is.character(unit_index) &&
      length(unit_index) == 2L && is.null(time_index)) {
    # Legacy positional call
    unit_var <- unit_index[1L]
    time_var <- unit_index[2L]
  } else if (!is.null(unit_index) && !is.null(time_index)) {
    unit_var <- unit_index
    time_var <- time_index
  } else if (!is.null(index) && is.character(index) && length(index) == 2L) {
    unit_var <- index[1L]
    time_var <- index[2L]
  } else {
    cli::cli_abort(
      "Must provide {.arg unit_index} and {.arg time_index}, or a length-2
       {.arg index} vector."
    )
  }

  if (!unit_var %in% names(data)) {
    cli::cli_abort(
      "Unit variable {.val {unit_var}} not found in {.arg data}."
    )
  }
  if (!time_var %in% names(data)) {
    cli::cli_abort(
      "Time variable {.val {time_var}} not found in {.arg data}."
    )
  }

  # Check for duplicates
  dup <- duplicated(data[, c(unit_var, time_var), drop = FALSE])
  if (any(dup)) {
    cli::cli_abort(
      "Duplicate (unit, time) observations found. Panel must have unique
       index pairs."
    )
  }

  # Sort by unit, then time
  ord <- order(data[[unit_var]], data[[time_var]])
  data <- data[ord, , drop = FALSE]
  rownames(data) <- NULL

  attr(data, "unit_var") <- unit_var
  attr(data, "time_var") <- time_var
  data
}


# ──────────────────────────────────────────────────────────────────────────────
# .check_balance — check if panel is balanced
# ──────────────────────────────────────────────────────────────────────────────

#' Check panel balance
#'
#' @param panel A panel data.frame from `.make_panel()`.
#' @return A list with elements `is_balanced`, `T_i` (named vector of per-unit
#'   time periods), `T_min`, `T_max`, `T_bar`, `N`.
#' @keywords internal
.check_balance <- function(panel) {
  unit_var <- attr(panel, "unit_var")
  T_i <- table(panel[[unit_var]])
  T_i <- as.integer(T_i)
  names(T_i) <- names(table(panel[[unit_var]]))
  list(
    is_balanced = (length(unique(T_i)) == 1L),
    T_i         = T_i,
    T_min       = min(T_i),
    T_max       = max(T_i),
    T_bar       = mean(T_i),
    N           = length(T_i)
  )
}


# ──────────────────────────────────────────────────────────────────────────────
# .lag_panel — panel-aware lag
# ──────────────────────────────────────────────────────────────────────────────

#' Panel-aware lag
#'
#' Computes lagged values within each cross-sectional unit. Returns `NA` for
#' the first `k` observations of each unit to avoid cross-unit contamination.
#'
#' @param panel A panel data.frame from `.make_panel()`.
#' @param vars Character vector of column names to lag.
#' @param k Integer lag order (positive = lag, negative = lead).
#' @return A data.frame (or vector if single var) of lagged values, same
#'   length as `nrow(panel)`.
#' @keywords internal
.lag_panel <- function(panel, vars, k = 1L) {
  unit_var <- attr(panel, "unit_var")
  time_var <- attr(panel, "time_var")
  units <- unique(panel[[unit_var]])

  single <- (length(vars) == 1L)
  result <- matrix(NA_real_, nrow = nrow(panel), ncol = length(vars))
  colnames(result) <- paste0(vars, "_L", k)

  for (u in units) {
    idx <- which(panel[[unit_var]] == u)
    Ti <- length(idx)
    if (abs(k) >= Ti) next
    for (j in seq_along(vars)) {
      x <- panel[[vars[j]]][idx]
      if (k > 0L) {
        result[idx[(k + 1L):Ti], j] <- x[1L:(Ti - k)]
      } else if (k < 0L) {
        ak <- abs(k)
        result[idx[1L:(Ti - ak)], j] <- x[(ak + 1L):Ti]
      } else {
        result[idx, j] <- x
      }
    }
  }

  if (single) {
    return(as.numeric(result[, 1L]))
  }
  as.data.frame(result)
}


# ──────────────────────────────────────────────────────────────────────────────
# .diff_panel — panel-aware differencing
# ──────────────────────────────────────────────────────────────────────────────

#' Panel-aware differencing
#'
#' Computes k-th differences within each cross-sectional unit.
#'
#' @param panel A panel data.frame from `.make_panel()`.
#' @param vars Character vector of column names to difference.
#' @param k Integer difference order (default 1).
#' @return A data.frame (or vector if single var) of differenced values.
#' @keywords internal
.diff_panel <- function(panel, vars, k = 1L) {
  unit_var <- attr(panel, "unit_var")
  units <- unique(panel[[unit_var]])

  single <- (length(vars) == 1L)
  result <- matrix(NA_real_, nrow = nrow(panel), ncol = length(vars))
  colnames(result) <- paste0("d.", vars)

  for (u in units) {
    idx <- which(panel[[unit_var]] == u)
    Ti <- length(idx)
    if (k >= Ti) next
    for (j in seq_along(vars)) {
      x <- panel[[vars[j]]][idx]
      result[idx[(k + 1L):Ti], j] <- x[(k + 1L):Ti] - x[1L:(Ti - k)]
    }
  }

  if (single) {
    return(as.numeric(result[, 1L]))
  }
  as.data.frame(result)
}


# ──────────────────────────────────────────────────────────────────────────────
# Exported formula helpers: L(), D(), Lrange()
# ──────────────────────────────────────────────────────────────────────────────

#' Lag operator for dcce formulas
#'
#' Creates a lagged version of a variable for use inside `dcce()` formulas.
#' This function is evaluated during formula processing by `dcce()` and should
#' not be called directly on raw vectors outside of a dcce formula context.
#'
#' @param x A numeric vector (column name evaluated within `dcce()`).
#' @param k Integer lag order. Default 1. Positive values lag, negative lead.
#' @return A numeric vector of the same length as `x` with leading `NA`s.
#' @export
#' @examples
#' x <- c(10, 20, 30, 40, 50)
#' L(x, 1)   # NA 10 20 30 40
#' L(x, 2)   # NA NA 10 20 30
L <- function(x, k = 1L) {
  k <- as.integer(k)
  n <- length(x)
  if (abs(k) >= n) return(rep(NA_real_, n))
  result <- rep(NA_real_, n)
  if (k > 0L) {
    result[(k + 1L):n] <- x[1L:(n - k)]
  } else if (k < 0L) {
    ak <- abs(k)
    result[1L:(n - ak)] <- x[(ak + 1L):n]
  } else {
    result <- x
  }
  result
}

#' Difference operator for dcce formulas
#'
#' Creates a differenced version of a variable for use inside `dcce()` formulas.
#'
#' @param x A numeric vector.
#' @param k Integer difference order. Default 1.
#' @return A numeric vector of the same length as `x` with leading `NA`s.
#' @export
#' @examples
#' x <- c(10, 20, 30, 40, 50)
#' D(x, 1)   # NA 10 10 10 10
#' D(x, 2)   # NA NA 20 20 20
D <- function(x, k = 1L) {
  k <- as.integer(k)
  n <- length(x)
  if (k >= n) return(rep(NA_real_, n))
  result <- rep(NA_real_, n)
  result[(k + 1L):n] <- x[(k + 1L):n] - x[1L:(n - k)]
  result
}

#' Lag range operator for dcce formulas
#'
#' Creates multiple lagged versions of a variable from lag `k0` to lag `k1`.
#' Useful for distributed lag specifications.
#'
#' @param x A numeric vector.
#' @param k0 Integer start lag (inclusive).
#' @param k1 Integer end lag (inclusive).
#' @return A matrix with columns `L[k0]` through `L[k1]`.
#' @export
#' @examples
#' x <- c(10, 20, 30, 40, 50)
#' Lrange(x, 0, 2)   # 3 columns: lag 0, lag 1, lag 2
Lrange <- function(x, k0, k1) {
  if (k0 > k1) {
    cli::cli_abort("{.arg k0} must be <= {.arg k1}.")
  }
  lags <- k0:k1
  out <- matrix(NA_real_, nrow = length(x), ncol = length(lags))
  colnames(out) <- paste0("L", lags)
  for (i in seq_along(lags)) {
    out[, i] <- L(x, lags[i])
  }
  out
}
