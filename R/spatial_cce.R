#' Spatial Common Correlated Effects
#'
#' Helpers to replace the global cross-sectional averages of classical
#' CCE with **spatially-weighted local averages**. Given a row-normalised
#' spatial weight matrix \eqn{W} (\eqn{N \times N}, zero diagonal, rows
#' summing to one), the cross-sectional average for unit \eqn{i} at
#' time \eqn{t} becomes
#' \deqn{\bar{y}_{it} = \sum_{j=1}^N w_{ij} y_{jt},}
#' i.e. a weighted average over unit \eqn{i}'s spatial neighbours rather
#' than a global mean across the whole panel. This allows the common
#' factor proxy to vary across units in a way that respects the spatial
#' topology of the data (geographical contiguity, trade links, input-
#' output connections, etc.).
#'
#' The weight matrix is supplied directly by the user; `dcce` makes no
#' assumptions about how it was constructed. It must be square with
#' dimensions equal to the number of units, and rows should ideally sum
#' to one (the helper will row-normalise silently if they do not).
#'
#' @name spatial_cce_internals
#' @keywords internal
NULL


#' Validate and normalise a spatial weight matrix
#'
#' @param W A numeric matrix (N x N).
#' @param units Character vector of unit identifiers in the panel, used
#'   both for dimension checks and for aligning rows/columns.
#' @return A row-normalised matrix with zero diagonal.
#' @keywords internal
.spatial_validate_W <- function(W, units) {
  if (!is.matrix(W)) W <- as.matrix(W)
  if (!is.numeric(W)) {
    cli::cli_abort("{.arg spatial_weights} must be a numeric matrix.")
  }
  N <- length(units)
  if (nrow(W) != N || ncol(W) != N) {
    cli::cli_abort(c(
      "{.arg spatial_weights} must be {N} x {N} (one row/column per unit).",
      "i" = "Got {nrow(W)} x {ncol(W)}."
    ))
  }

  # Align by row/column names if present; otherwise trust the order
  if (!is.null(rownames(W)) && !is.null(colnames(W))) {
    missing_rows <- setdiff(units, rownames(W))
    missing_cols <- setdiff(units, colnames(W))
    if (length(missing_rows) > 0L || length(missing_cols) > 0L) {
      cli::cli_abort(c(
        "{.arg spatial_weights} is missing entries for some units.",
        "i" = "Missing row{?s}: {.val {missing_rows}}",
        "i" = "Missing column{?s}: {.val {missing_cols}}"
      ))
    }
    W <- W[units, units, drop = FALSE]
  } else {
    rownames(W) <- colnames(W) <- units
  }

  # Zero out the diagonal
  diag(W) <- 0

  # Row-normalise (silent; if a row is all zero, leave it as zero)
  row_sums <- rowSums(W, na.rm = TRUE)
  nz <- row_sums > 0
  W[nz, ] <- W[nz, , drop = FALSE] / row_sums[nz]

  # Sanity check: no NAs left
  if (any(!is.finite(W))) {
    cli::cli_abort(
      "{.arg spatial_weights} contains non-finite values after normalisation."
    )
  }
  W
}


#' Build spatial (local) cross-sectional averages
#'
#' For each variable in \code{vars} and each unit \eqn{i} at time \eqn{t},
#' computes \eqn{\bar y^W_{i,t} = \sum_j w_{ij} y_{j,t}} where \eqn{W} is
#' the row-normalised spatial weight matrix. Appends these as new
#' columns on the panel, alongside lagged versions when requested.
#'
#' Compared to \code{.build_csa()} (which appends a single global series
#' per variable, broadcast to all units), the spatial variant produces
#' unit-specific values — each unit sees its own neighbourhood average.
#'
#' @param panel A panel data.frame from \code{.make_panel()}.
#' @param vars Character vector of variables to build spatial CSAs for.
#' @param W A validated spatial weight matrix from
#'   \code{.spatial_validate_W()}, with rows/columns labelled by unit id.
#' @param lags Integer scalar or named integer vector of lag orders.
#'   Default \code{0L} (contemporaneous only).
#' @return The panel with new \code{csa_*} columns appended.
#' @keywords internal
.build_spatial_csa <- function(panel, vars, W, lags = 0L) {
  unit_var <- attr(panel, "unit_var")
  time_var <- attr(panel, "time_var")
  units <- rownames(W)
  times <- sort(unique(panel[[time_var]]))
  N <- length(units)
  T_val <- length(times)

  if (length(lags) == 1L && is.null(names(lags))) {
    lag_vec <- stats::setNames(rep(as.integer(lags), length(vars)), vars)
  } else {
    lag_vec <- as.integer(lags[vars])
    names(lag_vec) <- vars
  }

  # For each variable, build the (N x T) matrix Y[i,t], then compute W %*% Y
  # to get the spatial lag matrix WY[i,t]. That matrix is then sliced back
  # into the long panel as a new column.
  for (v in vars) {
    Y <- matrix(NA_real_, N, T_val,
                dimnames = list(units, as.character(times)))
    for (i in seq_along(units)) {
      idx <- which(panel[[unit_var]] == units[i])
      t_pos <- match(panel[[time_var]][idx], times)
      Y[i, t_pos] <- panel[[v]][idx]
    }

    # Replace NA with 0 for the matrix multiply; the weight matrix
    # already handles the "no data" case by producing a weighted average
    # that ignores missing neighbours only if the user did that
    # construction upstream.
    Y_clean <- Y
    Y_clean[is.na(Y_clean)] <- 0
    WY <- W %*% Y_clean

    # Broadcast back to the long panel
    col_name <- paste0("csa_", v)
    panel[[col_name]] <- NA_real_
    for (i in seq_along(units)) {
      idx <- which(panel[[unit_var]] == units[i])
      t_pos <- match(panel[[time_var]][idx], times)
      panel[[col_name]][idx] <- WY[i, t_pos]
    }

    # Lagged spatial CSAs
    p_v <- lag_vec[v]
    if (p_v > 0L) {
      for (l in seq_len(p_v)) {
        lag_name <- paste0("csa_", v, "_L", l)
        panel[[lag_name]] <- .lag_panel(panel, col_name, l)
      }
    }
  }

  attr(panel, "unit_var") <- unit_var
  attr(panel, "time_var") <- time_var
  panel
}
