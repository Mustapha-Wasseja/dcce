#' Mean Group Estimator Internals
#'
#' Internal functions for the Pesaran & Smith (1995) Mean Group estimator.
#'
#' @name mg_estimator
#' @keywords internal
NULL

# ──────────────────────────────────────────────────────────────────────────────
# .unit_ols — OLS for a single cross-section
# ──────────────────────────────────────────────────────────────────────────────

#' Unit-level OLS
#'
#' Fast OLS for a single cross-sectional unit. Uses normal equations
#' rather than `lm()` for speed when called in a loop.
#'
#' @param y Numeric vector: dependent variable.
#' @param X Numeric matrix: regressors (including intercept if desired).
#' @return A list with elements:
#'   \describe{
#'     \item{b}{Coefficient vector.}
#'     \item{V}{Variance-covariance matrix of coefficients.}
#'     \item{e}{Residual vector.}
#'     \item{r2}{R-squared.}
#'     \item{df_resid}{Residual degrees of freedom.}
#'     \item{sigma2}{Residual variance (s^2).}
#'   }
#' @keywords internal
.unit_ols <- function(y, X) {
  y <- as.numeric(y)
  X <- as.matrix(X)
  n <- length(y)
  k <- ncol(X)
  df_resid <- n - k

  XtX <- crossprod(X)
  Xty <- crossprod(X, y)

  # Check condition
  rc <- rcond(XtX)
  if (is.na(rc) || rc < .Machine$double.eps^0.5) {
    # Use SVD-based pseudoinverse
    XtX_inv <- tryCatch(
      solve(XtX),
      error = function(e) .pinv(XtX)
    )
  } else {
    XtX_inv <- solve(XtX)
  }

  b <- as.numeric(XtX_inv %*% Xty)
  names(b) <- colnames(X)

  e <- as.numeric(y - X %*% b)
  sigma2 <- sum(e^2) / df_resid

  V <- sigma2 * XtX_inv
  colnames(V) <- rownames(V) <- colnames(X)

  SST <- sum((y - mean(y))^2)
  r2 <- if (SST > 0) 1 - sum(e^2) / SST else NA_real_


  list(
    b        = b,
    V        = V,
    e        = e,
    r2       = r2,
    df_resid = df_resid,
    sigma2   = sigma2
  )
}


# ──────────────────────────────────────────────────────────────────────────────
# .run_unit_loop — dispatcher with Rcpp fast path and optional parallelism
# ──────────────────────────────────────────────────────────────────────────────

#' Run unit-level OLS over a list of (y, X) pairs
#'
#' Dispatcher used by \code{dcce()}. Accepts a named list where each element
#' is a list with elements \code{y} and \code{X}. Returns a list of unit-level
#' OLS result lists with the same names. Internally dispatches to:
#' \itemize{
#'   \item the C++ batch routine (\code{.batch_ols_cpp}) when
#'         \code{fast = TRUE} and the compiled shared library is available;
#'   \item \code{parallel::mclapply()} when \code{n_cores > 1} and the
#'         platform is Unix/macOS;
#'   \item the pure-R \code{.unit_ols()} in all other cases.
#' }
#'
#' @param panel_list Named list of \code{list(y, X)} pairs.
#' @param fast Logical: use the compiled C++ fast path when available?
#'   Default \code{TRUE}.
#' @param n_cores Integer: number of cores for parallel unit estimation.
#'   Only effective on Unix/macOS. Default \code{1L} (no parallelism).
#' @return A named list of unit-level OLS results matching the structure
#'   returned by \code{.unit_ols()}.
#' @keywords internal
.run_unit_loop <- function(panel_list, fast = TRUE, n_cores = 1L) {

  # Empty panel list — return an empty named list rather than crashing
  # the C++ side. dcce() will then surface a clean diagnostic to the user.
  if (length(panel_list) == 0L) return(list())

  # Helper: post-process a single C++ result so its shape matches
  # pure-R .unit_ols(): b is a named numeric vector, V is a named matrix,
  # e is a plain numeric vector, and r2 / sigma2 are scalars.
  finalize_cpp <- function(res, u) {
    if (is.null(res) || is.null(res$b)) return(res)
    cn <- colnames(u$X)
    res$b <- as.numeric(res$b)
    if (!is.null(cn)) names(res$b) <- cn
    if (!is.null(res$V)) {
      res$V <- as.matrix(res$V)
      if (!is.null(cn) && nrow(res$V) == length(cn)) {
        rownames(res$V) <- cn
        colnames(res$V) <- cn
      }
    }
    if (!is.null(res$e)) res$e <- as.numeric(res$e)
    if (!is.null(res$r2))     res$r2 <- as.numeric(res$r2)
    if (!is.null(res$sigma2)) res$sigma2 <- as.numeric(res$sigma2)
    res
  }

  # Parallel path (Unix/macOS only; Windows silently falls back)
  if (n_cores > 1L && .Platform$OS.type == "unix" &&
      requireNamespace("parallel", quietly = TRUE)) {
    worker <- if (fast && .dcce_cpp_available()) {
      function(u) finalize_cpp(.unit_ols_cpp(u$X, u$y), u)
    } else {
      function(u) .unit_ols(u$y, u$X)
    }
    return(parallel::mclapply(panel_list, worker, mc.cores = n_cores))
  }

  # Sequential C++ fast path
  if (fast && .dcce_cpp_available()) {
    res <- tryCatch(
      .batch_ols_cpp(panel_list),
      error = function(e) {
        cli::cli_warn(c(
          "RcppArmadillo batch OLS failed; falling back to pure R.",
          "i" = "Error: {conditionMessage(e)}"
        ))
        NULL
      }
    )
    if (!is.null(res)) {
      # Attach column names from each unit's X back to b/V
      for (nm in names(res)) {
        res[[nm]] <- finalize_cpp(res[[nm]], panel_list[[nm]])
      }
      return(res)
    }
  }

  # Pure-R fallback (always available)
  lapply(panel_list, function(u) .unit_ols(u$y, u$X))
}


#' Check whether the compiled C++ unit-OLS routines are available
#' @keywords internal
.dcce_cpp_available <- function() {
  exists(".batch_ols_cpp", where = asNamespace("dcce"), inherits = FALSE) &&
    exists(".unit_ols_cpp",  where = asNamespace("dcce"), inherits = FALSE)
}


# ──────────────────────────────────────────────────────────────────────────────
# .mg_aggregate — Mean Group coefficient
# ──────────────────────────────────────────────────────────────────────────────

#' Mean Group aggregation
#'
#' Computes the unweighted average of unit-level coefficient vectors.
#'
#' @param coef_list A list of numeric coefficient vectors (one per unit).
#' @return Numeric vector: Mean Group coefficient.
#' @keywords internal
.mg_aggregate <- function(coef_list) {
  B <- do.call(rbind, coef_list)
  b_mg <- colMeans(B)
  b_mg
}


# ──────────────────────────────────────────────────────────────────────────────
# .mg_variance — Mean Group variance
# ──────────────────────────────────────────────────────────────────────────────

#' Mean Group variance
#'
#' Non-parametric variance estimator for the MG coefficient:
#' `V = (1/N^2) * sum_i (b_i - b_mg)(b_i - b_mg)'`
#'
#' @param coef_list A list of numeric coefficient vectors.
#' @param b_mg Numeric vector: the MG coefficient (from `.mg_aggregate()`).
#' @return A variance-covariance matrix.
#' @keywords internal
.mg_variance <- function(coef_list, b_mg) {
  N <- length(coef_list)
  K <- length(b_mg)
  B <- do.call(rbind, coef_list)

  # Center
  B_centered <- sweep(B, 2, b_mg)
  V <- crossprod(B_centered) / N^2

  colnames(V) <- rownames(V) <- names(b_mg)
  V
}
