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
