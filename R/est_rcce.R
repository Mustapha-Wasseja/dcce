#' Regularized CCE Estimator Internals
#'
#' Internal functions for the regularized CCE (rCCE) estimator of Juodis (2022).
#' Uses SVD of cross-sectional averages to extract principal components,
#' with the number of components selected by Ahn-Horenstein ER/GR criteria.
#'
#' @name rcce_estimator
#' @keywords internal
NULL


#' SVD of CSA matrix
#'
#' Computes thin SVD and returns right singular vectors (principal components
#' in the time dimension).
#'
#' @param Z Numeric matrix (T x K) of CSA variables.
#' @param npc Integer: number of principal components to retain. If `NULL`,
#'   determined by ER/GR criterion.
#' @param criterion Character: `"er"` (eigenvalue ratio) or `"gr"` (growth
#'   ratio). Default `"er"`.
#' @return A list with `V_k` (T x k matrix of PCs), `eigenvalues`, and `k`.
#' @keywords internal
.svd_csa <- function(Z, npc = NULL, criterion = c("er", "gr")) {
  criterion <- match.arg(criterion)
  Z <- as.matrix(Z)
  T_val <- nrow(Z)
  K <- ncol(Z)
  K_max <- floor(min(T_val, K) / 2)
  K_max <- max(K_max, 1L)

  s <- svd(Z, nu = min(T_val, K), nv = 0)
  eigenvalues <- s$d^2

  if (is.null(npc)) {
    if (criterion == "er") {
      npc <- .ahn_horenstein_er(eigenvalues, K_max)
    } else {
      npc <- .ahn_horenstein_gr(eigenvalues, K_max)
    }
  }

  npc <- min(npc, ncol(s$u))
  # Left singular vectors give T-dimensional PCs
  V_k <- s$u[, 1:npc, drop = FALSE] * sqrt(T_val)

  list(V_k = V_k, eigenvalues = eigenvalues, k = npc)
}


#' Ahn-Horenstein ER criterion
#'
#' Eigenvalue Ratio criterion: `k* = argmax_{k=1,...,K_max} lambda[k] / lambda[k+1]`
#'
#' @param eigenvalues Numeric vector of eigenvalues (decreasing order).
#' @param K_max Integer: maximum number of factors to consider.
#' @return Integer: estimated number of factors.
#' @keywords internal
.ahn_horenstein_er <- function(eigenvalues, K_max) {
  if (length(eigenvalues) < 2) return(1L)
  K_max <- min(K_max, length(eigenvalues) - 1L)
  if (K_max < 1L) return(1L)

  ratios <- eigenvalues[1:K_max] / eigenvalues[2:(K_max + 1)]
  as.integer(which.max(ratios))
}


#' Ahn-Horenstein GR criterion
#'
#' Growth Ratio criterion:
#' `GR(k) = log(1 + lambda[k]/V_k) / log(1 + lambda[k+1]/V_k)`
#' where `V_k = sum(lambda[k+1:end])`.
#'
#' @param eigenvalues Numeric vector of eigenvalues.
#' @param K_max Integer: maximum number of factors.
#' @return Integer: estimated number of factors.
#' @keywords internal
.ahn_horenstein_gr <- function(eigenvalues, K_max) {
  if (length(eigenvalues) < 3) return(1L)
  K_max <- min(K_max, length(eigenvalues) - 2L)
  if (K_max < 1L) return(1L)

  gr <- numeric(K_max)
  total <- sum(eigenvalues)
  for (k in 1:K_max) {
    V_k <- total - sum(eigenvalues[1:k])
    if (V_k <= 0) {
      gr[k] <- 0
      next
    }
    num <- log(1 + eigenvalues[k] / V_k)
    den <- log(1 + eigenvalues[k + 1] / V_k)
    gr[k] <- if (den > 0) num / den else 0
  }
  as.integer(which.max(gr))
}
