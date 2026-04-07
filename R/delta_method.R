#' Delta Method
#'
#' Computes standard errors for nonlinear transformations of parameters
#' using the delta method: `V(g(b)) = G' V(b) G` where `G = dg/db`.
#'
#' @param g A function of the parameter vector returning the transformed values.
#' @param b Numeric vector: parameter estimates.
#' @param V Numeric matrix: variance-covariance of `b`.
#' @param h Numeric: step size for numerical differentiation. Default `1e-6`.
#' @return A list with `estimate` (transformed values) and `vcov` (delta-method
#'   variance).
#' @keywords internal
.delta_method <- function(g, b, V, h = 1e-6) {
  est <- g(b)
  K <- length(b)
  M <- length(est)

  # Numerical Jacobian
  G <- matrix(0, M, K)
  for (j in seq_len(K)) {
    b_plus <- b
    b_minus <- b
    b_plus[j] <- b[j] + h
    b_minus[j] <- b[j] - h
    G[, j] <- (g(b_plus) - g(b_minus)) / (2 * h)
  }

  V_g <- G %*% V %*% t(G)

  list(estimate = est, vcov = V_g)
}
