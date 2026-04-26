#' Interactive Fixed Effects (IFE) Estimator
#'
#' Implements the iterative principal-components estimator of Bai (2009)
#' for panel data with interactive fixed effects:
#' \deqn{y_{it} = \alpha_i + \beta' x_{it} + \lambda_i' f_t + u_{it},}
#' where \eqn{f_t} is an \eqn{r}-vector of unobserved common factors and
#' \eqn{\lambda_i} is unit \eqn{i}'s vector of factor loadings. Unlike CCE
#' (which proxies \eqn{f_t} via cross-sectional averages), IFE estimates
#' \eqn{f_t} and \eqn{\lambda_i} directly via principal components on the
#' residual matrix.
#'
#' The algorithm iterates between:
#' \enumerate{
#'   \item Given \eqn{(\hat\beta, \hat\alpha_i)}, compute residuals
#'     \eqn{\hat e_{it}} and extract the first \eqn{r} principal components
#'     as \eqn{\hat f_t}, with loadings \eqn{\hat\lambda_i}.
#'   \item Given \eqn{(\hat f_t, \hat\lambda_i)}, re-estimate
#'     \eqn{\hat\beta} and \eqn{\hat\alpha_i} by pooled OLS on the
#'     defactored data.
#' }
#' until convergence. The number of factors \eqn{r} can be specified or
#' selected by information criteria (BIC3 of Bai & Ng 2002).
#'
#' **Important:** IFE estimates a **common** (pooled) \eqn{\beta}, not
#' heterogeneous unit-specific slopes. If you need heterogeneous slopes,
#' use CCE/DCCE instead.
#'
#' @name ife_estimator
#' @keywords internal
NULL


#' Post-process a dcce_fit into an IFE fit
#'
#' Called from \code{dcce()} when \code{model = "ife"}. Runs the Bai (2009)
#' iterative PC estimator on the panel stored inside the fit.
#'
#' @param panel Prepared panel data.frame.
#' @param y_name Character: dependent variable name.
#' @param x_names Character: regressor names.
#' @param n_factors Integer or \code{NULL}: number of factors. If
#'   \code{NULL}, the BIC3 criterion of Bai & Ng (2002) is used to
#'   select \eqn{r}.
#' @param max_iter Integer: maximum iterations. Default 100.
#' @param tol Numeric: convergence tolerance on beta. Default 1e-8.
#' @param include_constant Logical: include unit FE? Default TRUE.
#' @return A list with \code{coefficients}, \code{vcov}, \code{se},
#'   \code{factors}, \code{loadings}, \code{n_factors}, \code{residuals},
#'   \code{r2}, and related fields.
#' @keywords internal
.fit_ife <- function(panel, y_name, x_names,
                     n_factors = NULL,
                     max_iter  = 100L,
                     tol       = 1e-8,
                     include_constant = TRUE) {
  unit_var <- attr(panel, "unit_var")
  time_var <- attr(panel, "time_var")
  units <- unique(panel[[unit_var]])
  times <- sort(unique(panel[[time_var]]))
  N <- length(units)
  T_val <- length(times)
  K <- length(x_names)

  # Build (N*T) x K design matrix and y vector, with unit dummies
  # Filter to balanced sub-panel for the PC step
  Y_mat <- matrix(NA_real_, N, T_val)
  X_arr <- array(NA_real_, dim = c(N, T_val, K))

  for (i in seq_along(units)) {
    idx <- which(panel[[unit_var]] == units[i])
    t_pos <- match(panel[[time_var]][idx], times)
    Y_mat[i, t_pos] <- panel[[y_name]][idx]
    for (k in seq_len(K)) {
      X_arr[i, t_pos, k] <- panel[[x_names[k]]][idx]
    }
  }

  # Use only complete (balanced) rows/cols for the iterative step
  row_ok <- rowSums(is.na(Y_mat)) == 0 & apply(X_arr, 1, function(r) !any(is.na(r)))
  col_ok <- colSums(is.na(Y_mat[row_ok, , drop = FALSE])) == 0
  Y <- Y_mat[row_ok, col_ok, drop = FALSE]
  X <- X_arr[row_ok, col_ok, , drop = FALSE]
  N_eff <- nrow(Y)
  T_eff <- ncol(Y)

  if (N_eff < 3L || T_eff < 3L) {
    cli::cli_abort(
      "IFE: balanced sub-panel too small ({N_eff} x {T_eff}). Need at least 3 x 3."
    )
  }

  # Demean within unit (absorb unit FE)
  if (include_constant) {
    Y_dm <- Y - rowMeans(Y)
    X_dm <- X
    for (k in seq_len(K)) {
      X_dm[, , k] <- X[, , k] - rowMeans(X[, , k])
    }
  } else {
    Y_dm <- Y
    X_dm <- X
  }

  # Determine number of factors
  if (is.null(n_factors)) {
    n_factors <- .ife_select_factors(Y_dm, max_r = min(10L, T_eff - 1L, N_eff - 1L))
  }
  n_factors <- max(1L, as.integer(n_factors))

  # Initialise beta from pooled OLS (ignoring factors)
  # Stack: vec(Y_dm) ~ [X_dm[,,1] | ... | X_dm[,,K]] * beta
  y_vec <- as.numeric(t(Y_dm))
  X_pool <- matrix(NA_real_, N_eff * T_eff, K)
  for (k in seq_len(K)) {
    X_pool[, k] <- as.numeric(t(X_dm[, , k]))
  }
  colnames(X_pool) <- x_names
  beta <- solve(crossprod(X_pool), crossprod(X_pool, y_vec))
  beta <- as.numeric(beta)
  names(beta) <- x_names

  # Iterate
  for (it in seq_len(max_iter)) {
    beta_old <- beta

    # Residuals given current beta
    E <- Y_dm
    for (k in seq_len(K)) E <- E - beta[k] * X_dm[, , k]

    # PC on the residual matrix: SVD of E = U D V'
    # Factors = sqrt(T) * first r columns of V (T x r)
    # Loadings = E %*% F / T (N x r)
    sv <- svd(E, nu = n_factors, nv = n_factors)
    F_hat <- sv$v[, 1:n_factors, drop = FALSE] * sqrt(T_eff)
    L_hat <- E %*% F_hat / T_eff

    # Project out factors from Y_dm and X_dm
    P_F <- F_hat %*% solve(crossprod(F_hat), t(F_hat))   # T x T projection
    M_F <- diag(T_eff) - P_F

    Y_proj <- Y_dm %*% M_F
    X_proj <- array(NA_real_, dim = c(N_eff, T_eff, K))
    for (k in seq_len(K)) {
      X_proj[, , k] <- X_dm[, , k] %*% M_F
    }

    # Pooled OLS on factor-projected data
    y_proj_vec <- as.numeric(t(Y_proj))
    X_proj_pool <- matrix(NA_real_, N_eff * T_eff, K)
    for (k in seq_len(K)) {
      X_proj_pool[, k] <- as.numeric(t(X_proj[, , k]))
    }
    beta <- as.numeric(solve(crossprod(X_proj_pool), crossprod(X_proj_pool, y_proj_vec)))
    names(beta) <- x_names

    if (max(abs(beta - beta_old)) < tol) break
  }

  # Final residuals and variance
  E_final <- Y_dm
  for (k in seq_len(K)) E_final <- E_final - beta[k] * X_dm[, , k]
  E_defac <- E_final %*% M_F
  n_obs <- N_eff * T_eff
  sigma2 <- sum(E_defac^2) / (n_obs - N_eff - K - N_eff * n_factors)

  # Variance of beta (Bai 2009, simplified)
  V_beta <- sigma2 * solve(crossprod(X_proj_pool))
  rownames(V_beta) <- colnames(V_beta) <- x_names
  se_beta <- sqrt(diag(V_beta))

  # R-squared
  tss <- sum(Y_dm^2)
  rss <- sum(E_defac^2)
  r2 <- if (tss > 0) 1 - rss / tss else NA_real_

  list(
    coefficients = beta,
    vcov         = V_beta,
    se           = se_beta,
    factors      = F_hat,
    loadings     = L_hat,
    n_factors    = n_factors,
    residuals    = as.numeric(t(E_defac)),
    r2           = r2,
    rmse         = sqrt(mean(E_defac^2)),
    N            = N_eff,
    T_val        = T_eff,
    n_obs        = n_obs,
    iterations   = it
  )
}


#' Select number of factors via BIC3 (Bai & Ng 2002)
#' @keywords internal
.ife_select_factors <- function(Y_dm, max_r = 10L) {
  N <- nrow(Y_dm)
  T_val <- ncol(Y_dm)
  NT <- N * T_val
  bic <- rep(Inf, max_r)

  for (r in seq_len(max_r)) {
    sv <- svd(Y_dm, nu = r, nv = r)
    F_hat <- sv$v[, 1:r, drop = FALSE] * sqrt(T_val)
    L_hat <- Y_dm %*% F_hat / T_val
    E <- Y_dm - L_hat %*% t(F_hat)
    sigma2 <- sum(E^2) / NT
    pen <- r * ((N + T_val) / NT) * log(min(N, T_val))
    bic[r] <- log(sigma2) + pen
  }
  which.min(bic)
}
