#' Swamy Slope Heterogeneity Test
#'
#' Tests the null hypothesis that all slope coefficients are identical
#' across cross-sectional units, against the alternative of heterogeneous
#' slopes. Implements the Swamy (1970) chi-square test and the Pesaran &
#' Yamagata (2008) standardised dispersion statistic.
#'
#' The Swamy statistic is
#' \deqn{\tilde{S} = \sum_{i=1}^N (\hat\beta_i - \hat\beta^*)'
#'       \frac{X_i' M_\tau X_i}{\hat\sigma_i^2}
#'       (\hat\beta_i - \hat\beta^*),}
#' where \eqn{\hat\beta_i} is the unit-level OLS estimate,
#' \eqn{\hat\beta^*} is the weighted pooled estimate, and
#' \eqn{M_\tau = I - \tau(\tau'\tau)^{-1}\tau'} projects off the
#' intercept. Under \eqn{H_0} (homogeneous slopes), \eqn{\tilde{S}}
#' is asymptotically \eqn{\chi^2_{k(N-1)}}.
#'
#' Pesaran & Yamagata (2008) propose a standardised version that is
#' asymptotically standard normal:
#' \deqn{\tilde\Delta = \sqrt{N} \, \frac{N^{-1}\tilde{S} - k}{\sqrt{2k}}.}
#'
#' @param object A `dcce_fit` object (typically with `model = "mg"`,
#'   `"cce"`, or `"dcce"`).
#' @return An object of class `dcce_swamy` with elements `S_stat`,
#'   `delta_stat`, `df`, `p_swamy`, `p_delta`, `N`, `k`.
#'
#' @references
#' Swamy, P. A. V. B. (1970). Efficient inference in a random coefficient
#' regression model. *Econometrica*, 38(2), 311-323.
#'
#' Pesaran, M. H., & Yamagata, T. (2008). Testing slope homogeneity in
#' large panels. *Journal of Econometrics*, 142(1), 50-93.
#'
#' @export
swamy_test <- function(object) {
  if (!inherits(object, "dcce_fit")) {
    cli::cli_abort("{.arg object} must be a {.cls dcce_fit} object.")
  }

  coef_names <- object$x_names
  if (object$include_constant) coef_names <- c("(Intercept)", coef_names)

  # Only heterogeneous slopes (exclude intercept from the test)
  slope_names <- setdiff(coef_names, "(Intercept)")
  if (length(slope_names) == 0L) {
    cli::cli_abort("No slope coefficients available for the Swamy test.")
  }

  # Pooled estimate: weighted average using each unit's X'M_tau X
  # Here we use the simple weighted average with weights 1/sigma_i^2 * (X_i' M X_i)
  unit_results <- object$unit_results
  N <- length(unit_results)
  k <- length(slope_names)

  # Stack unit beta_i on slopes only
  B <- matrix(NA_real_, N, k)
  colnames(B) <- slope_names
  sigma2 <- numeric(N)
  XMX <- vector("list", N)

  for (i in seq_along(unit_results)) {
    res_i <- unit_results[[i]]
    if (!all(slope_names %in% names(res_i$b))) {
      B[i, ] <- NA_real_
      next
    }
    B[i, ] <- res_i$b[slope_names]
    sigma2[i] <- res_i$sigma2

    # X_i' M_tau X_i = V_i^{-1} * sigma_i^2, where V_i is the OLS vcov
    V_slopes <- res_i$V[slope_names, slope_names, drop = FALSE]
    # Invert (use pseudoinverse if singular)
    XMX_i <- tryCatch(solve(V_slopes) * res_i$sigma2,
                      error = function(e) .pinv(V_slopes) * res_i$sigma2)
    XMX[[i]] <- XMX_i
  }

  keep <- stats::complete.cases(B) & sigma2 > 0
  if (sum(keep) < 2L) {
    cli::cli_abort("Not enough valid units for the Swamy test.")
  }
  B <- B[keep, , drop = FALSE]
  sigma2 <- sigma2[keep]
  XMX <- XMX[keep]
  N <- nrow(B)

  # Weighted pooled estimate: beta_star = (sum XMX_i/sigma_i^2)^{-1} sum (XMX_i/sigma_i^2) b_i
  W_sum <- matrix(0, k, k)
  Wb_sum <- numeric(k)
  for (i in seq_len(N)) {
    Wi <- XMX[[i]] / sigma2[i]
    W_sum <- W_sum + Wi
    Wb_sum <- Wb_sum + Wi %*% B[i, ]
  }
  beta_star <- tryCatch(solve(W_sum, Wb_sum),
                        error = function(e) .pinv(W_sum) %*% Wb_sum)
  beta_star <- as.numeric(beta_star)
  names(beta_star) <- slope_names

  # Swamy statistic
  S_stat <- 0
  for (i in seq_len(N)) {
    diff_i <- B[i, ] - beta_star
    S_stat <- S_stat + as.numeric(t(diff_i) %*% (XMX[[i]] / sigma2[i]) %*% diff_i)
  }

  df <- k * (N - 1L)
  p_swamy <- stats::pchisq(S_stat, df = df, lower.tail = FALSE)

  # Pesaran-Yamagata delta
  delta_stat <- sqrt(N) * (S_stat / N - k) / sqrt(2 * k)
  p_delta <- 2 * stats::pnorm(-abs(delta_stat))

  out <- list(
    S_stat     = as.numeric(S_stat),
    delta_stat = as.numeric(delta_stat),
    df         = as.integer(df),
    p_swamy    = as.numeric(p_swamy),
    p_delta    = as.numeric(p_delta),
    N          = N,
    k          = k,
    beta_star  = beta_star
  )
  class(out) <- "dcce_swamy"
  out
}


#' Print a dcce_swamy object
#'
#' @param x A `dcce_swamy` object.
#' @param ... Ignored.
#' @return Invisibly returns `x`.
#' @export
print.dcce_swamy <- function(x, ...) {
  cat("\nSwamy / Pesaran-Yamagata Slope Homogeneity Test\n")
  cat("===============================================\n")
  cat(sprintf("N = %d, k = %d\n", x$N, x$k))
  cat("H0: slope coefficients are homogeneous across units\n\n")
  cat(sprintf("Swamy S           = %10.4f    df = %d    p-value = %.4f\n",
              x$S_stat, x$df, x$p_swamy))
  cat(sprintf("Pesaran-Yamagata  = %10.4f                 p-value = %.4f\n",
              x$delta_stat, x$p_delta))
  cat("\n")
  invisible(x)
}


#' Hausman-type Test: MG vs Pooled
#'
#' Tests the null that the pooled and MG estimators are both consistent
#' (i.e. slopes are homogeneous and the pooled estimator is efficient)
#' against the alternative that slopes are heterogeneous and only MG is
#' consistent. The test statistic is
#' \deqn{H = (\hat\beta_{MG} - \hat\beta_{pool})'
#'           [V_{MG} - V_{pool}]^{-1}
#'           (\hat\beta_{MG} - \hat\beta_{pool}) \sim \chi^2_k.}
#'
#' @param object A `dcce_fit` object.
#' @return An object of class `dcce_hausman` with the test statistic,
#'   degrees of freedom, and p-value.
#'
#' @export
hausman_test <- function(object) {
  if (!inherits(object, "dcce_fit")) {
    cli::cli_abort("{.arg object} must be a {.cls dcce_fit} object.")
  }

  coef_names <- object$x_names  # heterogeneous slopes
  slope_names <- coef_names  # exclude intercept (not in x_names)

  b_mg <- object$coefficients[slope_names]
  V_mg <- object$vcov[slope_names, slope_names, drop = FALSE]

  # Compute pooled estimate from the unit-level data
  unit_results <- object$unit_results
  N <- length(unit_results)

  # Weighted pooled: beta_pool = (sum (X_i' X_i / sigma_i^2))^{-1} sum(X_i'y_i/sigma_i^2)
  # Using the fact that for each unit:
  #   b_i = (X_i'X_i)^{-1} X_i'y_i, V_i = sigma_i^2 (X_i'X_i)^{-1}
  # So (X_i'X_i) / sigma_i^2 = V_i^{-1}, and (X_i'y_i)/sigma_i^2 = V_i^{-1} b_i
  Vinv_sum <- matrix(0, length(slope_names), length(slope_names))
  Vinv_b_sum <- numeric(length(slope_names))

  for (res_i in unit_results) {
    if (!all(slope_names %in% names(res_i$b))) next
    V_i <- res_i$V[slope_names, slope_names, drop = FALSE]
    b_i <- res_i$b[slope_names]
    Vinv_i <- tryCatch(solve(V_i), error = function(e) .pinv(V_i))
    Vinv_sum <- Vinv_sum + Vinv_i
    Vinv_b_sum <- Vinv_b_sum + Vinv_i %*% b_i
  }

  V_pool <- tryCatch(solve(Vinv_sum), error = function(e) .pinv(Vinv_sum))
  b_pool <- as.numeric(V_pool %*% Vinv_b_sum)
  names(b_pool) <- slope_names

  diff_vec <- b_mg - b_pool
  D <- V_mg - V_pool
  # Ensure D is non-negative definite
  D_inv <- tryCatch(solve(D), error = function(e) .pinv(D))
  H_stat <- as.numeric(t(diff_vec) %*% D_inv %*% diff_vec)
  df <- length(slope_names)
  p_val <- stats::pchisq(abs(H_stat), df = df, lower.tail = FALSE)

  out <- list(
    statistic = H_stat,
    df        = df,
    p_value   = p_val,
    b_mg      = b_mg,
    b_pool    = b_pool,
    diff      = diff_vec
  )
  class(out) <- "dcce_hausman"
  out
}


#' Print a dcce_hausman object
#'
#' @param x A `dcce_hausman` object.
#' @param ... Ignored.
#' @return Invisibly returns `x`.
#' @export
print.dcce_hausman <- function(x, ...) {
  cat("\nHausman-type Test: MG vs Pooled\n")
  cat("================================\n")
  cat("H0: slopes are homogeneous (pooled consistent)\n\n")
  cat(sprintf("H statistic = %10.4f    df = %d    p-value = %.4f\n",
              x$statistic, x$df, x$p_value))
  cat("\nMG vs Pooled coefficients:\n")
  tab <- data.frame(
    MG       = round(x$b_mg, 4),
    Pooled   = round(x$b_pool, 4),
    Diff     = round(x$diff, 4),
    row.names = names(x$b_mg)
  )
  print(tab)
  cat("\n")
  invisible(x)
}
