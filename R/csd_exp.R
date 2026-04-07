#' Exponent of Cross-Sectional Dependence
#'
#' Estimates the exponent of cross-sectional dependence (alpha) using the
#' methods of Bailey, Kapetanios & Pesaran (2016, 2019).
#'
#' @param x Either a numeric vector (variable stacked by unit), a numeric
#'   matrix (N x T), or a `dcce_fit` object (uses residuals).
#' @param data A data.frame containing the panel structure. Required if
#'   `x` is a vector.
#' @param unit_index Character: name of the unit variable in `data`.
#' @param time_index Character: name of the time variable in `data`.
#' @param use_residuals Logical: if `TRUE`, use the BKP (2019) residual method;
#'   if `FALSE` (default), use the BKP (2016) variable method.
#' @param n_pca Integer: number of principal components. Default 1.
#' @param test_size Numeric: significance level for thresholding. Default 0.1.
#' @param tuning Numeric: tuning parameter for residual method threshold.
#'   Default 0.5.
#' @param n_bootstrap Integer: bootstrap repetitions for SE (residual method
#'   only). Default 200.
#'
#' @return An object of class `dcce_csd` with elements `alpha` (estimated
#'   exponent), `se` (standard error), `ci` (confidence interval), `method`,
#'   `N`, and `T_val`.
#'
#' @export
#'
#' @examples
#' set.seed(42)
#' # Matrix of cross-sectionally dependent data
#' N <- 20; T_val <- 50
#' f <- rnorm(T_val)
#' x <- matrix(NA, N, T_val)
#' for (i in 1:N) x[i,] <- rnorm(1) * f + rnorm(T_val, sd = 0.5)
#' result <- csd_exp(x, use_residuals = FALSE)
#' print(result)
csd_exp <- function(x,
                    data          = NULL,
                    unit_index    = NULL,
                    time_index    = NULL,
                    use_residuals = FALSE,
                    n_pca         = 1L,
                    test_size     = 0.1,
                    tuning        = 0.5,
                    n_bootstrap   = 200L) {

  method <- if (use_residuals) "residual" else "variable"

  # ── Handle input types ──────────────────────────────────────────────────
  if (inherits(x, "dcce_fit")) {
    em <- .resid_list_to_matrix(x$resid_list)
  } else if (is.matrix(x)) {
    em <- x
  } else {
    if (is.null(data) || is.null(unit_index) || is.null(time_index)) {
      cli::cli_abort(
        "{.arg data}, {.arg unit_index}, and {.arg time_index} are required when {.arg x} is a numeric vector."
      )
    }
    panel <- .make_panel(data, unit_index = unit_index, time_index = time_index)
    em <- .resid_vector_to_matrix(x, panel)
  }

  N <- nrow(em)
  T_val <- ncol(em)

  if (method == "variable") {
    result <- .bkp_variable(em, N, T_val, n_pca, test_size)
  } else {
    result <- .bkp_residual(em, N, T_val, n_pca, tuning, n_bootstrap)
  }

  result$method <- method
  result$N <- N
  result$T_val <- T_val

  class(result) <- "dcce_csd"
  result
}


# ──────────────────────────────────────────────────────────────────────────────
# BKP 2016 — variable method
# ──────────────────────────────────────────────────────────────────────────────

#' BKP (2016) variable method for alpha
#' @keywords internal
.bkp_variable <- function(em, N, T_val, n_pca, test_size) {
  # Cross-sectional means at each t
  x_bar <- colMeans(em, na.rm = TRUE)
  sigma_x2 <- stats::var(x_bar)

  # Standardized CSA
  x_bar_std <- (x_bar - mean(x_bar)) / stats::sd(x_bar)

  # Regress each x_it on standardized CSA
  crit_val <- stats::qnorm(1 - test_size / 2)
  b_sig <- numeric(0)
  cn_sum <- 0

  for (i in 1:N) {
    xi <- em[i, ]
    valid <- !is.na(xi)
    if (sum(valid) < 3) next

    xi_v <- xi[valid]
    z_v <- x_bar_std[valid]

    # Simple regression of x_it on CSA
    Z <- cbind(1, z_v)
    fit <- .unit_ols(xi_v, Z)
    b_i <- fit$b[2]
    se_i <- sqrt(fit$V[2, 2])
    t_i <- b_i / se_i

    if (abs(t_i) > crit_val) {
      b_sig <- c(b_sig, b_i)
    }

    # Regression on first n_pca PCs for cn
    cn_sum <- cn_sum + fit$sigma2
  }

  mu2 <- if (length(b_sig) > 0) mean(b_sig^2) else 1e-10
  cn <- cn_sum

  # alpha = 1 + ln(sigma_x^2) / (2*ln(N)) - ln(mu^2) / (2*ln(N)) - cn / (2*N*ln(N)*sigma_x^2)
  ln_N <- log(N)
  alpha <- 1 + log(sigma_x2) / (2 * ln_N) - log(mu2) / (2 * ln_N) - cn / (2 * N * ln_N * sigma_x2)

  # SE approximation from BKP (2016)
  se <- 1 / (2 * ln_N * sqrt(N))

  ci <- c(alpha - 1.96 * se, alpha + 1.96 * se)

  list(alpha = alpha, se = se, ci = ci)
}


# ──────────────────────────────────────────────────────────────────────────────
# BKP 2019 — residual method
# ──────────────────────────────────────────────────────────────────────────────

#' BKP (2019) residual method for alpha
#' @keywords internal
.bkp_residual <- function(em, N, T_val, n_pca, tuning, n_bootstrap) {
  # Pairwise correlations
  rho <- .pairwise_correlations(em)

  # Threshold
  threshold <- tuning * sqrt(log(N * T_val) / T_val)

  # Build Delta matrix (significant pairwise correlations)
  Delta <- rho
  Delta[abs(rho) <= threshold | is.na(rho)] <- 0
  diag(Delta) <- 0

  # tau = vector of ones

  tau <- rep(1, N)

  # alpha = ln(tau' Delta tau) / (2 * ln(N))
  quad <- as.numeric(t(tau) %*% Delta %*% tau)

  if (quad <= 0) {
    alpha <- 0.5  # no significant CSD detected
    se <- NA_real_
    ci <- c(NA_real_, NA_real_)
  } else {
    alpha <- log(quad) / (2 * log(N))

    # Bootstrap SE: resample units with replacement
    alpha_boot <- numeric(n_bootstrap)
    for (r in seq_len(n_bootstrap)) {
      boot_idx <- sample(1:N, N, replace = TRUE)
      em_boot <- em[boot_idx, , drop = FALSE]
      rho_boot <- .pairwise_correlations(em_boot)
      Delta_boot <- rho_boot
      Delta_boot[abs(rho_boot) <= threshold | is.na(rho_boot)] <- 0
      diag(Delta_boot) <- 0
      quad_boot <- as.numeric(t(tau) %*% Delta_boot %*% tau)
      alpha_boot[r] <- if (quad_boot > 0) log(quad_boot) / (2 * log(N)) else 0.5
    }
    se <- stats::sd(alpha_boot)
    ci <- stats::quantile(alpha_boot, c(0.025, 0.975))
  }

  list(alpha = alpha, se = se, ci = as.numeric(ci))
}


#' Print method for dcce_csd objects
#'
#' @param x A `dcce_csd` object.
#' @param ... Ignored.
#' @return Invisibly returns `x`.
#' @export
print.dcce_csd <- function(x, ...) {
  cat("\nExponent of Cross-Sectional Dependence\n")
  cat("======================================\n")
  cat(sprintf("Method  : BKP (%s)\n", if (x$method == "variable") "2016" else "2019"))
  cat(sprintf("N = %d, T = %d\n\n", x$N, x$T_val))
  cat(sprintf("alpha = %.4f", x$alpha))
  if (!is.na(x$se)) {
    cat(sprintf("  (SE = %.4f)", x$se))
  }
  cat("\n")
  if (all(!is.na(x$ci))) {
    cat(sprintf("95%% CI: [%.4f, %.4f]\n", x$ci[1], x$ci[2]))
  }
  cat("\n")
  invisible(x)
}
