#' Bootstrap Inference for DCCE Models
#'
#' Computes bootstrap standard errors and confidence intervals for
#' `dcce_fit` objects using either cross-section or wild bootstrap.
#'
#' @param object A `dcce_fit` object.
#' @param type Character: `"crosssection"` (default) or `"wild"`.
#' @param reps Integer: number of bootstrap repetitions. Default 500.
#' @param percentile Logical: compute percentile CIs? Default `TRUE`.
#' @param cfresiduals Logical: for wild bootstrap, use common-factor
#'   residuals instead of defactored residuals? Default `FALSE`.
#' @param seed Integer: random seed for reproducibility. Default `NULL`.
#'
#' @return An object of class `dcce_boot` with elements:
#'   \describe{
#'     \item{se_boot}{Bootstrap standard errors.}
#'     \item{ci_lower}{Percentile CI lower bound (if `percentile=TRUE`).}
#'     \item{ci_upper}{Percentile CI upper bound.}
#'     \item{b_boot}{B x K matrix of bootstrap coefficient draws.}
#'     \item{reps}{Number of repetitions.}
#'     \item{type}{Bootstrap type.}
#'   }
#'
#' @export
#'
#' @examples
#' set.seed(42)
#' df <- data.frame(
#'   id = rep(1:10, each = 30),
#'   t  = rep(1:30, 10),
#'   y  = rnorm(300),
#'   x  = rnorm(300)
#' )
#' fit <- dcce(df, "id", "t", y ~ x, model = "mg", cross_section_vars = NULL)
#' boot_res <- bootstrap(fit, reps = 50)
#' print(boot_res)
bootstrap <- function(object,
                      type        = c("crosssection", "wild"),
                      reps        = 500L,
                      percentile  = TRUE,
                      cfresiduals = FALSE,
                      seed        = NULL) {

  if (!inherits(object, "dcce_fit")) {
    cli::cli_abort("{.arg object} must be a {.cls dcce_fit} object.")
  }
  type <- rlang::arg_match(type)

  if (!is.null(seed)) set.seed(seed)

  coef_names <- names(object$coefficients)
  K <- length(coef_names)

  b_boot <- matrix(NA_real_, reps, K)
  colnames(b_boot) <- coef_names

  if (type == "crosssection") {
    b_boot <- .bootstrap_crosssection(object, reps, coef_names)
  } else {
    b_boot <- .bootstrap_wild(object, reps, coef_names, cfresiduals)
  }

  # Compute bootstrap SE
  se_boot <- apply(b_boot, 2, stats::sd, na.rm = TRUE)
  names(se_boot) <- coef_names

  # Percentile CIs
  ci_lower <- ci_upper <- NULL
  if (percentile) {
    ci_lower <- apply(b_boot, 2, stats::quantile, probs = 0.025, na.rm = TRUE)
    ci_upper <- apply(b_boot, 2, stats::quantile, probs = 0.975, na.rm = TRUE)
    names(ci_lower) <- names(ci_upper) <- coef_names
  }

  result <- list(
    se_boot  = se_boot,
    ci_lower = ci_lower,
    ci_upper = ci_upper,
    b_boot   = b_boot,
    reps     = reps,
    type     = type
  )
  class(result) <- "dcce_boot"
  result
}


# ──────────────────────────────────────────────────────────────────────────────
# Cross-section bootstrap
# ──────────────────────────────────────────────────────────────────────────────

#' @keywords internal
.bootstrap_crosssection <- function(object, reps, coef_names) {
  K <- length(coef_names)
  b_boot <- matrix(NA_real_, reps, K)
  colnames(b_boot) <- coef_names

  unit_names <- names(object$unit_coefs)
  N <- length(unit_names)

  coef_mat <- do.call(rbind, lapply(object$unit_coefs, function(b) b[coef_names]))

  for (r in seq_len(reps)) {
    boot_idx <- sample(1:N, N, replace = TRUE)
    boot_coefs <- coef_mat[boot_idx, , drop = FALSE]
    b_boot[r, ] <- colMeans(boot_coefs)
  }

  b_boot
}


# ──────────────────────────────────────────────────────────────────────────────
# Wild bootstrap
# ──────────────────────────────────────────────────────────────────────────────

#' @keywords internal
.bootstrap_wild <- function(object, reps, coef_names, cfresiduals) {
  K <- length(coef_names)
  b_boot <- matrix(NA_real_, reps, K)
  colnames(b_boot) <- coef_names

  panel <- object$panel
  unit_var <- attr(panel, "unit_var")
  time_var <- attr(panel, "time_var")
  y_name <- object$y_name
  x_names <- object$x_names
  csa_colnames <- object$csa_colnames
  constant <- object$include_constant
  trend <- object$unit_trend

  unit_names <- names(object$unit_results)
  N <- length(unit_names)

  for (r in seq_len(reps)) {
    # Rademacher weights: one per unit
    w <- sample(c(-1, 1), N, replace = TRUE)

    coef_list_star <- list()
    for (i in seq_along(unit_names)) {
      u <- unit_names[i]
      idx <- which(panel[[unit_var]] == u)
      yi <- panel[[y_name]][idx]
      Xi <- .build_unit_X(panel, idx, x_names, csa_colnames, constant, trend)

      complete <- stats::complete.cases(cbind(yi, Xi))
      yi <- yi[complete]
      Xi <- Xi[complete, , drop = FALSE]

      if (length(yi) <= ncol(Xi)) next

      # Get fitted values and residuals from original fit
      fit_i <- object$unit_results[[u]]
      y_hat <- as.numeric(Xi %*% fit_i$b)
      e_i <- fit_i$e

      # Wild bootstrap: y* = y_hat + w_i * e_i
      y_star <- y_hat + w[i] * e_i

      # Re-estimate
      fit_star <- .unit_ols(y_star, Xi)
      coef_list_star[[u]] <- fit_star$b[coef_names]
    }

    if (length(coef_list_star) > 0) {
      b_boot[r, ] <- .mg_aggregate(coef_list_star)
    }
  }

  b_boot
}


#' Print method for dcce_boot objects
#'
#' @param x A `dcce_boot` object.
#' @param ... Ignored.
#' @return Invisibly returns `x`.
#' @export
print.dcce_boot <- function(x, ...) {
  cat("\nBootstrap Inference\n")
  cat("===================\n")
  cat(sprintf("Type: %s, Reps: %d\n\n", x$type, x$reps))

  tab <- data.frame(
    SE = x$se_boot,
    row.names = names(x$se_boot)
  )

  if (!is.null(x$ci_lower)) {
    tab$`CI_lo` <- x$ci_lower
    tab$`CI_hi` <- x$ci_upper
  }

  print(round(tab, 4))
  cat("\n")
  invisible(x)
}
