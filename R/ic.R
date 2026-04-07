#' Information Criteria for CSA Selection
#'
#' Implements the Margaritella & Westerlund (2023) information and panel
#' criteria for selecting the number of cross-sectional averages.
#'
#' @param object A `dcce_fit` object.
#' @param models Optional list of `dcce_fit` objects for PC criteria comparison.
#'
#' @return An object of class `dcce_ic` with `IC1`, `IC2`, `PC1`, `PC2`.
#' @export
#' @keywords internal
ic <- function(object, models = NULL) {
  if (!inherits(object, "dcce_fit")) {
    cli::cli_abort("{.arg object} must be a {.cls dcce_fit} object.")
  }

  N <- object$N
  T_bar <- round(object$T_bar)
  n_obs <- object$n_obs
  resids <- object$residuals
  S_M2 <- sum(resids^2) / n_obs

  m <- if (!is.null(object$csa_colnames)) length(object$csa_colnames) else 0L

  c_NT <- min(N, T_bar)
  pen1 <- m * (N + T_bar) / (N * T_bar) * log(N * T_bar / (N + T_bar))
  pen2 <- m * (N + T_bar) / (N * T_bar) * log(c_NT)

  IC1 <- log(S_M2) + pen1
  IC2 <- log(S_M2) + pen2

  # PC criteria need a reference model
  PC1 <- PC2 <- NA_real_
  if (!is.null(models) && length(models) > 0) {
    # Reference = model with most CSAs
    n_csa <- vapply(models, function(m) {
      if (!is.null(m$csa_colnames)) length(m$csa_colnames) else 0L
    }, integer(1))
    ref_idx <- which.max(n_csa)
    ref_resids <- models[[ref_idx]]$residuals
    S_Mb2 <- sum(ref_resids^2) / models[[ref_idx]]$n_obs

    PC1 <- S_M2 + m * S_Mb2 * (N + T_bar) / (N * T_bar) * log(N * T_bar / (N + T_bar))
    PC2 <- S_M2 + m * S_Mb2 * (N + T_bar) / (N * T_bar) * log(c_NT)
  }

  result <- list(
    IC1 = IC1, IC2 = IC2, PC1 = PC1, PC2 = PC2,
    S_M2 = S_M2, m = m, N = N, T_bar = T_bar
  )
  class(result) <- "dcce_ic"
  result
}

#' @export
print.dcce_ic <- function(x, ...) {
  cat("\nInformation Criteria for CSA Selection\n")
  cat("======================================\n")
  cat(sprintf("N = %d, T = %d, cross-section lags (m) = %d\n\n", x$N, round(x$T_bar), x$m))
  cat(sprintf("IC1 = %.4f\n", x$IC1))
  cat(sprintf("IC2 = %.4f\n", x$IC2))
  if (!is.na(x$PC1)) cat(sprintf("PC1 = %.4f\n", x$PC1))
  if (!is.na(x$PC2)) cat(sprintf("PC2 = %.4f\n", x$PC2))
  cat("\n")
  invisible(x)
}
