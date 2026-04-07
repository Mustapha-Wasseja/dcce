#' Pooled Mean Group (PMG) Estimator Internals
#'
#' Internal functions for the Pooled Mean Group (PMG) estimator of Pesaran,
#' Shin & Smith (1999). PMG imposes common long-run coefficients across
#' units while letting the speed of adjustment and short-run dynamics remain
#' heterogeneous.
#'
#' The model is the same ARDL\eqn{(p_y, p_x)} as CS-ARDL, but the long-run
#' coefficient vector \eqn{\theta} is pooled. This implementation uses a
#' two-step inverse-variance weighted pooling of the unit-level long-run
#' coefficients obtained from the CS-ARDL fit, a simplification of the
#' concentrated maximum-likelihood estimator of Pesaran, Shin & Smith (1999)
#' that is fast, consistent, and requires no numerical optimisation.
#'
#' @name pmg_estimator
#' @keywords internal
NULL


#' Pool unit-level long-run coefficients via inverse-variance weighting
#'
#' For each base regressor \eqn{k}, compute
#' \deqn{\hat\theta^{PMG}_k =
#'   \frac{\sum_i w_{ik} \hat\theta_{ik}}{\sum_i w_{ik}},
#'   \quad w_{ik} = 1 / \hat V(\hat\theta_{ik}),}
#' with pooled standard error \eqn{(\sum_i w_{ik})^{-1/2}}.
#'
#' @param unit_lr A list of unit-level results from \code{.csardl_unit_lr()}.
#' @return A list with \code{lr_coef}, \code{lr_se}, and the full pooled
#'   vcov (diagonal).
#' @keywords internal
.pmg_pool_lr <- function(unit_lr) {
  lr_names <- names(unit_lr[[1L]]$lr_coef)
  K <- length(lr_names)
  N <- length(unit_lr)

  Theta <- do.call(rbind, lapply(unit_lr, `[[`, "lr_coef"))
  colnames(Theta) <- lr_names

  pooled <- rep(NA_real_, K)
  pooled_se <- rep(NA_real_, K)
  names(pooled) <- lr_names
  names(pooled_se) <- lr_names

  for (k in seq_len(K)) {
    vars_k <- vapply(unit_lr, function(u) u$lr_vcov[k, k], numeric(1))
    vars_k[vars_k <= 0 | !is.finite(vars_k)] <- NA_real_

    if (all(is.na(vars_k))) next
    w <- 1 / vars_k
    w[is.na(w)] <- 0
    if (sum(w) <= 0) next

    pooled[k] <- sum(w * Theta[, k], na.rm = TRUE) / sum(w)
    pooled_se[k] <- sqrt(1 / sum(w))
  }

  list(
    lr_coef = pooled,
    lr_se   = pooled_se,
    lr_vcov = diag(pooled_se^2, nrow = K),
    Theta   = Theta
  )
}


#' Post-process a PMG fit
#'
#' Takes a fit that has already been run through \code{.csardl_postprocess()}
#' and overrides the long-run block with the pooled PMG estimates. The
#' adjustment speed remains the MG average.
#'
#' @param fit A \code{dcce_fit} object with CS-ARDL post-processing.
#' @return Augmented fit object.
#' @keywords internal
.pmg_postprocess <- function(fit) {
  if (is.null(fit$unit_lr) || length(fit$unit_lr) == 0L) return(fit)

  pooled <- .pmg_pool_lr(fit$unit_lr)
  lr_names <- names(pooled$lr_coef)

  fit$lr_coef <- pooled$lr_coef
  fit$lr_se   <- pooled$lr_se
  fit$lr_vcov <- pooled$lr_vcov
  rownames(fit$lr_vcov) <- colnames(fit$lr_vcov) <- lr_names
  fit$pmg_pooled <- TRUE
  fit
}
