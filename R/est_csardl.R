#' CS-ARDL Estimator Internals
#'
#' Internal functions for the Cross-Sectionally augmented ARDL (CS-ARDL)
#' estimator of Chudik, Mohaddes, Pesaran & Raissi (2016). Estimates an
#' ARDL\eqn{(p_y, p_x)} model with cross-sectional averages and recovers
#' long-run coefficients and the speed of adjustment via the delta method.
#'
#' The unit-level regression is
#' \deqn{y_{it} = \alpha_i + \sum_{p=1}^{P_y} \phi_{ip} y_{i,t-p}
#'                + \sum_{q=0}^{P_x} \beta'_{iq} x_{i,t-q}
#'                + \delta'_{i} \bar{z}_t + e_{it}.}
#'
#' The long-run coefficient on \eqn{x_k} is
#' \deqn{\theta_{ik} = \frac{\sum_{q=0}^{P_x} \beta_{ikq}}{1 - \sum_{p=1}^{P_y} \phi_{ip}}}
#' and the implied speed of adjustment is
#' \deqn{\varphi_i = -\left(1 - \sum_{p=1}^{P_y} \phi_{ip}\right).}
#'
#' @name csardl_estimator
#' @keywords internal
NULL


#' Classify ARDL regressor terms
#'
#' Given the dependent variable name, the list of regressor names, and the
#' panel, classify each regressor term into:
#'   - y-lag (e.g. \code{L(y,1)})
#'   - contemporaneous or lagged x (grouped by base variable name)
#'
#' @param y_name Character: dependent variable name.
#' @param x_names Character vector: regressor names from the parsed formula.
#' @return A list with
#'   \describe{
#'     \item{y_lag_terms}{Character vector of y-lag term names.}
#'     \item{y_lag_orders}{Integer vector of the lag orders matching
#'       \code{y_lag_terms}.}
#'     \item{x_groups}{Named list: each element is a character vector of
#'       terms (level + lags) belonging to a single base regressor.}
#'     \item{x_group_lags}{Named list: each element is an integer vector of
#'       lag orders matching \code{x_groups}.}
#'   }
#' @keywords internal
.csardl_classify_terms <- function(y_name, x_names) {
  # Pattern to match L(name, k)
  lag_pattern <- "^L\\(([^,]+),\\s*(-?\\d+)\\)$"

  y_lag_terms <- character(0)
  y_lag_orders <- integer(0)

  x_groups <- list()
  x_group_lags <- list()

  for (term in x_names) {
    m <- regmatches(term, regexec(lag_pattern, term))[[1L]]
    if (length(m) == 3L) {
      base <- m[2L]
      k <- as.integer(m[3L])
      if (base == y_name) {
        y_lag_terms <- c(y_lag_terms, term)
        y_lag_orders <- c(y_lag_orders, k)
      } else {
        if (is.null(x_groups[[base]])) {
          x_groups[[base]] <- character(0)
          x_group_lags[[base]] <- integer(0)
        }
        x_groups[[base]] <- c(x_groups[[base]], term)
        x_group_lags[[base]] <- c(x_group_lags[[base]], k)
      }
    } else {
      # Plain variable name (lag 0)
      base <- term
      if (is.null(x_groups[[base]])) {
        x_groups[[base]] <- character(0)
        x_group_lags[[base]] <- integer(0)
      }
      x_groups[[base]] <- c(x_groups[[base]], term)
      x_group_lags[[base]] <- c(x_group_lags[[base]], 0L)
    }
  }

  list(
    y_lag_terms  = y_lag_terms,
    y_lag_orders = y_lag_orders,
    x_groups     = x_groups,
    x_group_lags = x_group_lags
  )
}


#' Recover unit-level long-run coefficients and adjustment speed
#'
#' For a single unit, given the ARDL coefficient vector \code{b} and its
#' variance \code{V}, compute the long-run coefficient on each base regressor
#' and the speed of adjustment, with delta-method standard errors.
#'
#' @param b Numeric vector: unit-level coefficients from the ARDL regression
#'   (must have names matching the regressor columns).
#' @param V Numeric matrix: unit-level variance-covariance matrix.
#' @param classify List returned by \code{.csardl_classify_terms()}.
#' @return A list with
#'   \describe{
#'     \item{lr_coef}{Named numeric vector of long-run coefficients.}
#'     \item{lr_vcov}{Variance-covariance matrix of \code{lr_coef}.}
#'     \item{phi}{Adjustment speed \eqn{\varphi_i = -(1 - \sum \phi_p)}.}
#'     \item{phi_se}{Delta-method SE of \code{phi}.}
#'   }
#' @keywords internal
.csardl_unit_lr <- function(b, V, classify) {
  y_lag_terms <- classify$y_lag_terms
  x_groups <- classify$x_groups

  # Denominator: d = 1 - sum(phi)
  if (length(y_lag_terms) > 0L) {
    phi_sum <- sum(b[y_lag_terms])
  } else {
    phi_sum <- 0
  }
  d_val <- 1 - phi_sum

  # Long-run coefficient for each base regressor
  lr_coef <- numeric(length(x_groups))
  names(lr_coef) <- names(x_groups)
  for (k in names(x_groups)) {
    terms_k <- x_groups[[k]]
    lr_coef[k] <- sum(b[terms_k]) / d_val
  }

  # Joint delta method for [theta_1, ..., theta_K, phi]
  # Stack: full transformation g(b)
  g <- function(bvec) {
    if (length(y_lag_terms) > 0L) {
      phi_sum_local <- sum(bvec[y_lag_terms])
    } else {
      phi_sum_local <- 0
    }
    d_local <- 1 - phi_sum_local
    out <- numeric(length(x_groups) + 1L)
    for (j in seq_along(x_groups)) {
      out[j] <- sum(bvec[x_groups[[j]]]) / d_local
    }
    out[length(x_groups) + 1L] <- -d_local  # phi = -(1 - sum phi_p)
    out
  }

  # Subvector of b to feed delta method
  rel_names <- unique(c(y_lag_terms, unlist(x_groups)))
  b_sub <- b[rel_names]
  V_sub <- V[rel_names, rel_names, drop = FALSE]

  # Wrap g to operate only on rel_names
  g_sub <- function(bvec_sub) {
    full <- setNames(numeric(length(rel_names)), rel_names)
    full[names(bvec_sub)] <- bvec_sub
    g(full)
  }

  dm <- .delta_method(g_sub, b_sub, V_sub)

  K <- length(x_groups)
  lr_vcov <- dm$vcov[1:K, 1:K, drop = FALSE]
  rownames(lr_vcov) <- colnames(lr_vcov) <- names(x_groups)

  phi_val <- dm$estimate[K + 1L]
  phi_var <- dm$vcov[K + 1L, K + 1L]
  phi_se <- sqrt(max(phi_var, 0))

  list(
    lr_coef = lr_coef,
    lr_vcov = lr_vcov,
    phi     = phi_val,
    phi_se  = phi_se
  )
}


#' Aggregate unit-level long-run and adjustment results to the MG level
#'
#' @param unit_lr A list of unit-level results from \code{.csardl_unit_lr()}.
#' @return A list with MG long-run coefficients, MG variance, MG adjustment
#'   speed, and an inverse-variance weighted pooled long-run estimate (used
#'   by PMG).
#' @keywords internal
.csardl_mg_lr <- function(unit_lr) {
  N <- length(unit_lr)
  lr_names <- names(unit_lr[[1L]]$lr_coef)
  K <- length(lr_names)

  # Stack unit-level LR coefs
  Theta <- do.call(rbind, lapply(unit_lr, `[[`, "lr_coef"))
  colnames(Theta) <- lr_names

  # MG long-run
  lr_mg <- colMeans(Theta, na.rm = TRUE)

  # MG variance
  Theta_c <- sweep(Theta, 2, lr_mg)
  lr_vcov_mg <- crossprod(Theta_c) / (N * (N - 1L))
  rownames(lr_vcov_mg) <- colnames(lr_vcov_mg) <- lr_names

  # Adjustment speed (MG)
  phi_unit <- vapply(unit_lr, `[[`, numeric(1), "phi")
  phi_mg <- mean(phi_unit, na.rm = TRUE)
  phi_se_mg <- sqrt(var(phi_unit, na.rm = TRUE) / N)

  # Inverse-variance pooled LR (for PMG)
  pooled_lr <- rep(NA_real_, K)
  pooled_lr_se <- rep(NA_real_, K)
  names(pooled_lr) <- lr_names
  names(pooled_lr_se) <- lr_names

  for (k in seq_len(K)) {
    vars_k <- vapply(unit_lr, function(u) u$lr_vcov[k, k], numeric(1))
    vars_k[vars_k <= 0 | !is.finite(vars_k)] <- NA_real_
    if (any(!is.na(vars_k))) {
      w <- 1 / vars_k
      w[is.na(w)] <- 0
      if (sum(w) > 0) {
        pooled_lr[k] <- sum(w * Theta[, k], na.rm = TRUE) / sum(w)
        pooled_lr_se[k] <- sqrt(1 / sum(w))
      }
    }
  }

  list(
    lr_coef     = lr_mg,
    lr_vcov     = lr_vcov_mg,
    lr_se       = sqrt(diag(lr_vcov_mg)),
    phi         = phi_mg,
    phi_se      = phi_se_mg,
    pooled_lr   = pooled_lr,
    pooled_lr_se = pooled_lr_se,
    Theta       = Theta
  )
}


#' Post-process a CS-ARDL fit to attach LR / adjustment information
#'
#' Called from \code{dcce()} after the unit-level OLS and MG aggregation are
#' done. Classifies terms, computes unit-level LR and adjustment, aggregates
#' to MG, and returns the augmented fit object.
#'
#' @param fit A \code{dcce_fit} object from the main dcce() pipeline.
#' @return The fit object augmented with \code{lr_coef}, \code{lr_vcov},
#'   \code{lr_se}, \code{adjustment}, \code{adjustment_se}, and related fields.
#' @keywords internal
.csardl_postprocess <- function(fit) {
  y_name <- fit$y_name
  x_names <- fit$x_names

  classify <- .csardl_classify_terms(y_name, x_names)

  if (length(classify$x_groups) == 0L) {
    cli::cli_warn("CS-ARDL: no non-y regressors found. Returning MG fit only.")
    return(fit)
  }

  # Unit-level LR + adjustment
  unit_lr <- list()
  for (u in names(fit$unit_results)) {
    b_u <- fit$unit_results[[u]]$b
    V_u <- fit$unit_results[[u]]$V
    # Drop CSA/intercept rows not relevant here — the delta method only uses
    # names in y_lag_terms and x_groups.
    res <- tryCatch(
      .csardl_unit_lr(b_u, V_u, classify),
      error = function(e) NULL
    )
    if (!is.null(res)) unit_lr[[u]] <- res
  }

  if (length(unit_lr) == 0L) {
    cli::cli_warn("CS-ARDL: unit-level LR recovery failed for all units.")
    return(fit)
  }

  mg_res <- .csardl_mg_lr(unit_lr)

  fit$ardl_classify   <- classify
  fit$unit_lr         <- unit_lr
  fit$lr_coef         <- mg_res$lr_coef
  fit$lr_vcov         <- mg_res$lr_vcov
  fit$lr_se           <- mg_res$lr_se
  fit$adjustment      <- mg_res$phi
  fit$adjustment_se   <- mg_res$phi_se
  fit$pooled_lr       <- mg_res$pooled_lr
  fit$pooled_lr_se    <- mg_res$pooled_lr_se
  fit$unit_lr_matrix  <- mg_res$Theta

  fit
}
