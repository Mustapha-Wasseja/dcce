#' S3 Methods for dcce_fit Objects
#'
#' @name dcce_fit_methods
#' @keywords internal
NULL

#' Extract coefficients from a dcce_fit object
#'
#' @param object A `dcce_fit` object.
#' @param type Character: `"mg"` for Mean Group coefficients (default),
#'   `"unit"` for unit-level coefficients as a tibble.
#' @param ... Ignored.
#' @return A named numeric vector (for `"mg"`) or a tibble (for `"unit"`).
#' @export
coef.dcce_fit <- function(object, type = c("mg", "unit"), ...) {
  type <- rlang::arg_match(type)
  if (type == "mg") {
    return(object$coefficients)
  }
  # Unit-level coefficients as tibble
  B <- do.call(rbind, object$unit_coefs)
  unit_names <- rep(rownames(B), each = ncol(B))
  term_names <- rep(colnames(B), times = nrow(B))
  estimates <- as.numeric(t(B))
  tibble::tibble(
    unit     = unit_names,
    term     = term_names,
    estimate = estimates
  )
}

#' Extract variance-covariance matrix from a dcce_fit object
#'
#' @param object A `dcce_fit` object.
#' @param ... Ignored.
#' @return A variance-covariance matrix.
#' @export
vcov.dcce_fit <- function(object, ...) {
  object$vcov
}

#' Extract residuals from a dcce_fit object
#'
#' @param object A `dcce_fit` object.
#' @param ... Ignored.
#' @return A numeric vector.
#' @export
residuals.dcce_fit <- function(object, ...) {
  object$residuals
}

#' Print a dcce_fit object
#'
#' @param x A `dcce_fit` object.
#' @param ... Ignored.
#' @return Invisibly returns `x`.
#' @export
print.dcce_fit <- function(x, ...) {
  model_label <- switch(x$model,
    mg     = "Mean Group (MG)",
    cce    = "CCE (Mean Group)",
    dcce   = "DCCE (Mean Group)",
    pmg    = "Pooled Mean Group (PMG)",
    csdl   = "CS-DL (Long-run)",
    csardl = "CS-ARDL (Long-run)",
    rcce   = "Regularized CCE"
  )

  cat("\nDynamic Common Correlated Effects Estimation\n")
  cat("============================================\n")
  cat(sprintf("Estimator    : %s\n", model_label))
  cat(sprintf("Dep. variable: %s\n", x$y_name))
  cat(sprintf("Groups (N)   : %-10d Time periods (T): %d\n",
              x$N, round(x$T_bar)))
  cat(sprintf("Observations : %-10d", x$n_obs))
  if (!is.null(x$csa_colnames)) {
    cat(sprintf(" CSA lags        : %s",
                paste(x$cross_section_lags, collapse = ",")))
  }
  cat("\n\n")

  # Coefficient table
  b <- x$coefficients
  se <- x$se
  z <- b / se
  p <- 2 * stats::pnorm(-abs(z))
  ci_lo <- b - 1.96 * se
  ci_hi <- b + 1.96 * se

  tab <- cbind(
    Coef.     = b,
    Std.Err.  = se,
    z         = z,
    `P>|z|`   = p,
    `[95% lo` = ci_lo,
    `95% hi]` = ci_hi
  )

  cat(sprintf("%-25s %9s %9s %9s %9s %9s %9s\n",
              "", "Coef.", "Std. Err.", "z", "P>|z|", "[95% CI", "]"))
  for (i in seq_len(nrow(tab))) {
    cat(sprintf("%-25s %9.4f %9.4f %9.2f %9.3f %9.3f %9.3f\n",
                rownames(tab)[i],
                tab[i, 1], tab[i, 2], tab[i, 3], tab[i, 4],
                tab[i, 5], tab[i, 6]))
  }

  cat(sprintf("\nR-sq (MG): %.3f   RMSE: %.4f\n", x$r2, x$rmse))

  invisible(x)
}

#' Summary for a dcce_fit object
#'
#' @param object A `dcce_fit` object.
#' @param ... Ignored.
#' @return Invisibly returns `object`.
#' @export
summary.dcce_fit <- function(object, ...) {
  print(object)
  invisible(object)
}

#' Extract fitted values from a dcce_fit object
#'
#' @param object A `dcce_fit` object.
#' @param ... Ignored.
#' @return A numeric vector of fitted values (y-hat), or `NULL` if not
#'   available.
#' @export
fitted.dcce_fit <- function(object, ...) {
  object$fitted_values
}

#' Tidy a dcce_fit object
#'
#' Returns a tibble of MG coefficients with standard errors, test statistics,
#' p-values, and confidence intervals, compatible with the broom package.
#'
#' @param x A `dcce_fit` object.
#' @param ... Ignored.
#' @return A tibble with columns `term`, `estimate`, `std.error`,
#'   `statistic`, `p.value`, `conf.low`, `conf.high`, `type`.
#' @export
tidy.dcce_fit <- function(x, ...) {
  b <- x$coefficients
  se <- x$se
  z <- b / se
  p <- 2 * stats::pnorm(-abs(z))

  tibble::tibble(
    term      = names(b),
    estimate  = unname(b),
    std.error = unname(se),
    statistic = unname(z),
    p.value   = unname(p),
    conf.low  = unname(b - 1.96 * se),
    conf.high = unname(b + 1.96 * se),
    type      = "mg"
  )
}

#' Glance at a dcce_fit object
#'
#' Returns a single-row tibble of model summary statistics, compatible with
#' the broom package.
#'
#' @param x A `dcce_fit` object.
#' @param ... Ignored.
#' @return A single-row tibble.
#' @export
glance.dcce_fit <- function(x, ...) {
  tibble::tibble(
    nobs          = as.integer(x$n_obs),
    n_units       = as.integer(x$N),
    t_bar         = x$T_bar,
    r.squared     = x$r2,
    adj.r.squared = NA_real_,
    rmse          = x$rmse,
    cd_statistic  = NA_real_,
    cd_p.value    = NA_real_,
    estimator     = toupper(x$model),
    cr_lags       = if (is.null(x$cross_section_lags)) NA_integer_
                    else as.integer(x$cross_section_lags[1])
  )
}

#' Predict from a dcce_fit object
#'
#' @param object A `dcce_fit` object.
#' @param ... Ignored.
#' @return A numeric vector of in-sample predictions.
#' @export
predict.dcce_fit <- function(object, ...) {
  pred <- numeric(0)
  panel <- object$panel
  unit_var <- attr(panel, "unit_var")

  for (u in names(object$unit_results)) {
    idx <- which(panel[[unit_var]] == u)
    yi <- panel[[object$y_name]][idx]
    Xi <- .build_unit_X(panel, idx, object$x_names, object$csa_colnames,
                        object$include_constant, object$unit_trend)
    complete <- stats::complete.cases(cbind(yi, Xi))
    Xi_c <- Xi[complete, , drop = FALSE]

    fit_i <- object$unit_results[[u]]
    pred <- c(pred, as.numeric(Xi_c %*% fit_i$b))
  }
  pred
}
