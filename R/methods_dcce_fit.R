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
    amg    = "Augmented Mean Group (AMG)",
    pmg    = "Pooled Mean Group (PMG)",
    csdl   = "CS-DL (Long-run)",
    csardl = "CS-ARDL (Short + Long run)",
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

  has_lr  <- !is.null(x$lr_coef) && length(x$lr_coef) > 0L
  has_adj <- !is.null(x$adjustment)

  # --- Short-run (or main) coefficient table --------------------------------
  sr_label <- if (has_lr) "Short-run (unit-level MG) Estimates" else "Mean Group Coefficients"
  cat(sr_label, "\n")
  cat(strrep("-", nchar(sr_label)), "\n", sep = "")
  .print_coef_table(x$coefficients, x$se)

  # --- Adjustment block (CS-ARDL, PMG) --------------------------------------
  if (has_adj) {
    cat("\nAdjustment Term\n")
    cat("---------------\n")
    .print_coef_table(
      c(adjustment = x$adjustment),
      c(adjustment = x$adjustment_se)
    )
  }

  # --- Long-run block -------------------------------------------------------
  if (has_lr) {
    lr_label <- if (isTRUE(x$pmg_pooled)) {
      "Long-run Estimates (Pooled, PMG)"
    } else if (x$model == "csdl") {
      "Long-run Estimates (CS-DL, direct)"
    } else {
      "Long-run Estimates (MG)"
    }
    cat("\n", lr_label, "\n", sep = "")
    cat(strrep("-", nchar(lr_label)), "\n", sep = "")
    .print_coef_table(x$lr_coef, x$lr_se)
  }

  cat(sprintf("\nR-sq (MG): %.3f   RMSE: %.4f\n", x$r2, x$rmse))
  invisible(x)
}


#' Internal: print a single coefficient table
#' @keywords internal
.print_coef_table <- function(b, se) {
  z <- b / se
  p <- 2 * stats::pnorm(-abs(z))
  ci_lo <- b - 1.96 * se
  ci_hi <- b + 1.96 * se

  signif_codes <- ifelse(p < 0.001, "***",
                  ifelse(p < 0.01,  "**",
                  ifelse(p < 0.05,  "*",
                  ifelse(p < 0.1,   ".", ""))))

  cat(sprintf("%-22s %10s %10s %8s %9s %9s %9s %5s\n",
              "", "Coef.", "Std. Err.", "z", "P>|z|", "[95% CI", "]", "Sig"))
  for (i in seq_along(b)) {
    cat(sprintf("%-22s %10.4f %10.4f %8.2f %9.4f %9.4f %9.4f %5s\n",
                names(b)[i],
                b[i], se[i], z[i], p[i],
                ci_lo[i], ci_hi[i], signif_codes[i]))
  }
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
#' For long-run estimators (CS-ARDL, CS-DL, PMG) the tibble additionally
#' includes rows for long-run coefficients and the adjustment speed.
#'
#' @param x A `dcce_fit` object.
#' @param include_lr Logical: include long-run and adjustment rows when
#'   available? Default `TRUE`.
#' @param ... Ignored.
#' @return A tibble with columns `term`, `estimate`, `std.error`,
#'   `statistic`, `p.value`, `conf.low`, `conf.high`, `type`.
#' @export
tidy.dcce_fit <- function(x, include_lr = TRUE, ...) {
  make_row <- function(b, se, type_label) {
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
      type      = type_label
    )
  }

  out <- make_row(x$coefficients, x$se, "mg")

  if (isTRUE(include_lr)) {
    if (!is.null(x$adjustment)) {
      out <- rbind(out, make_row(
        c(adjustment = x$adjustment),
        c(adjustment = x$adjustment_se),
        "adjustment"
      ))
    }
    if (!is.null(x$lr_coef) && length(x$lr_coef) > 0L) {
      lr_type <- if (isTRUE(x$pmg_pooled)) "lr_pooled" else "lr"
      out <- rbind(out, make_row(x$lr_coef, x$lr_se, lr_type))
    }
  }

  out
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

#' Confidence intervals for a dcce_fit object
#'
#' @param object A `dcce_fit` object.
#' @param parm Character vector of parameter names (default: all).
#' @param level Confidence level. Default `0.95`.
#' @param type Character: `"mg"` (default) for the main MG coefficients,
#'   `"lr"` for long-run coefficients (CS-ARDL, CS-DL, PMG),
#'   `"adjustment"` for the speed of adjustment (CS-ARDL, PMG).
#' @param ... Ignored.
#' @return A matrix of confidence intervals with rows corresponding to
#'   parameters.
#' @export
confint.dcce_fit <- function(object, parm = NULL, level = 0.95,
                             type = c("mg", "lr", "adjustment"), ...) {
  type <- rlang::arg_match(type)
  z <- stats::qnorm(1 - (1 - level) / 2)

  if (type == "mg") {
    b <- object$coefficients
    se <- object$se
  } else if (type == "lr") {
    if (is.null(object$lr_coef)) {
      cli::cli_abort("No long-run coefficients found. Use {.code type = \"mg\"}.")
    }
    b <- object$lr_coef
    se <- object$lr_se
  } else {
    if (is.null(object$adjustment)) {
      cli::cli_abort("No adjustment coefficient found.")
    }
    b <- c(adjustment = object$adjustment)
    se <- c(adjustment = object$adjustment_se)
  }

  out <- cbind(
    lower = b - z * se,
    upper = b + z * se
  )
  colnames(out) <- c(
    sprintf("%.1f%%", 100 * (1 - level) / 2),
    sprintf("%.1f%%", 100 * (1 + level) / 2)
  )
  if (!is.null(parm)) {
    keep <- intersect(parm, rownames(out))
    out <- out[keep, , drop = FALSE]
  }
  out
}

#' Update a dcce_fit object
#'
#' @param object A `dcce_fit` object.
#' @param formula. Optional replacement formula (use \code{.} to keep
#'   existing parts).
#' @param ... Additional arguments passed to \code{dcce()} to override
#'   the original call.
#' @param evaluate Logical: evaluate the updated call? Default `TRUE`.
#' @return An updated `dcce_fit` object (if \code{evaluate = TRUE}) or
#'   an unevaluated call.
#' @export
update.dcce_fit <- function(object, formula. = NULL, ..., evaluate = TRUE) {
  call <- object$call
  extras <- match.call(expand.dots = FALSE)$`...`
  if (!is.null(formula.)) {
    call$formula <- stats::update.formula(stats::formula(call$formula), formula.)
  }
  if (length(extras) > 0L) {
    existing <- !is.na(match(names(extras), names(call)))
    for (a in names(extras)) call[[a]] <- extras[[a]]
  }
  if (evaluate) eval(call, parent.frame()) else call
}

#' Plot method for dcce_fit objects
#'
#' Produces coefficient distribution plots (one histogram per regressor)
#' showing the unit-level estimates and the MG mean.
#'
#' @param x A `dcce_fit` object.
#' @param which Character: `"coef"` (default) for coefficient histograms,
#'   or `"resid"` for a residuals-versus-time plot by unit.
#' @param ... Passed to the underlying graphics call.
#' @return Invisibly returns `x`.
#' @export
plot.dcce_fit <- function(x, which = c("coef", "resid"), ...) {
  which <- rlang::arg_match(which)

  if (which == "coef") {
    B <- do.call(rbind, x$unit_coefs)
    non_intercept <- setdiff(colnames(B), "(Intercept)")
    if (length(non_intercept) == 0L) {
      cli::cli_warn("No non-intercept coefficients to plot.")
      return(invisible(x))
    }
    old_par <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(old_par), add = TRUE)

    n_plots <- length(non_intercept)
    nr <- ceiling(sqrt(n_plots))
    nc <- ceiling(n_plots / nr)
    graphics::par(mfrow = c(nr, nc), mar = c(4, 4, 2, 1))

    for (v in non_intercept) {
      vals <- B[, v]
      graphics::hist(
        vals,
        main   = v,
        xlab   = "Unit-level estimate",
        col    = "lightsteelblue",
        border = "white",
        breaks = max(10, min(25, length(vals) / 3))
      )
      graphics::abline(v = mean(vals, na.rm = TRUE),
                       col = "firebrick", lwd = 2, lty = 2)
    }
  } else {
    # Residuals versus time
    resid_mat <- .resid_list_to_matrix(x$resid_list)
    graphics::matplot(
      t(resid_mat),
      type = "l",
      col  = grDevices::rgb(0, 0, 0, 0.2),
      lty  = 1,
      xlab = "Time index",
      ylab = "Residual",
      main = "Unit-level residuals"
    )
    graphics::abline(h = 0, col = "firebrick", lwd = 1.5)
  }

  invisible(x)
}


#' Predict from a dcce_fit object
#'
#' @param object A `dcce_fit` object.
#' @param newdata Optional data.frame with new observations. When supplied,
#'   predictions are computed using the Mean Group coefficients on the
#'   structural regressors (CSAs are not applied because they depend on the
#'   full cross-section). When \code{NULL}, in-sample unit-level fitted
#'   values are returned.
#' @param type Character: \code{"response"} (default) returns predicted y;
#'   \code{"xb"} is an alias for response (kept for compatibility with
#'   \pkg{marginaleffects}).
#' @param ... Ignored.
#' @return A numeric vector of predictions.
#' @export
predict.dcce_fit <- function(object, newdata = NULL,
                             type = c("response", "xb"), ...) {
  type <- rlang::arg_match(type)

  if (is.null(newdata)) {
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
    return(pred)
  }

  # newdata path: use MG coefficients on structural regressors only.
  # Useful for marginaleffects which needs a vectorised prediction method.
  b_mg <- object$coefficients
  cols <- names(b_mg)

  X <- matrix(0, nrow = nrow(newdata), ncol = length(cols))
  colnames(X) <- cols
  for (cn in cols) {
    if (cn == "(Intercept)") {
      X[, cn] <- 1
    } else if (cn == "(trend)") {
      X[, cn] <- seq_len(nrow(newdata))
    } else if (cn %in% names(newdata)) {
      X[, cn] <- as.numeric(newdata[[cn]])
    } else {
      cli::cli_abort("Variable {.val {cn}} not found in newdata.")
    }
  }
  as.numeric(X %*% b_mg)
}
