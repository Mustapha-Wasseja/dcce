#' Automatic Diagnostic Workflow for Panel Data with CSD
#'
#' Runs the recommended pre-estimation diagnostic sequence on a panel
#' regression specification and returns a structured report with a
#' suggested `dcce()` call. The workflow performs six steps:
#'
#' \enumerate{
#'   \item Panel summary (N, T, balance).
#'   \item Pesaran (2007) CIPS panel unit root test on each variable.
#'   \item Pesaran CD test on raw residuals (pooled OLS) to check for
#'         cross-sectional dependence.
#'   \item Westerlund (2007) cointegration test, if at least one variable
#'         is non-stationary.
#'   \item Rank condition classifier (De Vos et al. 2024) for a baseline
#'         static CCE fit.
#'   \item Information criterion selection of the optimal CSA lag order.
#' }
#'
#' Based on the results, the function chooses a recommended estimator
#' from \code{c("mg", "cce", "dcce", "pmg", "csardl")} and returns a
#' printable suggested `dcce()` call.
#'
#' @param data A panel data.frame.
#' @param unit_index Character: unit identifier column.
#' @param time_index Character: time identifier column.
#' @param formula Two-sided formula (levels, not differences).
#' @param max_cr_lags Integer: maximum CSA lag order to evaluate.
#'   Default \code{NULL} uses \eqn{\lfloor T^{1/3} \rfloor}.
#' @param significance Numeric: significance level used for decisions.
#'   Default \code{0.05}.
#' @param verbose Logical: print progress as each step runs. Default
#'   \code{TRUE}.
#' @param n_bootstrap Integer: bootstrap replications for the Westerlund
#'   p-values. Default \code{0L} (asymptotic).
#' @return An object of class \code{dcce_workflow} with elements
#'   \code{panel_summary}, \code{unit_root}, \code{csd_premodel},
#'   \code{cointegration}, \code{rank_condition}, \code{optimal_cr_lags},
#'   \code{recommendation}, and \code{call}.
#'
#' @export
#' @examples
#' data(pwt8)
#' wf <- dcce_workflow(
#'   data       = pwt8,
#'   unit_index = "country",
#'   time_index = "year",
#'   formula    = log_rgdpo ~ log_hc + log_ck + log_ngd,
#'   verbose    = FALSE
#' )
#' print(wf)
dcce_workflow <- function(data,
                          unit_index,
                          time_index,
                          formula,
                          max_cr_lags   = NULL,
                          significance  = 0.05,
                          verbose       = TRUE,
                          n_bootstrap   = 0L) {
  call <- match.call()

  # Prepare panel and parse formula once to get y_name and x_names
  panel <- .make_panel(data, unit_index = unit_index, time_index = time_index)
  parsed <- .parse_dcce_formula(formula, panel, unit_index, time_index)
  y_name <- parsed$y_name
  x_names <- parsed$x_names

  if (verbose) cli::cli_h1("dcce Diagnostic Workflow")

  # ---- Step 1: Panel summary ------------------------------------------------
  bal <- .check_balance(panel)
  panel_summary <- list(
    N           = bal$N,
    T_bar       = bal$T_bar,
    T_min       = bal$T_min,
    T_max       = bal$T_max,
    balanced    = bal$is_balanced
  )
  if (verbose) {
    cli::cli_h2("Step 1: Panel Summary")
    cli::cli_inform(c(
      "N = {panel_summary$N} units",
      "T (min/bar/max) = {panel_summary$T_min}/{round(panel_summary$T_bar,1)}/{panel_summary$T_max}",
      if (panel_summary$balanced) "Panel: balanced" else "Panel: unbalanced"
    ))
  }

  # Default max_cr_lags from Chudik-Pesaran rule
  if (is.null(max_cr_lags)) {
    max_cr_lags <- max(1L, floor(panel_summary$T_bar^(1/3)))
  }
  max_cr_lags <- as.integer(max_cr_lags)

  # ---- Step 2: Unit root tests on each variable ----------------------------
  if (verbose) cli::cli_h2("Step 2: CIPS Unit Root Tests")

  all_vars <- c(y_name, x_names)
  ur_results <- list()
  any_integrated <- FALSE

  for (v in all_vars) {
    mat <- tryCatch(
      .var_to_matrix(panel, v),
      error = function(e) NULL
    )
    if (is.null(mat)) {
      ur_results[[v]] <- list(
        statistic = NA_real_, p_value = NA_real_, integrated = NA
      )
      if (verbose) cli::cli_alert_warning("{v}: could not run CIPS")
      next
    }
    ct <- tryCatch(cips_test(mat, lags = 1L), error = function(e) NULL)
    if (is.null(ct)) {
      ur_results[[v]] <- list(
        statistic = NA_real_, p_value = NA_real_, integrated = NA
      )
      next
    }
    integrated <- ct$p_value > significance
    ur_results[[v]] <- list(
      statistic  = ct$statistic,
      p_value    = ct$p_value,
      integrated = integrated
    )
    any_integrated <- any_integrated || integrated
    if (verbose) {
      status <- if (integrated) "integrated (I(1))" else "stationary (I(0))"
      cli::cli_alert_info(
        "{v}: CIPS = {round(ct$statistic, 3)}, p = {round(ct$p_value, 4)}  -> {status}"
      )
    }
  }

  unit_root <- list(
    results_by_var = ur_results,
    any_integrated = any_integrated,
    decision       = if (any_integrated) "non-stationary variables present" else "all stationary"
  )

  # ---- Step 3: CSD on pooled residuals --------------------------------------
  if (verbose) cli::cli_h2("Step 3: Cross-Sectional Dependence (CD test)")

  lm_formula <- stats::as.formula(
    paste(y_name, "~", paste(x_names, collapse = " + "))
  )
  ols_fit <- tryCatch(
    stats::lm(lm_formula, data = as.data.frame(panel)),
    error = function(e) NULL
  )
  if (!is.null(ols_fit)) {
    cd_res <- tryCatch(
      pcd_test(stats::residuals(ols_fit),
               data = panel,
               unit_index = attr(panel, "unit_var"),
               time_index = attr(panel, "time_var"),
               test = "pesaran"),
      error = function(e) NULL
    )
  } else {
    cd_res <- NULL
  }

  if (!is.null(cd_res)) {
    cd_stat <- cd_res$statistics$statistic[1L]
    cd_p    <- cd_res$statistics$p_value[1L]
    csd_present <- cd_p < significance
  } else {
    cd_stat <- NA_real_; cd_p <- NA_real_; csd_present <- NA
  }

  csd_premodel <- list(
    cd_statistic = cd_stat,
    p_value      = cd_p,
    csd_present  = csd_present,
    decision     = if (isTRUE(csd_present)) "CSD present -> CCE/DCCE needed" else "no CSD -> MG sufficient"
  )
  if (verbose) {
    cli::cli_alert_info(
      "CD = {round(cd_stat, 3)}, p = {round(cd_p, 4)}  -> {csd_premodel$decision}"
    )
  }

  # ---- Step 4: Westerlund cointegration (if any integrated variables) ------
  cointegration <- NULL
  if (isTRUE(unit_root$any_integrated)) {
    if (verbose) cli::cli_h2("Step 4: Westerlund Cointegration Test")
    ct_res <- tryCatch(
      cointegration_test(
        data = panel,
        unit_index = attr(panel, "unit_var"),
        time_index = attr(panel, "time_var"),
        formula = formula,
        lags = 1L,
        test = c("ga", "gt"),
        n_bootstrap = n_bootstrap
      ),
      error = function(e) NULL
    )
    if (!is.null(ct_res)) {
      gt_row <- ct_res$statistics[ct_res$statistics$test == "Gt", ]
      ga_row <- ct_res$statistics[ct_res$statistics$test == "Ga", ]
      cointegrated <- (gt_row$p_value < significance) &&
                       (ga_row$p_value < significance)
      cointegration <- list(
        gt_statistic = gt_row$statistic,
        gt_p_value   = gt_row$p_value,
        ga_statistic = ga_row$statistic,
        ga_p_value   = ga_row$p_value,
        cointegrated = cointegrated,
        decision     = if (cointegrated) "cointegrated -> PMG / CS-ARDL recommended"
                       else "not cointegrated -> CCE in differences"
      )
      if (verbose) {
        cli::cli_alert_info(
          "Gt = {round(gt_row$statistic, 3)}, p = {round(gt_row$p_value, 4)}"
        )
        cli::cli_alert_info(
          "Ga = {round(ga_row$statistic, 3)}, p = {round(ga_row$p_value, 4)}"
        )
      }
    }
  }

  # ---- Step 5: Rank condition check ----------------------------------------
  if (verbose) cli::cli_h2("Step 5: Rank Condition Classifier")
  rank_res <- tryCatch({
    baseline_cce <- dcce(
      data = panel,
      unit_index = attr(panel, "unit_var"),
      time_index = attr(panel, "time_var"),
      formula = formula,
      model = "cce",
      cross_section_vars = ~ .,
      cross_section_lags = 0
    )
    rank_condition(baseline_cce)
  }, error = function(e) NULL)
  rank_condition_res <- if (!is.null(rank_res)) {
    list(RC = rank_res$RC, n_factors = rank_res$m, g = rank_res$g)
  } else {
    list(RC = NA_integer_, n_factors = NA_integer_, g = NA_integer_)
  }
  if (verbose && !is.null(rank_res)) {
    msg <- if (rank_condition_res$RC == 1L)
      "rank condition HOLDS" else "rank condition FAILS"
    cli::cli_alert_info(
      "Factors = {rank_condition_res$n_factors}, g = {rank_condition_res$g} -> {msg}"
    )
  }

  # ---- Step 6: IC for optimal CSA lags -------------------------------------
  if (verbose) cli::cli_h2("Step 6: CSA Lag Selection via IC")

  lags_try <- 0L:max_cr_lags
  models_ic <- lapply(lags_try, function(p) {
    tryCatch(
      dcce(
        data = panel,
        unit_index = attr(panel, "unit_var"),
        time_index = attr(panel, "time_var"),
        formula = formula,
        model = "cce",
        cross_section_vars = ~ .,
        cross_section_lags = p
      ),
      error = function(e) NULL
    )
  })
  valid <- !vapply(models_ic, is.null, logical(1))
  if (any(valid)) {
    models_ic_v <- models_ic[valid]
    lags_v <- lags_try[valid]
    ic_vals <- vapply(models_ic_v, function(m) {
      res <- ic(m, models = models_ic_v)
      c(IC1 = res$IC1, IC2 = res$IC2)
    }, numeric(2))
    ic_table <- data.frame(
      lags = lags_v,
      IC1  = ic_vals["IC1", ],
      IC2  = ic_vals["IC2", ]
    )
    best_lag <- lags_v[which.min(ic_table$IC1)]
  } else {
    ic_table <- data.frame(lags = integer(0), IC1 = numeric(0),
                           IC2 = numeric(0))
    best_lag <- 0L
  }
  optimal_cr_lags <- list(
    ic_table = ic_table,
    recommended_lags = as.integer(best_lag)
  )
  if (verbose) {
    cli::cli_alert_info("IC1 recommends cross_section_lags = {best_lag}")
  }

  # ---- Recommendation ------------------------------------------------------
  # Decision tree:
  #   - No CSD + all stationary -> "mg"
  #   - CSD + all stationary -> "cce" (lag 0) or "dcce"
  #   - CSD + integrated + cointegrated -> "pmg" (preferred) or "csardl"
  #   - CSD + integrated + not cointegrated -> "cce" in differences
  if (isTRUE(csd_premodel$csd_present)) {
    if (isTRUE(unit_root$any_integrated) &&
        !is.null(cointegration) && isTRUE(cointegration$cointegrated)) {
      rec_model <- "pmg"
    } else if (isTRUE(unit_root$any_integrated)) {
      rec_model <- "cce"  # use differences
    } else {
      rec_model <- "dcce"
    }
  } else {
    rec_model <- "mg"
  }

  rec_lags <- if (rec_model %in% c("dcce", "pmg", "csardl")) {
    max(best_lag, 1L)
  } else if (rec_model == "cce") {
    0L
  } else {
    0L
  }

  suggested_call <- .build_suggested_call(
    rec_model     = rec_model,
    unit_index    = unit_index,
    time_index    = time_index,
    y_name        = y_name,
    x_names       = x_names,
    rec_lags      = rec_lags,
    data_name     = as.character(substitute(data))
  )

  recommendation <- list(
    model              = rec_model,
    cross_section_lags = as.integer(rec_lags),
    suggested_call     = suggested_call
  )

  if (verbose) {
    cli::cli_h2("Recommendation")
    cli::cli_alert_success("Suggested model: {rec_model}")
    cli::cli_code(suggested_call)
  }

  out <- list(
    panel_summary     = panel_summary,
    unit_root         = unit_root,
    csd_premodel      = csd_premodel,
    cointegration     = cointegration,
    rank_condition    = rank_condition_res,
    optimal_cr_lags   = optimal_cr_lags,
    recommendation    = recommendation,
    call              = call
  )
  class(out) <- "dcce_workflow"
  out
}


#' Internal: build a suggested dcce() call string
#' @keywords internal
.build_suggested_call <- function(rec_model, unit_index, time_index,
                                   y_name, x_names, rec_lags, data_name) {
  rhs <- paste(x_names, collapse = " + ")
  form_str <- paste0(y_name, " ~ ", rhs)
  if (rec_model == "pmg") {
    # PMG needs explicit ARDL lags in the formula
    form_str <- paste0(y_name, " ~ L(", y_name, ", 1)")
    for (xv in x_names) {
      form_str <- paste0(form_str, " + ", xv, " + L(", xv, ", 1)")
    }
  }
  csa_part <- if (rec_model %in% c("mg")) {
    'cross_section_vars = NULL'
  } else {
    paste0('cross_section_vars = ~ ',
           paste(c(y_name, x_names), collapse = " + "))
  }
  lag_part <- if (rec_model %in% c("cce", "dcce", "pmg", "csardl")) {
    paste0(",\n  cross_section_lags = ", rec_lags)
  } else {
    ""
  }
  sprintf(
    paste0("dcce(\n",
           "  data       = %s,\n",
           "  unit_index = \"%s\",\n",
           "  time_index = \"%s\",\n",
           "  formula    = %s,\n",
           "  model      = \"%s\",\n",
           "  %s%s\n",
           ")"),
    data_name, unit_index, time_index, form_str, rec_model, csa_part, lag_part
  )
}


#' Internal: convert a single variable column of a panel to a (N x T) matrix
#' @keywords internal
.var_to_matrix <- function(panel, var) {
  unit_var <- attr(panel, "unit_var")
  time_var <- attr(panel, "time_var")
  units <- unique(panel[[unit_var]])
  times <- sort(unique(panel[[time_var]]))
  mat <- matrix(NA_real_, nrow = length(units), ncol = length(times))
  for (i in seq_along(units)) {
    idx <- which(panel[[unit_var]] == units[i])
    tpos <- match(panel[[time_var]][idx], times)
    mat[i, tpos] <- panel[[var]][idx]
  }
  # Drop columns (time periods) with all NA
  keep_cols <- colSums(is.na(mat)) < nrow(mat)
  mat <- mat[, keep_cols, drop = FALSE]
  # Keep only rows with enough data
  keep_rows <- rowSums(is.na(mat)) < ncol(mat)
  mat <- mat[keep_rows, , drop = FALSE]
  mat
}


#' Print a dcce_workflow object
#'
#' @param x A \code{dcce_workflow} object.
#' @param ... Ignored.
#' @return Invisibly returns \code{x}.
#' @export
print.dcce_workflow <- function(x, ...) {
  .wf_header <- function(txt) cat("\n== ", txt, " ==\n", sep = "")

  cat("\n===== dcce Diagnostic Workflow =====\n")

  .wf_header("Panel Summary")
  cat(sprintf("  N = %d units\n", x$panel_summary$N))
  cat(sprintf("  T (min/bar/max) = %d/%.1f/%d\n",
              x$panel_summary$T_min,
              x$panel_summary$T_bar,
              x$panel_summary$T_max))
  cat(sprintf("  Panel: %s\n",
              if (isTRUE(x$panel_summary$balanced)) "balanced" else "unbalanced"))

  .wf_header("Unit Root Tests (CIPS)")
  for (nm in names(x$unit_root$results_by_var)) {
    r <- x$unit_root$results_by_var[[nm]]
    status <- if (is.na(r$integrated)) "n/a" else
              if (r$integrated) "I(1)" else "I(0)"
    cat(sprintf("  %-20s  stat = %8.3f  p = %6.4f  (%s)\n",
                nm, r$statistic, r$p_value, status))
  }

  .wf_header("Cross-Sectional Dependence")
  cat(sprintf("  CD = %.3f, p = %.4f -> %s\n",
              x$csd_premodel$cd_statistic,
              x$csd_premodel$p_value,
              x$csd_premodel$decision))

  if (!is.null(x$cointegration)) {
    .wf_header("Westerlund Cointegration")
    cat(sprintf("  Gt = %.3f, p = %.4f\n",
                x$cointegration$gt_statistic,
                x$cointegration$gt_p_value))
    cat(sprintf("  Ga = %.3f, p = %.4f\n",
                x$cointegration$ga_statistic,
                x$cointegration$ga_p_value))
    cat(sprintf("  -> %s\n", x$cointegration$decision))
  }

  .wf_header("Rank Condition")
  rc <- x$rank_condition
  if (!is.na(rc$RC)) {
    msg <- if (rc$RC == 1L) "holds" else "FAILS"
    cat(sprintf("  m = %d, g = %d -> %s\n", rc$n_factors, rc$g, msg))
  } else {
    cat("  (could not be computed)\n")
  }

  .wf_header("Optimal CSA Lag")
  if (nrow(x$optimal_cr_lags$ic_table) > 0L) {
    cat(sprintf("  IC1 selects cross_section_lags = %d\n",
                x$optimal_cr_lags$recommended_lags))
  }

  .wf_header("Recommendation")
  cat(sprintf("  Suggested model: %s\n", x$recommendation$model))
  cat("\nSuggested dcce() call:\n\n")
  cat(x$recommendation$suggested_call, "\n")
  invisible(x)
}
