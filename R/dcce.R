#' Dynamic Common Correlated Effects Estimation
#'
#' Estimates heterogeneous coefficient panel data models with cross-sectional
#' dependence using Mean Group (MG), Common Correlated Effects (CCE), Dynamic
#' CCE (DCCE), and related estimators.
#'
#' @param data A data.frame containing the panel data.
#' @param unit_index Character: name of the unit (cross-section) variable.
#' @param time_index Character: name of the time variable.
#' @param formula A formula of the form `y ~ x1 + x2`. Supports `L()`, `D()`,
#'   and `Lrange()` operators for lags, differences, and lag ranges.
#' @param model Character: estimator to use. One of `"dcce"`, `"cce"`, `"mg"`,
#'   `"rcce"`, `"pmg"`, `"csdl"`, `"csardl"`. Default `"dcce"`.
#' @param cross_section_vars A one-sided formula specifying variables for
#'   cross-sectional averages, e.g. `~ x1 + x2`. Use `~ .` for all RHS
#'   variables plus the dependent variable. Use `NULL` for no CSAs (plain MG).
#' @param cross_section_lags Integer or named integer vector: number of lags of
#'   CSAs. Default 0. A single integer applies to all CSA variables.
#' @param pooled_vars A one-sided formula specifying which coefficients to
#'   constrain equal across units (pooled estimation). Default `NULL` (all
#'   heterogeneous).
#' @param include_constant Logical: include unit-specific intercepts? Default
#'   `TRUE`.
#' @param unit_trend Logical: include unit-specific linear trends? Default
#'   `FALSE`.
#' @param bias_correction Character: bias correction method. One of `"none"`,
#'   `"jackknife"`, `"recursive"`. Default `"none"`.
#' @param long_run_vars A one-sided formula specifying long-run variables
#'   (reserved for CSDL/CSARDL models). Default `NULL`.
#' @param long_run_model Character: long-run model specification (reserved).
#'   Default `NULL`.
#' @param csdl_xlags Integer: number of lags of \eqn{\Delta x} to include as
#'   short-run controls when `model = "csdl"`. Default 3.
#' @param run_cd_test Logical: run the Pesaran CD test on residuals? Default
#'   `FALSE`.
#' @param full_sample Logical: use the full (unbalanced) sample? Default
#'   `FALSE`.
#' @param verbose Logical: print progress messages? Default `FALSE`.
#' @param ... Additional arguments passed to model-specific estimators.
#'
#' @return An object of class `dcce_fit` (and a model-specific subclass).
#'
#' @export
#'
#' @examples
#' # Simple Mean Group estimation
#' df <- data.frame(
#'   id = rep(1:10, each = 20),
#'   t  = rep(1:20, 10),
#'   y  = rnorm(200),
#'   x  = rnorm(200)
#' )
#' fit <- dcce(df, unit_index = "id", time_index = "t",
#'             formula = y ~ x, model = "mg", cross_section_vars = NULL)
#' coef(fit)
dcce <- function(data, unit_index, time_index, formula,
                 model = c("dcce", "cce", "mg", "rcce", "pmg", "csdl", "csardl"),
                 cross_section_vars = ~ .,
                 cross_section_lags = 0L,
                 pooled_vars = NULL,
                 include_constant = TRUE,
                 unit_trend = FALSE,
                 bias_correction = c("none", "jackknife", "recursive"),
                 long_run_vars = NULL,
                 long_run_model = NULL,
                 csdl_xlags = 3L,
                 run_cd_test = FALSE,
                 full_sample = FALSE,
                 verbose = FALSE,
                 ...) {

  model <- rlang::arg_match(model)
  bias_correction <- rlang::arg_match(bias_correction)
  call <- match.call()

  # -- 1. Prepare panel -------------------------------------------------------
  panel <- .make_panel(data, unit_index = unit_index, time_index = time_index)
  bal <- .check_balance(panel)
  unit_var <- attr(panel, "unit_var")
  time_var <- attr(panel, "time_var")

  # -- 2. Parse formula -------------------------------------------------------
  parsed <- .parse_dcce_formula(formula, panel, unit_var, time_var)
  y_name <- parsed$y_name
  x_names <- parsed$x_names

  # -- 2b. CS-DL formula augmentation -----------------------------------------
  # For CS-DL, replace the LHS with Delta y and augment the RHS with
  # contemporaneous and lagged Delta x terms. The long-run coefficient is
  # the coefficient on the level of x (which is kept in x_names).
  csdl_lr_names <- NULL
  if (model == "csdl") {
    aug <- .csdl_augment(panel, y_name, x_names, p_x = as.integer(csdl_xlags))
    panel <- aug$panel
    y_name <- aug$y_diff_name
    x_names <- aug$rhs_terms
    csdl_lr_names <- aug$lr_names
  }

  # -- 3. Build CSAs if requested ---------------------------------------------
  csa_vars <- NULL
  csa_colnames <- NULL

  if (!is.null(cross_section_vars)) {
    csa_vars <- .resolve_cr_vars(cross_section_vars, y_name, x_names, panel)
    panel <- .build_csa(panel, csa_vars, lags = cross_section_lags)
    csa_colnames <- .get_csa_colnames(panel)
  }

  # -- 4. Build design matrices per unit --------------------------------------
  units <- unique(panel[[unit_var]])
  N <- length(units)

  # Determine regressors
  rhs_vars <- x_names
  if (include_constant) rhs_vars <- c("(Intercept)", rhs_vars)
  if (unit_trend) rhs_vars <- c(rhs_vars, "(trend)")

  unit_results <- list()
  dropped_units <- character(0)

  for (i in seq_along(units)) {
    u <- units[i]
    idx <- which(panel[[unit_var]] == u)
    yi <- panel[[y_name]][idx]
    Ti <- length(idx)

    # Build X matrix
    Xi <- .build_unit_X(panel, idx, x_names, csa_colnames,
                        include_constant, unit_trend)

    # Drop rows with NAs (from lags/diffs)
    complete <- stats::complete.cases(cbind(yi, Xi))
    yi <- yi[complete]
    Xi <- Xi[complete, , drop = FALSE]
    Ti_eff <- length(yi)

    # Check degrees of freedom
    K_total <- ncol(Xi)
    if (Ti_eff <= K_total) {
      dropped_units <- c(dropped_units, as.character(u))
      next
    }

    # Run OLS
    fit_i <- .unit_ols(yi, Xi)
    unit_results[[as.character(u)]] <- fit_i
  }

  if (length(dropped_units) > 0L) {
    cli::cli_warn(c(
      "Dropped {length(dropped_units)} unit{?s} due to insufficient degrees of freedom.",
      "i" = "Dropped: {paste(dropped_units, collapse = ', ')}"
    ))
  }

  N_used <- length(unit_results)
  if (N_used == 0L) {
    cli::cli_abort("No units remaining after dropping units with insufficient observations.")
  }

  # -- 5. Extract only the "interesting" coefficients (not CSA/intercept) ------
  # For MG aggregation, we want only the structural parameters
  coef_names_structural <- x_names
  if (include_constant) coef_names_structural <- c("(Intercept)", coef_names_structural)
  if (unit_trend) coef_names_structural <- c(coef_names_structural, "(trend)")

  coef_list_full <- lapply(unit_results, `[[`, "b")
  coef_list_structural <- lapply(coef_list_full, function(b) b[coef_names_structural])

  # -- 6. Aggregate -----------------------------------------------------------
  b_mg <- .mg_aggregate(coef_list_structural)
  V_mg <- .mg_variance(coef_list_structural, b_mg)

  # -- 6b. Bias correction (if requested) -------------------------------------
  jk_result <- NULL
  if (bias_correction == "jackknife" && model %in% c("dcce", "cce")) {
    jk_result <- .jackknife_bias_correction(
      panel, y_name, x_names, csa_colnames,
      unit_var, time_var, include_constant, unit_trend, b_mg
    )
    b_mg <- jk_result$b_jk
    names(b_mg) <- coef_names_structural
  } else if (bias_correction == "recursive") {
    # Stub: recursive mean adjustment not yet implemented
    cli::cli_warn("Recursive bias correction is not yet implemented; proceeding without correction.")
  }

  # -- 7. Compute summary statistics ------------------------------------------
  resid_list <- lapply(unit_results, `[[`, "e")
  all_resids <- unlist(resid_list)
  rmse <- sqrt(mean(all_resids^2))

  r2_units <- vapply(unit_results, `[[`, numeric(1), "r2")
  r2_mg <- mean(r2_units, na.rm = TRUE)

  # Total observations used
  n_obs <- sum(vapply(resid_list, length, integer(1)))

  # -- 8. Build output object -------------------------------------------------
  result <- list(
    call                  = call,
    coefficients          = b_mg,
    vcov                  = V_mg,
    se                    = sqrt(diag(V_mg)),
    unit_coefs            = coef_list_structural,
    unit_results          = unit_results,
    residuals             = all_resids,
    resid_list            = resid_list,
    fitted_values         = NULL,  # populated below
    r2                    = r2_mg,
    rmse                  = rmse,
    N                     = N_used,
    N_original            = N,
    T_bar                 = bal$T_bar,
    T_min                 = bal$T_min,
    T_max                 = bal$T_max,
    n_obs                 = n_obs,
    is_balanced           = bal$is_balanced,
    dropped_units         = dropped_units,
    y_name                = y_name,
    x_names               = x_names,
    cross_section_vars_used = csa_vars,
    csa_colnames          = csa_colnames,
    cross_section_lags    = cross_section_lags,
    model                 = model,
    formula               = formula,
    unit_index            = unit_index,
    time_index            = time_index,
    panel                 = panel,
    include_constant      = include_constant,
    unit_trend            = unit_trend,
    bias_correction       = bias_correction
  )

  # Assign S3 class
  model_class <- switch(model,
    mg    = "dcce_mg_fit",
    cce   = "dcce_cce_fit",
    dcce  = "dcce_dcce_fit",
    pmg   = "dcce_pmg_fit",
    csdl  = "dcce_csdl_fit",
    csardl = "dcce_csardl_fit",
    rcce  = "dcce_rcce_fit"
  )
  class(result) <- c(model_class, "dcce_fit")

  # -- 9. Long-run / adjustment post-processing -------------------------------
  if (model == "csardl") {
    result <- .csardl_postprocess(result)
  } else if (model == "pmg") {
    # PMG shares the CS-ARDL backbone but pools the long-run coefficients.
    result <- .csardl_postprocess(result)
    result <- .pmg_postprocess(result)
  } else if (model == "csdl") {
    result <- .csdl_postprocess(result, csdl_lr_names)
  }

  result
}


# ------------------------------------------------------------------------------
# Formula parsing helpers
# ------------------------------------------------------------------------------

#' Parse dcce formula
#'
#' Extracts the dependent variable name and regressor names from the formula.
#' Evaluates `L()`, `D()`, `Lrange()` calls by applying them panel-aware
#' and adding the resulting columns to the panel.
#'
#' @param formula A formula.
#' @param panel A panel data.frame.
#' @param unit_var Character: unit variable name.
#' @param time_var Character: time variable name.
#' @return A list with `y_name` and `x_names`.
#' @keywords internal
.parse_dcce_formula <- function(formula, panel, unit_var, time_var) {
  # Get the LHS
  lhs <- formula[[2L]]
  y_name <- .eval_formula_term(lhs, panel, unit_var, time_var)

  # Get the RHS terms
  rhs <- formula[[3L]]
  tt <- stats::terms(formula, data = panel)
  term_labels <- attr(tt, "term.labels")

  x_names <- character(0)
  for (tl in term_labels) {
    nm <- .eval_formula_term(parse(text = tl)[[1L]], panel, unit_var, time_var)
    x_names <- c(x_names, nm)
  }

  list(y_name = y_name, x_names = x_names)
}


#' Evaluate a formula term
#'
#' Handles plain variable names, `L()`, `D()`, and `Lrange()` calls.
#' Adds the computed column(s) to the parent panel in the calling
#' environment.
#'
#' @param expr An expression (from formula).
#' @param panel A panel data.frame.
#' @param unit_var Character: unit variable name.
#' @param time_var Character: time variable name.
#' @return Character vector of column name(s) created or referenced.
#' @keywords internal
.eval_formula_term <- function(expr, panel, unit_var, time_var) {
  if (is.name(expr)) {
    # Simple variable reference
    vname <- as.character(expr)
    if (!vname %in% names(panel)) {
      cli::cli_abort("Variable {.var {vname}} not found in data.")
    }
    return(vname)
  }

  if (is.call(expr)) {
    fn_name <- as.character(expr[[1L]])

    if (fn_name == "L") {
      var_expr <- expr[[2L]]
      k <- if (length(expr) >= 3L) eval(expr[[3L]]) else 1L
      var_name <- as.character(var_expr)
      col_name <- sprintf("L(%s,%d)", var_name, k)

      # Compute panel-aware lag and add to panel
      parent_env <- parent.frame(2L)  # .parse_dcce_formula's caller
      p <- get("panel", envir = parent_env)
      if (!col_name %in% names(p)) {
        p[[col_name]] <- .lag_panel(p, var_name, k)
        assign("panel", p, envir = parent_env)
      }
      return(col_name)
    }

    if (fn_name == "D") {
      var_expr <- expr[[2L]]
      k <- if (length(expr) >= 3L) eval(expr[[3L]]) else 1L
      var_name <- as.character(var_expr)
      col_name <- sprintf("D(%s,%d)", var_name, k)

      parent_env <- parent.frame(2L)
      p <- get("panel", envir = parent_env)
      if (!col_name %in% names(p)) {
        p[[col_name]] <- .diff_panel(p, var_name, k)
        assign("panel", p, envir = parent_env)
      }
      return(col_name)
    }

    if (fn_name == "Lrange") {
      var_expr <- expr[[2L]]
      k0 <- eval(expr[[3L]])
      k1 <- eval(expr[[4L]])
      var_name <- as.character(var_expr)
      col_names <- sprintf("L(%s,%d)", var_name, k0:k1)

      parent_env <- parent.frame(2L)
      p <- get("panel", envir = parent_env)
      for (ki in k0:k1) {
        cn <- sprintf("L(%s,%d)", var_name, ki)
        if (!cn %in% names(p)) {
          p[[cn]] <- .lag_panel(p, var_name, ki)
        }
      }
      assign("panel", p, envir = parent_env)
      return(col_names)
    }

    # Check for d.variable convention (common in panel data)
    if (fn_name == "d" || fn_name == "d.") {
      var_expr <- expr[[2L]]
      var_name <- as.character(var_expr)
      col_name <- paste0("d.", var_name)

      parent_env <- parent.frame(2L)
      p <- get("panel", envir = parent_env)
      if (!col_name %in% names(p)) {
        p[[col_name]] <- .diff_panel(p, var_name, 1L)
        assign("panel", p, envir = parent_env)
      }
      return(col_name)
    }
  }

  # Fallback: try to deparse the expression and check if it's a column name
  col_name <- deparse(expr)
  # Handle names with dots like d.log_rgdpo that are literal column names
  parent_env <- parent.frame(2L)
  p <- get("panel", envir = parent_env)
  if (col_name %in% names(p)) {
    return(col_name)
  }

  cli::cli_abort("Cannot parse formula term: {.code {col_name}}")
}


# ------------------------------------------------------------------------------
# CR variable resolution
# ------------------------------------------------------------------------------

#' Resolve cross-sectional average variables
#'
#' @param cr_formula A one-sided formula for CSA variables.
#' @param y_name Character: dependent variable name.
#' @param x_names Character: regressor names.
#' @param panel A panel data.frame.
#' @return Character vector of variable names to compute CSAs for.
#' @keywords internal
.resolve_cr_vars <- function(cr_formula, y_name, x_names, panel) {
  if (is.null(cr_formula)) return(NULL)

  # Check for ~ . pattern using rlang
  rhs <- rlang::f_rhs(cr_formula)
  if (identical(rhs, quote(.))) {
    # Use y + all x variables (in level form)
    # Extract base variable names from transformed names
    base_vars <- c(y_name, x_names)
    # Filter to only variables that exist in panel
    base_vars <- .extract_base_vars(base_vars, panel)
    return(unique(base_vars))
  }

  # Explicit variable list
  cr_terms <- attr(stats::terms(cr_formula), "term.labels")
  cr_terms
}


#' Extract base variable names
#'
#' Given potentially transformed names like "L(y,1)" or "D(x,1)",
#' extract the underlying variable names.
#'
#' @param var_names Character vector of variable names.
#' @param panel A panel data.frame.
#' @return Character vector of base variable names that exist in panel.
#' @keywords internal
.extract_base_vars <- function(var_names, panel) {
  base <- character(0)
  for (v in var_names) {
    # Check if it's a transformed name like L(x,1) or D(x,1)
    m <- regmatches(v, regexec("^[LD]\\(([^,]+),", v))
    if (length(m[[1L]]) > 1L) {
      base <- c(base, m[[1L]][2L])
    } else if (startsWith(v, "d.")) {
      # d.variable convention -- base is without d. prefix
      base_v <- sub("^d\\.", "", v)
      if (base_v %in% names(panel)) {
        base <- c(base, base_v)
      } else if (v %in% names(panel)) {
        base <- c(base, v)
      }
    } else if (v %in% names(panel)) {
      base <- c(base, v)
    }
  }
  unique(base)
}


# ------------------------------------------------------------------------------
# Unit X matrix builder
# ------------------------------------------------------------------------------

#' Build the regressor matrix for a single unit
#'
#' @param panel Panel data.frame.
#' @param idx Integer indices for this unit.
#' @param x_names Character: structural regressor names.
#' @param csa_colnames Character: CSA column names (or NULL).
#' @param include_constant Logical: include intercept?
#' @param unit_trend Logical: include time trend?
#' @return A numeric matrix.
#' @keywords internal
.build_unit_X <- function(panel, idx, x_names, csa_colnames,
                          include_constant, unit_trend) {
  Ti <- length(idx)
  cols <- list()

  if (include_constant) {
    cols[["(Intercept)"]] <- rep(1, Ti)
  }

  for (v in x_names) {
    cols[[v]] <- panel[[v]][idx]
  }

  if (unit_trend) {
    cols[["(trend)"]] <- seq_len(Ti)
  }

  # CSA columns
  if (!is.null(csa_colnames)) {
    for (v in csa_colnames) {
      cols[[v]] <- panel[[v]][idx]
    }
  }

  X <- do.call(cbind, cols)
  colnames(X) <- names(cols)
  X
}
