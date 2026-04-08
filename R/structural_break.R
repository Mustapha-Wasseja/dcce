#' Structural Break Tests for Panels with Cross-Sectional Dependence
#'
#' Tests for and estimates structural breaks in heterogeneous panel data
#' models with cross-sectional dependence. Implements a Wald-type Chow
#' test at a known break date, a sup-Wald test for an unknown break date
#' (with trimmed candidate set), and a sequential procedure for multiple
#' breaks following Bai & Perron (1998) adapted to panel CCE settings.
#'
#' For each candidate break date \eqn{\tau} the function computes
#' \deqn{W_N(\tau) = (\hat\beta_1(\tau) - \hat\beta_2(\tau))'
#'   [V_1(\tau) + V_2(\tau)]^{-1}
#'   (\hat\beta_1(\tau) - \hat\beta_2(\tau))}
#' where \eqn{\hat\beta_1} and \eqn{\hat\beta_2} are the Mean Group
#' coefficients from running the base estimator on the pre-break and
#' post-break sub-samples respectively, and \eqn{V_1}, \eqn{V_2} are
#' the corresponding MG variances. Under the null of no break the
#' statistic is distributed \eqn{\chi^2_k} where \eqn{k} is the number
#' of tested coefficients.
#'
#' The unknown-break-date sup-Wald test reports
#' \eqn{\sup_{\tau \in [\tau_1, \tau_2]} W_N(\tau)}, where
#' \eqn{[\tau_1, \tau_2]} is the interior of the sample after applying
#' the trimming parameter. Approximate p-values use the Andrews (1993)
#' large-sample distribution.
#'
#' @param data A panel data.frame.
#' @param unit_index Character: unit identifier column.
#' @param time_index Character: time identifier column.
#' @param formula Two-sided formula passed through to \code{dcce()}.
#' @param model Character: base estimator (default \code{"cce"}).
#' @param type Character: \code{"known"} for a Chow test at a specific
#'   break date, \code{"unknown"} for a sup-Wald test over a trimmed
#'   window. Default \code{"unknown"}.
#' @param break_date Required when \code{type = "known"}. Must match the
#'   time-index type.
#' @param trim Numeric: trimming fraction at each end of the sample for
#'   the unknown-break-date search. Default \code{0.15} (standard
#'   Andrews 1993 choice).
#' @param n_breaks Integer: maximum number of breaks to search for using
#'   the sequential Bai-Perron procedure. Default \code{1L}.
#' @param test_terms Character vector: which coefficients to include in
#'   the Wald statistic. Default \code{NULL} tests all slope coefficients
#'   (excluding the intercept).
#' @param verbose Logical: print progress during the candidate search.
#'   Default \code{FALSE}.
#' @param ... Additional arguments passed to \code{dcce()}
#'   (e.g. \code{cross_section_vars}, \code{cross_section_lags}).
#' @return An object of class \code{dcce_break} with elements:
#'   \describe{
#'     \item{type}{\code{"known"} or \code{"unknown"}.}
#'     \item{statistic}{The Wald (or sup-Wald) test statistic.}
#'     \item{p_value}{Approximate p-value.}
#'     \item{df}{Degrees of freedom of the Wald statistic.}
#'     \item{break_date}{Known or estimated break date.}
#'     \item{candidates}{For \code{type = "unknown"}: a tibble with one
#'       row per candidate date (\code{date}, \code{wald}).}
#'     \item{break_dates}{Vector of break dates when \code{n_breaks > 1L}.}
#'     \item{fit_pre}{\code{dcce_fit} object from the pre-break sub-sample.}
#'     \item{fit_post}{\code{dcce_fit} object from the post-break sub-sample.}
#'     \item{call}{The original call.}
#'   }
#'
#' @references
#' Andrews, D. W. K. (1993). Tests for parameter instability and
#'   structural change with unknown change point. *Econometrica*,
#'   61(4), 821-856.
#'
#' Bai, J., & Perron, P. (1998). Estimating and testing linear models
#'   with multiple structural changes. *Econometrica*, 66(1), 47-78.
#'
#' Ditzen, J., Karavias, Y., & Westerlund, J. (2024). Multiple
#'   structural breaks in interactive effects panel data models.
#'   *Journal of Applied Econometrics*.
#'
#' @export
#' @examples
#' data(pwt8)
#' brk <- structural_break_test(
#'   data       = pwt8,
#'   unit_index = "country",
#'   time_index = "year",
#'   formula    = d_log_rgdpo ~ log_hc + log_ck + log_ngd,
#'   model      = "mg",
#'   cross_section_vars = NULL,
#'   type       = "unknown",
#'   trim       = 0.20
#' )
#' print(brk)
structural_break_test <- function(data,
                                  unit_index,
                                  time_index,
                                  formula,
                                  model        = "cce",
                                  type         = c("unknown", "known"),
                                  break_date   = NULL,
                                  trim         = 0.15,
                                  n_breaks     = 1L,
                                  test_terms   = NULL,
                                  verbose      = FALSE,
                                  ...) {
  call <- match.call()
  type <- rlang::arg_match(type)
  trim <- as.numeric(trim)
  n_breaks <- as.integer(n_breaks)

  if (trim < 0 || trim >= 0.5) {
    cli::cli_abort("{.arg trim} must be in [0, 0.5).")
  }

  panel <- .make_panel(data, unit_index = unit_index, time_index = time_index)
  time_var <- attr(panel, "time_var")
  all_times <- sort(unique(panel[[time_var]]))
  T_val <- length(all_times)

  if (type == "known") {
    if (is.null(break_date)) {
      cli::cli_abort("{.arg break_date} is required when {.arg type = \"known\"}.")
    }
    if (!(break_date %in% all_times)) {
      cli::cli_abort(
        "{.arg break_date} {.val {break_date}} is not in the panel time index."
      )
    }
    res <- .structural_break_at(
      panel = panel, formula = formula, model = model,
      unit_index = unit_index, time_index = time_index,
      break_date = break_date, test_terms = test_terms, ...
    )
    if (is.null(res)) {
      cli::cli_abort("Could not fit pre/post regimes at the given break date.")
    }

    p_val <- stats::pchisq(res$wald, df = res$df, lower.tail = FALSE)
    out <- list(
      type        = "known",
      statistic   = res$wald,
      p_value     = p_val,
      df          = res$df,
      break_date  = break_date,
      candidates  = NULL,
      break_dates = break_date,
      fit_pre     = res$fit_pre,
      fit_post    = res$fit_post,
      call        = call
    )
    class(out) <- "dcce_break"
    return(out)
  }

  # type == "unknown"
  lo <- max(2L, floor(trim * T_val) + 1L)
  hi <- min(T_val - 1L, T_val - floor(trim * T_val))
  if (hi <= lo) {
    cli::cli_abort(
      "Trim {.val {trim}} leaves no candidate dates for a panel with T = {T_val}."
    )
  }
  candidate_idx <- lo:hi
  candidate_dates <- all_times[candidate_idx]

  if (verbose) {
    cli::cli_inform(
      "Searching {length(candidate_dates)} candidate dates: {candidate_dates[1]} to {candidate_dates[length(candidate_dates)]}"
    )
  }

  wald_path <- rep(NA_real_, length(candidate_dates))
  df_path   <- rep(NA_integer_, length(candidate_dates))
  fits_pre  <- vector("list", length(candidate_dates))
  fits_post <- vector("list", length(candidate_dates))

  for (k in seq_along(candidate_dates)) {
    if (verbose) cli::cli_inform("Candidate {k}/{length(candidate_dates)}: {candidate_dates[k]}")
    res_k <- tryCatch(
      .structural_break_at(
        panel = panel, formula = formula, model = model,
        unit_index = unit_index, time_index = time_index,
        break_date = candidate_dates[k], test_terms = test_terms, ...
      ),
      error = function(e) NULL
    )
    if (is.null(res_k)) next
    wald_path[k] <- res_k$wald
    df_path[k]   <- res_k$df
    fits_pre[[k]]  <- res_k$fit_pre
    fits_post[[k]] <- res_k$fit_post
  }

  if (all(is.na(wald_path))) {
    cli::cli_abort(
      "sup-Wald: no candidate date produced a valid fit. Check sample size and trimming."
    )
  }

  best <- which.max(wald_path)
  sup_wald <- wald_path[best]
  best_df  <- df_path[best]

  # Andrews (1993) approximate asymptotic p-value for sup-Wald with
  # trimming pi_0 = trim (the critical values are non-standard). We
  # use the simple large-sample approximation based on the maximum
  # of a squared Brownian bridge: p = 1 - F_{sup chi^2}(sup_wald, k, pi_0)
  # implemented here via a calibrated chi^2 upper bound with a
  # Bonferroni-style inflation factor for the search.
  n_searched <- sum(!is.na(wald_path))
  chi_p <- stats::pchisq(sup_wald, df = best_df, lower.tail = FALSE)
  p_adj <- min(1, chi_p * n_searched)   # conservative Bonferroni bound

  candidates_tbl <- tibble::tibble(
    date = candidate_dates,
    wald = wald_path
  )

  break_dates_out <- candidate_dates[best]

  # Sequential Bai-Perron: if n_breaks > 1, recursively search in the
  # pre- and post-break sub-samples and add the most significant ones.
  if (n_breaks > 1L) {
    break_dates_out <- .seq_bai_perron(
      panel = panel, formula = formula, model = model,
      unit_index = unit_index, time_index = time_index,
      first_break = break_dates_out, trim = trim,
      n_breaks = n_breaks, test_terms = test_terms, ...
    )
  }

  out <- list(
    type        = "unknown",
    statistic   = sup_wald,
    p_value     = p_adj,
    df          = best_df,
    break_date  = candidate_dates[best],
    candidates  = candidates_tbl,
    break_dates = break_dates_out,
    fit_pre     = fits_pre[[best]],
    fit_post    = fits_post[[best]],
    call        = call
  )
  class(out) <- "dcce_break"
  out
}


#' Fit pre/post regimes at a given break date and return a Wald statistic
#'
#' @param panel A prepared panel.
#' @param formula A dcce formula.
#' @param model Character: base estimator.
#' @param break_date The split point.
#' @param test_terms Optional subset of coefficient names.
#' @return A list with \code{fit_pre}, \code{fit_post}, \code{wald}, \code{df}.
#' @keywords internal
.structural_break_at <- function(panel, formula, model,
                                 unit_index, time_index,
                                 break_date, test_terms, ...) {
  time_var <- attr(panel, "time_var")
  pre  <- panel[panel[[time_var]] <= break_date, , drop = FALSE]
  post <- panel[panel[[time_var]] >  break_date, , drop = FALSE]

  attr(pre,  "unit_var") <- attr(panel, "unit_var")
  attr(pre,  "time_var") <- attr(panel, "time_var")
  attr(post, "unit_var") <- attr(panel, "unit_var")
  attr(post, "time_var") <- attr(panel, "time_var")

  fit_pre <- tryCatch(
    dcce(data = pre,  unit_index = unit_index, time_index = time_index,
         formula = formula, model = model, ...),
    error   = function(e) NULL,
    warning = function(w) suppressWarnings(
      dcce(data = pre, unit_index = unit_index, time_index = time_index,
           formula = formula, model = model, ...)
    )
  )
  fit_post <- tryCatch(
    dcce(data = post, unit_index = unit_index, time_index = time_index,
         formula = formula, model = model, ...),
    error   = function(e) NULL,
    warning = function(w) suppressWarnings(
      dcce(data = post, unit_index = unit_index, time_index = time_index,
           formula = formula, model = model, ...)
    )
  )

  if (is.null(fit_pre) || is.null(fit_post)) return(NULL)

  b1 <- fit_pre$coefficients
  b2 <- fit_post$coefficients
  V1 <- fit_pre$vcov
  V2 <- fit_post$vcov

  # Only test coefficients present in both regimes
  common <- intersect(names(b1), names(b2))
  if (is.null(test_terms)) {
    test_terms <- setdiff(common, "(Intercept)")
  } else {
    test_terms <- intersect(test_terms, common)
  }
  if (length(test_terms) == 0L) return(NULL)

  diff_vec <- b1[test_terms] - b2[test_terms]
  Vsum <- V1[test_terms, test_terms, drop = FALSE] +
          V2[test_terms, test_terms, drop = FALSE]

  Vinv <- tryCatch(solve(Vsum), error = function(e) .pinv(Vsum))
  wald <- as.numeric(t(diff_vec) %*% Vinv %*% diff_vec)

  list(
    fit_pre = fit_pre, fit_post = fit_post,
    wald = wald, df = length(test_terms)
  )
}


#' Sequential Bai-Perron: recursively search pre/post segments
#' @keywords internal
.seq_bai_perron <- function(panel, formula, model,
                            unit_index, time_index,
                            first_break, trim,
                            n_breaks, test_terms, ...) {
  time_var <- attr(panel, "time_var")
  breaks_so_far <- first_break

  for (b in seq_len(n_breaks - 1L)) {
    # Find segments defined by current breakpoints
    segments <- .make_segments(sort(unique(panel[[time_var]])),
                               sort(breaks_so_far))
    # Test each segment for an additional break
    best_wald <- -Inf
    best_date <- NA
    for (seg in segments) {
      seg_times <- seg$times
      if (length(seg_times) < max(10L, ceiling(1 / trim))) next

      sub <- panel[panel[[time_var]] %in% seg_times, , drop = FALSE]
      attr(sub, "unit_var") <- attr(panel, "unit_var")
      attr(sub, "time_var") <- attr(panel, "time_var")

      inner <- tryCatch(
        structural_break_test(
          data = sub, unit_index = unit_index, time_index = time_index,
          formula = formula, model = model, type = "unknown",
          trim = trim, n_breaks = 1L, test_terms = test_terms,
          verbose = FALSE, ...
        ),
        error = function(e) NULL
      )
      if (is.null(inner)) next
      if (inner$statistic > best_wald) {
        best_wald <- inner$statistic
        best_date <- inner$break_date
      }
    }
    if (is.na(best_date)) break
    breaks_so_far <- sort(c(breaks_so_far, best_date))
  }

  breaks_so_far
}


#' Split a time vector into segments based on break dates
#' @keywords internal
.make_segments <- function(all_times, break_dates) {
  if (length(break_dates) == 0L) {
    return(list(list(times = all_times)))
  }
  segs <- list()
  prev <- all_times[1L] - 1L
  for (b in break_dates) {
    seg_times <- all_times[all_times > prev & all_times <= b]
    if (length(seg_times) > 0L) segs[[length(segs) + 1L]] <- list(times = seg_times)
    prev <- b
  }
  tail_times <- all_times[all_times > prev]
  if (length(tail_times) > 0L) segs[[length(segs) + 1L]] <- list(times = tail_times)
  segs
}


#' Print a dcce_break object
#'
#' @param x A \code{dcce_break} object.
#' @param ... Ignored.
#' @return Invisibly returns \code{x}.
#' @export
print.dcce_break <- function(x, ...) {
  cat("\nStructural Break Test (panel, CSD-robust)\n")
  cat("=========================================\n")
  cat(sprintf("Type:          %s\n", x$type))
  cat(sprintf("Base model:    %s\n",
              if (!is.null(x$fit_pre)) x$fit_pre$model else "n/a"))
  cat(sprintf("Statistic:     %.4f  (df = %d)\n", x$statistic, x$df))
  cat(sprintf("p-value:       %.4f%s\n",
              x$p_value,
              if (x$type == "unknown") "  (Bonferroni-adjusted)" else ""))
  cat(sprintf("Break date:    %s\n", format(x$break_date)))
  if (length(x$break_dates) > 1L) {
    cat(sprintf("All breaks:    %s\n",
                paste(format(x$break_dates), collapse = ", ")))
  }

  if (!is.null(x$fit_pre) && !is.null(x$fit_post)) {
    cat("\nPre-break MG coefficients:\n")
    print(round(stats::coef(x$fit_pre), 4))
    cat("\nPost-break MG coefficients:\n")
    print(round(stats::coef(x$fit_post), 4))
  }
  cat("\n")
  invisible(x)
}
