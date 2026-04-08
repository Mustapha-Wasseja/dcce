#' Rolling-Window Panel Estimation
#'
#' Fits a sequence of \code{dcce()} models on overlapping time windows of
#' the panel, producing a time path of coefficient estimates. Useful for
#' detecting coefficient drift, parameter instability, or regime shifts
#' in long panels.
#'
#' At each window the function subsets the panel to observations whose
#' time index falls in \eqn{[t_k, t_k + \text{window} - 1]}, runs
#' \code{dcce()} with the user-supplied arguments, and collects the
#' Mean Group coefficient vector and its standard errors. The result is
#' returned as a \code{dcce_rolling} object containing the full list of
#' fits and a tidy tibble of coefficients indexed by window end-date.
#'
#' @param data A panel data.frame.
#' @param unit_index Character: unit identifier column.
#' @param time_index Character: time identifier column.
#' @param formula Two-sided model formula passed through to \code{dcce()}.
#' @param model Character: estimator to use (default \code{"cce"}). Same
#'   choices as \code{dcce()}.
#' @param window Integer: window length in time periods.
#' @param step Integer: number of time periods to advance between windows.
#'   Default \code{1L}.
#' @param min_units Integer: minimum number of cross-sectional units a
#'   window must retain after NA handling to be kept. Default \code{5L}.
#' @param verbose Logical: print progress? Default \code{FALSE}.
#' @param ... Additional arguments passed to \code{dcce()}
#'   (e.g. \code{cross_section_vars}, \code{cross_section_lags},
#'   \code{fast}).
#'
#' @return An object of class \code{dcce_rolling} containing:
#'   \describe{
#'     \item{fits}{List of \code{dcce_fit} objects, one per window
#'       (or \code{NULL} for windows that failed).}
#'     \item{coefficients}{A tibble with one row per
#'       (window end-date, term) combination, columns \code{window_end},
#'       \code{window_start}, \code{term}, \code{estimate}, \code{std.error},
#'       \code{conf.low}, \code{conf.high}.}
#'     \item{window}{The window length.}
#'     \item{step}{The step size.}
#'     \item{n_windows}{Number of successful windows.}
#'     \item{call}{The original call.}
#'   }
#'
#' @export
#' @examples
#' data(pwt8)
#' roll <- dcce_rolling(
#'   data       = pwt8,
#'   unit_index = "country",
#'   time_index = "year",
#'   formula    = d_log_rgdpo ~ log_hc + log_ck + log_ngd,
#'   model      = "cce",
#'   cross_section_vars = ~ .,
#'   window     = 20,
#'   step       = 2
#' )
#' print(roll)
dcce_rolling <- function(data,
                         unit_index,
                         time_index,
                         formula,
                         model     = "cce",
                         window    = 20L,
                         step      = 1L,
                         min_units = 5L,
                         verbose   = FALSE,
                         ...) {
  call <- match.call()
  window    <- as.integer(window)
  step      <- as.integer(step)
  min_units <- as.integer(min_units)

  if (window <= 1L) {
    cli::cli_abort("{.arg window} must be at least 2.")
  }
  if (step <= 0L) {
    cli::cli_abort("{.arg step} must be positive.")
  }

  panel <- .make_panel(data, unit_index = unit_index, time_index = time_index)
  time_var <- attr(panel, "time_var")
  all_times <- sort(unique(panel[[time_var]]))

  if (length(all_times) < window) {
    cli::cli_abort(
      "Panel has only {length(all_times)} time periods, but {.arg window} = {window}."
    )
  }

  # Enumerate window start positions
  starts <- seq(1L, length(all_times) - window + 1L, by = step)
  if (length(starts) == 0L) {
    cli::cli_abort("No valid windows for the chosen {.arg window} / {.arg step}.")
  }

  fits <- vector("list", length(starts))
  coef_rows <- list()

  for (k in seq_along(starts)) {
    t0 <- all_times[starts[k]]
    t1 <- all_times[starts[k] + window - 1L]

    if (verbose) {
      cli::cli_inform("Window {k}/{length(starts)}: {t0}-{t1}")
    }

    sub <- panel[panel[[time_var]] >= t0 & panel[[time_var]] <= t1, ,
                 drop = FALSE]
    # Restore panel attributes on the subset
    attr(sub, "unit_var") <- attr(panel, "unit_var")
    attr(sub, "time_var") <- attr(panel, "time_var")

    # Skip if too few units remain after subsetting
    n_units_window <- length(unique(sub[[attr(panel, "unit_var")]]))
    if (n_units_window < min_units) {
      fits[[k]] <- NULL
      next
    }

    fit_k <- tryCatch(
      dcce(
        data       = sub,
        unit_index = attr(panel, "unit_var"),
        time_index = attr(panel, "time_var"),
        formula    = formula,
        model      = model,
        ...
      ),
      error   = function(e) NULL,
      warning = function(w) {
        withCallingHandlers(
          dcce(
            data       = sub,
            unit_index = attr(panel, "unit_var"),
            time_index = attr(panel, "time_var"),
            formula    = formula,
            model      = model,
            ...
          ),
          warning = function(wi) invokeRestart("muffleWarning")
        )
      }
    )

    if (is.null(fit_k) || !inherits(fit_k, "dcce_fit")) {
      fits[[k]] <- NULL
      next
    }

    fits[[k]] <- fit_k

    b  <- fit_k$coefficients
    se <- fit_k$se
    if (length(b) == 0L) next

    coef_rows[[length(coef_rows) + 1L]] <- tibble::tibble(
      window_end   = t1,
      window_start = t0,
      term         = names(b),
      estimate     = unname(b),
      std.error    = unname(se),
      conf.low     = unname(b - 1.96 * se),
      conf.high    = unname(b + 1.96 * se)
    )
  }

  if (length(coef_rows) == 0L) {
    cli::cli_abort(
      "No windows produced a valid fit. Check {.arg window}, {.arg min_units}, or the model specification."
    )
  }

  coef_tbl <- do.call(rbind, coef_rows)

  out <- list(
    fits         = fits,
    coefficients = coef_tbl,
    window       = window,
    step         = step,
    n_windows    = sum(!vapply(fits, is.null, logical(1))),
    call         = call
  )
  class(out) <- "dcce_rolling"
  out
}


#' Print a dcce_rolling object
#'
#' @param x A \code{dcce_rolling} object.
#' @param ... Ignored.
#' @return Invisibly returns \code{x}.
#' @export
print.dcce_rolling <- function(x, ...) {
  cat("\nRolling-Window Panel Estimation\n")
  cat("===============================\n")
  cat(sprintf("Windows fitted:  %d\n", x$n_windows))
  cat(sprintf("Window length:   %d time periods\n", x$window))
  cat(sprintf("Step:            %d\n", x$step))

  terms <- unique(x$coefficients$term)
  cat("\nCoefficient path (per term):\n")
  for (tm in terms) {
    rows <- x$coefficients[x$coefficients$term == tm, , drop = FALSE]
    cat(sprintf("  %-20s  mean %8.4f   sd %7.4f   [%.3f, %.3f]\n",
                tm,
                mean(rows$estimate, na.rm = TRUE),
                stats::sd(rows$estimate, na.rm = TRUE),
                min(rows$estimate, na.rm = TRUE),
                max(rows$estimate, na.rm = TRUE)))
  }
  cat("\nUse as.data.frame(x$coefficients) for the full time path.\n")
  invisible(x)
}


#' Plot a dcce_rolling coefficient path
#'
#' Produces one line per regressor showing the rolling-window coefficient
#' against the window end-date, with a 95% confidence ribbon built from
#' the unit-loop standard errors.
#'
#' @param x A \code{dcce_rolling} object.
#' @param terms Character vector of terms to plot. Default \code{NULL}
#'   plots all non-intercept terms.
#' @param ... Passed to the underlying graphics call.
#' @return Invisibly returns \code{x}.
#' @export
plot.dcce_rolling <- function(x, terms = NULL, ...) {
  tbl <- x$coefficients
  if (is.null(terms)) {
    terms <- setdiff(unique(tbl$term), "(Intercept)")
  }
  tbl <- tbl[tbl$term %in% terms, , drop = FALSE]
  if (nrow(tbl) == 0L) {
    cli::cli_warn("No terms to plot.")
    return(invisible(x))
  }

  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par), add = TRUE)

  n_terms <- length(terms)
  nr <- ceiling(sqrt(n_terms))
  nc <- ceiling(n_terms / nr)
  graphics::par(mfrow = c(nr, nc), mar = c(4, 4, 2, 1))

  for (tm in terms) {
    sub <- tbl[tbl$term == tm, , drop = FALSE]
    ylim <- range(c(sub$conf.low, sub$conf.high), na.rm = TRUE)
    plot(
      sub$window_end, sub$estimate,
      type = "l",
      lwd  = 2,
      col  = "steelblue",
      xlab = "Window end",
      ylab = "Estimate",
      main = tm,
      ylim = ylim
    )
    graphics::polygon(
      c(sub$window_end, rev(sub$window_end)),
      c(sub$conf.low,  rev(sub$conf.high)),
      col    = grDevices::adjustcolor("steelblue", alpha.f = 0.2),
      border = NA
    )
    graphics::abline(h = 0, col = "gray60", lty = 2)
  }
  invisible(x)
}
