#' Impulse Response Functions for Dynamic Panel Models
#'
#' Computes impulse response functions (IRFs) from a fitted dynamic panel
#' model (DCCE, CS-ARDL, or PMG). The IRF traces the response of the
#' dependent variable to a one-unit shock in a specific regressor over
#' a given number of horizons, using the Mean Group ARDL coefficient
#' estimates and the autoregressive lag polynomial.
#'
#' Confidence bands are computed via the cross-section bootstrap: units
#' are resampled with replacement, the MG coefficients are recomputed,
#' and the IRF is re-traced. The 2.5 and 97.5 percent quantiles of the
#' bootstrap distribution form the 95 percent band.
#'
#' @param object A \code{dcce_fit} object from a dynamic model (must
#'   contain at least one \code{L(y, k)} term).
#' @param impulse Character: the regressor that receives the shock.
#' @param horizon Integer: number of periods to trace. Default 10.
#' @param boot_reps Integer: bootstrap replications for confidence bands.
#'   0 = no bands. Default 200.
#' @param seed Integer: random seed.
#'
#' @return An object of class \code{dcce_irf} containing the impulse
#'   response path (\code{$irf}), optional bootstrap lower/upper bands,
#'   and metadata (impulse variable name and horizon length).
#'
#' @export
#' @examples
#' data(dcce_sim)
#' fit <- dcce(
#'   data = dcce_sim, unit_index = "unit", time_index = "time",
#'   formula = y ~ L(y, 1) + x,
#'   model = "dcce", cross_section_vars = ~ ., cross_section_lags = 3
#' )
#' ir <- irf(fit, impulse = "x", horizon = 10, boot_reps = 0)
#' print(ir)
irf <- function(object,
                impulse,
                horizon   = 10L,
                boot_reps = 200L,
                seed      = NULL) {
  if (!inherits(object, "dcce_fit")) {
    cli::cli_abort("{.arg object} must be a {.cls dcce_fit} object.")
  }

  b  <- object$coefficients
  y_name  <- object$y_name
  x_names <- object$x_names

  # Find the AR lag coefficients
  lag_pattern <- paste0("^L\\(", y_name, ",(\\d+)\\)$")
  ar_terms <- grep(lag_pattern, names(b), value = TRUE)
  if (length(ar_terms) == 0L) {
    cli::cli_abort(
      "No lagged dependent variable found. IRFs require a dynamic model."
    )
  }
  ar_orders <- as.integer(sub(lag_pattern, "\\1", ar_terms))
  ar_coefs  <- b[ar_terms]
  max_p     <- max(ar_orders)

  if (!impulse %in% names(b)) {
    cli::cli_abort(
      "{.arg impulse} = {.val {impulse}} not found in the coefficient vector."
    )
  }
  beta_impulse <- b[impulse]

  # Compute IRF via recursion
  H <- as.integer(horizon)
  ir <- .compute_irf(ar_coefs, ar_orders, max_p, beta_impulse, H)

  # Bootstrap bands
  lower <- upper <- NULL
  if (boot_reps > 0L) {
    if (!is.null(seed)) set.seed(seed)

    unit_coefs <- object$unit_coefs
    N <- length(unit_coefs)
    ir_boot <- matrix(NA_real_, boot_reps, H + 1L)

    for (r in seq_len(boot_reps)) {
      idx <- sample.int(N, N, replace = TRUE)
      B_star <- do.call(rbind, unit_coefs[idx])
      b_star <- colMeans(B_star, na.rm = TRUE)
      ar_star <- b_star[ar_terms]
      beta_star <- b_star[impulse]
      ir_boot[r, ] <- .compute_irf(ar_star, ar_orders, max_p, beta_star, H)
    }
    lower <- apply(ir_boot, 2, stats::quantile, 0.025, na.rm = TRUE)
    upper <- apply(ir_boot, 2, stats::quantile, 0.975, na.rm = TRUE)
  }

  out <- list(
    irf     = ir,
    lower   = lower,
    upper   = upper,
    impulse = impulse,
    horizon = H
  )
  class(out) <- "dcce_irf"
  out
}


#' Compute IRF from AR coefficients and impulse coefficient
#' @keywords internal
.compute_irf <- function(ar_coefs, ar_orders, max_p, beta_impulse, H) {
  # y_h = phi_1 y_{h-1} + phi_2 y_{h-2} + ... + beta * shock_0
  # shock_0 = 1 at h=0, 0 afterwards
  ir <- numeric(H + 1L)
  ir[1] <- beta_impulse   # h = 0

  for (h in seq_len(H)) {
    val <- 0
    for (j in seq_along(ar_orders)) {
      lag_j <- ar_orders[j]
      if (h - lag_j >= 0 && (h - lag_j + 1L) <= length(ir)) {
        val <- val + ar_coefs[j] * ir[h - lag_j + 1L]
      }
    }
    ir[h + 1L] <- val
  }
  ir
}


#' Print a dcce_irf object
#'
#' @param x A \code{dcce_irf} object.
#' @param ... Ignored.
#' @return Invisibly returns \code{x}.
#' @export
print.dcce_irf <- function(x, ...) {
  cat("\nImpulse Response Function\n")
  cat("========================\n")
  cat(sprintf("Impulse: %s, Horizon: %d\n\n", x$impulse, x$horizon))
  tab <- data.frame(h = 0:x$horizon, irf = round(x$irf, 6))
  if (!is.null(x$lower)) {
    tab$lower <- round(x$lower, 6)
    tab$upper <- round(x$upper, 6)
  }
  print(tab, row.names = FALSE)
  cat("\n")
  invisible(x)
}


#' Plot a dcce_irf object
#'
#' @param x A \code{dcce_irf} object.
#' @param ... Passed to \code{plot()}.
#' @return Invisibly returns \code{x}.
#' @export
plot.dcce_irf <- function(x, ...) {
  h <- 0:x$horizon
  ylim <- range(c(x$irf, x$lower, x$upper), na.rm = TRUE)
  plot(h, x$irf, type = "l", lwd = 2, col = "steelblue",
       xlab = "Horizon", ylab = "Response",
       main = paste("IRF:", x$impulse),
       ylim = ylim, ...)
  if (!is.null(x$lower)) {
    graphics::polygon(
      c(h, rev(h)), c(x$lower, rev(x$upper)),
      col = grDevices::adjustcolor("steelblue", alpha.f = 0.2),
      border = NA
    )
  }
  graphics::abline(h = 0, col = "gray60", lty = 2)
  invisible(x)
}
