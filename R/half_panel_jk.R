#' Half-Panel Jackknife Bias Correction
#'
#' Implements the Chudik & Pesaran (2015) half-panel jackknife for
#' dynamic CCE/DCCE estimators. Each unit's time series is split into
#' two halves at the midpoint; MG estimates are computed on each half,
#' and the full-sample MG is bias-corrected as
#' \deqn{\hat\beta_{HPJ} = 2\hat\beta_{full}
#'   - \frac{1}{2}(\hat\beta_{half1} + \hat\beta_{half2}).}
#'
#' This targets the Nickell (1981) bias that arises in short-T dynamic
#' panels when the lagged dependent variable is correlated with the
#' unit fixed effect.
#'
#' @param panel_list Named list of \code{list(y, X)} pairs (from the
#'   main \code{dcce()} pipeline).
#' @param coef_names Character vector: structural coefficient names to
#'   extract from each half-sample fit.
#' @param fast Logical: use the C++ fast path.
#' @return A named numeric vector of the half-panel MG average (to be
#'   combined with the full-sample estimate in \code{dcce()}), or
#'   \code{NULL} if the correction cannot be applied.
#' @keywords internal
.half_panel_jackknife <- function(panel_list, coef_names, fast = TRUE) {
  N <- length(panel_list)
  if (N == 0L) return(NULL)

  # For each unit, split the time series at the midpoint
  b_half1 <- matrix(NA_real_, N, length(coef_names))
  b_half2 <- matrix(NA_real_, N, length(coef_names))
  colnames(b_half1) <- colnames(b_half2) <- coef_names

  for (i in seq_len(N)) {
    u <- panel_list[[i]]
    Ti <- length(u$y)
    mid <- floor(Ti / 2)
    if (mid < ncol(u$X) + 1L || (Ti - mid) < ncol(u$X) + 1L) next

    # First half
    fit1 <- tryCatch(
      .unit_ols(u$y[1:mid], u$X[1:mid, , drop = FALSE]),
      error = function(e) NULL
    )
    # Second half
    fit2 <- tryCatch(
      .unit_ols(u$y[(mid + 1):Ti], u$X[(mid + 1):Ti, , drop = FALSE]),
      error = function(e) NULL
    )

    if (!is.null(fit1) && all(coef_names %in% names(fit1$b))) {
      b_half1[i, ] <- fit1$b[coef_names]
    }
    if (!is.null(fit2) && all(coef_names %in% names(fit2$b))) {
      b_half2[i, ] <- fit2$b[coef_names]
    }
  }

  ok <- stats::complete.cases(b_half1) & stats::complete.cases(b_half2)
  if (sum(ok) < 2L) return(NULL)

  mg_half1 <- colMeans(b_half1[ok, , drop = FALSE])
  mg_half2 <- colMeans(b_half2[ok, , drop = FALSE])

  # Return the average of the two half-panel MG estimates
  0.5 * (mg_half1 + mg_half2)
}
