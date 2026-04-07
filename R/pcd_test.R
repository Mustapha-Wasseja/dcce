#' Cross-Sectional Dependence Tests
#'
#' Computes cross-sectional dependence test statistics for panel data residuals.
#' Implements the Pesaran (2015) CD test, CDw (Juodis & Reese 2022),
#' PEA (Fan et al. 2015), and CD* (Pesaran & Xie 2021).
#'
#' @param x Either a numeric vector of residuals (stacked by unit),
#'   a `dcce_fit` object, a numeric matrix (N x T) of residuals,
#'   or a data.frame containing the panel structure.
#' @param ... Arguments passed to methods.
#' @param data A data.frame containing the panel structure. Required if
#'   `x` is a vector.
#' @param unit_index Character: name of the unit variable in `data`.
#' @param time_index Character: name of the time variable in `data`.
#' @param test Character vector specifying which tests to compute.
#'   One or more of `"pesaran"`, `"cdw"`, `"pea"`, `"cdstar"`.
#'   Default uses all four.
#' @param n_reps Integer: number of Rademacher draws for CDw. Default 500.
#' @param n_pca Integer: number of principal components for CD* bias
#'   correction. Default 1.
#'
#' @return An object of class `dcce_cd` containing:
#'   \describe{
#'     \item{statistics}{A data.frame with columns `test`, `statistic`, `p_value`.}
#'     \item{N}{Number of cross-sectional units.}
#'     \item{T_bar}{Average time dimension.}
#'     \item{rho_ij}{Matrix of pairwise correlations (if retained).}
#'   }
#'
#' @export
#'
#' @examples
#' set.seed(42)
#' df <- data.frame(
#'   id = rep(1:10, each = 20),
#'   t  = rep(1:20, 10),
#'   e  = rnorm(200)
#' )
#' cd <- pcd_test(df$e, data = df, unit_index = "id", time_index = "t",
#'                test = "pesaran")
#' print(cd)
pcd_test <- function(x, ...) {
  UseMethod("pcd_test")
}


#' @rdname pcd_test
#' @export
pcd_test.dcce_fit <- function(x, ...,
                              test   = c("pesaran", "cdw", "pea", "cdstar"),
                              n_reps = 500L,
                              n_pca  = 1L) {
  test <- .validate_cd_tests(test)

  resid_list <- x$resid_list
  N <- length(resid_list)
  em <- .resid_list_to_matrix(resid_list)

  .pcd_test_core(em, N, test = test, n_reps = n_reps, n_pca = n_pca)
}


#' @rdname pcd_test
#' @export
pcd_test.data.frame <- function(x, ...,
                                unit_index = NULL,
                                time_index = NULL,
                                test   = c("pesaran", "cdw", "pea", "cdstar"),
                                n_reps = 500L,
                                n_pca  = 1L) {
  test <- .validate_cd_tests(test)

  if (is.null(unit_index) || is.null(time_index)) {
    cli::cli_abort(
      "{.arg unit_index} and {.arg time_index} are required when {.arg x} is a data.frame."
    )
  }

  panel <- .make_panel(x, unit_index = unit_index, time_index = time_index)
  # Assume last numeric column as residuals — caller should supply matrix or vector

  cli::cli_abort(
    "Use the vector or matrix method instead. Pass residuals directly, not a data.frame."
  )
}


#' @rdname pcd_test
#' @export
pcd_test.matrix <- function(x, ...,
                            test   = c("pesaran", "cdw", "pea", "cdstar"),
                            n_reps = 500L,
                            n_pca  = 1L) {
  test <- .validate_cd_tests(test)

  em <- x
  N <- nrow(em)

  .pcd_test_core(em, N, test = test, n_reps = n_reps, n_pca = n_pca)
}


#' @rdname pcd_test
#' @export
pcd_test.default <- function(x, ...,
                             data       = NULL,
                             unit_index = NULL,
                             time_index = NULL,
                             test   = c("pesaran", "cdw", "pea", "cdstar"),
                             n_reps = 500L,
                             n_pca  = 1L) {
  test <- .validate_cd_tests(test)

  if (!is.numeric(x)) {
    cli::cli_abort("{.arg x} must be numeric, a matrix, a data.frame, or a {.cls dcce_fit} object.")
  }

  if (is.null(data) || is.null(unit_index) || is.null(time_index)) {
    cli::cli_abort(
      "{.arg data}, {.arg unit_index}, and {.arg time_index} are required when {.arg x} is a numeric vector."
    )
  }

  panel <- .make_panel(data, unit_index = unit_index, time_index = time_index)
  em <- .resid_vector_to_matrix(x, panel)
  N <- nrow(em)

  .pcd_test_core(em, N, test = test, n_reps = n_reps, n_pca = n_pca)
}


#' Validate CD test names
#' @keywords internal
.validate_cd_tests <- function(test) {
  valid_tests <- c("pesaran", "cdw", "pea", "cdstar")
  bad <- setdiff(test, valid_tests)
  if (length(bad) > 0L) {
    cli::cli_abort(
      c("Unknown test{?s}: {.val {bad}}.",
        "i" = "Valid choices are {.val {valid_tests}}.")
    )

  }
  test
}


#' Core computation for pcd_test
#' @keywords internal
.pcd_test_core <- function(em, N, test, n_reps, n_pca) {
  T_val <- ncol(em)
  results <- data.frame(test = character(0), statistic = numeric(0),
                        p_value = numeric(0), stringsAsFactors = FALSE)

  # ── Compute pairwise correlations ──────────────────────────────────────────
  rho_ij <- .pairwise_correlations(em)

  # ── pesaran (Pesaran 2015) ────────────────────────────────────────────────
  if ("pesaran" %in% test) {
    cd_stat <- .cd_pesaran(rho_ij, N, T_val, em)
    p_val <- 2 * stats::pnorm(-abs(cd_stat))
    results <- rbind(results, data.frame(test = "pesaran", statistic = cd_stat,
                                         p_value = p_val, stringsAsFactors = FALSE))
  }

  # ── cdw (Juodis & Reese 2022) ────────────────────────────────────────────
  if ("cdw" %in% test) {
    cdw_stat <- .cd_weighted(em, n_reps)
    p_val <- 2 * stats::pnorm(-abs(cdw_stat))
    results <- rbind(results, data.frame(test = "cdw", statistic = cdw_stat,
                                         p_value = p_val, stringsAsFactors = FALSE))
  }

  # ── pea (Fan et al. 2015) ────────────────────────────────────────────────
  if ("pea" %in% test) {
    pea_stat <- .cd_pea(rho_ij, N, T_val, em)
    p_val <- 2 * stats::pnorm(-abs(pea_stat))
    results <- rbind(results, data.frame(test = "pea", statistic = pea_stat,
                                         p_value = p_val, stringsAsFactors = FALSE))
  }

  # ── cdstar (Pesaran & Xie 2021) ──────────────────────────────────────────
  if ("cdstar" %in% test) {
    cdstar_stat <- .cd_star(rho_ij, N, T_val, em, n_pca)
    p_val <- 2 * stats::pnorm(-abs(cdstar_stat))
    results <- rbind(results, data.frame(test = "cdstar", statistic = cdstar_stat,
                                         p_value = p_val, stringsAsFactors = FALSE))
  }

  out <- list(
    statistics = results,
    N          = N,
    T_bar      = T_val,
    rho_ij     = rho_ij
  )
  class(out) <- "dcce_cd"
  out
}


# ──────────────────────────────────────────────────────────────────────────────
# Internal helpers
# ──────────────────────────────────────────────────────────────────────────────

#' Compute pairwise correlations of residuals
#' @keywords internal
.pairwise_correlations <- function(em) {
  # em is N x T matrix; rows = units, cols = time
  N <- nrow(em)
  rho <- matrix(NA_real_, N, N)
  for (i in 1:(N - 1)) {
    for (j in (i + 1):N) {
      # Use common non-NA observations
      valid <- !is.na(em[i, ]) & !is.na(em[j, ])
      if (sum(valid) < 3L) next
      rho[i, j] <- stats::cor(em[i, valid], em[j, valid])
      rho[j, i] <- rho[i, j]
    }
  }
  diag(rho) <- 1
  rho
}

#' Pesaran (2015) CD statistic
#' @keywords internal
.cd_pesaran <- function(rho_ij, N, T_val, em) {
  # CD = sqrt(2T / (N*(N-1))) * sum_{i<j} rho_ij
  sum_rho <- 0
  count <- 0
  for (i in 1:(N - 1)) {
    for (j in (i + 1):N) {
      if (!is.na(rho_ij[i, j])) {
        # For unbalanced: weight by sqrt(T_ij)
        valid <- !is.na(em[i, ]) & !is.na(em[j, ])
        T_ij <- sum(valid)
        sum_rho <- sum_rho + sqrt(T_ij) * rho_ij[i, j]
        count <- count + 1
      }
    }
  }
  sqrt(2 / (N * (N - 1))) * sum_rho
}

#' CDw (Juodis & Reese 2022) weighted statistic
#' @keywords internal
.cd_weighted <- function(em, n_reps) {
  N <- nrow(em)
  T_val <- ncol(em)

  cdw_vals <- numeric(n_reps)
  for (r in seq_len(n_reps)) {
    # Rademacher weights
    w <- sample(c(-1, 1), N, replace = TRUE)
    # Weight residuals
    em_w <- em * w
    # Compute CD on weighted residuals
    sum_rho <- 0
    for (i in 1:(N - 1)) {
      for (j in (i + 1):N) {
        valid <- !is.na(em_w[i, ]) & !is.na(em_w[j, ])
        if (sum(valid) < 3L) next
        rho_w <- stats::cor(em_w[i, valid], em_w[j, valid])
        T_ij <- sum(valid)
        sum_rho <- sum_rho + sqrt(T_ij) * rho_w
      }
    }
    cdw_vals[r] <- sqrt(2 / (N * (N - 1))) * sum_rho
  }
  mean(cdw_vals)
}

#' PEA (Fan et al. 2015) power-enhanced statistic
#' @keywords internal
.cd_pea <- function(rho_ij, N, T_val, em) {
  cd <- .cd_pesaran(rho_ij, N, T_val, em)
  # Power enhancement term
  threshold <- 2 * sqrt(log(N) / T_val)
  pe_term <- 0
  for (i in 1:(N - 1)) {
    for (j in (i + 1):N) {
      if (!is.na(rho_ij[i, j]) && abs(rho_ij[i, j]) > threshold) {
        pe_term <- pe_term + abs(rho_ij[i, j])
      }
    }
  }
  cd + pe_term
}

#' CD* (Pesaran & Xie 2021) bias-corrected statistic
#' @keywords internal
.cd_star <- function(rho_ij, N, T_val, em, n_pca) {
  cd <- .cd_pesaran(rho_ij, N, T_val, em)

  # Estimate Theta from first n_pca principal components
  # Remove NAs for PCA: use complete cases
  em_complete <- em[, colSums(is.na(em)) == 0, drop = FALSE]
  if (ncol(em_complete) < 3) return(cd)

  pc <- tryCatch({
    svd_res <- svd(scale(t(em_complete)), nu = n_pca, nv = n_pca)
    svd_res$u[, 1:n_pca, drop = FALSE]
  }, error = function(e) return(NULL))

  if (is.null(pc)) return(cd)

  # Theta = average of squared loadings / N
  # Approximate from eigenvalues
  eigenvals <- svd(em_complete)$d^2 / (ncol(em_complete) - 1)
  Theta <- sum(eigenvals[1:min(n_pca, length(eigenvals))]) / sum(eigenvals)
  Theta <- min(Theta, 0.99)  # bound away from 1

  # CD* = (CD + sqrt(T/2) * Theta) / (1 - Theta)
  cd_star <- (cd + sqrt(T_val / 2) * Theta) / (1 - Theta)
  cd_star
}

#' Convert residual list to matrix
#' @keywords internal
.resid_list_to_matrix <- function(resid_list) {
  N <- length(resid_list)
  T_max <- max(vapply(resid_list, length, integer(1)))
  em <- matrix(NA_real_, N, T_max)
  for (i in seq_along(resid_list)) {
    Ti <- length(resid_list[[i]])
    em[i, 1:Ti] <- resid_list[[i]]
  }
  em
}

#' Convert stacked residual vector to matrix
#' @keywords internal
.resid_vector_to_matrix <- function(resids, panel) {
  unit_var <- attr(panel, "unit_var")
  time_var <- attr(panel, "time_var")
  units <- unique(panel[[unit_var]])
  times <- sort(unique(panel[[time_var]]))
  N <- length(units)
  T_val <- length(times)
  em <- matrix(NA_real_, N, T_val)

  for (i in seq_along(units)) {
    idx <- which(panel[[unit_var]] == units[i])
    t_idx <- match(panel[[time_var]][idx], times)
    em[i, t_idx] <- resids[idx]
  }
  em
}


#' Print method for dcce_cd objects
#'
#' @param x A `dcce_cd` object.
#' @param ... Ignored.
#' @return Invisibly returns `x`.
#' @export
print.dcce_cd <- function(x, ...) {
  # Map lowercase test names to display labels
  label_map <- c(
    pesaran = "CD",
    cdw     = "CDw",
    pea     = "PEA",
    cdstar  = "CD*"
  )

  cat("\nCross-Sectional Dependence Tests\n")
  cat("================================\n")
  cat(sprintf("N = %d, T = %d\n\n", x$N, round(x$T_bar)))

  for (i in seq_len(nrow(x$statistics))) {
    display_name <- label_map[x$statistics$test[i]]
    if (is.na(display_name)) display_name <- x$statistics$test[i]
    cat(sprintf("%-8s  statistic = %8.4f  p-value = %.4f\n",
                display_name,
                x$statistics$statistic[i],
                x$statistics$p_value[i]))
  }
  cat("\n")
  invisible(x)
}
