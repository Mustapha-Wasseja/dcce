#' Rank Condition Classifier
#'
#' Implements the De Vos et al. (2024) rank condition classifier for CCE
#' estimators. Tests whether the rank condition holds (RC=1) which ensures
#' CCE consistency.
#'
#' **Note:** Only valid for static panels. Emits a warning if called on
#' dynamic or PMG models.
#'
#' @param object A `dcce_fit` object.
#' @param criterion Character: `"er"` or `"gr"` for factor number selection.
#'   Default `"er"`.
#'
#' @return An object of class `dcce_rank` with elements `m` (number of
#'   common factors), `g` (rank of average factor loadings), `RC` (1 if
#'   rank condition holds, 0 otherwise).
#' @export
rank_condition <- function(object, criterion = c("er", "gr")) {
  if (!inherits(object, "dcce_fit")) {
    cli::cli_abort("{.arg object} must be a {.cls dcce_fit} object.")
  }
  criterion <- rlang::arg_match(criterion)

  # Warn for dynamic models

  if (inherits(object, c("dcce_dcce_fit", "dcce_pmg_fit", "dcce_csdl_fit",
                          "dcce_csardl_fit"))) {
    cli::cli_warn(
      "Rank condition classifier is only valid for static panels (CCE/MG)."
    )
  }

  panel <- object$panel
  unit_var <- attr(panel, "unit_var")
  x_names <- object$x_names

  # Build data matrix from observed variables
  units <- unique(panel[[unit_var]])
  N <- length(units)

  # Use all structural variables for factor estimation
  all_vars <- c(object$y_name, x_names)
  all_vars <- all_vars[all_vars %in% names(panel)]

  # Build N x T matrix from cross-sectional averages
  csa_colnames <- object$csa_colnames
  if (is.null(csa_colnames) || length(csa_colnames) == 0L) {
    result <- list(m = 0L, g = 0L, RC = 1L)
    class(result) <- "dcce_rank"
    return(result)
  }

  # Estimate m (number of common factors) using ER/GR on observed data
  time_var <- attr(panel, "time_var")
  times <- sort(unique(panel[[time_var]]))
  T_val <- length(times)

  # Build data matrix: stack all variables into N x T
  X_mat <- matrix(NA_real_, N, T_val)
  for (i in seq_along(units)) {
    idx <- which(panel[[unit_var]] == units[i])
    t_idx <- match(panel[[time_var]][idx], times)
    X_mat[i, t_idx] <- panel[[all_vars[1]]][idx]
  }

  # Remove NA columns
  valid_cols <- colSums(is.na(X_mat)) == 0
  X_clean <- X_mat[, valid_cols, drop = FALSE]

  if (ncol(X_clean) < 3 || nrow(X_clean) < 3) {
    result <- list(m = NA_integer_, g = NA_integer_, RC = NA_integer_)
    class(result) <- "dcce_rank"
    return(result)
  }

  eigenvalues <- svd(X_clean)$d^2
  K_max <- floor(min(N, ncol(X_clean)) / 2)

  if (criterion == "er") {
    m <- .ahn_horenstein_er(eigenvalues, K_max)
  } else {
    m <- .ahn_horenstein_gr(eigenvalues, K_max)
  }

  # Estimate g (rank of average factor loadings) from CSA matrix
  csa_mat <- matrix(NA_real_, T_val, length(csa_colnames))
  colnames(csa_mat) <- csa_colnames
  # Take CSA values from first unit (they're the same for all units at each t)
  first_unit <- units[1]
  idx1 <- which(panel[[unit_var]] == first_unit)
  t_idx1 <- match(panel[[time_var]][idx1], times)

  for (j in seq_along(csa_colnames)) {
    csa_mat[t_idx1, j] <- panel[[csa_colnames[j]]][idx1]
  }

  valid_rows <- rowSums(is.na(csa_mat)) == 0
  csa_clean <- csa_mat[valid_rows, , drop = FALSE]

  if (nrow(csa_clean) < 3 || ncol(csa_clean) < 1) {
    g <- 0L
  } else {
    g <- as.integer(Matrix::rankMatrix(csa_clean))
  }

  RC <- as.integer(!(g < m))

  result <- list(m = m, g = g, RC = RC)
  class(result) <- "dcce_rank"
  result
}

#' @export
print.dcce_rank <- function(x, ...) {
  cat("\nRank Condition Classifier\n")
  cat("========================\n")
  cat(sprintf("Common factors (m)     : %s\n", ifelse(is.na(x$m), "NA", as.character(x$m))))
  cat(sprintf("CSA rank (g)           : %s\n", ifelse(is.na(x$g), "NA", as.character(x$g))))
  cat(sprintf("Rank condition (RC)    : %s\n", ifelse(is.na(x$RC), "NA",
    ifelse(x$RC == 1L, "HOLDS (CCE consistent)", "FAILS (CCE may be inconsistent)"))))
  cat("\n")
  invisible(x)
}
