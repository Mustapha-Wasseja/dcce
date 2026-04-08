test_that(".spatial_validate_W row-normalises a raw contiguity matrix", {
  W <- matrix(c(0, 1, 1, 0,
                1, 0, 1, 0,
                1, 1, 0, 1,
                0, 0, 1, 0), 4, 4, byrow = TRUE)
  rownames(W) <- colnames(W) <- c("A", "B", "C", "D")
  W_norm <- .spatial_validate_W(W, units = c("A","B","C","D"))
  expect_true(all(abs(rowSums(W_norm) - 1) < 1e-12))
  expect_equal(diag(W_norm), c(A=0, B=0, C=0, D=0))
})


test_that(".spatial_validate_W errors on wrong dimension", {
  W <- matrix(0, 3, 3)
  expect_error(
    .spatial_validate_W(W, units = c("A","B","C","D")),
    regexp = "4 x 4"
  )
})


test_that(".spatial_validate_W aligns by rownames when present", {
  W <- matrix(c(0, 0.5, 0.5,
                0.3, 0, 0.7,
                1, 0, 0), 3, 3, byrow = TRUE)
  rownames(W) <- colnames(W) <- c("Z", "Y", "X")
  W_a <- .spatial_validate_W(W, units = c("X", "Y", "Z"))
  expect_equal(rownames(W_a), c("X", "Y", "Z"))
  expect_equal(colnames(W_a), c("X", "Y", "Z"))
})


test_that(".build_spatial_csa produces unit-specific CSAs", {
  set.seed(1)
  N <- 6; T_val <- 10
  units <- paste0("u", 1:N)
  df <- data.frame(
    id   = rep(units, each = T_val),
    t    = rep(1:T_val, N),
    y    = rnorm(N * T_val),
    x    = rnorm(N * T_val)
  )
  attr(df, "unit_var") <- "id"
  attr(df, "time_var") <- "t"

  # Nearest-neighbour ring
  W <- matrix(0, N, N)
  for (i in 1:N) {
    W[i, ifelse(i == 1, N, i - 1)] <- 1
    W[i, ifelse(i == N, 1, i + 1)] <- 1
  }
  rownames(W) <- colnames(W) <- units
  W_n <- .spatial_validate_W(W, units = units)

  out <- .build_spatial_csa(df, vars = c("y", "x"), W = W_n, lags = 0L)
  expect_true("csa_y" %in% names(out))
  expect_true("csa_x" %in% names(out))
  # Spatial CSA for a unit at time t should differ across units
  # (unlike the global mean, which is constant across units)
  t1_csa_y <- out$csa_y[out$t == 1]
  expect_equal(length(unique(round(t1_csa_y, 8))), N)
})


test_that("dcce() with spatial_weights runs and produces sensible output", {
  set.seed(2)
  N <- 15; T_val <- 25
  units <- paste0("u", 1:N)
  f <- rnorm(T_val)

  df <- data.frame(id = rep(units, each = T_val),
                   t  = rep(1:T_val, N))
  df$x <- NA_real_; df$y <- NA_real_
  for (i in 1:N) {
    idx <- ((i - 1) * T_val + 1):(i * T_val)
    lambda_i <- rnorm(1, 1, 0.3)
    df$x[idx] <- 0.4 * f + rnorm(T_val)
    df$y[idx] <- 1 + 0.5 * df$x[idx] + lambda_i * f + rnorm(T_val, sd = 0.3)
  }

  # Random row-stochastic weight matrix
  W <- matrix(runif(N * N), N, N)
  diag(W) <- 0
  W <- W / rowSums(W)
  rownames(W) <- colnames(W) <- units

  fit <- dcce(
    data       = df,
    unit_index = "id",
    time_index = "t",
    formula    = y ~ x,
    model      = "cce",
    cross_section_vars = ~ .,
    spatial_weights = W
  )
  expect_s3_class(fit, "dcce_cce_fit")
  expect_true(abs(coef(fit)["x"] - 0.5) < 0.25)
})


test_that("dcce() with spatial_weights matches classical CCE when W is uniform", {
  # If W[i,j] = 1/(N-1) for i != j and 0 on the diagonal, the spatial
  # CSA equals the leave-one-out mean which converges to the global
  # mean as N grows. Results should be numerically close to classical
  # CCE for a reasonably large N.
  set.seed(3)
  N <- 30; T_val <- 20
  units <- paste0("u", 1:N)
  f <- rnorm(T_val)
  df <- data.frame(id = rep(units, each = T_val),
                   t  = rep(1:T_val, N))
  df$x <- NA_real_; df$y <- NA_real_
  for (i in 1:N) {
    idx <- ((i - 1) * T_val + 1):(i * T_val)
    lambda_i <- rnorm(1, 1, 0.3)
    df$x[idx] <- 0.5 * f + rnorm(T_val)
    df$y[idx] <- 1 + 0.5 * df$x[idx] + lambda_i * f + rnorm(T_val, sd = 0.3)
  }

  W_unif <- matrix(1 / (N - 1), N, N)
  diag(W_unif) <- 0
  rownames(W_unif) <- colnames(W_unif) <- units

  fit_global  <- dcce(data = df, unit_index = "id", time_index = "t",
                      formula = y ~ x, model = "cce",
                      cross_section_vars = ~ .)
  fit_spatial <- dcce(data = df, unit_index = "id", time_index = "t",
                      formula = y ~ x, model = "cce",
                      cross_section_vars = ~ .,
                      spatial_weights = W_unif)

  # Difference should be O(1/N) — within 0.1 for N = 30
  expect_lt(abs(coef(fit_global)["x"] - coef(fit_spatial)["x"]), 0.1)
})
