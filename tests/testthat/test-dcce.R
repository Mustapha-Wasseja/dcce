test_that("dcce() with model='dcce' and cr_lags > 0 runs correctly", {
  set.seed(500)
  N <- 20
  T_val <- 40

  # DGP with common factor
  f <- cumsum(rnorm(T_val))

  df <- data.frame(
    id = rep(1:N, each = T_val),
    t  = rep(1:T_val, N)
  )
  df$x <- NA_real_
  df$y <- NA_real_

  for (i in 1:N) {
    idx <- ((i - 1) * T_val + 1):(i * T_val)
    lambda_i <- rnorm(1, 1, 0.5)
    beta_i <- rnorm(1, 0.5, 0.1)
    df$x[idx] <- 0.5 * f + rnorm(T_val)
    df$y[idx] <- 1 + beta_i * df$x[idx] + lambda_i * f + rnorm(T_val, sd = 0.3)
  }

  fit <- dcce(y ~ x, data = df, unit_index = "id", time_index = "t",
              model = "dcce", cross_section_vars = ~ ., cross_section_lags = 3)

  expect_s3_class(fit, "dcce_dcce_fit")
  expect_s3_class(fit, "dcce_fit")

  # Should have CSA and lagged CSA columns
  expect_true(any(grepl("_L1$", fit$csa_colnames)))
  expect_true(any(grepl("_L3$", fit$csa_colnames)))

  # Coefficient on x should be roughly 0.5
  expect_equal(unname(coef(fit)["x"]), 0.5, tolerance = 0.25)
})

test_that("dcce() with per-variable cr_lags", {
  set.seed(501)
  N <- 15
  T_val <- 35

  df <- data.frame(
    id = rep(1:N, each = T_val),
    t  = rep(1:T_val, N),
    x  = rnorm(N * T_val),
    y  = rnorm(N * T_val)
  )

  fit <- dcce(y ~ x, data = df, unit_index = "id", time_index = "t",
              model = "dcce", cross_section_vars = ~ x + y,
              cross_section_lags = c(x = 1L, y = 2L))

  expect_s3_class(fit, "dcce_dcce_fit")
  # Check CSA columns
  expect_true("csa_x_L1" %in% fit$csa_colnames)
  expect_false("csa_x_L2" %in% fit$csa_colnames)
  expect_true("csa_y_L2" %in% fit$csa_colnames)
})

test_that("dcce() with jackknife bias correction", {
  set.seed(502)
  N <- 20
  T_val <- 40

  df <- data.frame(
    id = rep(1:N, each = T_val),
    t  = rep(1:T_val, N),
    x  = rnorm(N * T_val),
    y  = rnorm(N * T_val)
  )

  fit_plain <- dcce(y ~ x, data = df, unit_index = "id", time_index = "t",
                    model = "dcce", cross_section_vars = ~ ., cross_section_lags = 2, bias_correction = "none")
  fit_jk    <- dcce(y ~ x, data = df, unit_index = "id", time_index = "t",
                    model = "dcce", cross_section_vars = ~ ., cross_section_lags = 2, bias_correction = "jackknife")

  expect_s3_class(fit_jk, "dcce_dcce_fit")

  # Jackknife should produce different coefficients
  expect_false(isTRUE(all.equal(coef(fit_plain), coef(fit_jk), tolerance = 1e-6)))
})

test_that("dcce() with model='dcce' and lagged dependent variable", {
  set.seed(503)
  N <- 20
  T_val <- 40
  f <- rnorm(T_val)

  df <- data.frame(
    id = rep(1:N, each = T_val),
    t  = rep(1:T_val, N)
  )
  df$y <- NA_real_
  df$x <- NA_real_

  for (i in 1:N) {
    idx <- ((i - 1) * T_val + 1):(i * T_val)
    lambda_i <- rnorm(1, 1, 0.3)
    x_i <- rnorm(T_val)
    y_i <- numeric(T_val)
    y_i[1] <- rnorm(1)
    for (tt in 2:T_val) {
      y_i[tt] <- 0.5 * y_i[tt - 1] + 0.3 * x_i[tt] + lambda_i * f[tt] + rnorm(1, sd = 0.3)
    }
    df$x[idx] <- x_i
    df$y[idx] <- y_i
  }

  fit <- dcce(y ~ L(y, 1) + x, data = df, unit_index = "id", time_index = "t",
              model = "dcce", cross_section_vars = ~ ., cross_section_lags = 3)

  expect_s3_class(fit, "dcce_dcce_fit")
  expect_true("L(y,1)" %in% names(coef(fit)))
})

test_that(".check_csa_rank warns on rank deficiency", {
  # Create a rank-deficient matrix
  Z <- cbind(1:10, 1:10, rnorm(10))
  expect_warning(.check_csa_rank(Z, "test_unit"), "rank deficient")

  # Full rank matrix should not warn
  Z_ok <- cbind(rnorm(10), rnorm(10))
  expect_silent(.check_csa_rank(Z_ok, "test_unit"))
})
