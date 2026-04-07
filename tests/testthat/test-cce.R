test_that(".partial_out correctly projects out Z", {
  set.seed(42)
  n <- 50
  Z <- cbind(1, rnorm(n))
  X <- cbind(rnorm(n), rnorm(n))

  MX <- .partial_out(X, Z)

  # Residuals should be orthogonal to Z
  expect_true(max(abs(crossprod(Z, MX))) < 1e-10)
})

test_that("dcce() with model='cce' runs correctly", {
  set.seed(100)
  N <- 20
  T_val <- 30
  # Common factor
  f <- rnorm(T_val)

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
              model = "cce", cross_section_vars = ~ ., cross_section_lags = 0)

  expect_s3_class(fit, "dcce_cce_fit")
  expect_s3_class(fit, "dcce_fit")

  # Should have CSA columns
  expect_true(length(fit$csa_colnames) > 0)

  # Coefficient on x should be roughly 0.5
  expect_equal(unname(coef(fit)["x"]), 0.5, tolerance = 0.2)
})

test_that("CCE with cr_lags=0 produces different results from MG without CSAs", {
  set.seed(200)
  N <- 15
  T_val <- 25
  f <- rnorm(T_val)

  df <- data.frame(
    id = rep(1:N, each = T_val),
    t  = rep(1:T_val, N)
  )
  for (i in 1:N) {
    idx <- ((i - 1) * T_val + 1):(i * T_val)
    lambda_i <- rnorm(1, 2, 1)
    df$x[idx] <- f + rnorm(T_val)
    df$y[idx] <- 1 + 0.5 * df$x[idx] + lambda_i * f + rnorm(T_val, sd = 0.5)
  }

  fit_mg  <- dcce(y ~ x, data = df, unit_index = "id", time_index = "t", model = "mg", cross_section_vars = NULL)
  fit_cce <- dcce(y ~ x, data = df, unit_index = "id", time_index = "t", model = "cce", cross_section_vars = ~ ., cross_section_lags = 0)

  # They should give different coefficients
  expect_false(isTRUE(all.equal(coef(fit_mg)["x"], coef(fit_cce)["x"],
                                tolerance = 0.001)))
})

test_that("dcce() with model='cce' and multiple regressors", {
  set.seed(300)
  N <- 15
  T_val <- 30
  f <- rnorm(T_val)

  df <- data.frame(
    id = rep(1:N, each = T_val),
    t  = rep(1:T_val, N)
  )
  df$x1 <- NA_real_
  df$x2 <- NA_real_
  df$y  <- NA_real_

  for (i in 1:N) {
    idx <- ((i - 1) * T_val + 1):(i * T_val)
    lambda_i <- rnorm(1, 1, 0.5)
    df$x1[idx] <- 0.3 * f + rnorm(T_val)
    df$x2[idx] <- 0.5 * f + rnorm(T_val)
    df$y[idx] <- 1 + 0.4 * df$x1[idx] + 0.6 * df$x2[idx] + lambda_i * f + rnorm(T_val, sd = 0.3)
  }

  fit <- dcce(y ~ x1 + x2, data = df, unit_index = "id", time_index = "t",
              model = "cce", cross_section_vars = ~ ., cross_section_lags = 0)

  expect_s3_class(fit, "dcce_cce_fit")
  expect_equal(unname(coef(fit)["x1"]), 0.4, tolerance = 0.2)
  expect_equal(unname(coef(fit)["x2"]), 0.6, tolerance = 0.2)
})

test_that(".pooled_vcov_pesaran computes correct variance", {
  coef_list <- list(
    c(x = 1),
    c(x = 2),
    c(x = 3)
  )
  b_pooled <- c(x = 2)
  V <- .pooled_vcov_pesaran(coef_list, b_pooled)

  # (1/9) * [(1-2)^2 + (2-2)^2 + (3-2)^2] = 2/9
  expect_equal(V[1, 1], 2 / 9, tolerance = 1e-10)
})
