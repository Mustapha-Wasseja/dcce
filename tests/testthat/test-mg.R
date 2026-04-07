test_that(".unit_ols returns correct OLS estimates", {
  set.seed(42)
  n <- 100
  x <- rnorm(n)
  y <- 2 + 3 * x + rnorm(n, sd = 0.5)
  X <- cbind(`(Intercept)` = 1, x = x)

  fit <- .unit_ols(y, X)

  expect_length(fit$b, 2)
  expect_equal(unname(fit$b["(Intercept)"]), 2, tolerance = 0.3)
  expect_equal(unname(fit$b["x"]), 3, tolerance = 0.3)
  expect_length(fit$e, n)
  expect_equal(fit$df_resid, n - 2)
  expect_true(fit$r2 > 0.8)
})

test_that(".mg_aggregate computes correct mean", {
  coef_list <- list(
    c(a = 1, b = 2),
    c(a = 3, b = 4),
    c(a = 5, b = 6)
  )
  b_mg <- .mg_aggregate(coef_list)
  expect_equal(b_mg, c(a = 3, b = 4))
})

test_that(".mg_variance computes correct MG variance", {
  coef_list <- list(
    c(a = 1, b = 2),
    c(a = 3, b = 4),
    c(a = 5, b = 6)
  )
  b_mg <- c(a = 3, b = 4)
  V <- .mg_variance(coef_list, b_mg)

  # V = (1/9) * [(1-3)(1-3)' + (3-3)(3-3)' + (5-3)(5-3)']
  # = (1/9) * [c(4,4;4,4) + c(0,0;0,0) + c(4,4;4,4)]
  # = (1/9) * c(8,8;8,8) = c(8/9, 8/9; 8/9, 8/9)
  expect_equal(V[1, 1], 8 / 9, tolerance = 1e-10)
  expect_equal(V[2, 2], 8 / 9, tolerance = 1e-10)
  expect_equal(V[1, 2], 8 / 9, tolerance = 1e-10)
})

test_that("dcce() with model='mg' runs on synthetic data", {
  set.seed(123)
  N <- 20
  T_val <- 30
  beta <- 0.5
  df <- data.frame(
    id = rep(1:N, each = T_val),
    t  = rep(1:T_val, N)
  )
  # Heterogeneous coefficients centered on beta
  beta_i <- rnorm(N, mean = beta, sd = 0.1)
  df$x <- rnorm(N * T_val)
  df$y <- NA_real_
  for (i in 1:N) {
    idx <- ((i - 1) * T_val + 1):(i * T_val)
    df$y[idx] <- 1 + beta_i[i] * df$x[idx] + rnorm(T_val, sd = 0.5)
  }

  fit <- dcce(y ~ x, data = df, unit_index = "id", time_index = "t", model = "mg", cross_section_vars = NULL)

  expect_s3_class(fit, "dcce_mg_fit")
  expect_s3_class(fit, "dcce_fit")
  expect_equal(fit$N, N)

  # MG coefficient should be close to 0.5
  expect_equal(unname(coef(fit)["x"]), beta, tolerance = 0.15)

  # Check that methods work
  expect_length(coef(fit), 2)  # intercept + x
  expect_true(is.matrix(vcov(fit)))
  expect_true(length(residuals(fit)) > 0)
})

test_that("dcce() with model='mg' handles formula with L()", {
  set.seed(456)
  N <- 15
  T_val <- 25
  df <- data.frame(
    id = rep(1:N, each = T_val),
    t  = rep(1:T_val, N),
    y  = rnorm(N * T_val),
    x  = rnorm(N * T_val)
  )

  fit <- dcce(y ~ L(x, 1), data = df, unit_index = "id", time_index = "t",
              model = "mg", cross_section_vars = NULL)

  expect_s3_class(fit, "dcce_mg_fit")
  expect_true("L(x,1)" %in% names(coef(fit)))
})

test_that("dcce() drops units with insufficient observations", {
  # Create panel with one unit having only 2 observations
  df <- data.frame(
    id = c(rep(1, 3), rep(2, 20), rep(3, 20)),
    t  = c(1:3, 1:20, 1:20),
    y  = rnorm(43),
    x1 = rnorm(43),
    x2 = rnorm(43),
    x3 = rnorm(43)
  )

  expect_warning(
    fit <- dcce(y ~ x1 + x2 + x3, data = df, unit_index = "id", time_index = "t",
                model = "mg", cross_section_vars = NULL),
    "Dropped"
  )
  expect_equal(fit$N, 2)
})

test_that("dcce() with model='mg' and trend", {
  set.seed(789)
  N <- 10
  T_val <- 30
  df <- data.frame(
    id = rep(1:N, each = T_val),
    t  = rep(1:T_val, N),
    y  = rnorm(N * T_val),
    x  = rnorm(N * T_val)
  )

  fit <- dcce(y ~ x, data = df, unit_index = "id", time_index = "t",
              model = "mg", cross_section_vars = NULL, unit_trend = TRUE)

  expect_true("(trend)" %in% names(coef(fit)))
})

test_that("print.dcce_fit produces output", {
  set.seed(111)
  N <- 10
  T_val <- 20
  df <- data.frame(
    id = rep(1:N, each = T_val),
    t  = rep(1:T_val, N),
    y  = rnorm(N * T_val),
    x  = rnorm(N * T_val)
  )
  fit <- dcce(y ~ x, data = df, unit_index = "id", time_index = "t", model = "mg", cross_section_vars = NULL)
  expect_output(print(fit), "Mean Group")
})
