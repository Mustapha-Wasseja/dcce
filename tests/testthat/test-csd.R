test_that("csd_exp with variable method on independent data", {
  set.seed(100)
  N <- 30
  T_val <- 50
  em <- matrix(rnorm(N * T_val), N, T_val)

  result <- csd_exp(em, use_residuals = FALSE)

  expect_s3_class(result, "dcce_csd")
  # For independent data, alpha should be close to 0.5
  expect_true(result$alpha < 1.0)
  expect_equal(result$method, "variable")
})

test_that("csd_exp with variable method on strongly dependent data", {
  set.seed(101)
  N <- 30
  T_val <- 100
  f <- rnorm(T_val)

  em <- matrix(NA, N, T_val)
  for (i in 1:N) {
    em[i, ] <- rnorm(1, 1, 0.2) * f + rnorm(T_val, sd = 0.3)
  }

  result <- csd_exp(em, use_residuals = FALSE)

  # With strong common factor, alpha should be close to 1
  expect_true(result$alpha > 0.7)
})

test_that("csd_exp with residual method", {
  set.seed(102)
  N <- 20
  T_val <- 50
  f <- rnorm(T_val)

  em <- matrix(NA, N, T_val)
  for (i in 1:N) {
    em[i, ] <- rnorm(1, 1, 0.5) * f + rnorm(T_val, sd = 0.5)
  }

  result <- csd_exp(em, use_residuals = TRUE, n_bootstrap = 50)

  expect_s3_class(result, "dcce_csd")
  expect_equal(result$method, "residual")
  expect_true(!is.na(result$se))
  expect_length(result$ci, 2)
})

test_that("csd_exp accepts dcce_fit object", {
  set.seed(103)
  N <- 15
  T_val <- 25
  df <- data.frame(
    id = rep(1:N, each = T_val),
    t  = rep(1:T_val, N),
    y  = rnorm(N * T_val),
    x  = rnorm(N * T_val)
  )

  fit <- dcce(y ~ x, data = df, unit_index = "id", time_index = "t", model = "mg", cross_section_vars = NULL)
  result <- csd_exp(fit, use_residuals = FALSE)

  expect_s3_class(result, "dcce_csd")
  expect_equal(result$N, N)
})

test_that("print.dcce_csd produces output", {
  set.seed(104)
  em <- matrix(rnorm(200), 10, 20)
  result <- csd_exp(em, use_residuals = FALSE)
  expect_output(print(result), "Exponent")
})
