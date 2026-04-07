test_that("ic() computes IC1 and IC2", {
  set.seed(300)
  N <- 15
  T_val <- 25
  df <- data.frame(
    id = rep(1:N, each = T_val),
    t  = rep(1:T_val, N),
    y  = rnorm(N * T_val),
    x  = rnorm(N * T_val)
  )
  fit <- dcce(y ~ x, data = df, unit_index = "id", time_index = "t",
              model = "cce", cross_section_vars = ~ ., cross_section_lags = 0)

  result <- ic(fit)
  expect_s3_class(result, "dcce_ic")
  expect_true(is.finite(result$IC1))
  expect_true(is.finite(result$IC2))
  expect_true(is.na(result$PC1))  # no reference models
})

test_that("ic() with reference models computes PC1 and PC2", {
  set.seed(301)
  N <- 15
  T_val <- 25
  df <- data.frame(
    id = rep(1:N, each = T_val),
    t  = rep(1:T_val, N),
    y  = rnorm(N * T_val),
    x  = rnorm(N * T_val)
  )
  fit0 <- dcce(y ~ x, data = df, unit_index = "id", time_index = "t", model = "mg", cross_section_vars = NULL)
  fit1 <- dcce(y ~ x, data = df, unit_index = "id", time_index = "t",
               model = "cce", cross_section_vars = ~ ., cross_section_lags = 0)
  fit2 <- dcce(y ~ x, data = df, unit_index = "id", time_index = "t",
               model = "dcce", cross_section_vars = ~ ., cross_section_lags = 2)

  result <- ic(fit1, models = list(fit0, fit1, fit2))
  expect_true(is.finite(result$PC1))
  expect_true(is.finite(result$PC2))
})

test_that("print.dcce_ic produces output", {
  set.seed(302)
  df <- data.frame(
    id = rep(1:10, each = 20),
    t  = rep(1:20, 10),
    y  = rnorm(200),
    x  = rnorm(200)
  )
  fit <- dcce(y ~ x, data = df, unit_index = "id", time_index = "t", model = "mg", cross_section_vars = NULL)
  result <- ic(fit)
  expect_output(print(result), "Information Criteria")
})
