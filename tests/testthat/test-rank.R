test_that("rank_condition runs on CCE model", {
  set.seed(400)
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

  result <- rank_condition(fit)
  expect_s3_class(result, "dcce_rank")
  expect_true(result$RC %in% c(0L, 1L))
})

test_that("rank_condition warns on dynamic models", {
  set.seed(401)
  N <- 15
  T_val <- 25
  df <- data.frame(
    id = rep(1:N, each = T_val),
    t  = rep(1:T_val, N),
    y  = rnorm(N * T_val),
    x  = rnorm(N * T_val)
  )
  fit <- dcce(y ~ x, data = df, unit_index = "id", time_index = "t",
              model = "dcce", cross_section_vars = ~ ., cross_section_lags = 2)

  expect_warning(rank_condition(fit), "only valid for static")
})

test_that("rank_condition with no CSAs returns RC=1", {
  set.seed(402)
  N <- 10
  T_val <- 20
  df <- data.frame(
    id = rep(1:N, each = T_val),
    t  = rep(1:T_val, N),
    y  = rnorm(N * T_val),
    x  = rnorm(N * T_val)
  )
  fit <- dcce(y ~ x, data = df, unit_index = "id", time_index = "t", model = "mg", cross_section_vars = NULL)

  result <- rank_condition(fit)
  expect_equal(result$RC, 1L)
})

test_that("print.dcce_rank produces output", {
  set.seed(403)
  N <- 10
  T_val <- 20
  df <- data.frame(
    id = rep(1:N, each = T_val),
    t  = rep(1:T_val, N),
    y  = rnorm(N * T_val),
    x  = rnorm(N * T_val)
  )
  fit <- dcce(y ~ x, data = df, unit_index = "id", time_index = "t",
              model = "cce", cross_section_vars = ~ ., cross_section_lags = 0)
  result <- rank_condition(fit)
  expect_output(print(result), "Rank Condition")
})
