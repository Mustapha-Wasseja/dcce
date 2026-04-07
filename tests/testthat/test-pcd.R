test_that("pcd_test with CD statistic on independent residuals", {
  set.seed(42)
  N <- 20
  T_val <- 30
  df <- data.frame(
    id = rep(1:N, each = T_val),
    t  = rep(1:T_val, N),
    e  = rnorm(N * T_val)
  )

  cd <- pcd_test(df$e, data = df, unit_index = "id", time_index = "t", test = "pesaran")

  expect_s3_class(cd, "dcce_cd")
  expect_equal(cd$N, N)
  # Under null (no CSD), CD should be small (|CD| < ~3)
  expect_true(abs(cd$statistics$statistic[1]) < 5)
  expect_true(cd$statistics$p_value[1] > 0.01)
})

test_that("pcd_test detects strong cross-sectional dependence", {
  set.seed(43)
  N <- 20
  T_val <- 50
  # Common factor creates strong CSD
  f <- rnorm(T_val)

  e <- matrix(NA_real_, N, T_val)
  for (i in 1:N) {
    lambda_i <- rnorm(1, 1, 0.3)
    e[i, ] <- lambda_i * f + rnorm(T_val, sd = 0.3)
  }

  cd <- pcd_test(e, test = "pesaran")

  expect_s3_class(cd, "dcce_cd")
  # Should detect CSD: large |CD| statistic
  expect_true(abs(cd$statistics$statistic[1]) > 5)
  expect_true(cd$statistics$p_value[1] < 0.01)
})

test_that("pcd_test accepts matrix input", {
  set.seed(44)
  em <- matrix(rnorm(100), 10, 10)
  cd <- pcd_test(em, test = "pesaran")
  expect_s3_class(cd, "dcce_cd")
  expect_equal(cd$N, 10)
})

test_that("pcd_test with CDw test", {
  set.seed(45)
  N <- 15
  T_val <- 25
  em <- matrix(rnorm(N * T_val), N, T_val)

  cd <- pcd_test(em, test = "cdw", n_reps = 100)

  expect_s3_class(cd, "dcce_cd")
  expect_equal(cd$statistics$test[1], "cdw")
  # Under null, CDw should be small
  expect_true(abs(cd$statistics$statistic[1]) < 5)
})

test_that("pcd_test with all tests", {
  set.seed(46)
  em <- matrix(rnorm(200), 10, 20)

  cd <- pcd_test(em, test = c("pesaran", "cdw", "pea", "cdstar"), n_reps = 50)

  expect_equal(nrow(cd$statistics), 4)
  expect_equal(cd$statistics$test, c("pesaran", "cdw", "pea", "cdstar"))
})

test_that("pcd_test accepts dcce_fit object", {
  set.seed(47)
  N <- 15
  T_val <- 20
  df <- data.frame(
    id = rep(1:N, each = T_val),
    t  = rep(1:T_val, N),
    y  = rnorm(N * T_val),
    x  = rnorm(N * T_val)
  )

  fit <- dcce(y ~ x, data = df, unit_index = "id", time_index = "t", model = "mg", cross_section_vars = NULL)
  cd <- pcd_test(fit, test = "pesaran")

  expect_s3_class(cd, "dcce_cd")
  expect_equal(cd$N, N)
})

test_that("print.dcce_cd produces output", {
  set.seed(48)
  em <- matrix(rnorm(100), 10, 10)
  cd <- pcd_test(em, test = "pesaran")
  expect_output(print(cd), "Cross-Sectional Dependence")
})
