test_that("dcce_bootstrap() alias matches bootstrap() exactly", {
  set.seed(42)
  df <- data.frame(
    id = rep(1:8, each = 25),
    t  = rep(1:25, 8),
    y  = rnorm(200),
    x  = rnorm(200)
  )
  fit <- dcce(data = df, unit_index = "id", time_index = "t",
              formula = y ~ x, model = "mg", cross_section_vars = NULL)

  b1 <- bootstrap(fit,       type = "crosssection", reps = 20, seed = 123)
  b2 <- dcce_bootstrap(fit,  type = "crosssection", reps = 20, seed = 123)

  expect_s3_class(b1, "dcce_boot")
  expect_s3_class(b2, "dcce_boot")
  expect_equal(b1$se_boot, b2$se_boot, tolerance = 1e-12)
})


test_that("dcce_bootstrap is exported from the dcce namespace", {
  expect_true("dcce_bootstrap" %in% getNamespaceExports("dcce"))
  expect_true(is.function(dcce::dcce_bootstrap))
})


test_that("bootstrap crosssection produces valid output", {
  set.seed(42)
  N <- 15
  T_val <- 25
  df <- data.frame(
    id = rep(1:N, each = T_val),
    t  = rep(1:T_val, N),
    y  = rnorm(N * T_val),
    x  = rnorm(N * T_val)
  )
  fit <- dcce(y ~ x, data = df, unit_index = "id", time_index = "t", model = "mg", cross_section_vars = NULL)
  boot_res <- bootstrap(fit, type = "crosssection", reps = 50, seed = 123)

  expect_s3_class(boot_res, "dcce_boot")
  expect_equal(boot_res$type, "crosssection")
  expect_equal(boot_res$reps, 50)
  expect_equal(nrow(boot_res$b_boot), 50)
  expect_equal(ncol(boot_res$b_boot), length(coef(fit)))
  expect_true(all(boot_res$se_boot > 0))
  expect_true(all(boot_res$ci_lower < boot_res$ci_upper))
})

test_that("bootstrap wild produces valid output", {
  set.seed(43)
  N <- 15
  T_val <- 25
  df <- data.frame(
    id = rep(1:N, each = T_val),
    t  = rep(1:T_val, N),
    y  = rnorm(N * T_val),
    x  = rnorm(N * T_val)
  )
  fit <- dcce(y ~ x, data = df, unit_index = "id", time_index = "t", model = "mg", cross_section_vars = NULL)
  boot_res <- bootstrap(fit, type = "wild", reps = 50, seed = 456)

  expect_s3_class(boot_res, "dcce_boot")
  expect_equal(boot_res$type, "wild")
  expect_true(all(boot_res$se_boot > 0))
})

test_that("bootstrap works with CCE model", {
  set.seed(44)
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
  boot_res <- bootstrap(fit, reps = 30, seed = 789)

  expect_s3_class(boot_res, "dcce_boot")
  expect_equal(nrow(boot_res$b_boot), 30)
})

test_that("bootstrap seed reproducibility", {
  set.seed(45)
  N <- 10
  T_val <- 20
  df <- data.frame(
    id = rep(1:N, each = T_val),
    t  = rep(1:T_val, N),
    y  = rnorm(N * T_val),
    x  = rnorm(N * T_val)
  )
  fit <- dcce(y ~ x, data = df, unit_index = "id", time_index = "t", model = "mg", cross_section_vars = NULL)

  b1 <- bootstrap(fit, reps = 20, seed = 999)
  b2 <- bootstrap(fit, reps = 20, seed = 999)

  expect_equal(b1$b_boot, b2$b_boot)
})

test_that("print.dcce_boot produces output", {
  set.seed(46)
  df <- data.frame(
    id = rep(1:10, each = 20),
    t  = rep(1:20, 10),
    y  = rnorm(200),
    x  = rnorm(200)
  )
  fit <- dcce(y ~ x, data = df, unit_index = "id", time_index = "t", model = "mg", cross_section_vars = NULL)
  boot_res <- bootstrap(fit, reps = 20, seed = 111)
  expect_output(print(boot_res), "Bootstrap")
})
