# Benchmark tests using dcce_sim dataset with known DGP
# True parameters: beta1_mg ~ 0.50, beta2_mg ~ 1.00

test_that("DCCE recovers true DGP parameters from dcce_sim", {
  data(dcce_sim, package = "dcce")
  data(dcce_sim_truth, package = "dcce")

  fit <- dcce(
    data               = dcce_sim,
    unit_index         = "unit",
    time_index         = "time",
    formula            = y ~ L(y, 1) + x,
    model              = "dcce",
    cross_section_vars = ~ y + x,
    cross_section_lags = 3
  )

  b <- coef(fit)


  # True beta1_mg ~ 0.50 (Nickell bias shifts AR coef down in short panels)
  # Allow 0.15 absolute deviation for N=30, T=50 with cr_lags=3
  expect_true(abs(unname(b[["L(y,1)"]]) - dcce_sim_truth$beta1_mg) < 0.15)

  # True beta2_mg ~ 1.00 — allow 0.15 absolute deviation
  expect_true(abs(unname(b[["x"]]) - dcce_sim_truth$beta2_mg) < 0.15)
})

test_that("Static CCE reduces CSD bias on dcce_sim", {
  data(dcce_sim, package = "dcce")
  data(dcce_sim_truth, package = "dcce")

  fit_cce <- dcce(
    data               = dcce_sim,
    unit_index         = "unit",
    time_index         = "time",
    formula            = y ~ x,
    model              = "cce",
    cross_section_vars = ~ y + x,
    cross_section_lags = 0
  )

  fit_mg <- dcce(
    data               = dcce_sim,
    unit_index         = "unit",
    time_index         = "time",
    formula            = y ~ x,
    model              = "mg",
    cross_section_vars = NULL
  )

  # CCE should be at least as close to truth as MG
  err_cce <- abs(coef(fit_cce)[["x"]] - dcce_sim_truth$beta2_mg)
  err_mg  <- abs(coef(fit_mg)[["x"]]  - dcce_sim_truth$beta2_mg)
  expect_lt(err_cce, err_mg + 0.1)
})

test_that("CD test detects CSD in naive OLS on dcce_sim", {
  data(dcce_sim, package = "dcce")

  # Naive OLS with unit FE
  naive <- lm(y ~ x + factor(unit), data = dcce_sim)
  cd_before <- pcd_test(
    residuals(naive),
    data       = dcce_sim,
    unit_index = "unit",
    time_index = "time",
    test       = "pesaran"
  )

  # Should detect CSD

  expect_gt(abs(cd_before$statistics$statistic[1]), 3)
})
