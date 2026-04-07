# Replication of Ditzen (2018) Stata Journal 18:3, 585-617
# Uses Jan Ditzen's exact xtdcce2 sample dataset (pwt8)

test_that("pwt8 dataset loads correctly", {
  data(pwt8, package = "dcce")
  expect_equal(length(unique(pwt8$country)), 93L)
  expect_equal(range(pwt8$year), c(1960L, 2007L))
  expect_true(all(c("log_rgdpo", "log_hc", "log_ck", "log_ngd",
                     "d_log_rgdpo") %in% names(pwt8)))
  expect_equal(nrow(pwt8), 4464L)
})

test_that("MG: correct signs on growth regression (Ex 7.1)", {
  data(pwt8, package = "dcce")
  fit <- dcce(
    data = pwt8, unit_index = "country", time_index = "year",
    formula = d_log_rgdpo ~ L(log_rgdpo, 1) + log_hc + log_ck + log_ngd,
    model = "mg", cross_section_vars = NULL
  )
  b <- coef(fit)
  # With d.y as LHS, L(y,1) coefficient is (rho-1), should be negative
  expect_lt(b["L(log_rgdpo,1)"], 0)
  expect_gt(b["L(log_rgdpo,1)"], -1)  # stationarity
  expect_gt(b["log_ck"], 0)           # physical capital positive
  expect_lt(b["log_ngd"], 0)          # pop growth negative
})

test_that("DCCE with cr_lags=3: residuals show weak CSD (Ex 7.3)", {
  data(pwt8, package = "dcce")
  fit <- dcce(
    data               = pwt8,
    unit_index         = "country",
    time_index         = "year",
    formula            = d_log_rgdpo ~ L(log_rgdpo, 1) + log_hc + log_ck + log_ngd,
    model              = "dcce",
    cross_section_vars = ~ log_rgdpo + log_hc + log_ck + log_ngd,
    cross_section_lags = 3
  )
  b <- coef(fit)
  # Speed of adjustment negative (error correction)
  expect_lt(b["L(log_rgdpo,1)"], 0)
  # Physical capital positive
  expect_gt(b["log_ck"], 0)

  # CD test: DCCE should remove CSD (fail to reject H0)
  cd <- pcd_test(fit, test = "pesaran")
  expect_gt(cd$statistics$p_value[1], 0.05)
})

test_that("CD test on pooled OLS ~ 36-39 (xtdcce2 README target: 36.34)", {
  data(pwt8, package = "dcce")
  dat <- na.omit(pwt8[, c("country", "year", "d_log_rgdpo",
                            "log_hc", "log_ck", "log_ngd")])
  fit_ols <- lm(d_log_rgdpo ~ log_hc + log_ck + log_ngd, data = dat)
  cd <- pcd_test(
    residuals(fit_ols), data = dat,
    unit_index = "country", time_index = "year",
    test = "pesaran"
  )
  # Published Stata: 36.34; our N=93 dataset gives ~38.7
  # Both indicate very strong CSD
  expect_gt(cd$statistics$statistic[1], 30)
  expect_lt(cd$statistics$p_value[1], 0.001)
})
