test_that("cointegration_test returns dcce_cointegration object", {
  data(pwt8)
  result <- cointegration_test(
    data       = pwt8,
    unit_index = "country",
    time_index = "year",
    formula    = log_rgdpo ~ log_hc + log_ck + log_ngd,
    test       = c("ga", "gt"),
    lags       = 1L
  )
  expect_s3_class(result, "dcce_cointegration")
  expect_true(
    all(c("test", "statistic", "p_value", "method") %in%
          names(result$statistics))
  )
  expect_equal(nrow(result$statistics), 2L)
  expect_equal(result$statistics$test, c("Gt", "Ga"))
})


test_that("cointegration_test runs all four statistics", {
  data(pwt8)
  result <- cointegration_test(
    data       = pwt8,
    unit_index = "country",
    time_index = "year",
    formula    = log_rgdpo ~ log_hc + log_ck,
    test       = c("ga", "gt", "pa", "pt"),
    lags       = 1L
  )
  expect_equal(nrow(result$statistics), 4L)
  expect_setequal(result$statistics$test, c("Gt", "Ga", "Pt", "Pa"))
  expect_true(all(is.finite(result$statistics$statistic)))
  expect_true(all(result$statistics$p_value >= 0 &
                    result$statistics$p_value <= 1))
})


test_that("cointegration_test Ga/Gt are negative for pwt8 growth data", {
  data(pwt8)
  result <- cointegration_test(
    data       = pwt8,
    unit_index = "country",
    time_index = "year",
    formula    = log_rgdpo ~ log_hc + log_ck,
    test       = c("ga", "gt"),
    lags       = 1L
  )
  expect_lt(result$statistics$statistic[result$statistics$test == "Gt"], 0)
  expect_lt(result$statistics$statistic[result$statistics$test == "Ga"], 0)
})


test_that("cointegration_test validates test argument", {
  data(pwt8)
  expect_error(
    cointegration_test(
      data = pwt8, unit_index = "country", time_index = "year",
      formula = log_rgdpo ~ log_hc, test = "bogus"
    ),
    regexp = "bogus"
  )
})


test_that("cointegration_test print method works", {
  data(pwt8)
  result <- cointegration_test(
    data = pwt8, unit_index = "country", time_index = "year",
    formula = log_rgdpo ~ log_hc, test = "gt", lags = 1L
  )
  expect_output(print(result), "Westerlund")
  expect_output(print(result), "Gt")
})


test_that("cointegration_test bootstrap runs and returns valid p-values", {
  data(pwt8)
  result <- cointegration_test(
    data        = pwt8,
    unit_index  = "country",
    time_index  = "year",
    formula     = log_rgdpo ~ log_hc,
    test        = "gt",
    lags        = 1L,
    n_bootstrap = 9L,   # tiny — just exercise the code path
    seed        = 42L
  )
  expect_equal(result$statistics$method, "bootstrap")
  expect_true(result$statistics$p_value >= 0 &&
                result$statistics$p_value <= 1)
})
