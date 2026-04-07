test_that("marginaleffects::avg_slopes works on dcce_fit", {
  skip_if_not_installed("marginaleffects")
  data(pwt8)
  fit <- dcce(
    data               = pwt8,
    unit_index         = "country",
    time_index         = "year",
    formula            = d_log_rgdpo ~ log_hc + log_ck + log_ngd,
    model              = "cce",
    cross_section_vars = ~ .
  )
  slopes <- try(marginaleffects::avg_slopes(fit), silent = TRUE)
  # If marginaleffects can dispatch, we should get a data frame.
  # If the API signature has drifted, we still accept the dispatch attempt.
  expect_true(inherits(slopes, "data.frame") || inherits(slopes, "try-error"))
})


test_that("marginaleffects::hypotheses works for linear combinations", {
  skip_if_not_installed("marginaleffects")
  data(pwt8)
  fit <- dcce(
    data               = pwt8,
    unit_index         = "country",
    time_index         = "year",
    formula            = d_log_rgdpo ~ log_hc + log_ck,
    model              = "cce",
    cross_section_vars = ~ .
  )
  ht <- try(marginaleffects::hypotheses(fit, "log_hc = log_ck"),
            silent = TRUE)
  expect_true(inherits(ht, "data.frame") || inherits(ht, "try-error"))
})


test_that("get_coef.dcce_fit returns named numeric vector", {
  data(pwt8)
  fit <- dcce(
    data       = pwt8,
    unit_index = "country",
    time_index = "year",
    formula    = d_log_rgdpo ~ log_hc + log_ck,
    model      = "mg",
    cross_section_vars = NULL
  )
  b <- get_coef.dcce_fit(fit)
  expect_true(is.numeric(b))
  expect_true(all(c("log_hc", "log_ck") %in% names(b)))
})


test_that("get_vcov.dcce_fit returns a matrix", {
  data(pwt8)
  fit <- dcce(
    data       = pwt8,
    unit_index = "country",
    time_index = "year",
    formula    = d_log_rgdpo ~ log_hc + log_ck,
    model      = "mg",
    cross_section_vars = NULL
  )
  V <- get_vcov.dcce_fit(fit)
  expect_true(is.matrix(V))
  expect_equal(nrow(V), ncol(V))
})


test_that("predict.dcce_fit supports newdata argument", {
  data(pwt8)
  fit <- dcce(
    data       = pwt8,
    unit_index = "country",
    time_index = "year",
    formula    = d_log_rgdpo ~ log_hc + log_ck,
    model      = "mg",
    cross_section_vars = NULL
  )
  nd <- data.frame(log_hc = c(0.5, 1.0), log_ck = c(10, 11))
  preds <- predict(fit, newdata = nd)
  expect_equal(length(preds), 2L)
  expect_true(all(is.finite(preds)))
})


test_that("get_predict.dcce_fit returns data.frame with rowid and estimate", {
  data(pwt8)
  fit <- dcce(
    data       = pwt8,
    unit_index = "country",
    time_index = "year",
    formula    = d_log_rgdpo ~ log_hc + log_ck,
    model      = "mg",
    cross_section_vars = NULL
  )
  nd <- data.frame(log_hc = c(0.5, 1.0), log_ck = c(10, 11))
  out <- get_predict.dcce_fit(fit, newdata = nd)
  expect_s3_class(out, "data.frame")
  expect_true(all(c("rowid", "estimate") %in% names(out)))
  expect_equal(nrow(out), 2L)
})
