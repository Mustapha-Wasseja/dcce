test_that("dcce_workflow returns dcce_workflow object", {
  data(pwt8)
  wf <- dcce_workflow(
    data       = pwt8,
    unit_index = "country",
    time_index = "year",
    formula    = log_rgdpo ~ log_hc + log_ck + log_ngd,
    verbose    = FALSE
  )
  expect_s3_class(wf, "dcce_workflow")
  expect_true(!is.null(wf$recommendation$model))
  expect_true(!is.null(wf$recommendation$suggested_call))
})


test_that("dcce_workflow recommendation is a valid model string", {
  data(pwt8)
  wf <- dcce_workflow(
    data = pwt8, unit_index = "country", time_index = "year",
    formula = log_rgdpo ~ log_hc + log_ck, verbose = FALSE
  )
  expect_true(wf$recommendation$model %in%
                c("mg", "cce", "dcce", "amg", "rcce", "pmg", "csdl", "csardl"))
})


test_that("dcce_workflow suggested_call parses as valid R", {
  data(pwt8)
  wf <- dcce_workflow(
    data = pwt8, unit_index = "country", time_index = "year",
    formula = log_rgdpo ~ log_hc + log_ck, verbose = FALSE
  )
  # The suggested call must parse without error
  expect_no_error(parse(text = wf$recommendation$suggested_call))
})


test_that("dcce_workflow verbose=TRUE prints progress messages", {
  data(pwt8)
  # cli::cli_h1() prints to stderr via message(), so use expect_message
  expect_message(
    dcce_workflow(data = pwt8, unit_index = "country", time_index = "year",
                  formula = log_rgdpo ~ log_hc, verbose = TRUE),
    regexp = "Diagnostic|Workflow|Panel"
  )
})


test_that("dcce_workflow fills all expected components", {
  data(pwt8)
  wf <- dcce_workflow(
    data = pwt8, unit_index = "country", time_index = "year",
    formula = log_rgdpo ~ log_hc + log_ck, verbose = FALSE
  )
  expect_true(!is.null(wf$panel_summary))
  expect_true(!is.null(wf$unit_root))
  expect_true(!is.null(wf$csd_premodel))
  expect_true(!is.null(wf$rank_condition))
  expect_true(!is.null(wf$optimal_cr_lags))
  expect_true(!is.null(wf$recommendation))
  # Panel summary
  expect_equal(wf$panel_summary$N, 93L)
})


test_that("dcce_workflow print method works", {
  data(pwt8)
  wf <- dcce_workflow(
    data = pwt8, unit_index = "country", time_index = "year",
    formula = log_rgdpo ~ log_hc, verbose = FALSE
  )
  expect_output(print(wf), "Workflow")
  expect_output(print(wf), "Suggested")
})
