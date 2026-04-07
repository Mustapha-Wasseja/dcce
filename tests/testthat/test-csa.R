test_that(".build_csa computes correct cross-sectional means", {
  # 3 units, 4 periods
  df <- data.frame(
    id = rep(1:3, each = 4),
    t  = rep(1:4, 3),
    y  = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12)
  )
  p <- .make_panel(df, unit_index = "id", time_index = "t")
  csa <- .build_csa(p, "y", lags = 0)

  # at t=1: mean(1,5,9)=5; at t=2: mean(2,6,10)=6; etc.
  expect_equal(csa$csa_y[csa$t == 1], rep(5, 3))
  expect_equal(csa$csa_y[csa$t == 2], rep(6, 3))
  expect_equal(csa$csa_y[csa$t == 3], rep(7, 3))
  expect_equal(csa$csa_y[csa$t == 4], rep(8, 3))
})

test_that(".build_csa with lags produces correct lagged CSAs", {
  df <- data.frame(
    id = rep(1:2, each = 5),
    t  = rep(1:5, 2),
    y  = c(2, 4, 6, 8, 10, 8, 6, 4, 2, 0)
  )
  p <- .make_panel(df, unit_index = "id", time_index = "t")
  csa <- .build_csa(p, "y", lags = 2)

  # Contemporaneous CSA: mean(2,8)=5, mean(4,6)=5, mean(6,4)=5, etc.
  expect_equal(csa$csa_y[csa$t == 1], rep(5, 2))
  expect_equal(csa$csa_y[csa$t == 2], rep(5, 2))

  # CSA lag 1: at t=1, L1 should be NA; at t=2, L1 = CSA at t=1 = 5
  expect_true(all(is.na(csa$csa_y_L1[csa$t == 1])))
  expect_equal(csa$csa_y_L1[csa$t == 2], rep(5, 2))

  # CSA lag 2: at t=1 and t=2 should be NA
  expect_true(all(is.na(csa$csa_y_L2[csa$t == 1])))
  expect_true(all(is.na(csa$csa_y_L2[csa$t == 2])))
  expect_equal(csa$csa_y_L2[csa$t == 3], rep(5, 2))
})

test_that(".build_csa with multiple variables", {
  df <- data.frame(
    id = rep(1:2, each = 3),
    t  = rep(1:3, 2),
    x  = c(10, 20, 30, 40, 50, 60),
    y  = c(1, 2, 3, 4, 5, 6)
  )
  p <- .make_panel(df, unit_index = "id", time_index = "t")
  csa <- .build_csa(p, c("x", "y"), lags = 0)

  expect_true("csa_x" %in% names(csa))
  expect_true("csa_y" %in% names(csa))
  # csa_x at t=1: mean(10,40)=25
  expect_equal(csa$csa_x[csa$t == 1], rep(25, 2))
  # csa_y at t=1: mean(1,4)=2.5
  expect_equal(csa$csa_y[csa$t == 1], rep(2.5, 2))
})

test_that(".build_csa handles unbalanced panels", {
  df <- data.frame(
    id = c(1, 1, 1, 2, 2),
    t  = c(1, 2, 3, 2, 3),
    y  = c(10, 20, 30, 40, 50)
  )
  p <- .make_panel(df, unit_index = "id", time_index = "t")
  csa <- .build_csa(p, "y", lags = 0)

  # t=1: only unit 1 observed -> mean = 10
  expect_equal(csa$csa_y[csa$t == 1], 10)
  # t=2: mean(20,40)=30
  expect_equal(csa$csa_y[csa$t == 2], rep(30, 2))
  # t=3: mean(30,50)=40
  expect_equal(csa$csa_y[csa$t == 3], rep(40, 2))
})

test_that(".build_csa per-variable lags", {
  df <- data.frame(
    id = rep(1:2, each = 4),
    t  = rep(1:4, 2),
    x  = c(1, 2, 3, 4, 5, 6, 7, 8),
    y  = c(10, 20, 30, 40, 50, 60, 70, 80)
  )
  p <- .make_panel(df, unit_index = "id", time_index = "t")
  csa <- .build_csa(p, c("x", "y"), lags = c(x = 1L, y = 2L))

  # x should have L1 but not L2

  expect_true("csa_x_L1" %in% names(csa))
  expect_false("csa_x_L2" %in% names(csa))

  # y should have L1 and L2
  expect_true("csa_y_L1" %in% names(csa))
  expect_true("csa_y_L2" %in% names(csa))
})

test_that(".build_csa preserves panel attributes", {
  df <- data.frame(id = rep(1:2, each = 3), t = rep(1:3, 2), y = 1:6)
  p <- .make_panel(df, unit_index = "id", time_index = "t")
  csa <- .build_csa(p, "y", lags = 0)
  expect_equal(attr(csa, "unit_var"), "id")
  expect_equal(attr(csa, "time_var"), "t")
})

test_that(".build_global_csa uses full_panel for means", {
  # Estimation sample: 2 units; full sample: 3 units
  df_est <- data.frame(
    id = rep(1:2, each = 3),
    t  = rep(1:3, 2),
    y  = c(10, 20, 30, 40, 50, 60)
  )
  df_full <- data.frame(
    id = rep(1:3, each = 3),
    t  = rep(1:3, 3),
    y  = c(10, 20, 30, 40, 50, 60, 70, 80, 90)
  )
  p_est <- .make_panel(df_est, unit_index = "id", time_index = "t")
  p_full <- .make_panel(df_full, unit_index = "id", time_index = "t")
  csa <- .build_global_csa(p_est, p_full, "y", lags = 0)

  # CSA should be from full panel: t=1: mean(10,40,70)=40
  expect_equal(csa$csa_y[csa$t == 1], rep(40, 2))
  # Not from est panel: mean(10,40)=25
  expect_false(all(csa$csa_y[csa$t == 1] == 25))
})

test_that(".build_cluster_csa computes group-specific CSAs", {
  df <- data.frame(
    id      = rep(1:4, each = 3),
    t       = rep(1:3, 4),
    y       = c(1, 2, 3, 4, 5, 6, 10, 20, 30, 40, 50, 60),
    cluster = rep(c("A", "A", "B", "B"), each = 3)
  )
  p <- .make_panel(df, unit_index = "id", time_index = "t")
  csa <- .build_cluster_csa(p, "y", lags = 0, by = "cluster")

  # Cluster A (ids 1,2): t=1 mean(1,4)=2.5
  idx_a_t1 <- which(csa$cluster == "A" & csa$t == 1)
  expect_equal(csa$csa_y[idx_a_t1], rep(2.5, 2))

  # Cluster B (ids 3,4): t=1 mean(10,40)=25
  idx_b_t1 <- which(csa$cluster == "B" & csa$t == 1)
  expect_equal(csa$csa_y[idx_b_t1], rep(25, 2))
})

test_that(".get_csa_colnames returns CSA columns", {
  df <- data.frame(id = rep(1:2, each = 3), t = rep(1:3, 2), y = 1:6)
  p <- .make_panel(df, unit_index = "id", time_index = "t")
  csa <- .build_csa(p, "y", lags = 1)
  cols <- .get_csa_colnames(csa)
  expect_true("csa_y" %in% cols)
  expect_true("csa_y_L1" %in% cols)
  expect_equal(length(cols), 2L)
})
