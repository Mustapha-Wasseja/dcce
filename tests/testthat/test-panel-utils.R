test_that(".make_panel validates input correctly", {
  df <- data.frame(id = rep(1:3, each = 4), t = rep(1:4, 3), y = 1:12)


  # Valid input
  p <- .make_panel(df, unit_index = "id", time_index = "t")
  expect_s3_class(p, "data.frame")
  expect_equal(attr(p, "unit_var"), "id")
  expect_equal(attr(p, "time_var"), "t")

  # Missing column

  expect_error(.make_panel(df, unit_index = "bad", time_index = "t"), "not found")

  # Duplicate observations
  df_dup <- rbind(df, df[1, ])
  expect_error(.make_panel(df_dup, unit_index = "id", time_index = "t"), "Duplicate")
})

test_that(".make_panel sorts by unit then time", {
  df <- data.frame(
    id = c(2, 2, 1, 1, 3, 3),
    t  = c(2, 1, 2, 1, 1, 2),
    y  = c(4, 3, 2, 1, 5, 6)
  )
  p <- .make_panel(df, unit_index = "id", time_index = "t")
  expect_equal(p$id, c(1, 1, 2, 2, 3, 3))
  expect_equal(p$t,  c(1, 2, 1, 2, 1, 2))
  expect_equal(p$y,  c(1, 2, 3, 4, 5, 6))
})

test_that(".check_balance works for balanced and unbalanced panels", {
  # Balanced
  df <- data.frame(id = rep(1:3, each = 4), t = rep(1:4, 3), y = 1:12)
  p <- .make_panel(df, unit_index = "id", time_index = "t")
  bal <- .check_balance(p)
  expect_true(bal$is_balanced)
  expect_equal(bal$N, 3L)
  expect_equal(bal$T_min, 4L)
  expect_equal(bal$T_max, 4L)

  # Unbalanced
  df2 <- df[-c(4, 8), ]
  p2 <- .make_panel(df2, unit_index = "id", time_index = "t")
  bal2 <- .check_balance(p2)
  expect_false(bal2$is_balanced)
  expect_equal(bal2$T_min, 3L)
  expect_equal(bal2$T_max, 4L)
})

test_that(".lag_panel — balanced panel, no cross-unit contamination", {
  df <- data.frame(id = rep(1:3, each = 4), t = rep(1:4, 3), y = 1:12)
  p <- .make_panel(df, unit_index = "id", time_index = "t")
  ly <- .lag_panel(p, "y", 1)

  # First obs of each unit is NA
  expect_true(is.na(ly[1]))   # unit 1, t=1
  expect_true(is.na(ly[5]))   # unit 2, t=1
  expect_true(is.na(ly[9]))   # unit 3, t=1

  # Second obs of unit 1 = y[t=1] = 1
  expect_equal(ly[2], 1)
  # Third obs of unit 1 = y[t=2] = 2
  expect_equal(ly[3], 2)
  # Second obs of unit 2 = y[t=1 of unit 2] = 5
  expect_equal(ly[6], 5)
})

test_that(".lag_panel — multiple lags", {
  df <- data.frame(id = rep(1:2, each = 5), t = rep(1:5, 2), y = 1:10)
  p <- .make_panel(df, unit_index = "id", time_index = "t")
  ly2 <- .lag_panel(p, "y", 2)

  expect_true(is.na(ly2[1]))   # unit 1, t=1
  expect_true(is.na(ly2[2]))   # unit 1, t=2
  expect_equal(ly2[3], 1)      # unit 1, t=3 -> y[t=1]=1
  expect_true(is.na(ly2[6]))   # unit 2, t=1
  expect_true(is.na(ly2[7]))   # unit 2, t=2
  expect_equal(ly2[8], 6)      # unit 2, t=3 -> y[t=1 of unit 2]=6
})

test_that(".lag_panel — multiple variables", {
  df <- data.frame(id = rep(1:2, each = 3), t = rep(1:3, 2),
                   x = c(10, 20, 30, 40, 50, 60),
                   y = c(1, 2, 3, 4, 5, 6))
  p <- .make_panel(df, unit_index = "id", time_index = "t")
  result <- .lag_panel(p, c("x", "y"), 1)
  expect_s3_class(result, "data.frame")
  expect_equal(ncol(result), 2L)
  expect_true(is.na(result[1, 1]))
  expect_equal(result[2, 1], 10)
  expect_equal(result[2, 2], 1)
})

test_that(".lag_panel — unbalanced panel", {
  df <- data.frame(
    id = c(1, 1, 1, 2, 2),
    t  = c(1, 2, 3, 2, 3),
    y  = c(10, 20, 30, 50, 60)
  )
  p <- .make_panel(df, unit_index = "id", time_index = "t")
  ly <- .lag_panel(p, "y", 1)

  expect_true(is.na(ly[1]))    # unit 1, first obs
  expect_equal(ly[2], 10)      # unit 1, t=2 -> y[t=1]=10
  expect_equal(ly[3], 20)      # unit 1, t=3 -> y[t=2]=20
  expect_true(is.na(ly[4]))    # unit 2, first obs
  expect_equal(ly[5], 50)      # unit 2, t=3 -> y[t=2]=50
})

test_that(".diff_panel works correctly", {
  df <- data.frame(id = rep(1:2, each = 4), t = rep(1:4, 2),
                   y = c(1, 3, 6, 10, 5, 7, 12, 20))
  p <- .make_panel(df, unit_index = "id", time_index = "t")
  dy <- .diff_panel(p, "y", 1)

  expect_true(is.na(dy[1]))    # first obs unit 1
  expect_equal(dy[2], 2)       # 3 - 1
  expect_equal(dy[3], 3)       # 6 - 3
  expect_equal(dy[4], 4)       # 10 - 6
  expect_true(is.na(dy[5]))    # first obs unit 2
  expect_equal(dy[6], 2)       # 7 - 5
})

test_that("L() works as standalone function", {
  x <- c(10, 20, 30, 40, 50)
  lx <- L(x, 1)
  expect_true(is.na(lx[1]))
  expect_equal(lx[2], 10)
  expect_equal(lx[5], 40)

  # Lag 0 = identity
  expect_equal(L(x, 0), x)

  # Lead (negative lag)
  fx <- L(x, -1)
  expect_equal(fx[1], 20)
  expect_true(is.na(fx[5]))
})

test_that("D() works as standalone function", {
  x <- c(1, 3, 6, 10)
  dx <- D(x, 1)
  expect_true(is.na(dx[1]))
  expect_equal(dx[2], 2)
  expect_equal(dx[3], 3)
  expect_equal(dx[4], 4)
})

test_that("Lrange() returns matrix of lags", {
  x <- c(10, 20, 30, 40, 50)
  m <- Lrange(x, 0, 2)
  expect_equal(ncol(m), 3L)
  expect_equal(colnames(m), c("L0", "L1", "L2"))
  # L0 = identity
  expect_equal(m[, 1], x)
  # L1
  expect_true(is.na(m[1, 2]))
  expect_equal(unname(m[2, 2]), 10)
})
