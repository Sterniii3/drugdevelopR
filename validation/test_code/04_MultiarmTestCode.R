#' @editor Lukas D Sauer
#' @editDate 2022-12-23
test_that("04.01", {
  # Multiarm time-to-event endpoints
  res <- optimal_multiarm(alpha = 0.025,
                          beta = 0.1,
                          hr1 = 0.75, hr2 = 0.85,
                          ec = 0.6,
                          n2min = 10, n2max = 200, stepn2 = 1,
                          hrgomin = 0.71, hrgomax = 0.9, stephrgo = 0.01,
                          steps1 = 1, stepm1 = 0.95, stepl1 = 0.85,
                          b1 = 1000, b2 = 2000, b3 = 3000,
                          num_cl = 3,
                          c02 = 100, c03 = 150,
                          c2 = 0.75, c3 = 1,
                          strategy = 1
               )
  expect_equal(res$u, 8.56, tolerance = 0.005)
  expect_equal(res$n2, 141)
  expect_equal(res$n3, 391)
  expect_equal(res$n, 532)
  expect_equal(res$HRgo, 0.82)
  expect_equal(res$pgo, 0.72)
})
#' @editor Lukas D Sauer
#' @editDate 2022-12-23
test_that("04.02", {
  # Multiarm time-to-event endpoints
  res <- optimal_multiarm(alpha = 0.025,
                          beta = 0.1,
                          hr1 = 0.75, hr2 = 0.85,
                          ec = 0.6,
                          n2min = 10, n2max = 200, stepn2 = 1,
                          hrgomin = 0.71, hrgomax = 0.9, stephrgo = 0.01,
                          steps1 = 1, stepm1 = 0.95, stepl1 = 0.85,
                          b1 = 1000, b2 = 2000, b3 = 3000,
                          num_cl = 3,
                          c02 = 100, c03 = 150,
                          c2 = 0.75, c3 = 1,
                          strategy = 2
  )
  expect_equal(res$u, 4.36, tolerance = 0.005)
  expect_equal(res$n2, 79)
  expect_equal(res$n3, 426)
  expect_equal(res$n, 505)
  expect_equal(res$HRgo, 0.78)
  expect_equal(res$pgo, 0.65)
})
#' @editor Lukas D Sauer
#' @editDate 2022-12-23
test_that("04.03", {
  # Multiarm time-to-event endpoints
  res <- optimal_multiarm(alpha = 0.025,
                          beta = 0.1,
                          hr1 = 0.75, hr2 = 0.85,
                          ec = 0.6,
                          n2min = 10, n2max = 200, stepn2 = 1,
                          hrgomin = 0.71, hrgomax = 0.9, stephrgo = 0.01,
                          steps1 = 1, stepm1 = 0.95, stepl1 = 0.85,
                          b1 = 1000, b2 = 2000, b3 = 3000,
                          num_cl = 3,
                          c02 = 100, c03 = 150,
                          c2 = 0.75, c3 = 1,
                          strategy = 3,
                          N = 520
  )
  # Strategy 1 = best promising
  expect_equal(res[1,]$u, 8.34, tolerance = 0.005)
  expect_equal(res[1,]$n2, 133)
  expect_equal(res[1,]$n3, 383)
  expect_equal(res[1,]$n, 516)
  expect_equal(res[1,]$HRgo, 0.82)
  # Strategy 2 = all promising -- shouldn't change compared to 04.02
  expect_equal(res[2,]$u, 4.36, tolerance = 0.005)
  expect_equal(res[2,]$n2, 79)
  expect_equal(res[2,]$n3, 426)
  expect_equal(res[2,]$n, 505)
  expect_equal(res[2,]$HRgo, 0.78)
  expect_equal(res[2,]$pgo, 0.65)
})
#' @editor Lukas D Sauer
#' @editDate 2022-12-23
test_that("04.04", {
  # Multiarm time-to-event endpoints
  start_time_3 = Sys.time()
  res <- optimal_multiarm(alpha = 0.025,
                          beta = 0.1,
                          hr1 = 0.75, hr2 = 0.85,
                          ec = 0.6,
                          n2min = 10, n2max = 200, stepn2 = 1,
                          hrgomin = 0.71, hrgomax = 0.9, stephrgo = 0.01,
                          steps1 = 1, stepm1 = 0.95, stepl1 = 0.85,
                          b1 = 1000, b2 = 2000, b3 = 3000,
                          num_cl = 3,
                          c02 = 100, c03 = 150,
                          c2 = 0.75, c3 = 1,
                          strategy = 2
  )
  end_time_3 = Sys.time()
  time_elapsed_num_cl_3 = end_time_3 - start_time_3
  start_time_1 = Sys.time()
  res <- optimal_multiarm(alpha = 0.025,
                          beta = 0.1,
                          hr1 = 0.75, hr2 = 0.85,
                          ec = 0.6,
                          n2min = 10, n2max = 200, stepn2 = 1,
                          hrgomin = 0.71, hrgomax = 0.9, stephrgo = 0.01,
                          steps1 = 1, stepm1 = 0.95, stepl1 = 0.85,
                          b1 = 1000, b2 = 2000, b3 = 3000,
                          num_cl = 1,
                          c02 = 100, c03 = 150,
                          c2 = 0.75, c3 = 1,
                          strategy = 2
  )
  end_time_1 = Sys.time()
  time_elapsed_num_cl_1 = end_time_1 - start_time_1
  expect_true(time_elapsed_num_cl_1 > time_elapsed_num_cl_3)
})
#' @editor Lukas D Sauer
#' @editDate 2022-12-23
test_that("04.05", {
  # Multiarm binary endpoints
  res <- optimal_multiarm_binary(alpha = 0.025,
                          beta = 0.1,
                          p0 = 0.5, p11 = 0.3, p12 = 0.4,
                          n2min = 10, n2max = 400, stepn2 = 2,
                          rrgomin = 0.71, rrgomax = 0.9, steprrgo = 0.01,
                          steps1 = 1, stepm1 = 0.95, stepl1 = 0.85,
                          b1 = 1000, b2 = 2000, b3 = 3000,
                          num_cl = 3,
                          c02 = 100, c03 = 150,
                          c2 = 0.75, c3 = 1,
                          strategy = 1
  )
  expect_equal(res$u, 1398.06, tolerance = 0.005)
  expect_equal(res$n2, 370)
  expect_equal(res$n3, 289)
  expect_equal(res$n, 659)
  expect_equal(res$HRgo, 0.87)
})