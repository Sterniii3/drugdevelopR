#' @editor Lukas D Sauer
#' @editDate 2022-12-29
test_that("05.01", {
  # Multiple time-to-event endpoints
  res <- optimal_multiple_tte(alpha = 0.025,
                          beta = 0.1,
                          hr1 = 0.75, hr2 = 0.85,
                          id1 = 210, id2 = 420,
                          n2min = 20, n2max = 200, stepn2 = 4,
                          hrgomin = 0.70, hrgomax = 0.86, stephrgo = 0.02,
                          steps1 = 1, stepm1 = 0.95, stepl1 = 0.85,
                          b11 = 1000, b21 = 1500, b31 = 2000,
                          b12 = 1000, b22 = 2000, b32 = 3000,
                          num_cl = 12,
                          c02 = 100, c03 = 150,
                          c2 = 0.75, c3 = 1,
                          rho = 0.6,
                          fixed = FALSE,
               )
  expect_equal(res$u, 0, tolerance = 0.005)
  expect_equal(res$n2, 0)
  expect_equal(res$n3, 0)
  expect_equal(res$n, 0)
  expect_equal(res$HRgo, 0)
  expect_equal(res$pgo, 0)
})
#' @editor Lukas D Sauer
#' @editDate 2022-12-29
test_that("05.02", {
  # Multiple time-to-event endpoints with sample size constraint
  res <- optimal_multiple_tte(alpha = 0.025,
                              beta = 0.1,
                              hr1 = 0.75, hr2 = 0.85,
                              id1 = 210, id2 = 420,
                              n2min = 20, n2max = 200, stepn2 = 4,
                              hrgomin = 0.70, hrgomax = 0.86, stephrgo = 0.02,
                              steps1 = 1, stepm1 = 0.95, stepl1 = 0.85,
                              b11 = 1000, b21 = 1500, b31 = 2000,
                              b12 = 1000, b22 = 2000, b32 = 3000,
                              num_cl = 12,
                              c02 = 100, c03 = 150,
                              c2 = 0.75, c3 = 1,
                              rho = 0.6,
                              fixed = FALSE,
                              N = 300
                              
  )
  expect_equal(res$u, 0, tolerance = 0.005)
  expect_equal(res$n2, 0)
  expect_equal(res$n3, 0)
  expect_equal(res$n, 0)
  expect_equal(res$HRgo, 0)
  expect_equal(res$pgo, 0)
})
#' @editor Lukas D Sauer
#' @editDate 2022-12-29
test_that("05.03", {
  # Multiple time-to-event endpoints -- no cost limit
  res_nolim <- optimal_multiple_tte(alpha = 0.025,
                              beta = 0.1,
                              hr1 = 0.75, hr2 = 0.85,
                              id1 = 210, id2 = 420,
                              n2min = 20, n2max = 200, stepn2 = 4,
                              hrgomin = 0.70, hrgomax = 0.86, stephrgo = 0.02,
                              steps1 = 1, stepm1 = 0.95, stepl1 = 0.85,
                              b11 = 1000, b21 = 1500, b31 = 2000,
                              b12 = 1000, b22 = 2000, b32 = 3000,
                              num_cl = 12,
                              c02 = 100, c03 = 150,
                              c2 = 0.75, c3 = 1,
                              rho = 0.6,
                              fixed = TRUE,
  )
  # With cost limit
  res_lim <- optimal_multiple_tte(alpha = 0.025,
                              beta = 0.1,
                              hr1 = 0.75, hr2 = 0.85,
                              id1 = 210, id2 = 420,
                              n2min = 20, n2max = 200, stepn2 = 4,
                              hrgomin = 0.70, hrgomax = 0.86, stephrgo = 0.02,
                              steps1 = 1, stepm1 = 0.95, stepl1 = 0.85,
                              b11 = 1000, b21 = 1500, b31 = 2000,
                              b12 = 1000, b22 = 2000, b32 = 3000,
                              num_cl = 12,
                              c02 = 100, c03 = 150,
                              c2 = 0.75, c3 = 1,
                              rho = 0.6,
                              fixed = TRUE,
                              K = 400
  )
  expect_equal(res_nolim$u, 0, tolerance = 0.005)
  expect_equal(res_nolim$n2, 0)
  expect_equal(res_nolim$n3, 0)
  expect_equal(res_nolim$n, 0)
  expect_equal(res_nolim$HRgo, 0)
  expect_equal(res_nolim$pgo, 0)
  expect_equal(res_lim$u, 0, tolerance = 0.005)
  expect_equal(res_lim$n2, 0)
  expect_equal(res_lim$n3, 0)
  expect_equal(res_lim$n, 0)
  expect_equal(res_lim$HRgo, 0)
  expect_equal(res_lim$pgo, 0)
})
#' @editor Lukas D Sauer
#' @editDate 2022-12-29
test_that("05.04", {
  # Multiple time-to-event endpoints
  res <- optimal_multiple_tte(alpha = 0.025,
                              beta = 0.1,
                              hr1 = 0.75, hr2 = 0.85,
                              id1 = 210, id2 = 420,
                              n2min = 20, n2max = 200, stepn2 = 4,
                              hrgomin = 0.70, hrgomax = 0.86, stephrgo = 0.02,
                              steps1 = 1, stepm1 = 0.95, stepl1 = 0.85,
                              b11 = 1000, b21 = 1500, b31 = 2000,
                              b12 = 1000, b22 = 2000, b32 = 3000,
                              num_cl = 12,
                              c02 = 100, c03 = 150,
                              c2 = 0.75, c3 = 1,
                              rho = 0.6,
                              fixed = TRUE,
                              S = 0.7
  )
  expect_equal(res$u, 0, tolerance = 0.005)
  expect_equal(res$n2, 0)
  expect_equal(res$n3, 0)
  expect_equal(res$n, 0)
  expect_equal(res$HRgo, 0)
  expect_equal(res$pgo, 0)
})
#' @editor Lukas D Sauer
#' @editDate 2022-12-29
test_that("05.05", {
  # Multiple time-to-event endpoints
  start_time_3 = Sys.time()
  optimal_multiple_tte(alpha = 0.025,
                              beta = 0.1,
                              hr1 = 0.75, hr2 = 0.85,
                              id1 = 210, id2 = 420,
                              n2min = 20, n2max = 200, stepn2 = 4,
                              hrgomin = 0.70, hrgomax = 0.86, stephrgo = 0.02,
                              steps1 = 1, stepm1 = 0.95, stepl1 = 0.85,
                              b11 = 1000, b21 = 1500, b31 = 2000,
                              b12 = 1000, b22 = 2000, b32 = 3000,
                              num_cl = 12,
                              c02 = 100, c03 = 150,
                              c2 = 0.75, c3 = 1,
                              rho = 0.6,
                              fixed = FALSE,
  )
  end_time_3 = Sys.time()
  time_elapsed_num_cl_3 = end_time_3 - start_time_3
  start_time_1 = Sys.time()
  optimal_multiple_tte(alpha = 0.025,
                              beta = 0.1,
                              hr1 = 0.75, hr2 = 0.85,
                              id1 = 210, id2 = 420,
                              n2min = 20, n2max = 200, stepn2 = 4,
                              hrgomin = 0.70, hrgomax = 0.86, stephrgo = 0.02,
                              steps1 = 1, stepm1 = 0.95, stepl1 = 0.85,
                              b11 = 1000, b21 = 1500, b31 = 2000,
                              b12 = 1000, b22 = 2000, b32 = 3000,
                              num_cl = 6,
                              c02 = 100, c03 = 150,
                              c2 = 0.75, c3 = 1,
                              rho = 0.6,
                              fixed = FALSE,
  )
  end_time_1 = Sys.time()
  time_elapsed_num_cl_1 = end_time_1 - start_time_1
  expect_true(time_elapsed_num_cl_1 > time_elapsed_num_cl_3)
})
#' @editor Lukas D Sauer
#' @editDate 2022-12-29
test_that("05.06", {
  # Multiple normally distributed endpoints -- no sample size constraint
  res_nolim <- optimal_multiple_normal(alpha = 0.05,
                              beta = 0.1,
                              Delta1 = 0.75, Delta2 = 0.85,
                              n2min = 20, n2max = 200, stepn2 = 4,
                              kappamin = 0.02, kappamax = 0.2, kappago = 0.02,
                              b1 = 1000, b2 = 2000, b3 = 3000,
                              num_cl = 12,
                              c02 = 15, c03 = 20,
                              c2 = 0.675, c3 = 0.72,
                              fixed = FALSE,
                              rho = 0.6,
                              sigma1 = 8, sigma2 = 12,
                              in1 = 210, in2 = 420,
                              relaxed = TRUE
  )
  expect_equal(res_nolim$u, 0, tolerance = 0.005)
  expect_equal(res_nolim$n2, 0)
  expect_equal(res_nolim$n3, 0)
  expect_equal(res_nolim$n, 0)
  expect_equal(res_nolim$HRgo, 0)
  expect_equal(res_nolim$pgo, 0)
  # Multiple normally distributed endpoints -- with sample size constraint
  res_lim <- optimal_multiple_normal(alpha = 0.05,
                                     beta = 0.1,
                                     Delta1 = 0.75, Delta2 = 0.85,
                                     n2min = 20, n2max = 200, stepn2 = 4,
                                     kappamin = 0.02, kappamax = 0.2, kappago = 0.02,
                                     b1 = 1000, b2 = 2000, b3 = 3000,
                                     num_cl = 12,
                                     c02 = 15, c03 = 20,
                                     c2 = 0.675, c3 = 0.72,
                                     fixed = FALSE,
                                     rho = 0.6,
                                     sigma1 = 8, sigma2 = 12,
                                     in1 = 210, in2 = 420,
                                     relaxed = TRUE,
                                     N = 300,
  )
  expect_equal(res_lim$u, 0, tolerance = 0.005)
  expect_equal(res_lim$n2, 0)
  expect_equal(res_lim$n3, 0)
  expect_equal(res_lim$n, 0)
  expect_equal(res_lim$HRgo, 0)
  expect_equal(res_lim$pgo, 0)
})
#' @editor Lukas D Sauer
#' @editDate 2022-12-29
test_that("05.07", {
  # Multiple normally distributed endpoints -- parallel computing
  start_time_3 = Sys.time()
  optimal_multiple_normal(alpha = 0.05,
                                       beta = 0.1,
                                       Delta1 = 0.75, Delta2 = 0.85,
                                       n2min = 20, n2max = 200, stepn2 = 4,
                                       kappamin = 0.02, kappamax = 0.2, kappago = 0.02,
                                       b1 = 1000, b2 = 2000, b3 = 3000,
                                       num_cl = 12,
                                       c02 = 15, c03 = 20,
                                       c2 = 0.675, c3 = 0.72,
                                       fixed = FALSE,
                                       rho = 0.6,
                                       sigma1 = 8, sigma2 = 12,
                                       in1 = 210, in2 = 420,
                                       relaxed = TRUE
  )
  end_time_3 = Sys.time()
  time_elapsed_num_cl_3 = end_time_3 - start_time_3
  start_time_1 = Sys.time()
  optimal_multiple_normal(alpha = 0.05,
                          beta = 0.1,
                          Delta1 = 0.75, Delta2 = 0.85,
                          n2min = 20, n2max = 200, stepn2 = 4,
                          kappamin = 0.02, kappamax = 0.2, kappago = 0.02,
                          b1 = 1000, b2 = 2000, b3 = 3000,
                          num_cl = 6,
                          c02 = 15, c03 = 20,
                          c2 = 0.675, c3 = 0.72,
                          fixed = FALSE,
                          rho = 0.6,
                          sigma1 = 8, sigma2 = 12,
                          in1 = 210, in2 = 420,
                          relaxed = TRUE
  )
  end_time_1 = Sys.time()
  time_elapsed_num_cl_1 = end_time_1 - start_time_1
  expect_true(time_elapsed_num_cl_1 > time_elapsed_num_cl_3)
})
#' @editor Lukas D Sauer
#' @editDate 2022-12-29
test_that("05.08", {
  # Multiple normally distributed endpoints -- with and without cost limit
  res_nolim <- optimal_multiple_normal(alpha = 0.05,
                                       beta = 0.1,
                                       Delta1 = 0.75, Delta2 = 0.85,
                                       n2min = 20, n2max = 200, stepn2 = 4,
                                       kappamin = 0.02, kappamax = 0.2, kappago = 0.02,
                                       b1 = 1000, b2 = 2000, b3 = 3000,
                                       num_cl = 12,
                                       c02 = 15, c03 = 20,
                                       c2 = 0.675, c3 = 0.72,
                                       fixed = TRUE,
                                       rho = 0.6,
                                       sigma1 = 8, sigma2 = 12,
                                       in1 = 210, in2 = 420,
                                       relaxed = TRUE
  )
  res_lim <- optimal_multiple_normal(alpha = 0.05,
                                       beta = 0.1,
                                       Delta1 = 0.75, Delta2 = 0.85,
                                       n2min = 20, n2max = 200, stepn2 = 4,
                                       kappamin = 0.02, kappamax = 0.2, kappago = 0.02,
                                       b1 = 1000, b2 = 2000, b3 = 3000,
                                       num_cl = 12,
                                       c02 = 15, c03 = 20,
                                       c2 = 0.675, c3 = 0.72,
                                       fixed = TRUE,
                                       rho = 0.6,
                                       sigma1 = 8, sigma2 = 12,
                                       in1 = 210, in2 = 420,
                                       relaxed = TRUE,
                                       K = 400
  )
  expect_equal(res_nolim$u, 0, tolerance = 0.005)
  expect_equal(res_nolim$n2, 0)
  expect_equal(res_nolim$n3, 0)
  expect_equal(res_nolim$n, 0)
  expect_equal(res_nolim$HRgo, 0)
  expect_equal(res_nolim$pgo, 0)
  expect_equal(res_lim$u, 0, tolerance = 0.005)
  expect_equal(res_lim$n2, 0)
  expect_equal(res_lim$n3, 0)
  expect_equal(res_lim$n, 0)
  expect_equal(res_lim$HRgo, 0)
  expect_equal(res_lim$pgo, 0)
})
#' @editor Lukas D Sauer
#' @editDate 2022-12-29
test_that("05.09", {
  # Multiple normally distributed endpoints
  # -- strict effect size combination rule, constraint on success probability
  res_nolim <- optimal_multiple_normal(alpha = 0.05,
                                       beta = 0.1,
                                       Delta1 = 0.75, Delta2 = 0.85,
                                       n2min = 20, n2max = 200, stepn2 = 4,
                                       kappamin = 0.02, kappamax = 0.2, kappago = 0.02,
                                       b1 = 1000, b2 = 2000, b3 = 3000,
                                       num_cl = 12,
                                       c02 = 15, c03 = 20,
                                       c2 = 0.675, c3 = 0.72,
                                       fixed = TRUE,
                                       rho = 0.6,
                                       sigma1 = 8, sigma2 = 12,
                                       in1 = 210, in2 = 420,
                                       relaxed = FALSE
  )
  res_lim <- optimal_multiple_normal(alpha = 0.05,
                                       beta = 0.1,
                                       Delta1 = 0.75, Delta2 = 0.85,
                                       n2min = 20, n2max = 200, stepn2 = 4,
                                       kappamin = 0.02, kappamax = 0.2, kappago = 0.02,
                                       b1 = 1000, b2 = 2000, b3 = 3000,
                                       num_cl = 12,
                                       c02 = 15, c03 = 20,
                                       c2 = 0.675, c3 = 0.72,
                                       fixed = TRUE,
                                       rho = 0.6,
                                       sigma1 = 8, sigma2 = 12,
                                       in1 = 210, in2 = 420,
                                       relaxed = FALSE,
                                     S = 0.7
  )
  expect_equal(res_nolim$u, 0, tolerance = 0.005)
  expect_equal(res_nolim$n2, 0)
  expect_equal(res_nolim$n3, 0)
  expect_equal(res_nolim$n, 0)
  expect_equal(res_nolim$HRgo, 0)
  expect_equal(res_nolim$pgo, 0)
  expect_equal(res_lim$u, 0, tolerance = 0.005)
  expect_equal(res_lim$n2, 0)
  expect_equal(res_lim$n3, 0)
  expect_equal(res_lim$n, 0)
  expect_equal(res_lim$HRgo, 0)
  expect_equal(res_lim$pgo, 0)
})
#' @editor Lukas D Sauer
#' @editDate 2022-12-29
test_that("05.10", {
  # Multiple normally distributed endpoints
  # -- with and without
  res <- optimal_multiple_normal(alpha = 0.05,
                                       beta = 0.1,
                                       Delta1 = 0.75, Delta2 = 0.85,
                                       n2min = 20, n2max = 200, stepn2 = 4,
                                       kappamin = 0.02, kappamax = 0.2, kappago = 0.02,
                                       b1 = 1000, b2 = 2000, b3 = 3000,
                                       num_cl = 12,
                                       c02 = 15, c03 = 20,
                                       c2 = 0.675, c3 = 0.72,
                                       fixed = FALSE,
                                       rho = 0.6,
                                       sigma1 = 8, sigma2 = 12,
                                       in1 = 210, in2 = 420,
                                       relaxed = FALSE
  )
  expect_equal(res$u, 0, tolerance = 0.005)
  expect_equal(res$n2, 0)
  expect_equal(res$n3, 0)
  expect_equal(res$n, 0)
  expect_equal(res$HRgo, 0)
  expect_equal(res$pgo, 0)
})