#' @editor Lukas D Sauer
#' @editDate 2022-12-29
test_that("05.01", {
  # Multiple time-to-event endpoints (run time approx. 1 day)
  set.seed(123)
  res <- optimal_multiple_tte(alpha = 0.025,
                          beta = 0.1,
                          hr1 = 0.75, hr2 = 0.85,
                          id1 = 210, id2 = 420,
                          n2min = 100, n2max = 300, stepn2 = 4,
                          hrgomin = 0.80, hrgomax = 0.90, stephrgo = 0.02,
                          steps1 = 1, stepm1 = 0.95, stepl1 = 0.85,
                          b11 = 1000, b21 = 1500, b31 = 2000,
                          b12 = 1000, b22 = 2000, b32 = 3000,
                          num_cl = 12,
                          c02 = 100, c03 = 150,
                          c2 = 0.75, c3 = 1,
                          rho = 0.6,
                          fixed = FALSE,
               )
  
  expect_equal(res$u, -212.24, tolerance = 0.005)
  expect_equal(res$n2, 100)
  expect_equal(res$n3, 262)
  expect_equal(res$n, 362)
  expect_equal(res$HRgo, 0.88)
  expect_equal(res$pgo, 0.74)
})
#' @editor Lukas D Sauer
#' @editDate 2022-12-29
test_that("05.02", {
  # Multiple time-to-event endpoints with sample size constraint
  # (run time approx. 1 day)
  set.seed(123)
  res <- optimal_multiple_tte(alpha = 0.025,
                              beta = 0.1,
                              hr1 = 0.75, hr2 = 0.85,
                              id1 = 210, id2 = 420,
                              n2min = 100, n2max = 300, stepn2 = 4,
                              hrgomin = 0.80, hrgomax = 0.9, stephrgo = 0.02,
                              steps1 = 1, stepm1 = 0.95, stepl1 = 0.85,
                              b11 = 1000, b21 = 1500, b31 = 2000,
                              b12 = 1000, b22 = 2000, b32 = 3000,
                              num_cl = 12,
                              c02 = 100, c03 = 150,
                              c2 = 0.75, c3 = 1,
                              rho = 0.6,
                              fixed = FALSE,
                              N = 330
                              
  )
  
  expect_equal(res$u, -218.92, tolerance = 0.005)
  expect_equal(res$n2, 100)
  expect_equal(res$n3, 219)
  expect_equal(res$n, 319)
  expect_equal(res$HRgo, 0.86)
})
#' @editor Lukas D Sauer
#' @editDate 2022-12-29
test_that("05.03", {
  # Multiple time-to-event endpoints -- no cost limit (duration: several seconds)
  set.seed(123)
  res_nolim <- optimal_multiple_tte(alpha = 0.025,
                              beta = 0.1,
                              hr1 = 0.75, hr2 = 0.85,
                              id1 = 210, id2 = 420,
                              n2min = 100, n2max = 300, stepn2 = 4,
                              hrgomin = 0.80, hrgomax = 0.90, stephrgo = 0.02,
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
                              n2min = 100, n2max = 300, stepn2 = 4,
                              hrgomin = 0.80, hrgomax = 0.90, stephrgo = 0.02,
                              steps1 = 1, stepm1 = 0.95, stepl1 = 0.85,
                              b11 = 1000, b21 = 1500, b31 = 2000,
                              b12 = 1000, b22 = 2000, b32 = 3000,
                              num_cl = 12,
                              c02 = 100, c03 = 150,
                              c2 = 0.75, c3 = 1,
                              rho = 0.6,
                              fixed = TRUE,
                              K = 600
  )
  
  expect_equal(res_nolim$u, 144.96, tolerance = 0.005)
  expect_equal(res_nolim$n2, 172)
  expect_equal(res_nolim$n3, 408)
  expect_equal(res_nolim$n, 580)
  expect_equal(res_nolim$HRgo, 0.86)
  expect_equal(res_nolim$K2, 229)
  expect_equal(res_nolim$K3, 532)

  expect_equal(res_lim$u, 130.67, tolerance = 0.005)
  expect_equal(res_lim$n2, 112)
  expect_equal(res_lim$n3, 301)
  expect_equal(res_lim$n, 413)
  expect_equal(res_lim$HRgo, 0.84)
  expect_equal(res_lim$K2, 184)
  expect_equal(res_lim$K3, 414)
})
#' @editor Lukas D Sauer
#' @editDate 2022-12-29
test_that("05.04", {
  # Multiple time-to-event endpoints (duration: several seconds)
  set.seed(123)
  res <- optimal_multiple_tte(alpha = 0.025,
                              beta = 0.1,
                              hr1 = 0.75, hr2 = 0.85,
                              id1 = 210, id2 = 420,
                              n2min = 100, n2max = 300, stepn2 = 4,
                              hrgomin = 0.80, hrgomax = 0.90, stephrgo = 0.02,
                              steps1 = 1, stepm1 = 0.95, stepl1 = 0.85,
                              b11 = 1000, b21 = 1500, b31 = 2000,
                              b12 = 1000, b22 = 2000, b32 = 3000,
                              num_cl = 12,
                              c02 = 100, c03 = 150,
                              c2 = 0.75, c3 = 1,
                              rho = 0.6,
                              fixed = TRUE,
                              S = 0.6
  )
  
  expect_equal(res$u, 132.6, tolerance = 0.005)
  expect_equal(res$n2, 280)
  expect_equal(res$n3, 467)
  expect_equal(res$n, 747)
  expect_equal(res$HRgo, 0.86)
  expect_equal(res$sProg, 0.6)
  expect_equal(res$OS, 0.54)
})
#' @editor Lukas D Sauer
#' @editDate 2022-12-29
test_that("05.05", {
  # Multiple time-to-event endpoints (duration: approx. one day)
  start_time_3 = Sys.time()
  set.seed(123)
  optimal_multiple_tte(alpha = 0.025,
                              beta = 0.1,
                              hr1 = 0.75, hr2 = 0.85,
                              id1 = 210, id2 = 420,
                              n2min = 200, n2max = 300, stepn2 = 4,
                              hrgomin = 0.86, hrgomax = 0.90, stephrgo = 0.02,
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
                              n2min = 200, n2max = 300, stepn2 = 4,
                              hrgomin = 0.86, hrgomax = 0.90, stephrgo = 0.02,
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
  # (Duration: approx. one day)
  set.seed(123)
  res_nolim <- optimal_multiple_normal(alpha = 0.05,
                              beta = 0.1,
                              Delta1 = 0.75, Delta2 = 0.8,
                              n2min = 80, n2max = 160, stepn2 = 4,
                              kappamin = 0.02, kappamax = 0.1, stepkappa = 0.02,
                              b1 = 1000, b2 = 2000, b3 = 3000,
                              num_cl = 12,
                              c02 = 100, c03 = 150,
                              c2 = 0.75, c3 = 1,
                              fixed = FALSE,
                              rho = 0.5,
                              sigma1 = 2, sigma2 = 1,
                              in1 = 300, in2 = 600,
                              relaxed = TRUE
  )
  
  expect_equal(res_nolim$u, 1765, tolerance = 0.005)
  expect_equal(res_nolim$n2, 160)
  expect_equal(res_nolim$n3, 74)
  expect_equal(res_nolim$n, 234)
  expect_equal(res_nolim$Kappa, 0.06)
  # Multiple normally distributed endpoints -- with sample size constraint
  
  res_lim <- optimal_multiple_normal(alpha = 0.05,
                                     beta = 0.1,
                                     Delta1 = 0.75, Delta2 = 0.8,
                                     n2min = 80, n2max = 160, stepn2 = 4,
                                     kappamin = 0.02, kappamax = 0.1, stepkappa = 0.02,
                                     b1 = 1000, b2 = 2000, b3 = 3000,
                                     num_cl = 12,
                                     c02 = 100, c03 = 150,
                                     c2 = 0.75, c3 = 1,
                                     fixed = FALSE,
                                     rho = 0.5,
                                     sigma1 = 2, sigma2 = 1,
                                     in1 = 300, in2 = 600,
                                     relaxed = TRUE,
                                     N = 190,
  )
  
  expect_equal(res_lim$u, 1384, tolerance = 0.005)
  expect_equal(res_lim$n2, 80)
  expect_equal(res_lim$n3, 109)
  expect_equal(res_lim$n, 189)
  expect_equal(res_lim$Kappa, 0.08)
})
#' @editor Lukas D Sauer
#' @editDate 2022-12-29
test_that("05.07", {
  # Multiple normally distributed endpoints -- parallel computing
  # (Duration: approx. one day)
  start_time_3 = Sys.time()
  set.seed(123)
  optimal_multiple_normal(alpha = 0.05,
                                       beta = 0.1,
                                       Delta1 = 0.75, Delta2 = 0.8,
                                       n2min = 80, n2max = 160, stepn2 = 4,
                                       kappamin = 0.02, kappamax = 0.04, stepkappa = 0.02,
                                       b1 = 1000, b2 = 2000, b3 = 3000,
                                       num_cl = 12,
                                       c02 = 100, c03 = 150,
                                       c2 = 0.75, c3 = 1,
                                       fixed = FALSE,
                                       rho = 0.5,
                                       sigma1 = 2, sigma2 = 1,
                                       in1 = 300, in2 = 600,
                                       relaxed = TRUE
  )
  end_time_3 = Sys.time()
  time_elapsed_num_cl_3 = end_time_3 - start_time_3
  start_time_1 = Sys.time()
  optimal_multiple_normal(alpha = 0.05,
                          beta = 0.1,
                          Delta1 = 0.75, Delta2 = 0.8,
                          n2min = 80, n2max = 160, stepn2 = 4,
                          kappamin = 0.02, kappamax = 0.04, stepkappa = 0.02,
                          b1 = 1000, b2 = 2000, b3 = 3000,
                          num_cl = 6,
                          c02 = 100, c03 = 150,
                          c2 = 0.75, c3 = 1,
                          fixed = FALSE,
                          rho = 0.5,
                          sigma1 = 2, sigma2 = 1,
                          in1 = 300, in2 = 600,
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
  # Duration approx. 1 hour
  set.seed(123)
  res_nolim <- optimal_multiple_normal(alpha = 0.05,
                                       beta = 0.1,
                                       Delta1 = 0.75, Delta2 = 0.8,
                                       n2min = 80, n2max = 160, stepn2 = 4,
                                       kappamin = 0.02, kappamax = 0.1, stepkappa = 0.02,
                                       b1 = 1000, b2 = 2000, b3 = 3000,
                                       num_cl = 12,
                                       c02 = 100, c03 = 150,
                                       c2 = 0.75, c3 = 1,
                                       fixed = TRUE,
                                       rho = 0.5,
                                       sigma1 = 2, sigma2 = 1,
                                       in1 = 300, in2 = 600,
                                       relaxed = TRUE
  )
  res_lim <- optimal_multiple_normal(alpha = 0.05,
                                       beta = 0.1,
                                       Delta1 = 0.75, Delta2 = 0.8,
                                       n2min = 80, n2max = 160, stepn2 = 4,
                                       kappamin = 0.02, kappamax = 0.1, stepkappa = 0.02,
                                       b1 = 1000, b2 = 2000, b3 = 3000,
                                       num_cl = 12,
                                       c02 = 100, c03 = 150,
                                       c2 = 0.75, c3 = 1,
                                       fixed = TRUE,
                                       rho = 0.5,
                                       sigma1 = 2, sigma2 = 1,
                                       in1 = 300, in2 = 600,
                                       relaxed = TRUE,
                                       K = 400
  )
  
  expect_equal(res_nolim$u, 1133, tolerance = 0.005)
  expect_equal(res_nolim$n2, 160)
  expect_equal(res_nolim$K2, 220)
  expect_equal(res_nolim$K3, 215)
  expect_equal(res_nolim$Kappa, 0.02)
  expect_equal(res_lim$u, 1099, tolerance = 0.005)
  expect_equal(res_lim$n2, 104)
  expect_equal(res_lim$K2, 178)
  expect_equal(res_lim$K3, 220)
  expect_equal(res_lim$Kappa, 0.02)
})
#' @editor Lukas D Sauer
#' @editDate 2022-12-29
test_that("05.09", {
  # Multiple normally distributed endpoints
  # -- constraint on success probability
  # Duration several minutes
  set.seed(123)
  res_nolim <- optimal_multiple_normal(alpha = 0.05,
                                       beta = 0.1,
                                       Delta1 = 0.75, Delta2 = 0.8,
                                       n2min = 80, n2max = 160, stepn2 = 4,
                                       kappamin = 0.02, kappamax = 0.1, stepkappa = 0.02,
                                       b1 = 1000, b2 = 2000, b3 = 3000,
                                       num_cl = 12,
                                       c02 = 100, c03 = 150,
                                       c2 = 0.75, c3 = 1,
                                       fixed = TRUE,
                                       rho = 0.5,
                                       sigma1 = 2, sigma2 = 1,
                                       in1 = 300, in2 = 600,
                                       relaxed = TRUE
  )
  res_lim <- optimal_multiple_normal(alpha = 0.05,
                                       beta = 0.1,
                                       Delta1 = 0.75, Delta2 = 0.8,
                                       n2min = 80, n2max = 160, stepn2 = 4,
                                       kappamin = 0.02, kappamax = 0.1, stepkappa = 0.02,
                                       b1 = 1000, b2 = 2000, b3 = 3000,
                                       num_cl = 12,
                                       c02 = 100, c03 = 150,
                                       c2 = 0.75, c3 = 1,
                                       fixed = TRUE,
                                       rho = 0.5,
                                       sigma1 = 2, sigma2 = 1,
                                       in1 = 300, in2 = 600,
                                       relaxed = TRUE,
                                       S = 0.85
  )
  
  expect_equal(res_nolim$u, 1133, tolerance = 0.005)
  expect_equal(res_nolim$n2, 160)
  expect_equal(res_nolim$sProg, 0.83)
  expect_equal(res_nolim$sProg1, 0.15)
  expect_equal(res_nolim$sProg2, 0.64)
  expect_equal(res_nolim$sProg3, 0.048)
  expect_equal(res_nolim$Kappa, 0.02)
  
  expect_equal(res_lim$u, -9999, tolerance = 0.005)
})
#' @editor Lukas D Sauer
#' @editDate 2022-12-29
test_that("05.10", {
  # Multiple normally distributed endpoints
  # -- relaxed and strict rule
  # Duration several minutes
  set.seed(123)
  res_relax <- optimal_multiple_normal(alpha = 0.05,

                                       beta = 0.1,
                                       Delta1 = 0.75, Delta2 = 0.8,
                                       n2min = 80, n2max = 160, stepn2 = 4,
                                       kappamin = 0.02, kappamax = 0.1, stepkappa = 0.02,
                                       b1 = 1000, b2 = 2000, b3 = 3000,
                                       num_cl = 12,
                                       c02 = 100, c03 = 150,
                                       c2 = 0.75, c3 = 1,
                                       fixed = TRUE,
                                       rho = 0.5,
                                       sigma1 = 2, sigma2 = 1,
                                       in1 = 300, in2 = 600,
                                       relaxed = TRUE
  )
  res_strict <- optimal_multiple_normal(alpha = 0.05,
                                     beta = 0.1,
                                     Delta1 = 0.75, Delta2 = 0.8,
                                     n2min = 80, n2max = 160, stepn2 = 4,
                                     kappamin = 0.02, kappamax = 0.1, stepkappa = 0.02,
                                     b1 = 1000, b2 = 2000, b3 = 3000,
                                     num_cl = 12,
                                     c02 = 100, c03 = 150,
                                     c2 = 0.75, c3 = 1,
                                     fixed = TRUE,
                                     rho = 0.5,
                                     sigma1 = 2, sigma2 = 1,
                                     in1 = 300, in2 = 600,
                                     relaxed = FALSE,
  )
  
  expect_equal(res_relax$u, 1133, tolerance = 0.005)
  expect_equal(res_relax$n2, 160)
  expect_equal(res_relax$Kappa, 0.02)
  expect_equal(res_relax$pgo, 0.99)

  expect_equal(res_strict$u, 410, tolerance = 0.005)
  expect_equal(res_strict$n2, 132)
  expect_equal(res_strict$Kappa, 0.02)
  expect_equal(res_strict$pgo, 0.98)
})


