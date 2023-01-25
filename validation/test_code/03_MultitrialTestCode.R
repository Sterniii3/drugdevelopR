#' @editor Lukas D Sauer
#' @editDate 2022-12-23
test_that("03.01", {
  # Multi-trial time-to-event endpoints
  res <- optimal_multitrial(alpha = 0.025,
                            beta = 0.1,
                            hr1 = 0.69, hr2 = 0.88,
                            xi2 = 0.7, xi3 = 0.7,
                            d2min = 10, d2max = 400, stepd2 = 2,
                            hrgomin = 0.65, hrgomax = 0.95, stephrgo = 0.01,
                            b1 = 1000, b2 = 3000, b3 = 5000,
                            num_cl = 3,
                            c02 = 100, c03 = 150,
                            c2 = 0.75, c3 = 1,
                            fixed = FALSE,
                            w = 0.3,
                            id1 = 210, id2 = 420,
                            case = 2, strategy = TRUE
               )
  # Strategies
  expect_equal(res$Strategy, c(1, 2, 3, 23))
  # Strategy 1
  expect_equal(res[1,]$u, -1.55, tolerance = 0.005)
  expect_equal(res[1,]$n2, 180)
  expect_equal(res[1,]$n3, 238)
  expect_equal(res[1,]$n, 418)
  expect_equal(res[1,]$d2, 126)
  expect_equal(res[1,]$d3, 167)
  expect_equal(res[1,]$d, 293)
  expect_equal(res[1,]$HRgo, 0.75)
  # Strategy 2
  expect_equal(res[2,]$u, -94.13, tolerance = 0.005)
  expect_equal(res[2,]$n2, 172)
  expect_equal(res[2,]$n3, 172)
  expect_equal(res[2,]$n, 344)
  expect_equal(res[2,]$d2, 120)
  expect_equal(res[2,]$d3, 122)
  expect_equal(res[2,]$d, 242)
  expect_equal(res[2,]$HRgo, 0.72)
  # Strategy 3
  expect_equal(res[3,]$u, -11.67, tolerance = 0.005)
  expect_equal(res[3,]$n2, 220)
  expect_equal(res[3,]$n3, 252)
  expect_equal(res[3,]$n, 472)
  expect_equal(res[3,]$d2, 154)
  expect_equal(res[3,]$d3, 177)
  expect_equal(res[3,]$d, 331)
  expect_equal(res[3,]$HRgo, 0.72)
  # Strategy 23
  expect_equal(res[4,]$u, 45.84, tolerance = 0.005)
  expect_equal(res[4,]$n2, 180)
  expect_equal(res[4,]$n3, 194)
  expect_equal(res[4,]$n, 374)
  expect_equal(res[4,]$d2, 126)
  expect_equal(res[4,]$d3, 136)
  expect_equal(res[4,]$d, 262)
  expect_equal(res[4,]$HRgo, 0.73)
})
#' @editor Lukas D Sauer
#' @editDate 2022-12-23
test_that("03.02", {
  # Multi-trial time-to-event endpoints
  res <- optimal_multitrial(alpha = 0.025,
                            beta = 0.1,
                            hr1 = 0.69, hr2 = 0.88,
                            xi2 = 0.7, xi3 = 0.7,
                            d2min = 10, d2max = 400, stepd2 = 2,
                            hrgomin = 0.65, hrgomax = 0.95, stephrgo = 0.01,
                            b1 = 1000, b2 = 3000, b3 = 5000,
                            num_cl = 3,
                            c02 = 100, c03 = 150,
                            c2 = 0.75, c3 = 1,
                            fixed = FALSE,
                            w = 0.3,
                            id1 = 210, id2 = 420,
                            case = 3, strategy = 1
  )
  # Strategies
  expect_equal(res$Strategy, 1)
  expect_equal(res$u, -150.95, tolerance = 0.005)
  expect_equal(res$n2, 176)
  expect_equal(res$n3, 146)
  expect_equal(res$n, 322)
  expect_equal(res$HRgo, 0.68)
  expect_equal(res$pgo, 0.22)
})
#' @editor Lukas D Sauer
#' @editDate 2022-12-23
test_that("03.03", {
  # Multi-trial time-to-event endpoints
  res <- optimal_multitrial(alpha = 0.025,
                            beta = 0.1,
                            hr1 = 0.69, hr2 = 0.88,
                            xi2 = 0.7, xi3 = 0.7,
                            d2min = 10, d2max = 400, stepd2 = 2,
                            hrgomin = 0.65, hrgomax = 0.95, stephrgo = 0.01,
                            b1 = 1000, b2 = 3000, b3 = 5000,
                            num_cl = 3,
                            c02 = 100, c03 = 150,
                            c2 = 0.75, c3 = 1,
                            fixed = FALSE,
                            w = 0.3,
                            id1 = 210, id2 = 420,
                            case = 2, strategy = TRUE,
                            K = 500
  )
  # Strategies
  expect_equal(res$Strategy, c(1, 2, 3, 23))
  # Strategy 1
  expect_equal(res[1,]$u, -3.75, tolerance = 0.005)
  expect_equal(res[1,]$n2, 180)
  expect_equal(res[1,]$n3, 212)
  expect_equal(res[1,]$n, 392)
  expect_equal(res[1,]$K2, 235)
  expect_equal(res[1,]$K3, 261)
  expect_equal(res[1,]$HRgo, 0.74)
  # Strategy 2
  expect_equal(res[2,]$u, -94.13, tolerance = 0.005)
  expect_equal(res[2,]$n2, 172)
  expect_equal(res[2,]$n3, 172)
  expect_equal(res[2,]$n, 344)
  expect_equal(res[2,]$d2, 120)
  expect_equal(res[2,]$d3, 122)
  expect_equal(res[2,]$d, 242)
  expect_equal(res[2,]$HRgo, 0.72)
  # Strategy 3
  expect_equal(res[3,]$u, -28.05, tolerance = 0.005)
  expect_equal(res[3,]$n2, 206)
  expect_equal(res[3,]$n3, 150)
  expect_equal(res[3,]$n, 356)
  expect_equal(res[3,]$K2, 254)
  expect_equal(res[3,]$K3, 243)
  expect_equal(res[3,]$HRgo, 0.68)
  # Strategy 23
  expect_equal(res[4,]$u, 45.84, tolerance = 0.005)
  expect_equal(res[4,]$n2, 180)
  expect_equal(res[4,]$n3, 194)
  expect_equal(res[4,]$n, 374)
  expect_equal(res[4,]$d2, 126)
  expect_equal(res[4,]$d3, 136)
  expect_equal(res[4,]$d, 262)
  expect_equal(res[4,]$HRgo, 0.73)
})
#' @editor Lukas D Sauer
#' @editDate 2022-12-23
test_that("03.04", {
  # Multi-trial time-to-event endpoints
  expect_error(optimal_multitrial(alpha = 0.025,
                            beta = 0.1,
                            hr1 = 0.69, hr2 = 0.88,
                            xi2 = 0.7, xi3 = 0.7,
                            d2min = 10, d2max = 400, stepd2 = 2,
                            hrgomin = 0.65, hrgomax = 0.95, stephrgo = 0.01,
                            b1 = 1000, b2 = 3000, b3 = 5000,
                            num_cl = 3,
                            c02 = 100, c03 = 150,
                            c2 = 0.75, c3 = 1,
                            fixed = FALSE,
                            w = 0.3,
                            id1 = 210, id2 = 420,
                            case = 3, strategy = 2
  ))
})
#' @editor Lukas D Sauer
#' @editDate 2022-12-23
test_that("03.05", {
  # Multi-trial time-to-event endpoints
  res <- optimal_multitrial(alpha = 0.025,
                            beta = 0.1,
                            hr1 = 0.69, hr2 = 0.88,
                            xi2 = 0.7, xi3 = 0.7,
                            d2min = 10, d2max = 400, stepd2 = 2,
                            hrgomin = 0.65, hrgomax = 0.95, stephrgo = 0.01,
                            b1 = 1000, b2 = 3000, b3 = 5000,
                            num_cl = 3,
                            c02 = 100, c03 = 150,
                            c2 = 0.75, c3 = 1,
                            fixed = TRUE,
                            w = 0.3,
                            id1 = 210, id2 = 420,
                            case = 2, strategy = TRUE
  )
  # Strategies
  expect_equal(res$Strategy, c(1, 2, 3, 23))
  # Strategy 1
  expect_equal(res[1,]$u, 1166.41, tolerance = 0.005)
  expect_equal(res[1,]$n2, 420)
  expect_equal(res[1,]$n3, 1044)
  expect_equal(res[1,]$n, 1464)
  expect_equal(res[1,]$d2, 294)
  expect_equal(res[1,]$d3, 731)
  expect_equal(res[1,]$d, 1025)
  expect_equal(res[1,]$HRgo, 0.86)
  # Strategy 2
  expect_equal(res[2,]$u, 811.67, tolerance = 0.005)
  expect_equal(res[2,]$n2, 386)
  expect_equal(res[2,]$n3, 1040)
  expect_equal(res[2,]$n, 1426)
  expect_equal(res[2,]$d2, 270)
  expect_equal(res[2,]$d3, 728)
  expect_equal(res[2,]$d, 998)
  expect_equal(res[2,]$HRgo, 0.85)
  # Strategy 3
  expect_equal(res[3,]$u, 1045.41, tolerance = 0.005)
  expect_equal(res[3,]$n2, 420)
  expect_equal(res[3,]$n3, 1386)
  expect_equal(res[3,]$n, 1806)
  expect_equal(res[3,]$d2, 294)
  expect_equal(res[3,]$d3, 972)
  expect_equal(res[3,]$d, 1266)
  expect_equal(res[3,]$HRgo, 0.82)
  # Strategy 23
  expect_equal(res[4,]$u, 45.84, tolerance = 0.005)
  expect_equal(res[4,]$n2, 180)
  expect_equal(res[4,]$n3, 194)
  expect_equal(res[4,]$n, 374)
  expect_equal(res[4,]$d2, 126)
  expect_equal(res[4,]$d3, 136)
  expect_equal(res[4,]$d, 262)
  expect_equal(res[4,]$HRgo, 0.73)
  expect_equal(res[4,]$pgo, 0.09)
})
#' @editor Lukas D Sauer
#' @editDate 2022-12-23
test_that("03.06", {
  # Multi-trial time-to-event endpoints
  res <- optimal_multitrial(alpha = 0.025,
                            beta = 0.1,
                            hr1 = 0.69, hr2 = 0.88,
                            xi2 = 0.7, xi3 = 0.7,
                            d2min = 10, d2max = 400, stepd2 = 2,
                            hrgomin = 0.65, hrgomax = 0.95, stephrgo = 0.01,
                            b1 = 1000, b2 = 3000, b3 = 5000,
                            num_cl = 3,
                            c02 = 100, c03 = 150,
                            c2 = 0.75, c3 = 1,
                            fixed = FALSE,
                            w = 0.3,
                            id1 = 210, id2 = 420,
                            case = 1, strategy = 1
  )
  # Strategies
  expect_equal(res$Strategy, 1)
  # Strategy 1
  expect_equal(res$n2, 206) # optimal sample size in phase II
  expect_equal(res$n3, 354) # resulting sample size in phase III
  expect_equal(res$n, 560) # resulting total sample size
  expect_equal(res$u, 432, tolerance = 0.05) # expected utility
  expect_equal(res$HRgo, 0.84) # threshold for proceeding to phase III
  expect_equal(res$d2, 144) # expected number of events in phase II
  expect_equal(res$d3, 248) # expected number of events in phase III
  expect_equal(res$d, 392) # total expected number of events
})
#' @editor Lukas D Sauer
#' @editDate 2022-12-23
test_that("03.07", {
  # Multi-trial binary endpoints
  res <- optimal_multitrial_binary(alpha = 0.025,
                            beta = 0.1,
                            p0 = 0.6, p11 = 0.3, p12 = 0.5,
                            n2min = 10, n2max = 400, stepn2 = 2,
                            rrgomin = 0.70, rrgomax = 0.95, steprrgo = 0.01,
                            b1 = 1000, b2 = 3000, b3 = 5000,
                            num_cl = 3,
                            c02 = 100, c03 = 150,
                            c2 = 0.75, c3 = 1,
                            fixed = FALSE,
                            w = 0.3,
                            in1 = 30, in2 = 60,
                            case = 3, strategy = TRUE,
  )
  # Strategies
  expect_equal(res$Strategy, c(1, 3, 4))
  # Strategy 1
  expect_equal(res[1,]$u, 71.66, tolerance = 0.005)
  expect_equal(res[1,]$RRgo, 0.71)
  expect_equal(res[1,]$n2, 130)
  expect_equal(res[1,]$n3, 178)
  expect_equal(res[1,]$n, 308)
  # Strategy 3
  expect_equal(res[2,]$u, 585.59, tolerance = 0.005)
  expect_equal(res[2,]$RRgo, 0.75)
  expect_equal(res[2,]$n2, 276)
  expect_equal(res[2,]$n3, 276)
  expect_equal(res[2,]$n, 552)
  # Strategy 4
  expect_equal(res[3,]$u, 718.18, tolerance = 0.005)
  expect_equal(res[3,]$RRgo, 0.73)
  expect_equal(res[3,]$n2, 316)
  expect_equal(res[3,]$n3, 304)
  expect_equal(res[3,]$n, 552)
})
#' @editor Lukas D Sauer
#' @editDate 2022-12-23
test_that("03.08", {
  # Multi-trial binary endpoints
  res <- optimal_multitrial_binary(alpha = 0.025,
                                   beta = 0.1,
                                   p0 = 0.6, p11 = 0.3, p12 = 0.5,
                                   n2min = 10, n2max = 400, stepn2 = 2,
                                   rrgomin = 0.70, rrgomax = 0.95, steprrgo = 0.01,
                                   b1 = 1000, b2 = 3000, b3 = 5000,
                                   num_cl = 3,
                                   c02 = 100, c03 = 150,
                                   c2 = 0.75, c3 = 1,
                                   fixed = FALSE,
                                   w = 0.3,
                                   in1 = 30, in2 = 60,
                                   case = 2, strategy = 23,
  )
  # Strategies
  expect_equal(res$Strategy, c(23))
  expect_equal(res$u, 810.94, tolerance = 0.005)
  expect_equal(res$n2, 220)
  expect_equal(res$n3, 184)
})
#' @editor Lukas D Sauer
#' @editDate 2022-12-23
test_that("03.09", {
  # Multi-trial binary endpoints
  res <- optimal_multitrial_binary(alpha = 0.025,
                                   beta = 0.1,
                                   p0 = 0.6, p11 = 0.3, p12 = 0.5,
                                   n2min = 10, n2max = 400, stepn2 = 2,
                                   rrgomin = 0.70, rrgomax = 0.95, steprrgo = 0.01,
                                   b1 = 1000, b2 = 3000, b3 = 5000,
                                   num_cl = 3,
                                   c02 = 100, c03 = 150,
                                   c2 = 0.75, c3 = 1,
                                   fixed = TRUE,
                                   w = 0.3,
                                   in1 = 30, in2 = 60,
                                   case = 1, strategy = TRUE,
  )
  # Strategies
  expect_equal(res$Strategy, c(1, 2))
  # Strategy 1
  expect_equal(res[1,]$u, 1742.40, tolerance = 0.005)
  expect_equal(res[1,]$RRgo, 0.88)
  expect_equal(res[1,]$n2, 162)
  expect_equal(res[1,]$n3, 160)
  expect_equal(res[1,]$n, 322)
  # Strategy 2
  expect_equal(res[2,]$u, 1878.56, tolerance = 0.005)
  expect_equal(res[2,]$RRgo, 0.82)
  expect_equal(res[2,]$n2, 200)
  expect_equal(res[2,]$n3, 292)
  expect_equal(res[2,]$n, 492)
})
#' @editor Lukas D Sauer
#' @editDate 2022-12-23
test_that("03.10", {
  # Multi-trial binary endpoints
  res <- optimal_multitrial_binary(alpha = 0.025,
                                   beta = 0.1,
                                   p0 = 0.6, p11 = 0.3, p12 = 0.5,
                                   n2min = 10, n2max = 400, stepn2 = 2,
                                   rrgomin = 0.70, rrgomax = 0.95, steprrgo = 0.01,
                                   b1 = 1000, b2 = 3000, b3 = 5000,
                                   num_cl = 3,
                                   c02 = 100, c03 = 150,
                                   c2 = 0.75, c3 = 1,
                                   fixed = TRUE,
                                   w = 0.3,
                                   in1 = 30, in2 = 60,
                                   case = 2, strategy = 3,
  )
  res_constrained <- optimal_multitrial_binary(alpha = 0.025,
                                   beta = 0.1,
                                   p0 = 0.6, p11 = 0.3, p12 = 0.5,
                                   n2min = 10, n2max = 400, stepn2 = 2,
                                   rrgomin = 0.70, rrgomax = 0.95, steprrgo = 0.01,
                                   b1 = 1000, b2 = 3000, b3 = 5000,
                                   num_cl = 3,
                                   c02 = 100, c03 = 150,
                                   c2 = 0.75, c3 = 1,
                                   fixed = TRUE,
                                   w = 0.3,
                                   in1 = 30, in2 = 60,
                                   case = 2, strategy = 3,
                                   N = 600
  )
  # Unconstrained
  expect_equal(res$Strategy, 3)
  expect_equal(res$u, 1332.94, tolerance = 0.005)
  expect_equal(res$n2, 242)
  expect_equal(res$n3, 414)
  expect_equal(res$n, 656)
  # Constrained
  expect_equal(res_constrained$u, 1313.11, tolerance = 0.005)
  expect_equal(res_constrained$n2, 186)
  expect_equal(res_constrained$n3, 414)
  expect_equal(res_constrained$n, 600)
})
#' @editor Lukas D Sauer
#' @editDate 2022-12-23
test_that("03.11", {
  # Multi-trial normal endpoints
  res <- optimal_multitrial_normal(alpha = 0.05,
                                   beta = 0.1,
                                   Delta1 = 0.375, Delta = 0.5,
                                   n2min = 20, n2max = 500, stepn2 = 4,
                                   kappamin = 0.02, kappamax = 0.4, stepkappa = 0.01,
                                   b1 = 3000, b2 = 8000, b3 = 10000,
                                   num_cl = 3,
                                   c02 = 15, c03 = 20,
                                   c2 = 0.675, c3 = 0.72,
                                   fixed = FALSE,
                                   w = 0.5,
                                   in1 = 300, in2 = 600,
                                   a = 0.25, b = 0.75,
                                   case = 3, strategy = TRUE,
  )
  # Strategies
  expect_equal(res$Strategy, c(1, 3, 4))
  # Strategy 1
  expect_equal(res[1,]$u, 1654.08, tolerance = 0.005)
  expect_equal(res[1,]$Kappa, 0.16)
  expect_equal(res[1,]$n2, 364)
  expect_equal(res[1,]$n3, 708)
  expect_equal(res[1,]$n, 1004)
  # Strategy 3
  expect_equal(res[2,]$u, 1308.56, tolerance = 0.005)
  expect_equal(res[2,]$Kappa, 0.16)
  expect_equal(res[2,]$n2, 296)
  expect_equal(res[2,]$n3, 708)
  expect_equal(res[2,]$n, 1004)
  expect_equal(res[2,]$sProg, 0.68)
  # Strategy 4
  expect_equal(res[3,]$u, 1843.1, tolerance = 0.005)
  expect_equal(res[3,]$Kappa, 0.18)
  expect_equal(res[3,]$n2, 342)
  expect_equal(res[3,]$n3, 888)
  expect_equal(res[3,]$n, 1230)
  expect_equal(res[2,]$sProg, 0.86)
})
#' @editor Lukas D Sauer
#' @editDate 2022-12-23
test_that("03.12", {
  # Testing that parallel computing has an effect
  start_time_3 = Sys.time()
  res <- optimal_multitrial_normal(alpha = 0.05,
                                   beta = 0.1,
                                   Delta1 = 0.375, Delta = 0.5,
                                   n2min = 10, n2max = 500, stepn2 = 2,
                                   kappamin = 0.01, kappamax = 0.5, stepkappa = 0.01,
                                   b1 = 3000, b2 = 8000, b3 = 10000,
                                   num_cl = 3,
                                   c02 = 15, c03 = 20,
                                   c2 = 0.675, c3 = 0.72,
                                   fixed = FALSE,
                                   w = 0.5,
                                   in1 = 300, in2 = 600,
                                   a = 0, b = 0.75,
                                   case = 3, strategy = TRUE,
  )
  end_time_3 = Sys.time()
  time_elapsed_num_cl_3 = end_time_3 - start_time_3
  start_time_1 = Sys.time()
  res <- optimal_multitrial_normal(alpha = 0.05,
                                   beta = 0.1,
                                   Delta1 = 0.375, Delta = 0.5,
                                   n2min = 10, n2max = 500, stepn2 = 2,
                                   kappamin = 0.01, kappamax = 0.5, stepkappa = 0.01,
                                   b1 = 3000, b2 = 8000, b3 = 10000,
                                   num_cl = 1,
                                   c02 = 15, c03 = 20,
                                   c2 = 0.675, c3 = 0.72,
                                   fixed = FALSE,
                                   w = 0.5,
                                   in1 = 300, in2 = 600,
                                   a = 0, b = 0.75,
                                   case = 3, strategy = TRUE,
  )
  end_time_1 = Sys.time()
  time_elapsed_num_cl_1 = end_time_1 - start_time_1
  expect_true(time_elapsed_num_cl_1 > time_elapsed_num_cl_3)
})
#' @editor Lukas D Sauer
#' @editDate 2022-12-23
test_that("03.13", {
  # Multi-trial normal endpoints
  res <- optimal_multitrial_normal(alpha = 0.05,
                                   beta = 0.1,
                                   Delta1 = 0.375, Delta = 0.5,
                                   n2min = 10, n2max = 500, stepn2 = 2,
                                   kappamin = 0.01, kappamax = 0.5, stepkappa = 0.01,
                                   b1 = 3000, b2 = 8000, b3 = 10000,
                                   num_cl = 3,
                                   c02 = 15, c03 = 20,
                                   c2 = 0.675, c3 = 0.72,
                                   fixed = TRUE,
                                   w = 0.5,
                                   in1 = 300, in2 = 600,
                                   a = 0, b = 0.75,
                                   case = 3, strategy = TRUE,
  )
  res_prob_constraint <- optimal_multitrial_normal(alpha = 0.05,
                                   beta = 0.1,
                                   Delta1 = 0.375, Delta = 0.5,
                                   n2min = 10, n2max = 500, stepn2 = 2,
                                   kappamin = 0.01, kappamax = 0.5, stepkappa = 0.01,
                                   b1 = 3000, b2 = 8000, b3 = 10000,
                                   num_cl = 3,
                                   c02 = 15, c03 = 20,
                                   c2 = 0.675, c3 = 0.72,
                                   fixed = TRUE,
                                   w = 0.5,
                                   in1 = 300, in2 = 600,
                                   a = 0, b = 0.75,
                                   case = 3, strategy = TRUE,
                                   S = 0.82
  )
  # Unconstrained
  expect_equal(res$Strategy, c(1, 3, 4))
  # Strategy 1
  expect_equal(res[1, ]$u, 1514.35, tolerance = 0.005)
  expect_equal(res[1, ]$Kappa, 0.16)
  expect_equal(res[1, ]$n2, 440)
  expect_equal(res[1, ]$n3, 830)
  expect_equal(res[1, ]$n, 1270)
  expect_equal(res[1, ]$sProg, 0.82)
  # Strategy 3
  expect_equal(res[2, ]$u, 1116.60, tolerance = 0.005)
  expect_equal(res[2, ]$Kappa, 0.16)
  expect_equal(res[2, ]$n2, 364)
  expect_equal(res[2, ]$n3, 882)
  expect_equal(res[2, ]$n, 1246)
  expect_equal(res[2, ]$sProg, 0.68)
  # Strategy 4
  expect_equal(res[3, ]$u, 1395.35, tolerance = 0.005)
  expect_equal(res[3, ]$Kappa, 0.18)
  expect_equal(res[3, ]$n2, 424)
  expect_equal(res[3, ]$n3, 1128)
  expect_equal(res[3, ]$n, 1552)
  expect_equal(res[3, ]$sProg, 0.86)
  # Constrained
  expect_equal(res_constrained$Strategy, c(1, 3, 4))
  # Strategy 1
  expect_equal(res_constrained[1, ]$u, 1514.35, tolerance = 0.005)
  expect_equal(res_constrained[1, ]$Kappa, 0.16)
  expect_equal(res_constrained[1, ]$n2, 440)
  expect_equal(res_constrained[1, ]$n3, 830)
  expect_equal(res_constrained[1, ]$n, 1270)
  expect_equal(res_constrained[1, ]$sProg, 0.82)
  # Strategy 3
  expect_equal(res_constrained[2, ]$u, -9999)
  # Strategy 4
  expect_equal(res_constrained[3, ]$u, 1395.35, tolerance = 0.005)
  expect_equal(res_constrained[3, ]$Kappa, 0.18)
  expect_equal(res_constrained[3, ]$n2, 424)
  expect_equal(res_constrained[3, ]$n3, 1128)
  expect_equal(res_constrained[3, ]$n, 1552)
  expect_equal(res_constrained[3, ]$sProg, 0.86)

})
#' @editor Lukas D Sauer
#' @editDate 2022-12-23
test_that("03.14", {
  # Multi-trial normal endpoints
  res <- optimal_multitrial_normal(alpha = 0.05,
                                   beta = 0.1,
                                   Delta1 = 0.375, Delta = 0.5,
                                   n2min = 10, n2max = 500, stepn2 = 2,
                                   kappamin = 0.01, kappamax = 0.5, stepkappa = 0.01,
                                   b1 = 3000, b2 = 8000, b3 = 10000,
                                   num_cl = 3,
                                   c02 = 15, c03 = 20,
                                   c2 = 0.675, c3 = 0.72,
                                   fixed = TRUE,
                                   w = 0.5,
                                   in1 = 300, in2 = 600,
                                   a = 0, b = 0.75,
                                   case = 2, strategy = 3,
  )
  # Strategy 1
  expect_equal(res$u, 1749.97, tolerance = 0.005)
  expect_equal(res$Kappa, 0.16)
  expect_equal(res$n2, 416)
  expect_equal(res$n3, 876)
  expect_equal(res$n, 1292)
})