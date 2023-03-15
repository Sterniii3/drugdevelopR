#' @editor Lukas D Sauer
#' @editDate 2022-12-07
test_that("02.01", {
  # Adjustment method "additive" for time-to-event endpoints
  res <- optimal_bias(alpha = 0.025,
               beta = 0.1,
               hr1 = 0.69, hr2 = 0.88,
               xi2 = 0.7, xi3 = 0.7,
               d2min = 20, d2max = 100, stepd2 = 5,
               hrgomin = 0.7, hrgomax = 0.9, stephrgo = 0.02,
               steps1 = 1, stepm1 = 0.95, stepl1 = 0.85,
               b1 = 1000, b2 = 2000, b3 = 3000,
               num_cl = 12,
               c02 = 100, c03 = 150, # fixed cost
               c2 = 0.75, c3 = 1, # variable per patient cost
               fixed = FALSE,
               w = 0.3,
               id1 = 210, id2 = 420,
               adj = "additive",
               alphaCImin = 0.3, alphaCImax = 0.5, stepalphaCI = 0.025,
               lambdamin = NULL, lambdamax = NULL, steplambda = NULL
               )
  expect_equal(res$n2, 122)
  expect_equal(res$n3, 200)
  expect_equal(res$n, 322)
  expect_equal(res$u, 78, tolerance = 0.005)
  expect_equal(res$HRgo, 0.78)
  expect_equal(res$d2, 85)
  expect_equal(res$d3, 140)
  expect_equal(res$d, 225)
  
})
#' @editor Lukas D Sauer
#' @editDate 2022-12-07
test_that("02.02", {
  # Adjustment method "multiplicative" for time-to-event endpoints
  res <- optimal_bias(alpha = 0.025,
                      beta = 0.1,
                      hr1 = 0.69, hr2 = 0.88,
                      xi2 = 0.7, xi3 = 0.7,
                      d2min = 20, d2max = 100, stepd2 = 5,
                      hrgomin = 0.7, hrgomax = 0.9, stephrgo = 0.02,
                      steps1 = 1, stepm1 = 0.95, stepl1 = 0.85,
                      b1 = 1000, b2 = 2000, b3 = 3000,
                      num_cl = 12,
                      c02 = 100, c03 = 150, # fixed cost
                      c2 = 0.75, c3 = 1, # variable per patient cost
                      fixed = FALSE,
                      w = 0.3,
                      id1 = 210, id2 = 420,
                      adj = "multiplicative",
                      alphaCImin = NULL, alphaCImax = NULL, stepalphaCI = NULL,
                      lambdamin = 0.5, lambdamax = 1, steplambda = 0.05
  )
  expect_equal(res$n2, 136)
  expect_equal(res$n3, 244)
  expect_equal(res$n, 380)
  expect_equal(res$u, 99, tolerance = 0.005)
  expect_equal(res$HRgo, 0.76)
  expect_equal(res$pgo, 0.38)
  
})

#' @editor Lukas D Sauer
#' @editDate 2022-12-07
test_that("02.03", {
  # Adjustment method "both" for time-to-event endpoints
  res <- optimal_bias(alpha = 0.025,
                      beta = 0.1,
                      hr1 = 0.69, hr2 = 0.88,
                      xi2 = 0.7, xi3 = 0.7,
                      d2min = 20, d2max = 100, stepd2 = 5,
                      hrgomin = 0.7, hrgomax = 0.9, stephrgo = 0.02,
                      steps1 = 1, stepm1 = 0.95, stepl1 = 0.85,
                      b1 = 1000, b2 = 2000, b3 = 3000,
                      num_cl = 12,
                      c02 = 100, c03 = 150, # fixed cost
                      c2 = 0.75, c3 = 1, # variable per patient cost
                      fixed = FALSE,
                      w = 0.3,
                      id1 = 210, id2 = 420,
                      adj = "both",
                      alphaCImin = 0.3, alphaCImax = 0.5, stepalphaCI = 0.025,
                      lambdamin = 0.5, lambdamax = 1, steplambda = 0.05,
                      N = 350
  )
  expect_equal(res$Method, c("multipl.", "add."))
  # Additive method
  expect_equal(res[2,]$n2, 122)
  expect_equal(res[2,]$n3, 200)
  expect_equal(res[2,]$n, 322)
  expect_equal(res[2,]$u, 78, tolerance = 0.005)
  expect_equal(res[2,]$HRgo, 0.78)
  expect_equal(res[2,]$d2, 85)
  expect_equal(res[2,]$d3, 140)
  expect_equal(res[2,]$d, 225)
  # Multiplicative method
  expect_equal(res[1,]$n2, 100)
  expect_equal(res[1,]$n3, 240)
  expect_equal(res[1,]$n, 340)
  expect_equal(res[1,]$u, 98, tolerance = 0.0001)
  
})

#' @editor Lukas D Sauer
#' @editDate 2022-12-07
test_that("02.04", {
  # Adjustment method "all" for time-to-event endpoints
  res <- optimal_bias(alpha = 0.025,
                      beta = 0.1,
                      hr1 = 0.69, hr2 = 0.88,
                      xi2 = 0.7, xi3 = 0.7,
                      d2min = 20, d2max = 100, stepd2 = 5,
                      hrgomin = 0.7, hrgomax = 0.9, stephrgo = 0.02,
                      steps1 = 1, stepm1 = 0.95, stepl1 = 0.85,
                      b1 = 1000, b2 = 2000, b3 = 3000,
                      num_cl = 12,
                      c02 = 100, c03 = 150, # fixed cost
                      c2 = 0.75, c3 = 1, # variable per patient cost
                      fixed = FALSE,
                      w = 0.3,
                      id1 = 210, id2 = 420,
                      adj = "all",
                      alphaCImin = 0.3, alphaCImax = 0.5, stepalphaCI = 0.025,
                      lambdamin = 0.5, lambdamax = 1, steplambda = 0.05,
                      N = 350
  )
  expect_equal(res$Method, c("multipl.", "add.", "multipl2.", "add2."))
  # Additive
  expect_equal(res[2,]$n2, 122)
  expect_equal(res[2,]$n3, 200)
  expect_equal(res[2,]$n, 322)
  expect_equal(res[2,]$u, 78, tolerance = 0.005)
  expect_equal(res[2,]$HRgo, 0.78)
  expect_equal(res[2,]$d2, 85)
  expect_equal(res[2,]$d3, 140)
  expect_equal(res[2,]$d, 225)
  # Multiplicative
  expect_equal(res[1,]$n2, 100)
  expect_equal(res[1,]$n3, 240)
  expect_equal(res[1,]$n, 340)
  expect_equal(res[1,]$u, 98, tolerance = 0.0001)
  # Multiplicative 2
  expect_equal(res[3,]$n2, 108)
  expect_equal(res[3,]$n3, 200)
  expect_equal(res[3,]$n, 308)
  expect_equal(res[3,]$u, 96.69, tolerance = 0.0001)
  expect_equal(res[3,]$Adj, 0.75)
  # Additive 2
  expect_equal(res[4,]$n2, 136)
  expect_equal(res[4,]$n3, 206)
  expect_equal(res[4,]$n, 342)
  expect_equal(res[4,]$u, 77.40, tolerance = 0.0001)
  expect_equal(res[4,]$Adj, 0.475)
})
#' @editor Lukas D Sauer
#' @editDate 2022-12-13
test_that("02.05", {
  # Without and with cost constraint
  res <- optimal_bias(alpha = 0.025,
                      beta = 0.1,
                      hr1 = 0.69, hr2 = 0.88,
                      xi2 = 0.7, xi3 = 0.7,
                      d2min = 20, d2max = 100, stepd2 = 5,
                      hrgomin = 0.7, hrgomax = 0.9, stephrgo = 0.02,
                      steps1 = 1, stepm1 = 0.95, stepl1 = 0.85,
                      b1 = 1000, b2 = 2000, b3 = 3000,
                      num_cl = 12,
                      c02 = 100, c03 = 150, # fixed cost
                      c2 = 0.75, c3 = 1, # variable per patient cost
                      fixed = TRUE,
                      w = 0.3,
                      id1 = 210, id2 = 420,
                      adj = "additive",
                      alphaCImin = 0.3, alphaCImax = 0.5, stepalphaCI = 0.025,
                      lambdamin = NULL, lambdamax = NULL, steplambda = NULL
  )
  expect_equal(res$n2, 144)
  expect_equal(res$n3, 478)
  expect_equal(res$u, 865.02)
  expect_equal(res$Adj, 0.5)
  expect_equal(res$K2, 208)
  expect_equal(res$K3, 608)
  res <- optimal_bias(alpha = 0.025,
                      beta = 0.1,
                      hr1 = 0.69, hr2 = 0.88,
                      xi2 = 0.7, xi3 = 0.7,
                      d2min = 20, d2max = 100, stepd2 = 5,
                      hrgomin = 0.7, hrgomax = 0.9, stephrgo = 0.02,
                      steps1 = 1, stepm1 = 0.95, stepl1 = 0.85,
                      b1 = 1000, b2 = 2000, b3 = 3000,
                      num_cl = 12,
                      c02 = 100, c03 = 150, # fixed cost
                      c2 = 0.75, c3 = 1, # variable per patient cost
                      fixed = TRUE,
                      w = 0.3,
                      id1 = 210, id2 = 420,
                      adj = "additive",
                      alphaCImin = 0.3, alphaCImax = 0.5, stepalphaCI = 0.025,
                      lambdamin = NULL, lambdamax = NULL, steplambda = NULL,
                      K = 400
  )
  expect_equal(res$n2, 44)
  expect_equal(res$n3, 172)
  expect_equal(res$u, 474.18)
  expect_equal(res$Adj, 0.475)
  expect_equal(res$K2, 133)
  expect_equal(res$K3, 263)
  
})
#' @editor Lukas D Sauer
#' @editDate 2022-12-14
test_that("02.06", {
  # With and without probability constraint
  # for time-to-event endpoints
  res <- optimal_bias(alpha = 0.025,
                      beta = 0.1,
                      hr1 = 0.69, hr2 = 0.88,
                      xi2 = 0.7, xi3 = 0.7,
                      d2min = 20, d2max = 100, stepd2 = 5,
                      hrgomin = 0.7, hrgomax = 0.9, stephrgo = 0.02,
                      steps1 = 1, stepm1 = 0.95, stepl1 = 0.85,
                      b1 = 1000, b2 = 2000, b3 = 3000,
                      num_cl = 12,
                      c02 = 100, c03 = 150, # fixed cost
                      c2 = 0.75, c3 = 1, # variable per patient cost
                      fixed = TRUE,
                      w = 0.3,
                      id1 = 210, id2 = 420,
                      adj = "additive",
                      alphaCImin = 0.3, alphaCImax = 0.5, stepalphaCI = 0.025,
                      lambdamin = NULL, lambdamax = NULL, steplambda = NULL
  )
  expect_equal(res$n2, 144)
  expect_equal(res$n3, 478)
  expect_equal(res$u, 865.02)
  expect_equal(res$Adj, 0.5)
  expect_equal(res$sProg, 0.68)
  expect_equal(res$sProg1, 0.07)
  expect_equal(res$sProg2, 0.21)
  expect_equal(res$sProg3, 0.39)
  res <- optimal_bias(alpha = 0.025,
                      beta = 0.1,
                      hr1 = 0.69, hr2 = 0.88,
                      xi2 = 0.7, xi3 = 0.7,
                      d2min = 20, d2max = 100, stepd2 = 5,
                      hrgomin = 0.7, hrgomax = 0.9, stephrgo = 0.02,
                      steps1 = 1, stepm1 = 0.95, stepl1 = 0.85,
                      b1 = 1000, b2 = 2000, b3 = 3000,
                      num_cl = 12,
                      c02 = 100, c03 = 150, # fixed cost
                      c2 = 0.75, c3 = 1, # variable per patient cost
                      fixed = TRUE,
                      w = 0.3,
                      id1 = 210, id2 = 420,
                      adj = "additive",
                      alphaCImin = 0.3, alphaCImax = 0.5, stepalphaCI = 0.025,
                      lambdamin = NULL, lambdamax = NULL, steplambda = NULL,
                      S = 0.7
  )
  expect_equal(res$n2, 144)
  expect_equal(res$n3, 710)
  expect_equal(res$u, 769.80)
  expect_equal(res$Adj, 0.5)
  expect_equal(res$sProg, 0.7)
  expect_equal(res$sProg1, 0.07)
  expect_equal(res$sProg2, 0.21)
  expect_equal(res$sProg3, 0.42)
})
#' @editor Lukas D Sauer
#' @editDate 2022-12-14
test_that("02.07", {
  # No bias adjustment is equal to the basic setting
  res <- optimal_bias(alpha = 0.025,
                      beta = 0.1,
                      hr1 = 0.69, hr2 = 0.88,
                      xi2 = 0.7, xi3 = 0.7,
                      d2min = 20, d2max = 100, stepd2 = 5,
                      hrgomin = 0.7, hrgomax = 0.9, stephrgo = 0.02,
                      steps1 = 1, stepm1 = 0.95, stepl1 = 0.85,
                      b1 = 1000, b2 = 2000, b3 = 3000,
                      num_cl = 12,
                      c02 = 100, c03 = 150, # fixed cost
                      c2 = 0.75, c3 = 1, # variable per patient cost
                      fixed = FALSE,
                      w = 0.3,
                      id1 = 210, id2 = 420,
                      adj = "both",
                      alphaCImin = 0.5, alphaCImax = 0.5, stepalphaCI = 1,
                      lambdamin = 1, lambdamax = 1, steplambda = 1
  )
  expect_equal(res[1,]$n2, 122)
  expect_equal(res[1,]$n3, 210)
  expect_equal(res[1,]$n, 332)
  expect_equal(res[1,]$u, 75.8, tolerance = 0.005)
  expect_equal(res[2,]$n2, 122)
  expect_equal(res[2,]$n3, 210)
  expect_equal(res[2,]$n, 332)
  expect_equal(res[2,]$u, 75.8, tolerance = 0.005)
  res <- optimal_tte(alpha = 0.025,
                      beta = 0.1,
                      hr1 = 0.69, hr2 = 0.88,
                      xi2 = 0.7, xi3 = 0.7,
                      d2min = 20, d2max = 100, stepd2 = 5,
                      hrgomin = 0.7, hrgomax = 0.9, stephrgo = 0.02,
                      steps1 = 1, stepm1 = 0.95, stepl1 = 0.85,
                      b1 = 1000, b2 = 2000, b3 = 3000,
                      num_cl = 12,
                      c02 = 100, c03 = 150, # fixed cost
                      c2 = 0.75, c3 = 1, # variable per patient cost
                      fixed = FALSE,
                      w = 0.3,
                      id1 = 210, id2 = 420,
  )
  expect_equal(res$n2, 122)
  expect_equal(res$n3, 210)
  expect_equal(res$n, 332)
  expect_equal(res$u, 75.8, tolerance = 0.005)
})
#' @editor Lukas D Sauer
#' @editDate 2022-12-14
test_that("02.08", {
  # Parallel computing works
  start_time_3 = Sys.time()
  res <- optimal_bias(alpha = 0.025,
                      beta = 0.1,
                      hr1 = 0.69, hr2 = 0.88,
                      xi2 = 0.7, xi3 = 0.7,
                      d2min = 20, d2max = 100, stepd2 = 5,
                      hrgomin = 0.7, hrgomax = 0.9, stephrgo = 0.02,
                      steps1 = 1, stepm1 = 0.95, stepl1 = 0.85,
                      b1 = 1000, b2 = 2000, b3 = 3000,
                      num_cl = 12,
                      c02 = 100, c03 = 150, # fixed cost
                      c2 = 0.75, c3 = 1, # variable per patient cost
                      fixed = FALSE,
                      w = 0.3,
                      id1 = 210, id2 = 420,
                      adj = "additive",
                      alphaCImin = 0.3, alphaCImax = 0.5, stepalphaCI = 0.025,
                      lambdamin = NULL, lambdamax = NULL, steplambda = NULL
  )
  end_time_3 = Sys.time()
  time_elapsed_num_cl_3 = end_time_3 - start_time_3
  start_time_1 = Sys.time()
  res <- optimal_bias(alpha = 0.025,
                      beta = 0.1,
                      hr1 = 0.69, hr2 = 0.88,
                      xi2 = 0.7, xi3 = 0.7,
                      d2min = 20, d2max = 100, stepd2 = 5,
                      hrgomin = 0.7, hrgomax = 0.9, stephrgo = 0.02,
                      steps1 = 1, stepm1 = 0.95, stepl1 = 0.85,
                      b1 = 1000, b2 = 2000, b3 = 3000,
                      num_cl = 6,
                      c02 = 100, c03 = 150, # fixed cost
                      c2 = 0.75, c3 = 1, # variable per patient cost
                      fixed = FALSE,
                      w = 0.3,
                      id1 = 210, id2 = 420,
                      adj = "additive",
                      alphaCImin = 0.3, alphaCImax = 0.5, stepalphaCI = 0.025,
                      lambdamin = NULL, lambdamax = NULL, steplambda = NULL
  )
  end_time_1 = Sys.time()
  time_elapsed_num_cl_1 = end_time_1 - start_time_1
  expect_true(time_elapsed_num_cl_1 > time_elapsed_num_cl_3)
})
#' @editor Lukas D Sauer
#' @editDate 2022-12-14
test_that("02.09", {
  # Bias adjustment for normally distributed endpoints
  res <- optimal_bias_normal(
    alpha = 0.05,
    beta = 0.1,
    Delta1 = 0.625, Delta2 = 0.325,
    n2min = 20, n2max = 400, stepn2 = 4,
    kappamin = 0.02, kappamax = 0.4, stepkappa = 0.02,
    steps1 = 0, stepm1 = 0.5, stepl1 = 0.8,
    b1 = 3000, b2 = 8000, b3 = 10000,
    num_cl = 12,
    c02 = 15, c03 = 20,
    c2 = 0.675, c3 = 0.72,
    fixed = FALSE,
    w = 0.5,
    in1 = 300, in2 = 600,
    a = 0.25, b = 0.75,
    adj = "multiplicative",
    lambdamin = 0.7, lambdamax = 0.9, steplambda = 0.01
  )
  expect_equal(res$n2, 192)
  expect_equal(res$n3, 474)
  expect_equal(res$n, 666)
  expect_equal(res$Kappa, 0.12)
  expect_equal(res$u, 2899.32, tolerance = 0.005)
})
#' @editor Lukas D Sauer
#' @editDate 2022-12-14
test_that("02.10", {
  # Bias adjustment for normally distributed endpoints
  # with fixed prior treatment effect and all adjustment
  # methods
  res <- optimal_bias_normal(
    alpha = 0.05,
    beta = 0.1,
    Delta1 = 0.625, Delta2 = 0.325,
    n2min = 20, n2max = 400, stepn2 = 4,
    kappamin = 0.02, kappamax = 0.4, stepkappa = 0.02,
    steps1 = 0, stepm1 = 0.5, stepl1 = 0.8,
    b1 = 3000, b2 = 8000, b3 = 10000,
    num_cl = 12,
    c02 = 15, c03 = 20,
    c2 = 0.675, c3 = 0.72,
    fixed = TRUE,
    w = 0.5,
    in1 = 300, in2 = 600,
    a = 0.25, b = 0.75,
    adj = "all",
    alphaCImin = 0.25, alphaCImax = 0.5, stepalphaCI = 0.025,
    lambdamin = 0.7, lambdamax = 0.9, steplambda = 0.01
  )
  expect_equal(res$Method, c("multipl.", "add.", "multipl2.", "add2."))
  # Multiplicative
  expect_equal(res[1,]$n2, 88)
  expect_equal(res[1,]$n3, 310)
  expect_equal(res[1,]$n, 398)
  expect_equal(res[1,]$u, 3861.76, tolerance = 0.0001)
  expect_equal(res[1,]$Adj, 0.7, tolerance = 0.0001)
  # Additive
  expect_equal(res[2,]$n2, 96)
  expect_equal(res[2,]$n3, 306)
  expect_equal(res[2,]$n, 402)
  expect_equal(res[2,]$u, 3631.51, tolerance = 0.0001)
  expect_equal(res[2,]$Adj, 0.25, tolerance = 0.0001)
  # Multiplicative 2
  expect_equal(res[3,]$n2, 88)
  expect_equal(res[3,]$n3, 306)
  expect_equal(res[3,]$n, 394)
  expect_equal(res[3,]$u, 3860.12, tolerance = 0.0001)
  expect_equal(res[3,]$Adj, 0.7)
  # Additive 2
  expect_equal(res[4,]$n2, 96)
  expect_equal(res[4,]$n3, 312)
  expect_equal(res[4,]$n, 408)
  expect_equal(res[4,]$u, 3631.28, tolerance = 0.0001)
  expect_equal(res[4,]$Adj, 0.25)
})
#' @editor Lukas D Sauer
#' @editDate 2022-12-20
test_that("02.11", {
  res <- optimal_bias_binary(
    alpha = 0.025,
    beta = 0.1,
    p0 = 0.6, p11 = 0.3, p12= 0.5,
    n2min = 10, n2max = 500, stepn2 = 2,
    rrgomin = 0.7, rrgomax = 0.9, steprrgo = 0.01,
    steps1 = 1, stepm1 = 0.95, stepl1 = 0.85,
    b1 = 1000, b2 = 2000, b3 = 3000,
    num_cl = 12,
    c02 = 100, c03 = 150,
    c2 = 0.75, c3 = 1,
    fixed = FALSE,
    w = 0.3,
    in1 = 30, in2 = 60,
    adj = "additive",
    alphaCImin = 0.1, alphaCImax = 0.5, stepalphaCI = 0.025,
  )
  expect_equal(res$n2, 166)
  expect_equal(res$n3, 264)
  expect_equal(res$n, 430)
  expect_equal(res$RRgo, 0.82)
  expect_equal(res$u, 605.91, tolerance = 0.005)
  expect_equal(res$Adj, 0.275)
})
#' @editor Lukas D Sauer
#' @editDate 2022-12-20
test_that("02.12", {
  res <- optimal_bias_binary(
    alpha = 0.025,
    beta = 0.1,
    p0 = 0.6, p11 = 0.3, p12= 0.5,
    n2min = 10, n2max = 500, stepn2 = 2,
    rrgomin = 0.7, rrgomax = 0.9, steprrgo = 0.01,
    steps1 = 1, stepm1 = 0.95, stepl1 = 0.85,
    b1 = 1000, b2 = 2000, b3 = 3000,
    num_cl = 12,
    c02 = 100, c03 = 150,
    c2 = 0.75, c3 = 1,
    fixed = TRUE,
    w = 0.3,
    in1 = 30, in2 = 60,
    adj = "both",
    alphaCImin = 0.1, alphaCImax = 0.5, stepalphaCI = 0.025,
    lambdamin = 0.5, lambdamax = 1, steplambda = 0.05
  )
  expect_equal(res[1,]$Method, "multipl.")
  expect_equal(res[1,]$n2, 206)
  expect_equal(res[1,]$n3, 340)
  expect_equal(res[1,]$n, 546)
  expect_equal(res[1,]$u, 2116.67, tolerance = 0.005)
  expect_equal(res[1,]$RRgo, 0.8)
  expect_equal(res[1,]$Adj, 0.65)
  expect_equal(res[2,]$Method, "add.")
  expect_equal(res[2,]$n2, 190)
  expect_equal(res[2,]$n3, 226)
  expect_equal(res[2,]$n, 416)
  expect_equal(res[2,]$u, 1996.10, tolerance = 0.005)
  expect_equal(res[2,]$RRgo, 0.77)
  expect_equal(res[2,]$Adj, 0.1)
})
