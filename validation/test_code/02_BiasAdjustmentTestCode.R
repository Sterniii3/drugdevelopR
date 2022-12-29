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
               num_cl = 3,
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
                      num_cl = 3,
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
                      num_cl = 3,
                      c02 = 100, c03 = 150, # fixed cost
                      c2 = 0.75, c3 = 1, # variable per patient cost
                      fixed = FALSE,
                      w = 0.3,
                      id1 = 210, id2 = 420,
                      adj = "both",
                      alphaCImin = 0.3, alphaCImax = 0.5, stepalphaCI = 0.025,
                      lambdamin = 0.5, lambdamax = 1, steplambda = 0.02,
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
                      num_cl = 3,
                      c02 = 100, c03 = 150, # fixed cost
                      c2 = 0.75, c3 = 1, # variable per patient cost
                      fixed = FALSE,
                      w = 0.3,
                      id1 = 210, id2 = 420,
                      adj = "all",
                      alphaCImin = 0.3, alphaCImax = 0.5, stepalphaCI = 0.025,
                      lambdamin = 0.5, lambdamax = 1, steplambda = 0.02,
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
  expect_equal(res[3,]$n2, 122)
  expect_equal(res[3,]$n3, 204)
  expect_equal(res[3,]$n, 326)
  expect_equal(res[3,]$u, 97.00, tolerance = 0.0001)
  expect_equal(res[3,]$Adj, 0.78)
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
                      num_cl = 3,
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
  expect_equal(res$n3, 456)
  expect_equal(res$u, 853.02)
  expect_equal(res$Adj, 0.475)
  expect_equal(res$K2, 208)
  expect_equal(res$K3, 582)
  res <- optimal_bias(alpha = 0.025,
                      beta = 0.1,
                      hr1 = 0.69, hr2 = 0.88,
                      xi2 = 0.7, xi3 = 0.7,
                      d2min = 20, d2max = 100, stepd2 = 5,
                      hrgomin = 0.7, hrgomax = 0.9, stephrgo = 0.02,
                      steps1 = 1, stepm1 = 0.95, stepl1 = 0.85,
                      b1 = 1000, b2 = 2000, b3 = 3000,
                      num_cl = 3,
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
                      num_cl = 3,
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
  expect_equal(res$n3, 456)
  expect_equal(res$u, 853.02)
  expect_equal(res$Adj, 0.475)
  expect_equal(res$sProg, 0.66)
  expect_equal(res$sProg1, 0.07)
  expect_equal(res$sProg2, 0.21)
  expect_equal(res$sProg3, 0.38)
  res <- optimal_bias(alpha = 0.025,
                      beta = 0.1,
                      hr1 = 0.69, hr2 = 0.88,
                      xi2 = 0.7, xi3 = 0.7,
                      d2min = 20, d2max = 100, stepd2 = 5,
                      hrgomin = 0.7, hrgomax = 0.9, stephrgo = 0.02,
                      steps1 = 1, stepm1 = 0.95, stepl1 = 0.85,
                      b1 = 1000, b2 = 2000, b3 = 3000,
                      num_cl = 3,
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
  expect_equal(res$Adj, 0.375)
  expect_equal(res$sProg, 0.71)
  expect_equal(res$sProg1, 0.06)
  expect_equal(res$sProg2, 0.19)
  expect_equal(res$sProg3, 0.46)
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
                      num_cl = 3,
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
                      num_cl = 3,
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
                      num_cl = 3,
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
                      num_cl = 1,
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
    num_cl = 3,
    c02 = 15, c03 = 20,
    c2 = 0.675, c3 = 0.72,
    fixed = FALSE,
    w = 0.5,
    in1 = 300, in2 = 600,
    a = 0.25, b = 0.75,
    adj = "multiplicative",
    lambdamin = 0.7, lambdamax = 0.9, steplambda = 0.01
  )
  expect_equal(res$n2, 236)
  expect_equal(res$n3, 548)
  expect_equal(res$n, 784)
  expect_equal(res$Kappa, 0.12)
  expect_equal(res$u, 2779.11, tolerance = 0.005)
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
    num_cl = 3,
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
  expect_equal(res[1,]$n2, 102)
  expect_equal(res[1,]$n3, 544)
  expect_equal(res[1,]$n, 646)
  expect_equal(res[1,]$u, 4336.42, tolerance = 0.0001)
  expect_equal(res[1,]$Adj, 0.5, tolerance = 0.0001)
  # Additive
  expect_equal(res[2,]$n2, 130)
  expect_equal(res[2,]$n3, 456)
  expect_equal(res[2,]$n, 586)
  expect_equal(res[2,]$u, 3871.26, tolerance = 0.0001)
  expect_equal(res[2,]$Adj, 0.1, tolerance = 0.0001)
  # Multiplicative 2
  expect_equal(res[3,]$n2, 98)
  expect_equal(res[3,]$n3, 536)
  expect_equal(res[3,]$n, 634)
  expect_equal(res[3,]$u, 4336.37, tolerance = 0.0001)
  expect_equal(res[3,]$Adj, 0.5)
  # Additive 2
  expect_equal(res[4,]$n2, 134)
  expect_equal(res[4,]$n3, 426)
  expect_equal(res[4,]$n, 560)
  expect_equal(res[4,]$u, 3870.71, tolerance = 0.0001)
  expect_equal(res[4,]$Adj, 0.1)
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
    num_cl = 3,
    c02 = 100, c03 = 150,
    c2 = 0.75, c3 = 1,
    fixed = FALSE,
    w = 0.3,
    in1 = 30, in2 = 60,
    adj = "additive",
    alphaCImin = 0.1, alphaCImax = 0.5, stepalphaCI = 0.025,
  )
  expect_equal(res$n2, 158)
  expect_equal(res$n3, 262)
  expect_equal(res$n, 420)
  expect_equal(res$RRgo, 0.86)
  expect_equal(res$u, 708.24, tolerance = 0.005)
  expect_equal(res$Adj, 0.4)
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
    num_cl = 3,
    c02 = 100, c03 = 150,
    c2 = 0.75, c3 = 1,
    fixed = TRUE,
    w = 0.3,
    in1 = 30, in2 = 60,
    adj = "both",
    alphaCImin = 0.1, alphaCImax = 0.5, stepalphaCI = 0.025,
    lambdamin = 0.5, lambdamax = 1, steplambda = 0.02
  )
  expect_equal(res[1,]$Method, "multipl.")
  expect_equal(res[1,]$n2, 198)
  expect_equal(res[1,]$n3, 294)
  expect_equal(res[1,]$n, 492)
  expect_equal(res[1,]$u, 2180.86, tolerance = 0.005)
  expect_equal(res[1,]$RRgo, 0.82)
  expect_equal(res[1,]$Adj, 0.64)
  expect_equal(res[2,]$Method, "add.")
  expect_equal(res[2,]$n2, 178)
  expect_equal(res[2,]$n3, 196)
  expect_equal(res[2,]$n, 374)
  expect_equal(res[2,]$u, 2062.42, tolerance = 0.005)
  expect_equal(res[2,]$RRgo, 0.82)
  expect_equal(res[2,]$Adj, 0.1)
})
