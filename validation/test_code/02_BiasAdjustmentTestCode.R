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
  expect_equal(res$Method, c("multipl.", "add.", "multipl2", "add2"))
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

