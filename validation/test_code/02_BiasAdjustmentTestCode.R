#' @editor Lukas D Sauer
#' @editDate 2022-12-07
test_that("02.01", {
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



