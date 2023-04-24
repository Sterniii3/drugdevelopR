test_that("optimal_bias_normal works for true treatment effects", {
  skip_on_cran()
  expect_equal(optimal_bias_normal(w=0.3, Delta1 = 0.375, Delta2 = 0.625,
                                   in1=300, in2=600, a = 0.25, b = 0.75,
                                   n2min = 20, n2max = 100, stepn2 = 4,
                                   kappamin = 0.02, kappamax = 0.2, stepkappa = 0.02,
                                   adj = "both", alpha = 0.05, beta = 0.1,
                                   lambdamin = 0.6, lambdamax = 1, steplambda = 0.05, 
                                   alphaCImin = 0.25, alphaCImax = 0.5, stepalphaCI = 0.025,
                                   c2 = 0.675, c3 = 0.72, c02 = 15, c03 = 20,
                                   K = Inf, N = Inf, S = -Inf,
                                   steps1 = 0, stepm1 = 0.5, stepl1 = 0.8,
                                   b1 = 3000, b2 = 8000, b3 = 10000,
                                   fixed = TRUE,num_cl = 2)$u, c(1968.77,1890.91))
})

test_that("optimal_bias_normal works when using a prior distribution", {
  skip_on_cran()
  expect_equal(optimal_bias_normal(w=0.3, Delta1 = 0.375, Delta2 = 0.625,
                                   in1=300, in2=600, a = 0.25, b = 0.75,
                                   n2min = 80, n2max = 100, stepn2 = 4,
                                   kappamin = 0.02, kappamax = 0.2, stepkappa = 0.02,
                                   adj = "additive", alpha = 0.05, beta = 0.1,
                                   lambdamin = 0.6, lambdamax = 1, steplambda = 0.05, 
                                   alphaCImin = 0.25, alphaCImax = 0.5, stepalphaCI = 0.025,
                                   c2 = 0.675, c3 = 0.72, c02 = 15, c03 = 20,
                                   K = Inf, N = Inf, S = -Inf,
                                   steps1 = 0, stepm1 = 0.5, stepl1 = 0.8,
                                   b1 = 3000, b2 = 8000, b3 = 10000,
                                   fixed = FALSE,num_cl = 2)$u, 2274.03)
})

test_that("optimal_bias_normal works for method all", {
  skip_on_cran()
  expect_equal(optimal_bias_normal(w=0.3, Delta1 = 0.375, Delta2 = 0.625,
                                   in1=300, in2=600, a = 0.25, b = 0.75,
                                   n2min = 80, n2max = 100, stepn2 = 4,
                                   kappamin = 0.02, kappamax = 0.2, stepkappa = 0.02,
                                   adj = "all", alpha = 0.05, beta = 0.1,
                                   lambdamin = 0.6, lambdamax = 1, steplambda = 0.05, 
                                   alphaCImin = 0.25, alphaCImax = 0.5, stepalphaCI = 0.025,
                                   c2 = 0.675, c3 = 0.72, c02 = 15, c03 = 20,
                                   K = Inf, N = Inf, S = -Inf,
                                   steps1 = 0, stepm1 = 0.5, stepl1 = 0.8,
                                   b1 = 3000, b2 = 8000, b3 = 10000,
                                   fixed = TRUE,num_cl = 2)$u, c(1968.77,1890.91,1966.07,1890.91))
})