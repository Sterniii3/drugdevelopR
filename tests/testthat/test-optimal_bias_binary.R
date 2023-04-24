test_that("optimal_bias_binary works for fixed treatment effects", {
  skip_on_cran()
  expect_equal(optimal_bias_binary(w = 0.3, p0 = 0.6, p11 = 0.3, p12 = 0.5,
                                   in1 = 30, in2 = 60,
                                   n2min = 20, n2max = 100, stepn2 = 4, 
                                   rrgomin = 0.7, rrgomax = 0.9, steprrgo = 0.05,
                                   adj = "both", alpha = 0.05, beta = 0.1, 
                                   lambdamin = 0.5, lambdamax = 1, steplambda = 0.05,
                                   alphaCImin = 0.025, alphaCImax = 0.5, stepalphaCI = 0.025,
                                   c2 = 0.75, c3 = 1, c02 = 100, c03 = 150,
                                   K = Inf, N = Inf, S = -Inf,
                                   steps1 = 1, stepm1 = 0.95, stepl1 = 0.85,
                                   b1 = 1000, b2 = 2000, b3 = 3000,
                                   fixed = TRUE, num_cl = 2)$u, c(2110.63,1971.77))
})

test_that("optimal_bias_binary works when using a prior distribution", {
  skip_on_cran()
  expect_equal(optimal_bias_binary(w = 0.3, p0 = 0.6, p11 = 0.3, p12 = 0.5,
                                   in1 = 30, in2 = 60,
                                   n2min = 20, n2max = 100, stepn2 = 4, 
                                   rrgomin = 0.7, rrgomax = 0.9, steprrgo = 0.05,
                                   adj = "additive", alpha = 0.05, beta = 0.1, 
                                   lambdamin = 0.5, lambdamax = 1, steplambda = 0.05,
                                   alphaCImin = 0.025, alphaCImax = 0.5, stepalphaCI = 0.025,
                                   c2 = 0.75, c3 = 1, c02 = 100, c03 = 150,
                                   K = Inf, N = Inf, S = -Inf,
                                   steps1 = 1, stepm1 = 0.95, stepl1 = 0.85,
                                   b1 = 1000, b2 = 2000, b3 = 3000,
                                   fixed = FALSE, num_cl = 2)$u, 688.38)
})

test_that("optimal_bias_binary works for fixed treatment effects with method all", {
  skip_on_cran()
  expect_equal(optimal_bias_binary(w = 0.3, p0 = 0.6, p11 = 0.3, p12 = 0.5,
                                   in1 = 30, in2 = 60,
                                   n2min = 20, n2max = 100, stepn2 = 4, 
                                   rrgomin = 0.7, rrgomax = 0.9, steprrgo = 0.05,
                                   adj = "all", alpha = 0.05, beta = 0.1, 
                                   lambdamin = 0.5, lambdamax = 1, steplambda = 0.05,
                                   alphaCImin = 0.025, alphaCImax = 0.5, stepalphaCI = 0.025,
                                   c2 = 0.75, c3 = 1, c02 = 100, c03 = 150,
                                   K = Inf, N = Inf, S = -Inf,
                                   steps1 = 1, stepm1 = 0.95, stepl1 = 0.85,
                                   b1 = 1000, b2 = 2000, b3 = 3000,
                                   fixed = TRUE, num_cl = 2)$u, c(2110.63,1971.77,2110.05,1889.68))
})