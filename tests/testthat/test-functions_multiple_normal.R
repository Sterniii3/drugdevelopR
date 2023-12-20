test_that("Utility increases with higher treatment effect", {
  skip_on_cran()
  set.seed(61216)
  expect_lte(utility_multiple_normal(kappa = 0.1, n2 = 50, alpha = 0.025, beta = 0.1,
                                       Delta1 = 0.375, Delta2 = 0.625, in1 = 300, in2 = 600,
                                       sigma1 = 2, sigma2 = 1, 
                                       c2 = 0.75, c3 = 1, c02 = 100, c03 = 150, 
                                       K = Inf, N = Inf, S = -Inf,
                                       steps1 = 0, stepm1 = 0.5, stepl1 = 0.8,  
                                       b1 = 1000, b2 = 2000, b3 = 3000, 
                                       fixed = TRUE, rho = 0.5, relaxed = "TRUE")[1], 
               utility_multiple_normal(kappa = 0.1, n2 = 50, alpha = 0.025, beta = 0.1,
                                       Delta1 = 0.625, Delta2 = 0.825, in1 = 300, in2 = 600,
                                       sigma1 = 2, sigma2 = 1, 
                                       c2 = 0.75, c3 = 1, c02 = 100, c03 = 150, 
                                       K = Inf, N = Inf, S = -Inf,
                                       steps1 = 0, stepm1 = 0.5, stepl1 = 0.8,  
                                       b1 = 1000, b2 = 2000, b3 = 3000, 
                                       fixed = TRUE, rho = 0.5, relaxed = "TRUE")[1])
})

test_that("Ess_multiple_normal works", {
  skip_on_cran()
  set.seed(61216)
  expect_equal(Ess_multiple_normal(kappa = 0.1, n2 = 50, alpha = 0.025, beta = 0.1,
                                 Delta1 = 0.375, Delta2 = 0.625, in1 = 300, in2 = 600,
                                 sigma1 = 2, sigma2 = 1, fixed = FALSE, rho = 0.3,
                                 rsamp = get_sample_multiple_normal(Delta1 = 0.375,
                                                                    Delta2 = 0.625,
                                                                    in1 = 300,
                                                                    in2 = 600,
                                                                    rho = 0.3)),
               375,
               tolerance = 1e-03)
})

