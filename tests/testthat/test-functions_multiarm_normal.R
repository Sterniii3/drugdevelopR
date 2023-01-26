test_that("Utility increases with higher treatment effect", {
  expect_equal(utility_multiarm_normal(n2 = 50, kappa = 0.8, alpha = 0.05, beta = 0.1,
                                       Delta1 = 0.375, Delta2 = 0.625, strategy = 1,
                                       c2 = 0.75, c3 = 1, c02 = 100, c03 = 150,
                                       K = Inf, N = Inf, S = -Inf, 
                                       steps1 = 0, stepm1 = 0.5,  stepl1 = 0.8,
                                       b1 = 1000, b2 = 2000, b3 = 3000)[1], 
               utility_multiarm_normal(n2 = 50, kappa = 0.8, alpha = 0.05, beta = 0.1,
                                       Delta1 = 0.375, Delta2 = 0.625, strategy = 1,
                                       c2 = 0.75, c3 = 1, c02 = 100, c03 = 150,
                                       K = Inf, N = Inf, S = -Inf, 
                                       steps1 = 0, stepm1 = 0.5,  stepl1 = 0.8,
                                       b1 = 1000, b2 = 2000, b3 = 3000)[1])
})

test_that("Utility increases with higher treatment effect", {
  expect_equal(utility_multiarm_normal(n2 = 50, kappa = 0.8, alpha = 0.05, beta = 0.1,
                                       Delta1 = 0.375, Delta2 = 0.625, strategy = 2,
                                       c2 = 0.75, c3 = 1, c02 = 100, c03 = 150,
                                       K = Inf, N = Inf, S = -Inf, 
                                       steps1 = 0, stepm1 = 0.5,  stepl1 = 0.8,
                                       b1 = 1000, b2 = 2000, b3 = 3000)[1], 
               utility_multiarm_normal(n2 = 50, kappa = 0.8, alpha = 0.05, beta = 0.1,
                                       Delta1 = 0.375, Delta2 = 0.625, strategy = 2,
                                       c2 = 0.75, c3 = 1, c02 = 100, c03 = 150,
                                       K = Inf, N = Inf, S = -Inf, 
                                       steps1 = 0, stepm1 = 0.5,  stepl1 = 0.8,
                                       b1 = 1000, b2 = 2000, b3 = 3000)[1])
})
