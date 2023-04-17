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

test_that("pgo_normal: Setting 1.21", {
  expect_equal(pgo_normal(kappa = 0.1, n2 = 50, Delta1 = 0.375, Delta2 = 0.625,
                          strategy = 1, case = 21), 0.22276717)
})

test_that("pgo_normal: Setting 1.22", {
  expect_equal(pgo_normal(kappa = 0.1, n2 = 50, Delta1 = 0.375, Delta2 = 0.625,
                          strategy = 1, case = 22), 0.73960867)
})

test_that("pgo_normal: Setting 2.31", {
  expect_equal(pgo_normal(kappa = 0.1, n2 = 50, Delta1 = 0.375, Delta2 = 0.625,
                          strategy = 2, case = 31), 0.19557404)
})

test_that("pgo_normal: Setting 2.32", {
  expect_equal(pgo_normal(kappa = 0.1, n2 = 50, Delta1 = 0.375, Delta2 = 0.625,
                          strategy = 2, case = 32), 0.56359338)
})
