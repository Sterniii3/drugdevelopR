test_that("Utility increases for higher treatment effect", {
  expect_lte(utility_normal_L(kappa = 0.1, n2 = 50, Adj = 0.8, 
                                alpha = 0.025, beta = 0.1, w = 0.3,
                                Delta1 = 0.375, Delta2 = 0.625, in1 = 300, in2 = 600, 
                                a = 0.25, b = 0.75, 
                                c2 = 0.75, c3 = 1, c02 = 100, c03 = 150,
                                K = Inf, N = Inf, S = -Inf, 
                                steps1 = 0, stepm1 = 0.5, stepl1 = 0.8,
                                b1 = 3000, b2 = 8000, b3 = 10000,fixed = TRUE)[1],
             utility_normal_L(kappa = 0.1, n2 = 50, Adj = 0.8, 
                              alpha = 0.025, beta = 0.1, w = 0.3,
                              Delta1 = 0.575, Delta2 = 0.825, in1 = 300, in2 = 600, 
                              a = 0.25, b = 0.75, 
                              c2 = 0.75, c3 = 1, c02 = 100, c03 = 150,
                              K = Inf, N = Inf, S = -Inf, 
                              steps1 = 0, stepm1 = 0.5, stepl1 = 0.8,
                              b1 = 3000, b2 = 8000, b3 = 10000,fixed = TRUE)[1])
})

test_that("Utility increases for higher treatment effect", {
  expect_lte(utility_normal_L(kappa = 0.1, n2 = 50, Adj = 0.8, 
                              alpha = 0.025, beta = 0.1, w = 0.3,
                              Delta1 = 0.375, Delta2 = 0.625, in1 = 300, in2 = 600, 
                              a = 0.25, b = 0.75, 
                              c2 = 0.75, c3 = 1, c02 = 100, c03 = 150,
                              K = Inf, N = Inf, S = -Inf, 
                              steps1 = 0, stepm1 = 0.5, stepl1 = 0.8,
                              b1 = 3000, b2 = 8000, b3 = 10000,fixed = FALSE)[1],
             utility_normal_L(kappa = 0.1, n2 = 50, Adj = 0.8, 
                              alpha = 0.025, beta = 0.1, w = 0.3,
                              Delta1 = 0.575, Delta2 = 0.825, in1 = 300, in2 = 600, 
                              a = 0.25, b = 0.75, 
                              c2 = 0.75, c3 = 1, c02 = 100, c03 = 150,
                              K = Inf, N = Inf, S = -Inf, 
                              steps1 = 0, stepm1 = 0.5, stepl1 = 0.8,
                              b1 = 3000, b2 = 8000, b3 = 10000,fixed = FALSE)[1])
})

test_that("Utility increases for higher treatment effect", {
  expect_lte(utility_normal_R(kappa = 0.1, n2 = 50, Adj = 0.4, 
                              alpha = 0.025, beta = 0.1, w = 0.3,
                              Delta1 = 0.375, Delta2 = 0.625, in1 = 300, in2 = 600, 
                              a = 0.25, b = 0.75, 
                              c2 = 0.75, c3 = 1, c02 = 100, c03 = 150,
                              K = Inf, N = Inf, S = -Inf, 
                              steps1 = 0, stepm1 = 0.5, stepl1 = 0.8,
                              b1 = 3000, b2 = 8000, b3 = 10000,fixed = TRUE)[1],
             utility_normal_R(kappa = 0.1, n2 = 50, Adj = 0.4, 
                              alpha = 0.025, beta = 0.1, w = 0.3,
                              Delta1 = 0.575, Delta2 = 0.825, in1 = 300, in2 = 600, 
                              a = 0.25, b = 0.75, 
                              c2 = 0.75, c3 = 1, c02 = 100, c03 = 150,
                              K = Inf, N = Inf, S = -Inf, 
                              steps1 = 0, stepm1 = 0.5, stepl1 = 0.8,
                              b1 = 3000, b2 = 8000, b3 = 10000,fixed = TRUE)[1])
})

test_that("Utility increases for higher treatment effect", {
  expect_lte(utility_normal_R(kappa = 0.1, n2 = 50, Adj = 0.4, 
                              alpha = 0.025, beta = 0.1, w = 0.3,
                              Delta1 = 0.375, Delta2 = 0.625, in1 = 300, in2 = 600, 
                              a = 0.25, b = 0.75, 
                              c2 = 0.75, c3 = 1, c02 = 100, c03 = 150,
                              K = Inf, N = Inf, S = -Inf, 
                              steps1 = 0, stepm1 = 0.5, stepl1 = 0.8,
                              b1 = 3000, b2 = 8000, b3 = 10000,fixed = FALSE)[1],
             utility_normal_R(kappa = 0.1, n2 = 50, Adj = 0.4, 
                              alpha = 0.025, beta = 0.1, w = 0.3,
                              Delta1 = 0.575, Delta2 = 0.825, in1 = 300, in2 = 600, 
                              a = 0.25, b = 0.75, 
                              c2 = 0.75, c3 = 1, c02 = 100, c03 = 150,
                              K = Inf, N = Inf, S = -Inf, 
                              steps1 = 0, stepm1 = 0.5, stepl1 = 0.8,
                              b1 = 3000, b2 = 8000, b3 = 10000,fixed = FALSE)[1])
})

test_that("Utility increases for higher treatment effect", {
  expect_lte(utility_normal_L2(kappa = 0.1, n2 = 50, Adj = 0.8, 
                              alpha = 0.025, beta = 0.1, w = 0.3,
                              Delta1 = 0.375, Delta2 = 0.625, in1 = 300, in2 = 600, 
                              a = 0.25, b = 0.75, 
                              c2 = 0.75, c3 = 1, c02 = 100, c03 = 150,
                              K = Inf, N = Inf, S = -Inf, 
                              steps1 = 0, stepm1 = 0.5, stepl1 = 0.8,
                              b1 = 3000, b2 = 8000, b3 = 10000,fixed = TRUE)[1],
             utility_normal_L2(kappa = 0.1, n2 = 50, Adj = 0.8, 
                              alpha = 0.025, beta = 0.1, w = 0.3,
                              Delta1 = 0.575, Delta2 = 0.825, in1 = 300, in2 = 600, 
                              a = 0.25, b = 0.75, 
                              c2 = 0.75, c3 = 1, c02 = 100, c03 = 150,
                              K = Inf, N = Inf, S = -Inf, 
                              steps1 = 0, stepm1 = 0.5, stepl1 = 0.8,
                              b1 = 3000, b2 = 8000, b3 = 10000,fixed = TRUE)[1])
})

test_that("Utility increases for higher treatment effect", {
  expect_lte(utility_normal_L2(kappa = 0.1, n2 = 50, Adj = 0.8, 
                              alpha = 0.025, beta = 0.1, w = 0.3,
                              Delta1 = 0.375, Delta2 = 0.625, in1 = 300, in2 = 600, 
                              a = 0.25, b = 0.75, 
                              c2 = 0.75, c3 = 1, c02 = 100, c03 = 150,
                              K = Inf, N = Inf, S = -Inf, 
                              steps1 = 0, stepm1 = 0.5, stepl1 = 0.8,
                              b1 = 3000, b2 = 8000, b3 = 10000,fixed = FALSE)[1],
             utility_normal_L2(kappa = 0.1, n2 = 50, Adj = 0.8, 
                              alpha = 0.025, beta = 0.1, w = 0.3,
                              Delta1 = 0.575, Delta2 = 0.825, in1 = 300, in2 = 600, 
                              a = 0.25, b = 0.75, 
                              c2 = 0.75, c3 = 1, c02 = 100, c03 = 150,
                              K = Inf, N = Inf, S = -Inf, 
                              steps1 = 0, stepm1 = 0.5, stepl1 = 0.8,
                              b1 = 3000, b2 = 8000, b3 = 10000,fixed = FALSE)[1])
})

test_that("Utility increases for higher treatment effect", {
  expect_lte(utility_normal_R2(kappa = 0.1, n2 = 50, Adj = 0.4, 
                              alpha = 0.025, beta = 0.1, w = 0.3,
                              Delta1 = 0.375, Delta2 = 0.625, in1 = 300, in2 = 600, 
                              a = 0.25, b = 0.75, 
                              c2 = 0.75, c3 = 1, c02 = 100, c03 = 150,
                              K = Inf, N = Inf, S = -Inf, 
                              steps1 = 0, stepm1 = 0.5, stepl1 = 0.8,
                              b1 = 3000, b2 = 8000, b3 = 10000,fixed = TRUE)[1],
             utility_normal_R2(kappa = 0.1, n2 = 50, Adj = 0.4, 
                              alpha = 0.025, beta = 0.1, w = 0.3,
                              Delta1 = 0.575, Delta2 = 0.825, in1 = 300, in2 = 600, 
                              a = 0.25, b = 0.75, 
                              c2 = 0.75, c3 = 1, c02 = 100, c03 = 150,
                              K = Inf, N = Inf, S = -Inf, 
                              steps1 = 0, stepm1 = 0.5, stepl1 = 0.8,
                              b1 = 3000, b2 = 8000, b3 = 10000,fixed = TRUE)[1])
})

test_that("Utility increases for higher treatment effect", {
  expect_lte(utility_normal_R2(kappa = 0.1, n2 = 50, Adj = 0.4, 
                              alpha = 0.025, beta = 0.1, w = 0.3,
                              Delta1 = 0.375, Delta2 = 0.625, in1 = 300, in2 = 600, 
                              a = 0.25, b = 0.75, 
                              c2 = 0.75, c3 = 1, c02 = 100, c03 = 150,
                              K = Inf, N = Inf, S = -Inf, 
                              steps1 = 0, stepm1 = 0.5, stepl1 = 0.8,
                              b1 = 3000, b2 = 8000, b3 = 10000,fixed = FALSE)[1],
             utility_normal_R2(kappa = 0.1, n2 = 50, Adj = 0.4, 
                              alpha = 0.025, beta = 0.1, w = 0.3,
                              Delta1 = 0.575, Delta2 = 0.825, in1 = 300, in2 = 600, 
                              a = 0.25, b = 0.75, 
                              c2 = 0.75, c3 = 1, c02 = 100, c03 = 150,
                              K = Inf, N = Inf, S = -Inf, 
                              steps1 = 0, stepm1 = 0.5, stepl1 = 0.8,
                              b1 = 3000, b2 = 8000, b3 = 10000,fixed = FALSE)[1])
})

