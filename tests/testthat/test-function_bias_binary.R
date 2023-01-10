test_that("Utility increases with smaller risk ratio", {
  expect_lte(utility_binary_L(n2 = 50, RRgo = 0.8, Adj = 0.1, w = 0.3,
                              p0 = 0.6, p11 =  0.3, p12 = 0.5, 
                              in1 = 300, in2 = 600, alpha = 0.025, beta = 0.1,
                              c2 = 0.75, c3 = 1, c02 = 100, c03 = 150,
                              K = Inf, N = Inf, S = -Inf,
                              steps1 = 1, stepm1 = 0.95, stepl1 = 0.85,
                              b1 = 1000, b2 = 2000, b3 = 3000,fixed = TRUE)[1], 
             utility_binary_L(n2 = 50, RRgo = 0.8, Adj = 0.1, w = 0.3,
                              p0 = 0.8, p11 =  0.2, p12 = 0.4, 
                              in1 = 300, in2 = 600, alpha = 0.025, beta = 0.1,
                              c2 = 0.75, c3 = 1, c02 = 100, c03 = 150,
                              K = Inf, N = Inf, S = -Inf,
                              steps1 = 1, stepm1 = 0.95, stepl1 = 0.85,
                              b1 = 1000, b2 = 2000, b3 = 3000,fixed = TRUE)[1])
})

test_that("Utility increases with smaller risk ratio", {
  expect_lte(utility_binary_L(n2 = 50, RRgo = 0.8, Adj = 0.1, w = 0.3,
                              p0 = 0.6, p11 =  0.3, p12 = 0.5, 
                              in1 = 300, in2 = 600, alpha = 0.025, beta = 0.1,
                              c2 = 0.75, c3 = 1, c02 = 100, c03 = 150,
                              K = Inf, N = Inf, S = -Inf,
                              steps1 = 1, stepm1 = 0.95, stepl1 = 0.85,
                              b1 = 1000, b2 = 2000, b3 = 3000,fixed = FALSE)[1], 
             utility_binary_L(n2 = 50, RRgo = 0.8, Adj = 0.1, w = 0.3,
                              p0 = 0.8, p11 =  0.2, p12 = 0.4, 
                              in1 = 300, in2 = 600, alpha = 0.025, beta = 0.1,
                              c2 = 0.75, c3 = 1, c02 = 100, c03 = 150,
                              K = Inf, N = Inf, S = -Inf,
                              steps1 = 1, stepm1 = 0.95, stepl1 = 0.85,
                              b1 = 1000, b2 = 2000, b3 = 3000,fixed = FALSE)[1])
})

test_that("Utility increases with smaller risk ratio", {
  expect_lte(utility_binary_R(n2 = 50, RRgo = 0.8, Adj = 0.8, w = 0.3,
                              p0 = 0.6, p11 =  0.3, p12 = 0.5, 
                              in1 = 300, in2 = 600, alpha = 0.025, beta = 0.1,
                              c2 = 0.75, c3 = 1, c02 = 100, c03 = 150,
                              K = Inf, N = Inf, S = -Inf,
                              steps1 = 1, stepm1 = 0.95, stepl1 = 0.85,
                              b1 = 1000, b2 = 2000, b3 = 3000,fixed = TRUE)[1], 
             utility_binary_R(n2 = 50, RRgo = 0.8, Adj = 0.8, w = 0.3,
                              p0 = 0.8, p11 =  0.2, p12 = 0.4, 
                              in1 = 300, in2 = 600, alpha = 0.025, beta = 0.1,
                              c2 = 0.75, c3 = 1, c02 = 100, c03 = 150,
                              K = Inf, N = Inf, S = -Inf,
                              steps1 = 1, stepm1 = 0.95, stepl1 = 0.85,
                              b1 = 1000, b2 = 2000, b3 = 3000,fixed = TRUE)[1])
})

test_that("Utility increases with smaller risk ratio", {
  expect_lte(utility_binary_R(n2 = 50, RRgo = 0.8, Adj = 0.8, w = 0.3,
                              p0 = 0.6, p11 =  0.3, p12 = 0.5, 
                              in1 = 300, in2 = 600, alpha = 0.025, beta = 0.1,
                              c2 = 0.75, c3 = 1, c02 = 100, c03 = 150,
                              K = Inf, N = Inf, S = -Inf,
                              steps1 = 1, stepm1 = 0.95, stepl1 = 0.85,
                              b1 = 1000, b2 = 2000, b3 = 3000,fixed = FALSE)[1], 
             utility_binary_R(n2 = 50, RRgo = 0.8, Adj = 0.8, w = 0.3,
                              p0 = 0.8, p11 =  0.2, p12 = 0.4, 
                              in1 = 300, in2 = 600, alpha = 0.025, beta = 0.1,
                              c2 = 0.75, c3 = 1, c02 = 100, c03 = 150,
                              K = Inf, N = Inf, S = -Inf,
                              steps1 = 1, stepm1 = 0.95, stepl1 = 0.85,
                              b1 = 1000, b2 = 2000, b3 = 3000,fixed = FALSE)[1])
})


test_that("Utility increases with smaller risk ratio", {
  expect_lte(utility_binary_L2(n2 = 50, RRgo = 0.8, Adj = 0.1, w = 0.3,
                              p0 = 0.6, p11 =  0.3, p12 = 0.5, 
                              in1 = 300, in2 = 600, alpha = 0.025, beta = 0.1,
                              c2 = 0.75, c3 = 1, c02 = 100, c03 = 150,
                              K = Inf, N = Inf, S = -Inf,
                              steps1 = 1, stepm1 = 0.95, stepl1 = 0.85,
                              b1 = 1000, b2 = 2000, b3 = 3000,fixed = TRUE)[1], 
             utility_binary_L2(n2 = 50, RRgo = 0.8, Adj = 0.1, w = 0.3,
                              p0 = 0.8, p11 =  0.2, p12 = 0.4, 
                              in1 = 300, in2 = 600, alpha = 0.025, beta = 0.1,
                              c2 = 0.75, c3 = 1, c02 = 100, c03 = 150,
                              K = Inf, N = Inf, S = -Inf,
                              steps1 = 1, stepm1 = 0.95, stepl1 = 0.85,
                              b1 = 1000, b2 = 2000, b3 = 3000,fixed = TRUE)[1])
})

test_that("Utility increases with smaller risk ratio", {
  expect_lte(utility_binary_L2(n2 = 50, RRgo = 0.8, Adj = 0.1, w = 0.3,
                              p0 = 0.6, p11 =  0.3, p12 = 0.5, 
                              in1 = 300, in2 = 600, alpha = 0.025, beta = 0.1,
                              c2 = 0.75, c3 = 1, c02 = 100, c03 = 150,
                              K = Inf, N = Inf, S = -Inf,
                              steps1 = 1, stepm1 = 0.95, stepl1 = 0.85,
                              b1 = 1000, b2 = 2000, b3 = 3000,fixed = FALSE)[1], 
             utility_binary_L2(n2 = 50, RRgo = 0.8, Adj = 0.1, w = 0.3,
                              p0 = 0.8, p11 =  0.2, p12 = 0.4, 
                              in1 = 300, in2 = 600, alpha = 0.025, beta = 0.1,
                              c2 = 0.75, c3 = 1, c02 = 100, c03 = 150,
                              K = Inf, N = Inf, S = -Inf,
                              steps1 = 1, stepm1 = 0.95, stepl1 = 0.85,
                              b1 = 1000, b2 = 2000, b3 = 3000,fixed = FALSE)[1])
})

test_that("Utility increases with smaller risk ratio", {
  expect_lte(utility_binary_R2(n2 = 50, RRgo = 0.8, Adj = 0.8, w = 0.3,
                              p0 = 0.6, p11 =  0.3, p12 = 0.5, 
                              in1 = 300, in2 = 600, alpha = 0.025, beta = 0.1,
                              c2 = 0.75, c3 = 1, c02 = 100, c03 = 150,
                              K = Inf, N = Inf, S = -Inf,
                              steps1 = 1, stepm1 = 0.95, stepl1 = 0.85,
                              b1 = 1000, b2 = 2000, b3 = 3000,fixed = TRUE)[1], 
             utility_binary_R2(n2 = 50, RRgo = 0.8, Adj = 0.8, w = 0.3,
                              p0 = 0.8, p11 =  0.2, p12 = 0.4, 
                              in1 = 300, in2 = 600, alpha = 0.025, beta = 0.1,
                              c2 = 0.75, c3 = 1, c02 = 100, c03 = 150,
                              K = Inf, N = Inf, S = -Inf,
                              steps1 = 1, stepm1 = 0.95, stepl1 = 0.85,
                              b1 = 1000, b2 = 2000, b3 = 3000,fixed = TRUE)[1])
})

test_that("Utility increases with smaller risk ratio", {
  expect_lte(utility_binary_R2(n2 = 50, RRgo = 0.8, Adj = 0.8, w = 0.3,
                              p0 = 0.6, p11 =  0.3, p12 = 0.5, 
                              in1 = 300, in2 = 600, alpha = 0.025, beta = 0.1,
                              c2 = 0.75, c3 = 1, c02 = 100, c03 = 150,
                              K = Inf, N = Inf, S = -Inf,
                              steps1 = 1, stepm1 = 0.95, stepl1 = 0.85,
                              b1 = 1000, b2 = 2000, b3 = 3000,fixed = FALSE)[1], 
             utility_binary_R2(n2 = 50, RRgo = 0.8, Adj = 0.8, w = 0.3,
                              p0 = 0.8, p11 =  0.2, p12 = 0.4, 
                              in1 = 300, in2 = 600, alpha = 0.025, beta = 0.1,
                              c2 = 0.75, c3 = 1, c02 = 100, c03 = 150,
                              K = Inf, N = Inf, S = -Inf,
                              steps1 = 1, stepm1 = 0.95, stepl1 = 0.85,
                              b1 = 1000, b2 = 2000, b3 = 3000,fixed = FALSE)[1])
})
