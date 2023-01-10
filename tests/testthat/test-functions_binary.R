test_that("Utility increases with lower risk ratio", {
  expect_lte(utility_binary(n2 = 100, RRgo = 0.8, w = 0.3, 
                            p0 = 0.5, p11 = 0.4, p12 = 0.6, in1 = 30, in2 = 60,
                            alpha = 0.05, beta = 0.1, 
                            c2 = 0.75, c3 = 1, c02 = 100, c03 = 150, 
                            K = Inf, N = Inf, S = -Inf,
                            steps1 = 1, stepm1 = 0.95, stepl1 = 0.85,
                            b1 = 1000, b2 = 2000, b3 = 3000,
                            gamma = 0, fixed = TRUE)[1], 
             utility_binary(n2 = 100, RRgo = 0.8, w = 0.3, 
                            p0 = 0.8, p11 = 0.2, p12 = 0.6, in1 = 30, in2 = 60,
                            alpha = 0.05, beta = 0.1, 
                            c2 = 0.75, c3 = 1, c02 = 100, c03 = 150, 
                            K = Inf, N = Inf, S = -Inf,
                            steps1 = 1, stepm1 = 0.95, stepl1 = 0.85,
                            b1 = 1000, b2 = 2000, b3 = 3000,
                            gamma = 0, fixed = TRUE)[1])
})

test_that("Utility is -9999 if sample size constraint is not satisfied", {
  expect_lte(utility_binary(n2 = 100, RRgo = 0.8, w = 0.3, 
                            p0 = 0.5, p11 = 0.4, p12 = 0.6, in1 = 30, in2 = 60,
                            alpha = 0.05, beta = 0.1, 
                            c2 = 0.75, c3 = 1, c02 = 100, c03 = 150, 
                            K = Inf, N = 100, S = -Inf,
                            steps1 = 1, stepm1 = 0.95, stepl1 = 0.85,
                            b1 = 1000, b2 = 2000, b3 = 3000,
                            gamma = 0, fixed = TRUE)[1], -9999)
})


test_that("Utility is -9999 if minimul success probability constraint is not satisfied", {
  expect_lte(utility_binary(n2 = 100, RRgo = 0.8, w = 0.3, 
                            p0 = 0.5, p11 = 0.4, p12 = 0.6, in1 = 30, in2 = 60,
                            alpha = 0.05, beta = 0.1, 
                            c2 = 0.75, c3 = 1, c02 = 100, c03 = 150, 
                            K = Inf, N = Inf, S = 0.9,
                            steps1 = 1, stepm1 = 0.95, stepl1 = 0.85,
                            b1 = 1000, b2 = 2000, b3 = 3000,
                            gamma = 0, fixed = TRUE)[1], -9999)
})
