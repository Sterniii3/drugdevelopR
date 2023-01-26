test_that("Utility increases with smaller hazard ratio", {
  expect_lte(utility2(d2 = 50, HRgo = 0.8,  w = 0.3, 
                        hr1 =  0.69, hr2 = 0.81, id1 = 210, id2 = 420, 
                        alpha = 0.025, beta = 0.1, xi2 = 0.7, xi3 = 0.7,
                        c2 = 0.75, c3 = 1, c02 = 100, c03 = 150,
                        K = Inf, N = Inf, S = -Inf,
                        b1 = 1000, b2 = 2000, b3 = 3000, case = 2, fixed = TRUE)[1], 
             utility2(d2 = 50, HRgo = 0.8,  w = 0.3, 
                      hr1 =  0.5, hr2 = 0.7, id1 = 210, id2 = 420, 
                      alpha = 0.025, beta = 0.1, xi2 = 0.7, xi3 = 0.7,
                      c2 = 0.75, c3 = 1, c02 = 100, c03 = 150,
                      K = Inf, N = Inf, S = -Inf,
                      b1 = 1000, b2 = 2000, b3 = 3000, case = 2, fixed = TRUE)[1])
})

test_that("Utility increases with smaller hazard ratio", {
  expect_lte(utility2(d2 = 50, HRgo = 0.8,  w = 0.3, 
                      hr1 =  0.69, hr2 = 0.81, id1 = 210, id2 = 420, 
                      alpha = 0.025, beta = 0.1, xi2 = 0.7, xi3 = 0.7,
                      c2 = 0.75, c3 = 1, c02 = 100, c03 = 150,
                      K = Inf, N = Inf, S = -Inf,
                      b1 = 1000, b2 = 2000, b3 = 3000, case = 1, fixed = FALSE)[1], 
             utility2(d2 = 50, HRgo = 0.8,  w = 0.3, 
                      hr1 =  0.5, hr2 = 0.7, id1 = 210, id2 = 420, 
                      alpha = 0.025, beta = 0.1, xi2 = 0.7, xi3 = 0.7,
                      c2 = 0.75, c3 = 1, c02 = 100, c03 = 150,
                      K = Inf, N = Inf, S = -Inf,
                      b1 = 1000, b2 = 2000, b3 = 3000, case = 1, fixed = FALSE)[1])
})

test_that("Utility increases with smaller hazard ratio", {
  expect_lte(utility3(d2 = 50, HRgo = 0.8,  w = 0.3, 
                      hr1 =  0.69, hr2 = 0.81, id1 = 210, id2 = 420, 
                      alpha = 0.025, beta = 0.1, xi2 = 0.7, xi3 = 0.7,
                      c2 = 0.75, c3 = 1, c02 = 100, c03 = 150,
                      K = Inf, N = Inf, S = -Inf,
                      b1 = 1000, b2 = 2000, b3 = 3000, case = 3, fixed = TRUE)[1], 
             utility3(d2 = 50, HRgo = 0.8,  w = 0.3, 
                      hr1 =  0.5, hr2 = 0.7, id1 = 210, id2 = 420, 
                      alpha = 0.025, beta = 0.1, xi2 = 0.7, xi3 = 0.7,
                      c2 = 0.75, c3 = 1, c02 = 100, c03 = 150,
                      K = Inf, N = Inf, S = -Inf,
                      b1 = 1000, b2 = 2000, b3 = 3000, case = 3, fixed = TRUE)[1])
})

test_that("Utility increases with smaller hazard ratio", {
  expect_lte(utility3(d2 = 50, HRgo = 0.8,  w = 0.3, 
                      hr1 =  0.69, hr2 = 0.81, id1 = 210, id2 = 420, 
                      alpha = 0.025, beta = 0.1, xi2 = 0.7, xi3 = 0.7,
                      c2 = 0.75, c3 = 1, c02 = 100, c03 = 150,
                      K = Inf, N = Inf, S = -Inf,
                      b1 = 1000, b2 = 2000, b3 = 3000, case = 2, fixed = FALSE)[1], 
             utility3(d2 = 50, HRgo = 0.8,  w = 0.3, 
                      hr1 =  0.5, hr2 = 0.7, id1 = 210, id2 = 420, 
                      alpha = 0.025, beta = 0.1, xi2 = 0.7, xi3 = 0.7,
                      c2 = 0.75, c3 = 1, c02 = 100, c03 = 150,
                      K = Inf, N = Inf, S = -Inf,
                      b1 = 1000, b2 = 2000, b3 = 3000, case = 2, fixed = FALSE)[1])
})

test_that("Utility increases with smaller hazard ratio", {
  expect_lte(utility4(d2 = 50, HRgo = 0.8,  w = 0.3, 
                      hr1 =  0.69, hr2 = 0.81, id1 = 210, id2 = 420, 
                      alpha = 0.025, beta = 0.1, xi2 = 0.7, xi3 = 0.7,
                      c2 = 0.75, c3 = 1, c02 = 100, c03 = 150,
                      K = Inf, N = Inf, S = -Inf,
                      b1 = 1000, b2 = 2000, b3 = 3000, case = 4, fixed = TRUE)[1], 
             utility4(d2 = 50, HRgo = 0.8,  w = 0.3, 
                      hr1 =  0.5, hr2 = 0.7, id1 = 210, id2 = 420, 
                      alpha = 0.025, beta = 0.1, xi2 = 0.7, xi3 = 0.7,
                      c2 = 0.75, c3 = 1, c02 = 100, c03 = 150,
                      K = Inf, N = Inf, S = -Inf,
                      b1 = 1000, b2 = 2000, b3 = 3000, case = 4, fixed = TRUE)[1])
})

test_that("Utility increases with smaller hazard ratio", {
  expect_lte(utility4(d2 = 50, HRgo = 0.8,  w = 0.3, 
                      hr1 =  0.69, hr2 = 0.81, id1 = 210, id2 = 420, 
                      alpha = 0.025, beta = 0.1, xi2 = 0.7, xi3 = 0.7,
                      c2 = 0.75, c3 = 1, c02 = 100, c03 = 150,
                      K = Inf, N = Inf, S = -Inf,
                      b1 = 1000, b2 = 2000, b3 = 3000, case = 3, fixed = FALSE)[1], 
             utility4(d2 = 50, HRgo = 0.8,  w = 0.3, 
                      hr1 =  0.5, hr2 = 0.7, id1 = 210, id2 = 420, 
                      alpha = 0.025, beta = 0.1, xi2 = 0.7, xi3 = 0.7,
                      c2 = 0.75, c3 = 1, c02 = 100, c03 = 150,
                      K = Inf, N = Inf, S = -Inf,
                      b1 = 1000, b2 = 2000, b3 = 3000, case = 3, fixed = FALSE)[1])
})

