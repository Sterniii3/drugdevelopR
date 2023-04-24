## tests for two phase III trials ##

test_that("Higher treatment effect leads to higher utility", {
  skip_on_cran()
  expect_lte(utility2_normal(kappa = 0.1, n2 = 50,  alpha = 0.025, beta = 0.1, w = 0.3,
                               Delta1 = 0.375, Delta2 = 0.625, in1 = 300, in2 = 600, 
                               a = 0.25, b = 0.75, 
                               c2 = 0.675, c3 = 0.72, c02 = 15, c03 = 20,
                               K = Inf, N = Inf, S = -Inf, 
                               b1 = 3000, b2 = 8000, b3 = 10000, 
                               case = 2, fixed = TRUE)[1], 
               utility2_normal(kappa = 0.1, n2 = 50,  alpha = 0.025, beta = 0.1, w = 0.3,
                               Delta1 = 0.5, Delta2 = 0.8, in1 = 300, in2 = 600, 
                               a = 0.25, b = 0.75, 
                               c2 = 0.675, c3 = 0.72, c02 = 15, c03 = 20,
                               K = Inf, N = Inf, S = -Inf, 
                               b1 = 3000, b2 = 8000, b3 = 10000, 
                               case = 2, fixed = TRUE)[1])
})

test_that("Higher treatment effect leads to higher utility", {
  skip_on_cran()
  expect_lte(utility2_normal(kappa = 0.1, n2 = 50,  alpha = 0.025, beta = 0.1, w = 0.3,
                             Delta1 = 0.375, Delta2 = 0.625, in1 = 300, in2 = 600, 
                             a = 0.25, b = 0.75, 
                             c2 = 0.675, c3 = 0.72, c02 = 15, c03 = 20,
                             K = Inf, N = Inf, S = -Inf, 
                             b1 = 3000, b2 = 8000, b3 = 10000, 
                             case = 1, fixed = FALSE)[1], 
             utility2_normal(kappa = 0.1, n2 = 50,  alpha = 0.025, beta = 0.1, w = 0.3,
                             Delta1 = 0.5, Delta2 = 0.8, in1 = 300, in2 = 600, 
                             a = 0.25, b = 0.75, 
                             c2 = 0.675, c3 = 0.72, c02 = 15, c03 = 20,
                             K = Inf, N = Inf, S = -Inf, 
                             b1 = 3000, b2 = 8000, b3 = 10000, 
                             case = 1, fixed = FALSE)[1])
})

test_that("Function returns -9999 if constraint can not be satisfied", {
  skip_on_cran()
  expect_equal(utility2_normal(kappa = 0.1, n2 = 50,  alpha = 0.025, beta = 0.1, w = 0.3,
                             Delta1 = 0.375, Delta2 = 0.625, in1 = 300, in2 = 600, 
                             a = 0.25, b = 0.75, 
                             c2 = 0.675, c3 = 0.72, c02 = 15, c03 = 20,
                             K = Inf, N = Inf, S = 0.8, 
                             b1 = 3000, b2 = 8000, b3 = 10000, 
                             case = 2, fixed = FALSE)[1], -9999)
})

## tests for three phase III trials ##

test_that("Higher treatment effect leads to higher utility", {
  skip_on_cran()
  expect_lte(utility3_normal(kappa = 0.1, n2 = 50,  alpha = 0.025, beta = 0.1, w = 0.3,
                             Delta1 = 0.375, Delta2 = 0.625, in1 = 300, in2 = 600, 
                             a = 0.25, b = 0.75, 
                             c2 = 0.675, c3 = 0.72, c02 = 15, c03 = 20,
                             K = Inf, N = Inf, S = -Inf, 
                             b1 = 3000, b2 = 8000, b3 = 10000, 
                             case = 3, fixed = TRUE)[1], 
             utility3_normal(kappa = 0.1, n2 = 50,  alpha = 0.025, beta = 0.1, w = 0.3,
                             Delta1 = 0.5, Delta2 = 0.8, in1 = 300, in2 = 600, 
                             a = 0.25, b = 0.75, 
                             c2 = 0.675, c3 = 0.72, c02 = 15, c03 = 20,
                             K = Inf, N = Inf, S = -Inf, 
                             b1 = 3000, b2 = 8000, b3 = 10000, 
                             case = 3, fixed = TRUE)[1])
})

test_that("Higher treatment effect leads to higher utility", {
  skip_on_cran()
  expect_lte(utility3_normal(kappa = 0.1, n2 = 50,  alpha = 0.025, beta = 0.1, w = 0.3,
                             Delta1 = 0.375, Delta2 = 0.625, in1 = 300, in2 = 600, 
                             a = 0.25, b = 0.75, 
                             c2 = 0.675, c3 = 0.72, c02 = 15, c03 = 20,
                             K = Inf, N = Inf, S = -Inf, 
                             b1 = 3000, b2 = 8000, b3 = 10000, 
                             case = 2, fixed = FALSE)[1], 
             utility3_normal(kappa = 0.1, n2 = 50,  alpha = 0.025, beta = 0.1, w = 0.3,
                             Delta1 = 0.5, Delta2 = 0.8, in1 = 300, in2 = 600, 
                             a = 0.25, b = 0.75, 
                             c2 = 0.675, c3 = 0.72, c02 = 15, c03 = 20,
                             K = Inf, N = Inf, S = -Inf, 
                             b1 = 3000, b2 = 8000, b3 = 10000, 
                             case = 2, fixed = FALSE)[1])
})

test_that("Function returns -9999 if constraint can not be satisfied", {
  skip_on_cran()
  expect_equal(utility3_normal(kappa = 0.1, n2 = 50,  alpha = 0.025, beta = 0.1, w = 0.3,
                               Delta1 = 0.375, Delta2 = 0.625, in1 = 300, in2 = 600, 
                               a = 0.25, b = 0.75, 
                               c2 = 0.675, c3 = 0.72, c02 = 15, c03 = 20,
                               K = Inf, N = Inf, S = 0.8, 
                               b1 = 3000, b2 = 8000, b3 = 10000, 
                               case = 2, fixed = TRUE)[1], -9999)
})


test_that("Probability for small is smaller than all", {
  skip_on_cran()
  expect_lte(EPsProg3_normal(kappa = 0.1, n2 = 50, alpha = 0.025, beta = 0.1,
                             w = 0.3, Delta1 = 0.375, Delta2 = 0.625,
                             in1 = 300, in2 = 600, a = 0.25, b = 0.75, 
                             case = 3, size = "small", fixed = FALSE), 
             EPsProg3_normal(kappa = 0.1, n2 = 50, alpha = 0.025, beta = 0.1,
                             w = 0.3, Delta1 = 0.375, Delta2 = 0.625,
                             in1 = 300, in2 = 600, a = 0.25, b = 0.75, 
                             case = 3, size = "all", fixed = FALSE))
})

test_that("Probability for large is smaller than all", {
  skip_on_cran()
  expect_lte(EPsProg3_normal(kappa = 0.1, n2 = 50, alpha = 0.025, beta = 0.1,
                             w = 0.3, Delta1 = 0.375, Delta2 = 0.625,
                             in1 = 300, in2 = 600, a = 0.25, b = 0.75, 
                             case = 3, size = "large", fixed = FALSE), 
             EPsProg3_normal(kappa = 0.1, n2 = 50, alpha = 0.025, beta = 0.1,
                             w = 0.3, Delta1 = 0.375, Delta2 = 0.625,
                             in1 = 300, in2 = 600, a = 0.25, b = 0.75, 
                             case = 3, size = "all", fixed = FALSE))
})





## tests for four phase III trials ##


test_that("Higher treatment effect leads to higher utility", {
  expect_lte(utility4_normal(kappa = 0.1, n2 = 50,  alpha = 0.025, beta = 0.1, w = 0.3,
                             Delta1 = 0.375, Delta2 = 0.625, in1 = 300, in2 = 600, 
                             a = 0.25, b = 0.75, 
                             c2 = 0.675, c3 = 0.72, c02 = 15, c03 = 20,
                             K = Inf, N = Inf, S = -Inf, 
                             b1 = 3000, b2 = 8000, b3 = 10000, 
                             case = 4, fixed = TRUE)[1], 
             utility4_normal(kappa = 0.1, n2 = 50,  alpha = 0.025, beta = 0.1, w = 0.3,
                             Delta1 = 0.5, Delta2 = 0.8, in1 = 300, in2 = 600, 
                             a = 0.25, b = 0.75, 
                             c2 = 0.675, c3 = 0.72, c02 = 15, c03 = 20,
                             K = Inf, N = Inf, S = -Inf, 
                             b1 = 3000, b2 = 8000, b3 = 10000, 
                             case = 4, fixed = TRUE)[1])
})

test_that("Higher treatment effect leads to higher utility", {
  skip_on_cran()
  expect_lte(utility4_normal(kappa = 0.1, n2 = 50,  alpha = 0.025, beta = 0.1, w = 0.3,
                             Delta1 = 0.375, Delta2 = 0.625, in1 = 300, in2 = 600, 
                             a = 0.25, b = 0.75, 
                             c2 = 0.675, c3 = 0.72, c02 = 15, c03 = 20,
                             K = Inf, N = Inf, S = -Inf, 
                             b1 = 3000, b2 = 8000, b3 = 10000, 
                             case = 3, fixed = FALSE)[1], 
             utility4_normal(kappa = 0.1, n2 = 50,  alpha = 0.025, beta = 0.1, w = 0.3,
                             Delta1 = 0.5, Delta2 = 0.8, in1 = 300, in2 = 600, 
                             a = 0.25, b = 0.75, 
                             c2 = 0.675, c3 = 0.72, c02 = 15, c03 = 20,
                             K = Inf, N = Inf, S = -Inf, 
                             b1 = 3000, b2 = 8000, b3 = 10000, 
                             case = 3, fixed = FALSE)[1])
})

## tests for case 23 ##

test_that("utility23_normal works",{
  skip_on_cran()
  expect_equal(utility23_normal(n2 = 50, kappa = 0.2, w = 0.3,
                                alpha = 0.025, beta = 0.1, Delta1 = 0.375, Delta2 = 0.625,
                                in1 = 300, in2 = 600, a = 0.25, b = 0.75,
                                c2 = 0.675, c3 = 0.72, c02 = 15, c03 = 20,
                                b1 = 3000, b2 = 8000, b3 = 10000)[1], 1844.2313)
})
