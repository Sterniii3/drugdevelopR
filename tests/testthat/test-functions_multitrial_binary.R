## Tests for two trials in phase III ##

test_that("Utility increases for lower risk ratio", {
  skip_on_cran()
  expect_lte(utility2_binary(n2 = 50, RRgo = 0.8, w = 0.3,
                             p0 = 0.6, p11 =  0.3, p12 = 0.5,
                             in1 = 300, in2 = 600, alpha = 0.025, beta = 0.1,
                             c2 = 0.75, c3 = 1, c02 = 100, c03 = 150,
                             K = Inf, N = Inf, S = -Inf,
                             b1 = 1000, b2 = 2000, b3 = 3000,
                             case = 2, fixed = TRUE)[1], 
             utility2_binary(n2 = 50, RRgo = 0.8, w = 0.3,
                             p0 = 0.8, p11 =  0.2, p12 = 0.5,
                             in1 = 300, in2 = 600, alpha = 0.025, beta = 0.1,
                             c2 = 0.75, c3 = 1, c02 = 100, c03 = 150,
                             K = Inf, N = Inf, S = -Inf,
                             b1 = 1000, b2 = 2000, b3 = 3000,
                             case = 2, fixed = TRUE)[1])
})

test_that("Utility increases for lower risk ratio", {
  skip_on_cran()
  expect_lte(utility2_binary(n2 = 50, RRgo = 0.8, w = 0.3,
                             p0 = 0.6, p11 =  0.3, p12 = 0.5,
                             in1 = 300, in2 = 600, alpha = 0.025, beta = 0.1,
                             c2 = 0.75, c3 = 1, c02 = 100, c03 = 150,
                             K = Inf, N = Inf, S = -Inf,
                             b1 = 1000, b2 = 2000, b3 = 3000,
                             case = 1, fixed = TRUE)[1], 
             utility2_binary(n2 = 50, RRgo = 0.8, w = 0.3,
                             p0 = 0.8, p11 =  0.2, p12 = 0.5,
                             in1 = 300, in2 = 600, alpha = 0.025, beta = 0.1,
                             c2 = 0.75, c3 = 1, c02 = 100, c03 = 150,
                             K = Inf, N = Inf, S = -Inf,
                             b1 = 1000, b2 = 2000, b3 = 3000,
                             case = 1, fixed = TRUE)[1])
})

test_that("Utility2_binary works for case 1", {
  skip_on_cran()
  expect_lte(utility2_binary(n2 = 50, RRgo = 0.8, w = 0.3,
                             p0 = 0.6, p11 =  0.3, p12 = 0.5,
                             in1 = 300, in2 = 600, alpha = 0.025, beta = 0.1,
                             c2 = 0.75, c3 = 1, c02 = 100, c03 = 150,
                             K = Inf, N = Inf, S = -Inf,
                             b1 = 1000, b2 = 2000, b3 = 3000,
                             case = 1, fixed = FALSE)[1],1196.528)
})

## Tests for three trials in phase III ##

test_that("Utility increases for lower risk ratio", {
  skip_on_cran()
  expect_lte(utility3_binary(n2 = 50, RRgo = 0.8, w = 0.3,
                             p0 = 0.6, p11 =  0.3, p12 = 0.5,
                             in1 = 300, in2 = 600, alpha = 0.025, beta = 0.1,
                             c2 = 0.75, c3 = 1, c02 = 100, c03 = 150,
                             K = Inf, N = Inf, S = -Inf,
                             b1 = 1000, b2 = 2000, b3 = 3000,
                             case = 3, fixed = TRUE)[1], 
             utility3_binary(n2 = 50, RRgo = 0.8, w = 0.3,
                             p0 = 0.8, p11 =  0.2, p12 = 0.5,
                             in1 = 300, in2 = 600, alpha = 0.025, beta = 0.1,
                             c2 = 0.75, c3 = 1, c02 = 100, c03 = 150,
                             K = Inf, N = Inf, S = -Inf,
                             b1 = 1000, b2 = 2000, b3 = 3000,
                             case = 3, fixed = TRUE)[1])
})

test_that("Utility increases for lower risk ratio", {
  skip_on_cran()
  expect_lte(utility3_binary(n2 = 100, RRgo = 0.8, w = 0.3,
                             p0 = 0.6, p11 =  0.3, p12 = 0.5,
                             in1 = 300, in2 = 600, alpha = 0.025, beta = 0.1,
                             c2 = 0.75, c3 = 1, c02 = 100, c03 = 150,
                             K = Inf, N = Inf, S = -Inf,
                             b1 = 1000, b2 = 2000, b3 = 3000,
                             case = 2, fixed = TRUE)[1], 
             utility3_binary(n2 = 100, RRgo = 0.8, w = 0.3,
                             p0 = 0.8, p11 =  0.2, p12 = 0.5,
                             in1 = 300, in2 = 600, alpha = 0.025, beta = 0.1,
                             c2 = 0.75, c3 = 1, c02 = 100, c03 = 150,
                             K = Inf, N = Inf, S = -Inf,
                             b1 = 1000, b2 = 2000, b3 = 3000,
                             case = 2, fixed = TRUE)[1])
})

test_that("EPsProg3_binary works for case 3", {
  skip_on_cran()
  expect_equal(EPsProg3_binary(RRgo = 0.8, n2 = 50,  alpha = 0.025, beta = 0.1, 
                               p0 = 0.6,  w = 0.3, p11 =  0.3, p12 = 0.5, 
                               in1 = 300, in2 = 600, case = 3, size = "small",
                               fixed = FALSE),0.004263322)
})

test_that("EPsProg3_binary works for case 3", {
  skip_on_cran()
  expect_equal(EPsProg3_binary(RRgo = 0.8, n2 = 50,  alpha = 0.025, beta = 0.1, 
                               p0 = 0.6,  w = 0.3, p11 =  0.3, p12 = 0.5, 
                               in1 = 300, in2 = 600, case = 3, size = "large",
                               fixed = FALSE),0.29460057)
})

test_that("EPsProg3_binary works for case 3", {
  skip_on_cran()
  expect_equal(EPsProg3_binary(RRgo = 0.8, n2 = 50,  alpha = 0.025, beta = 0.1, 
                               p0 = 0.6,  w = 0.3, p11 =  0.3, p12 = 0.5, 
                               in1 = 300, in2 = 600, case = 3, size = "all",
                               fixed = FALSE),0.38905358)
})

## Tests for four trials in phase III ##

test_that("Utility increases for lower risk ratio", {
  skip_on_cran()
  expect_lte(utility4_binary(n2 = 50, RRgo = 0.8, w = 0.3,
                             p0 = 0.6, p11 =  0.3, p12 = 0.5,
                             in1 = 300, in2 = 600, alpha = 0.025, beta = 0.1,
                             c2 = 0.75, c3 = 1, c02 = 100, c03 = 150,
                             K = Inf, N = Inf, S = -Inf,
                             b1 = 1000, b2 = 2000, b3 = 3000,
                             case = 3, fixed = TRUE)[1], 
             utility4_binary(n2 = 50, RRgo = 0.8, w = 0.3,
                             p0 = 0.8, p11 =  0.2, p12 = 0.5,
                             in1 = 300, in2 = 600, alpha = 0.025, beta = 0.1,
                             c2 = 0.75, c3 = 1, c02 = 100, c03 = 150,
                             K = Inf, N = Inf, S = -Inf,
                             b1 = 1000, b2 = 2000, b3 = 3000,
                             case = 3, fixed = TRUE)[1])
})

test_that("Utility is -9999 if constraint is not met", {
  skip_on_cran()
  expect_lte(utility4_binary(n2 = 50, RRgo = 0.8, w = 0.3,
                             p0 = 0.6, p11 =  0.3, p12 = 0.5,
                             in1 = 300, in2 = 600, alpha = 0.025, beta = 0.1,
                             c2 = 0.75, c3 = 1, c02 = 100, c03 = 150,
                             K = 200, N = Inf, S = -Inf,
                             b1 = 1000, b2 = 2000, b3 = 3000,
                             case = 4, fixed = FALSE)[1], -9999)
})

test_that("Probability for large is smaller than all", {
  skip_on_cran()
  expect_lte(EPsProg4_binary(RRgo = 0.8, n2 = 50,  alpha = 0.025, beta = 0.1, 
                             p0 = 0.6,  w = 0.3, p11 =  0.3, p12 = 0.5, 
                             in1 = 300, in2 = 600, case = 3, size = "large",
                             fixed = FALSE), 
             EPsProg4_binary(RRgo = 0.8, n2 = 50,  alpha = 0.025, beta = 0.1, 
                             p0 = 0.6,  w = 0.3, p11 =  0.3, p12 = 0.5, 
                             in1 = 300, in2 = 600, case = 3, size = "all",
                             fixed = FALSE))
})

test_that("Probability for small is smaller than all", {
  skip_on_cran()
  expect_lte(EPsProg4_binary(RRgo = 0.8, n2 = 50,  alpha = 0.025, beta = 0.1, 
                             p0 = 0.6,  w = 0.3, p11 =  0.3, p12 = 0.5, 
                             in1 = 300, in2 = 600, case = 3, size = "small",
                             fixed = FALSE), 
             EPsProg4_binary(RRgo = 0.8, n2 = 50,  alpha = 0.025, beta = 0.1, 
                             p0 = 0.6,  w = 0.3, p11 =  0.3, p12 = 0.5, 
                             in1 = 300, in2 = 600, case = 3, size = "all",
                             fixed = FALSE))
})

## Tests for case 23 ##

test_that("Utility function works for case 23", {
  skip_on_cran()
  expect_equal(utility23_binary(n2 = 50, RRgo = 0.8,  w = 0.3,
                              alpha = 0.05, beta = 0.1,
                              p0 = 0.6, p11 =  0.3, p12 = 0.5, 
                              in1 = 300, in2 = 600,
                              c2 = 0.75, c3 = 1, c02 = 100, c03 = 150,
                              b1 = 1000, b2 = 2000, b3 = 3000)[1], 740.60422)
})

