## Tests for two trials in phase III ##

test_that("Utility increases with smaller hazard ratio", {
  skip_on_cran()
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
  skip_on_cran()
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

test_that("Probability increases for lower hazard ratio", {
  skip_on_cran()
  expect_lte(EPsProg2(HRgo = 0.8, d2 = 50,  alpha = 0.025, beta = 0.1, 
                      w = 0.3, hr1 =  0.69, hr2 = 0.81,id1 = 210, id2 = 420,
                      case = 2, size = "all", fixed = FALSE), 
             EPsProg2(HRgo = 0.8, d2 = 50,  alpha = 0.025, beta = 0.1, 
                      w = 0.3, hr1 =  0.6, hr2 = 0.7,id1 = 210, id2 = 420,
                      case = 2, size = "all", fixed = FALSE))
})

test_that("EPsProg2 works", {
  skip_on_cran()
  expect_equal(EPsProg2(HRgo = 0.8, d2 = 50,  alpha = 0.025, beta = 0.1, 
                        w = 0.3, hr1 =  0.69, hr2 = 0.81,id1 = 210, id2 = 420,
                        case = 2, size = "small", fixed = FALSE), 0.0360759794)
})

test_that("Probability for large is smaller than all", {
  skip_on_cran()
  expect_lte(EPsProg2(HRgo = 0.8, d2 = 50,  alpha = 0.025, beta = 0.1, 
                      w = 0.3, hr1 =  0.69, hr2 = 0.81,id1 = 210, id2 = 420,
                      case = 2, size = "large", fixed = FALSE), 
             EPsProg2(HRgo = 0.8, d2 = 50,  alpha = 0.025, beta = 0.1, 
                      w = 0.3, hr1 =  0.69, hr2 = 0.81,id1 = 210, id2 = 420,
                      case = 2, size = "all", fixed = FALSE))
})


## Tests for three trials in phase III ##

test_that("Utility increases with smaller hazard ratio", {
  skip_on_cran()
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
  skip_on_cran()
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

test_that("Probability for large is smaller than all", {
  skip_on_cran()
  expect_lte(EPsProg3(HRgo = 0.8, d2 = 50,  alpha = 0.025, beta = 0.1, 
                      w = 0.3, hr1 =  0.69, hr2 = 0.81,id1 = 210, id2 = 420,
                      case = 2, size = "large", fixed = TRUE), 
             EPsProg3(HRgo = 0.8, d2 = 50,  alpha = 0.025, beta = 0.1, 
                      w = 0.3, hr1 =  0.69, hr2 = 0.81,id1 = 210, id2 = 420,
                      case = 2, size = "all", fixed = TRUE))
})

test_that("Probability for small is smaller than all", {
  skip_on_cran()
  expect_lte(EPsProg3(HRgo = 0.8, d2 = 50,  alpha = 0.025, beta = 0.1, 
                      w = 0.3, hr1 =  0.69, hr2 = 0.81,id1 = 210, id2 = 420,
                      case = 2, size = "small", fixed = TRUE), 
             EPsProg3(HRgo = 0.8, d2 = 50,  alpha = 0.025, beta = 0.1, 
                      w = 0.3, hr1 =  0.69, hr2 = 0.81,id1 = 210, id2 = 420,
                      case = 2, size = "all", fixed = TRUE))
})

test_that("EPsProg3 works for case 3", {
  skip_on_cran()
  expect_equal(EPsProg3(HRgo = 0.8, d2 = 50,  alpha = 0.025, beta = 0.1, 
                        w = 0.3, hr1 =  0.69, hr2 = 0.81, id1 = 210, id2 = 420,
                        case = 3, size = "small", fixed = FALSE),0.0077498414)
})

test_that("EPsProg3 works for case 3", {
  skip_on_cran()
  expect_equal(EPsProg3(HRgo = 0.8, d2 = 50,  alpha = 0.025, beta = 0.1, 
                        w = 0.3, hr1 =  0.69, hr2 = 0.81, id1 = 210, id2 = 420,
                        case = 3, size = "large", fixed = FALSE),0.083166677)
})

test_that("EPsProg3 works for case 3", {
  skip_on_cran()
  expect_equal(EPsProg3(HRgo = 0.8, d2 = 50,  alpha = 0.025, beta = 0.1, 
                        w = 0.3, hr1 =  0.69, hr2 = 0.81, id1 = 210, id2 = 420,
                        case = 3, size = "all", fixed = FALSE),0.159433104)
})

## Tests for four trials in phase III ##

test_that("Utility increases with smaller hazard ratio", {
  skip_on_cran()
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
  skip_on_cran()
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


## Tests for setting 23 ##

test_that("Epgo23 works", {
  skip_on_cran()
  expect_equal(Epgo23(HRgo = 0.8, d2 = 50,  w = 0.3, alpha = 0.025, beta = 0.1,
                       hr1 =  0.69, hr2 = 0.81, id1 = 280, id2 = 420), 0.165258433)
})

test_that("EPsProg23 works for case 3", {
  skip_on_cran()
  expect_equal(EPsProg23(HRgo = 0.8, d2 = 50,  alpha = 0.025, beta = 0.1, w = 0.3,
                         hr1 =  0.69, hr2 = 0.81, id1 = 280, id2 = 420, 
                         case = 2, size = "large", ymin = 0.5), 0.057610215)
})

test_that("Probability increases for lower hazard ratio", {
  skip_on_cran()
  expect_lte(EPsProg23(HRgo = 0.8, d2 = 50,  alpha = 0.025, beta = 0.1, w = 0.3,
                         hr1 =  0.69, hr2 = 0.81, id1 = 280, id2 = 420, 
                         case = 2, size = "all", ymin = 0.5), 
             EPsProg23(HRgo = 0.8, d2 = 50,  alpha = 0.025, beta = 0.1, w = 0.3,
                       hr1 =  0.6, hr2 = 0.8, id1 = 280, id2 = 420, 
                       case = 2, size = "all", ymin = 0.5))
})

test_that("utility23 works", {
  skip_on_cran()
  expect_equal(utility23(d2 = 50, HRgo = 0.8,  w = 0.3,
                       hr1 =  0.69, hr2 = 0.81, id1 = 280, id2 = 420, 
                       alpha = 0.025, beta = 0.1, xi2 = 0.7, xi3 = 0.7,
                       c2 = 0.75, c3 = 1, c02 = 100, c03 = 150,
                       b1 = 1000, b2 = 2000, b3 = 3000)[1], -247.53105)
})

