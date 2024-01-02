

test_that("Utility increases with lower hazard ratio", {
  skip_on_cran()
  set.seed(61216)
  expect_lte(utility_multiple_tte(n2 = 50, HRgo = 0.8, alpha = 0.025, beta = 0.1,
                                  hr1 = 0.75, hr2 = 0.80, id1 = 300, id2 = 600,
                                  c2 = 0.75, c3 = 1, c02 = 100, c03 = 150,
                                  K = Inf, N = Inf, S = -Inf, 
                                  steps1 = 1, stepm1 = 0.95, stepl1 = 0.85,
                                  b11 = 1000, b21 = 2000, b31 = 3000,
                                  b12 = 1000, b22 = 1500, b32 = 2000,
                                  fixed = TRUE, rho = 0.3)[1], 
             utility_multiple_tte(n2 = 50, HRgo = 0.8, alpha = 0.025, beta = 0.1,
                                  hr1 = 0.5, hr2 = 0.6, id1 = 300, id2 = 600,
                                  c2 = 0.75, c3 = 1, c02 = 100, c03 = 150,
                                  K = Inf, N = Inf, S = -Inf, 
                                  steps1 = 1, stepm1 = 0.95, stepl1 = 0.85,
                                  b11 = 1000, b21 = 2000, b31 = 3000,
                                  b12 = 1000, b22 = 1500, b32 = 2000,
                                  fixed = TRUE, rho = 0.3)[1])
})

test_that("probability to go to phase III increases with lower hazard ratio", {
  skip_on_cran()
  set.seed(61216)
  expect_lte(pgo_multiple_tte(HRgo = 0.8, n2 = 50,  hr1 = 0.75, hr2 = 0.80,
                              id1 = 300, id2 = 600, fixed = FALSE, rho = 0.3), 
             pgo_multiple_tte(HRgo = 0.8, n2 = 50,  hr1 = 0.7, hr2 = 0.70,
                              id1 = 300, id2 = 600, fixed = FALSE, rho = 0.3))
})

test_that("probability to go to phase III increases with lower hazard ratio", {
  skip_on_cran()
  set.seed(61216)
  expect_equal(os_tte(HRgo = 0.8, n2 = 50, alpha = 0.05, beta = 0.1,
                    hr1 = 0.75, hr2 = 0.80,id1 = 300, id2 = 600,
                    fixed = FALSE, rho = 0.3,
                    rsamp = get_sample_multiple_tte(hr1 = 0.75,
                                                    hr2 = 0.80,
                                                    id1 = 300,
                                                    id2 = 600,
                                                    rho = 0.3)),
               0.3391463, tolerance = 1e-05)
})


test_that("pw works", {
  skip_on_cran()
  set.seed(61216)
  expect_equal(pw(n2 = 50,hr1 = 0.75, hr2 = 0.80,
                  id1 = 300, id2 = 600, fixed = FALSE, rho = 0.3), 0.54823305)
})

test_that("os_tte works", {
  skip_on_cran()
  set.seed(61216)
  expect_equal(os_tte(HRgo = 0.8, n2 = 50, alpha = 0.05, beta = 0.1,
                      hr1 = 0.75, hr2 = 0.80,id1 = 300, id2 = 600,
                      fixed = FALSE, rho = 0.3,
                      rsamp = get_sample_multiple_tte(hr1 = 0.75,
                                                      hr2 = 0.80,
                                                      id1 = 300,
                                                      id2 = 600,
                                                      rho = 0.3)),
               0.3391463, tolerance = 1e-05)
})
