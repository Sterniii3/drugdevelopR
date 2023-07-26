test_that("optimal_multiple_normal works for fixed = TRUE", {
  skip_on_cran()
  expect_equal(optimal_multiple_normal(Delta1 = 0.75, Delta2 = 0.80,
                                       in1=300, in2=600, sigma1 = 1, sigma2= 1,
                                       n2min = 30, n2max = 90, stepn2 = 10,
                                       kappamin = 0.05, kappamax = 0.2, stepkappa = 0.05,
                                       alpha = 0.05, beta = 0.1,
                                       c2 = 0.75, c3 = 1, c02 = 100, c03 = 150,
                                       K = Inf, N = Inf, S = -Inf,
                                       steps1 = 0,  stepm1 = 0.5, stepl1 = 0.8,
                                       b1 = 1000, b2 = 2000, b3 = 3000, 
                                       rho = 0.5, relaxed = TRUE,
                                       fixed = TRUE,  num_cl = 2)$u, 1096.15)
})
