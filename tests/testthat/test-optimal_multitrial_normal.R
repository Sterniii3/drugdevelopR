test_that("optimal_multitrial_normal works for case 1", {
  skip_on_cran()
  expect_equal(optimal_multitrial_normal(w=0.3, Delta1 = 0.375, Delta2 = 0.625,
                                         in1=300, in2=600, a = 0.25, b = 0.75,
                                         n2min = 20, n2max = 100,stepn2 = 4,
                                         kappamin = 0.02, kappamax = 0.2,stepkappa = 0.02,
                                         alpha = 0.05, beta = 0.1,
                                         c2 = 0.675, c3 = 0.72, c02 = 15, c03 = 20,
                                         K = Inf, N = Inf, S = -Inf,
                                         b1 = 3000, b2 = 8000, b3 = 10000,
                                         case = 1, strategy = TRUE,
                                         fixed = TRUE,num_cl = 2)$u, c(1890.94,1928.79))
})
