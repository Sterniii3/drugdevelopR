
# Test setup


#' @editor Lukas D Sauer
#' @editDate 2022-09-13
test_that("01.01", {
  res = optimal_tte(
    alpha = 0.025, # significance level
    beta = 0.1, # 1- power
    hr1 = 0.8, hr2 = 0.65, # assumed treatment effects
    xi2 = 0.7, xi3 = 0.7, # event rates
    d2min = 10, d2max = 400, stepd2 = 1, # optimization region for the number of events
    hrgomin = -log(0.95), hrgomax = -log(0.7), stephrgo = 0.1, # optimization
    # region for the threshold values
    steps1 = 1, stepm1 = 0.95, stepl1 = 0.85, # boundaries for the effect size
    # categories small, medium and large
    b1 = 1000, b2 = 3000, b3 = 5000, # expected gains for each effect size
    num_cl = 3, # number of clusters
    c02 = 100, c03 = 150, # fixed cost for phase II and phase III
    c2 = 7.5, c3 = 10, # variable per-patient cost in phase II and phase III
    fixed = FALSE, # use a prior distribution
    w = 0.3, # weight for the prior distribution
    id1 = 210, id2 = 420, # amount of information (number of events) for prior
    # true treatment effect
  )
  expect_equal(res$n2, 206) # optimal sample size in phase II
  expect_equal(res$n3, 354) # resulting sample size in phase III
  expect_equal(res$n, 560) # resulting total sample size
  expect_equal(res$u, 432) # expected utility
  expect_equal(res$d2, 144) # expected number of events in phase II
  expect_equal(res$d3, 248) # expected number of events in phase III
  expect_equal(res$d, 392) # total expected number of events
})



