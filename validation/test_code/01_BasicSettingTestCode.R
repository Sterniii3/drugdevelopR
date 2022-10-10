#' @editor Lukas D Sauer
#' @editDate 2022-09-13
test_that("01.01", {
  # Testing time to event with prior distribution
  res = optimal_tte(
    alpha = 0.025, # significance level
    beta = 0.1, # 1 - power
    hr1 = 0.69, hr2 = 0.88, # assumed treatment effects
    xi2 = 0.7, xi3 = 0.7, # event rates
    d2min = 10, d2max = 400, stepd2 = 1, # optimization region for the number of events
    hrgomin = 0.71, hrgomax = 0.95, stephrgo = 0.01, # optimization
    # region for the threshold values
    steps1 = 1, stepm1 = 0.95, stepl1 = 0.85, # boundaries for the effect size
    # categories small, medium and large
    b1 = 1000, b2 = 3000, b3 = 5000, # expected gains for each effect size
    num_cl = 3, # number of clusters
    c02 = 100, c03 = 150, # fixed cost for phase II and phase III
    c2 = 0.75, c3 = 1, # variable per-patient cost in phase II and phase III
    fixed = FALSE, # use a prior distribution
    w = 0.3, # weight for the prior distribution
    id1 = 210, id2 = 420, # amount of information (number of events) for prior
    # true treatment effect
  )
  expect_equal(res$n2, 206) # optimal sample size in phase II
  expect_equal(res$n3, 354) # resulting sample size in phase III
  expect_equal(res$n, 560) # resulting total sample size
  expect_equal(res$u, 432, tolerance = 0.05) # expected utility
  expect_equal(res$HRgo, 0.84) # threshold for proceeding to phase III
  expect_equal(res$d2, 144) # expected number of events in phase II
  expect_equal(res$d3, 248) # expected number of events in phase III
  expect_equal(res$d, 392) # total expected number of events
})

#' @editor Lukas D Sauer
#' @editDate 2022-10-10
test_that("01.02", {
  # Time to event with fixed treatment effects
  res = optimal_tte(
    alpha = 0.025, # significance level
    beta = 0.1, # 1 - power
    hr1 = 0.69, hr2 = 0.80, # assumed treatment effects
    xi2 = 0.7, xi3 = 0.7, # event rates
    d2min = 10, d2max = 400, stepd2 = 1, # optimization region for the number of events
    hrgomin = 0.71, hrgomax = 0.95, stephrgo = 0.01, # optimization
    # region for the threshold values
    steps1 = 1, stepm1 = 0.95, stepl1 = 0.85, # boundaries for the effect size
    # categories small, medium and large
    b1 = 1000, b2 = 3000, b3 = 5000, # expected gains for each effect size
    num_cl = 3, # number of clusters
    c02 = 100, c03 = 150, # fixed cost for phase II and phase III
    c2 = 0.75, c3 = 1, # variable per-patient cost in phase II and phase III
    fixed = TRUE, # use a prior distribution
    w = NULL, # weight for the prior distribution
    id1 = NULL, id2 = NULL, # amount of information (number of events) for prior
    # true treatment effect
  )
  expect_equal(res$n2, 240) # optimal sample size in phase II
  expect_equal(res$u, 352) # expected utility
  expect_equal(res$HRgo, 0.88) # threshold for proceeding to phase III
  expect_equal(res$pgo, 0.73) # probability to proceed to phase III
  expect_equal(res$d2, 168) # expected number of events in phase II
  expect_equal(res$d3, 546) # expected number of events in phase III
  expect_equal(res$d, 714) # total expected number of events
})

#' @editor Lukas D Sauer
#' @editDate 2022-10-10
test_that("01.03", {
  # Testing time to event with cost constraints
  res = optimal_tte(
    alpha = 0.025, # significance level
    beta = 0.1, # 1 - power
    hr1 = 0.69, hr2 = 0.88, # assumed treatment effects
    xi2 = 0.7, xi3 = 0.7, # event rates
    d2min = 10, d2max = 400, stepd2 = 1, # optimization region for the number of events
    hrgomin = 0.71, hrgomax = 0.95, stephrgo = 0.01, # optimization
    # region for the threshold values
    steps1 = 1, stepm1 = 0.95, stepl1 = 0.85, # boundaries for the effect size
    # categories small, medium and large
    b1 = 1000, b2 = 3000, b3 = 5000, # expected gains for each effect size
    num_cl = 3, # number of clusters
    c02 = 100, c03 = 150, # fixed cost for phase II and phase III
    c2 = 0.75, c3 = 1, # variable per-patient cost in phase II and phase III
    fixed = FALSE, # use a prior distribution
    w = 0.6, # weight for the prior distribution
    id1 = 210, id2 = 420, # amount of information (number of events) for prior
    # true treatment effect
    K = 750 # cost constraint
  )
  expect_equal(res$n2, 228) # optimal sample size in phase II
  expect_equal(res$u, 996, tolerance = 0.05) # expected utility
  expect_equal(res$HRgo, 0.84)
  expect_equal(res$K, 750) # expected number of events in phase II
  expect_equal(res$K2, 271) # expected number of events in phase III
  expect_equal(res$K3, 530) # total expected number of events
})

#' @editor Lukas D Sauer
#' @editDate 2022-10-10
test_that("01.04", {
  # Testing time to event with sample size constraint
  res = optimal_tte(
    alpha = 0.025, # significance level
    beta = 0.1, # 1 - power
    hr1 = 0.69, hr2 = 0.88, # assumed treatment effects
    xi2 = 0.7, xi3 = 0.7, # event rates
    d2min = 10, d2max = 400, stepd2 = 1, # optimization region for the number of events
    hrgomin = 0.71, hrgomax = 0.95, stephrgo = 0.01, # optimization
    # region for the threshold values
    steps1 = 1, stepm1 = 0.95, stepl1 = 0.85, # boundaries for the effect size
    # categories small, medium and large
    b1 = 1000, b2 = 3000, b3 = 5000, # expected gains for each effect size
    num_cl = 3, # number of clusters
    c02 = 100, c03 = 150, # fixed cost for phase II and phase III
    c2 = 0.75, c3 = 1, # variable per-patient cost in phase II and phase III
    fixed = FALSE, # use a prior distribution
    w = 0.6, # weight for the prior distribution
    id1 = 210, id2 = 420, # amount of information (number of events) for prior
    # true treatment effect
    N = 500 # sample size constraint
  )
  expect_equal(res$n2, 170) # optimal sample size in phase II
  expect_equal(res$n3, 328) # resulting sample size in phase III
  expect_equal(res$n, 498) # resulting total sample size
  expect_equal(res$u, 956, tolerance = 0.05) # expected utility
  expect_equal(res$HRgo, 0.83) # threshold for proceeding to phase III
})



#' @editor Lukas D Sauer
#' @editDate 2022-10-10
test_that("01.05", {
  # Testing time to event with success probability constraint
  res = optimal_tte(
    alpha = 0.025, # significance level
    beta = 0.1, # 1 - power
    hr1 = 0.69, hr2 = 0.88, # assumed treatment effects
    xi2 = 0.7, xi3 = 0.7, # event rates
    d2min = 10, d2max = 400, stepd2 = 1, # optimization region for the number of events
    hrgomin = 0.71, hrgomax = 0.95, stephrgo = 0.01, # optimization
    # region for the threshold values
    steps1 = 1, stepm1 = 0.95, stepl1 = 0.85, # boundaries for the effect size
    # categories small, medium and large
    b1 = 1000, b2 = 3000, b3 = 5000, # expected gains for each effect size
    num_cl = 3, # number of clusters
    c02 = 100, c03 = 150, # fixed cost for phase II and phase III
    c2 = 0.75, c3 = 1, # variable per-patient cost in phase II and phase III
    fixed = FALSE, # use a prior distribution
    w = 0.6, # weight for the prior distribution
    id1 = 210, id2 = 420, # amount of information (number of events) for prior
    # true treatment effect
    S = 0.6 # constraint on success probability
  )
  expect_equal(res$n2, 470) # optimal sample size in phase II
  expect_equal(res$u, 899, tolerance = 0.01) # expected utility
  expect_equal(res$HRgo, 0.89) # threshold for proceeding to phase III
  expect_equal(res$pgo, 0.77) # probability to proceed to phase III
  expect_equal(res$sProg, 0.6) # probability of a successful program
})

#' @editor Lukas D Sauer
#' @editDate 2022-10-10
test_that("01.06", {
  # Testing time to event while to skipping phase II
  res = optimal_tte(
    alpha = 0.025, # significance level
    beta = 0.1, # 1 - power
    hr1 = 0.69, hr2 = 0.88, # assumed treatment effects
    xi2 = 0.7, xi3 = 0.7, # event rates
    d2min = 10, d2max = 400, stepd2 = 1, # optimization region for the number of events
    hrgomin = 0.71, hrgomax = 0.95, stephrgo = 0.01, # optimization
    # region for the threshold values
    steps1 = 1, stepm1 = 0.95, stepl1 = 0.85, # boundaries for the effect size
    # categories small, medium and large
    b1 = 1000, b2 = 3000, b3 = 5000, # expected gains for each effect size
    num_cl = 3, # number of clusters
    c02 = 100, c03 = 150, # fixed cost for phase II and phase III
    c2 = 0.75, c3 = 1, # variable per-patient cost in phase II and phase III
    fixed = FALSE, # use a prior distribution
    w = 0.6, # weight for the prior distribution
    id1 = 210, id2 = 420, # amount of information (number of events) for prior
    # true treatment effect
    skipII = TRUE # skip phase II
  )
  expect_equal(res$n2, 824) # optimal sample size in phase II
  expect_equal(res$u, 1706) # expected utility
  expect_equal(res$HRgo, 0.76) # threshold for proceeding to phase III
})

#' @editor Lukas D Sauer
#' @editDate 2022-10-10