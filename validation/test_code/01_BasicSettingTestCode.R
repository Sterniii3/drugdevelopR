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
    num_cl = 12, # number of clusters
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
    hr1 = 0.8, hr2 = NULL, # assumed treatment effects
    xi2 = 0.7, xi3 = 0.7, # event rates
    d2min = 10, d2max = 400, stepd2 = 1, # optimization region for the number of events
    hrgomin = 0.71, hrgomax = 0.95, stephrgo = 0.01, # optimization
    # region for the threshold values
    steps1 = 1, stepm1 = 0.95, stepl1 = 0.85, # boundaries for the effect size
    # categories small, medium and large
    b1 = 1000, b2 = 3000, b3 = 5000, # expected gains for each effect size
    num_cl = 12, # number of clusters
    c02 = 100, c03 = 150, # fixed cost for phase II and phase III
    c2 = 0.75, c3 = 1, # variable per-patient cost in phase II and phase III
    fixed = TRUE, # use a prior distribution
    w = NULL, # weight for the prior distribution
    id1 = NULL, id2 = NULL, # amount of information (number of events) for prior
    # true treatment effect
  )
  expect_equal(res$n2, 240) # optimal sample size in phase II
  expect_equal(res$u, 352, tolerance = 0.001) # expected utility
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
    num_cl = 12, # number of clusters
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
  expect_equal(res$K2, 271)
  expect_equal(res$K3, 478)
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
    num_cl = 12, # number of clusters
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
    num_cl = 12, # number of clusters
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
    num_cl = 12, # number of clusters
    c02 = 100, c03 = 150, # fixed cost for phase II and phase III
    c2 = 0.75, c3 = 1, # variable per-patient cost in phase II and phase III
    fixed = FALSE, # use a prior distribution
    w = 0.6, # weight for the prior distribution
    id1 = 210, id2 = 420, # amount of information (number of events) for prior
    # true treatment effect
    skipII = TRUE # skip phase II
  )
  expect_equal(res[[2]]$n3, 824) # optimal sample size in phase III after skipping phase II
  expect_equal(res[[2]]$u, 1706, tolerance = 0.0001) # expected utility
  expect_equal(res[[2]]$HR, 0.76) # threshold for proceeding to phase III
})

#' @editor Lukas D Sauer
#' @editDate 2022-10-11
test_that("01.07", {
  # Testing time to event while modeling differing population structures in phase II and III
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
    num_cl = 12, # number of clusters
    c02 = 100, c03 = 150, # fixed cost for phase II and phase III
    c2 = 0.75, c3 = 1, # variable per-patient cost in phase II and phase III
    fixed = FALSE, # use a prior distribution
    w = 0.6, # weight for the prior distribution
    id1 = 210, id2 = 420, # amount of information (number of events) for prior
    # true treatment effect
    gamma = 0.025
  )
  expect_equal(res$n2, 310) # optimal sample size in phase II
  expect_equal(res$u, 1207, tolerance = 0.0005) # expected utility
  expect_equal(res$HRgo, 0.86) # threshold for proceeding to phase III
})

#' @editor Lukas D Sauer
#' @editDate 2022-10-12
test_that("01.08", {
  # Testing binary endpoints with fixed effects
  res = optimal_binary(alpha = 0.025,
                       beta = 0.1,
                       p0 = 0.6, p11 = 0.5, p12 = NULL,
                       n2min = 10, n2max = 500, stepn2 = 2,
                       rrgomin = 0.7, rrgomax = 0.9, steprrgo = 0.01,
                       steps1 = 1, stepm1 = 0.95, stepl1 = 0.85,
                       b1 = 1000, b2 = 3000, b3 = 5000,
                       num_cl = 12,
                       c02 = 100, c03 = 150,
                       c2 = 0.75, c3 = 1,
                       fixed = TRUE,
                       w = NULL,
                       in1 = NULL, in2 = NULL
                       )
  expect_equal(res$n2, 204)
  expect_equal(res$u, 299, tolerance = 0.0005)
  expect_equal(res$RRgo, 0.90)
  expect_equal(res$K, Inf)
  expect_equal(res$K2, 253)
  expect_equal(res$K3, 810)
})

#' @editor Lukas D Sauer
#' @editDate 2022-10-12
test_that("01.09", {
  # Testing binary endpoints with treatment effects modelled on a prior distribution
  res = optimal_binary(alpha = 0.025,
                       beta = 0.1,
                       p0 = 0.6, p11 = 0.3, p12 = 0.5,
                       n2min = 10, n2max = 500, stepn2 = 2,
                       rrgomin = 0.7, rrgomax = 0.9, steprrgo = 0.01,
                       steps1 = 1, stepm1 = 0.95, stepl1 = 0.85,
                       b1 = 1000, b2 = 3000, b3 = 5000,
                       num_cl = 12,
                       c02 = 100, c03 = 150,
                       c2 = 0.75, c3 = 1,
                       fixed = FALSE,
                       w = 0.4,
                       in1 = 30, in2 = 60
  )
  expect_equal(res$n2, 224)
  expect_equal(res$u, 1542, tolerance = 0.0005)
  expect_equal(res$RRgo, 0.89)
})

#' @editor Lukas D Sauer
#' @editDate 2022-10-11
test_that("01.10", {
  # Testing normally distributed endpoints with prior distribution on effects
  res = optimal_normal(
    alpha = 0.025,
    beta = 0.1,
    Delta1 = 0.375, Delta2 = 0.5,
    n2min = 10, n2max = 500, stepn2 = 2,
    kappamin = 0.01, kappamax = 0.5, stepkappa = 0.01,
    steps1 = 0, stepm1 = 0.375, stepl1 = 0.625,
    b1 = 625, b2 = 2000, b3 = 10000,
    num_cl = 12,
    c02 = 15, c03 = 20,
    c2 = 0.675, c3 = 0.72,
    fixed = FALSE,
    w = 0.5,
    in1 = 300, in2 = 600,
    a = 0, b = 0.75
  )
  expect_equal(res$n2, 86) # optimal sample size in phase II
  expect_equal(res$u, 337, tolerance = 0.001) # expected utility
  expect_equal(res$Kappa, 0.19) # threshold for proceeding to phase III
})

#' @editor Lukas D Sauer
#' @editDate 2022-10-11
test_that("01.11", {
  # Testing normally distributed endpoints with fixed effects
  res = optimal_normal(
    alpha = 0.025,
    beta = 0.1,
    Delta1 = 0.625, Delta2 = NULL,
    n2min = 10, n2max = 500, stepn2 = 2,
    kappamin = 0.01, kappamax = 0.5, stepkappa = 0.01,
    steps1 = 0, stepm1 = 0.375, stepl1 = 0.625,
    b1 = 625, b2 = 2000, b3 = 10000,
    num_cl = 12,
    c02 = 15, c03 = 20,
    c2 = 0.675, c3 = 0.72,
    fixed = TRUE,
    w = NULL,
    # 300 events in phase II and 600 events in phase III, Wo finde ich diese Parameter?
    # Vielleicht ist folgendes gemeint?
    in1 = NULL, in2 = NULL,
    a = NULL, b = NULL
  )
  expect_equal(res$n2, 78) # optimal sample size in phase II
  expect_equal(res$u, 944, tolerance = 0.0001) # expected utility
  expect_equal(res$Kappa, 0.12) # threshold for proceeding to phase III
  expect_equal(res$sProg, 0.83) # probability of a successful program
  expect_equal(res$sProg1, 0.51) # probability of a successful program with small treatment effect
  expect_equal(res$sProg2, 0.30) # probability of a successful program with medium treatment effect
  expect_equal(res$sProg3, 0.02) # probability of a successful program with large treatment effect
  }
)



#' @editor Lukas D Sauer
#' @editDate 2022-10-11
test_that("01.12", {
  # Testing that parallel computing has an effect
  start_time_3 = Sys.time()
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
    num_cl = 12, # number of clusters
    c02 = 100, c03 = 150, # fixed cost for phase II and phase III
    c2 = 0.75, c3 = 1, # variable per-patient cost in phase II and phase III
    fixed = FALSE, # use a prior distribution
    w = 0.3, # weight for the prior distribution
    id1 = 210, id2 = 420, # amount of information (number of events) for prior
    # true treatment effect
  )
  end_time_3 = Sys.time()
  time_elapsed_01_02_num_cl_3 = end_time_3 - start_time_3
  start_time_1 = Sys.time()
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
    num_cl = 1, # number of clusters
    c02 = 100, c03 = 150, # fixed cost for phase II and phase III
    c2 = 0.75, c3 = 1, # variable per-patient cost in phase II and phase III
    fixed = FALSE, # use a prior distribution
    w = 0.3, # weight for the prior distribution
    id1 = 210, id2 = 420, # amount of information (number of events) for prior
    # true treatment effect
  )
  end_time_1 = Sys.time()
  time_elapsed_01_02_num_cl_1 = end_time_1 - start_time_1
  expect_true(time_elapsed_01_02_num_cl_1 > time_elapsed_01_02_num_cl_3)
}
)
#' @editor Lukas D Sauer
#' @editDate 2023-03-09
test_that("01.14", {
  # Comparing results of Epgo_binary() and En3_binary() to SAS results
  data <- haven::read_sas("./validation/ref/valref_binary.sas7bdat")
  res_r <- rep(0, nrow(data))
  for(i in (1:nrow(data))){
    res_r[i] <- En3_binary(RRgo = data$RRgo_vec[i], n2 = data$n2_vec[i],
                           alpha = data$alpha_vec[i], beta = data$beta_vec[i],
                           p0 = data$p0_vec[i], w = NULL, p11 = data$p11_vec[i],
                           p12 = NULL, in1 = NULL, in2 = NULL, fixed = TRUE)
  }
  expect_equal(data$res_sas, res_r)
})
