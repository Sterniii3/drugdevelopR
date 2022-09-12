test_that("pgo_tte: Setting 1.21", {
  expect_equal(pgo_tte(HRgo = 0.8, n2 = 48 ,
                       ec = 0.6, hr1 = 0.7, hr2 = 0.8, 
                       strategy = 1, case = 21), 
               0.4524518)
})

test_that("pgo_tte works: Setting 1.22", {
  expect_equal(pgo_tte(HRgo = 0.8, n2 = 48 ,
                       ec = 0.6, hr1 = 0.7, hr2 = 0.8, 
                       strategy = 1, case = 22), 
               0.27173794)
})

test_that("pgo_tte works: Setting 2.31", {
  expect_equal(pgo_tte(HRgo = 0.8, n2 = 48 ,
                       ec = 0.6, hr1 = 0.7, hr2 = 0.8, 
                       strategy = 2, case = 31), 
               0.22826206)
})

test_that("pgo_tte works: Setting 2.32", {
  expect_equal(pgo_tte(HRgo = 0.8, n2 = 48 ,
                       ec = 0.6, hr1 = 0.7, hr2 = 0.8, 
                       strategy = 2, case = 32), 
               0.15580054)
})


test_that("ss_tte works: Setting l=1", {
  expect_equal(ss_tte(alpha = 0.05, beta = 0.1,
                      ec = 0.6, ek = 0.8, y = 0.5, l=1), 
               199.823105)
})

test_that("ss_tte works: Setting l=2", {
  expect_equal(ss_tte(alpha = 0.05, beta = 0.1, 
                      ec = 0.6, ek = 0.8, y = 0.5, l=2), 
               357.94104)
})

test_that("Ess_tte works: Setting: Case 1", {
  expect_equal(Ess_tte(HRgo = 0.8 ,n2 = 50 ,alpha = 0.05, beta = 0.1,
                       ec = 0.6, hr1 = 0.7, hr2 = 0.8, strategy = NULL, case = 1), 
               0)
})

test_that("Ess_tte works: Setting: 1.21", {
  expect_equal(Ess_tte(HRgo = 0.8 ,n2 = 50 ,alpha = 0.05, beta = 0.1,
                       ec = 0.6, hr1 = 0.7, hr2 = 0.8, strategy = 1, case = 21), 
               116.791988)
})

test_that("Ess_tte works: Setting: 1.22", {
  expect_equal(Ess_tte(HRgo = 0.8 ,n2 = 50 ,alpha = 0.05, beta = 0.1,
                       ec = 0.6, hr1 = 0.7, hr2 = 0.8, strategy = 1, case = 22), 
               72.736693)
})

test_that("Ess_tte works: Setting: 2.21", {
  expect_equal(Ess_tte(HRgo = 0.8 ,n2 = 50 ,alpha = 0.05, beta = 0.1,
                       ec = 0.6, hr1 = 0.7, hr2 = 0.8, strategy = 2, case = 21), 
               83.737493)
})

test_that("Ess_tte works: Setting: 2.22", {
  expect_equal(Ess_tte(HRgo = 0.8 ,n2 = 50 ,alpha = 0.05, beta = 0.1,
                       ec = 0.6, hr1 = 0.7, hr2 = 0.8, strategy = 2, case = 22), 
               47.967458)
})

test_that("Ess_tte works: Setting: 2.31", {
  expect_equal(Ess_tte(HRgo = 0.8 ,n2 = 50 ,alpha = 0.05, beta = 0.1,
                       ec = 0.6, hr1 = 0.7, hr2 = 0.8, strategy = 2, case = 31), 
               157.001055)
})

test_that("Ess_tte works: Setting: 2.32", {
  expect_equal(Ess_tte(HRgo = 0.8 ,n2 = 50 ,alpha = 0.05, beta = 0.1,
                       ec = 0.6, hr1 = 0.7, hr2 = 0.8, strategy = 2, case = 32), 
               109.85048)
})

test_that("PsProg_tte works: Setting: Case 1", {
  expect_equal(PsProg_tte(HRgo = 0.8 ,n2 = 50 ,alpha = 0.05, beta = 0.1,
                          ec = 0.6, hr1 = 0.7, hr2 = 0.8, step1 = 1, step2 = 0.95,
                          strategy = NULL, case = 1), 
               0)
})

test_that("PsProg_tte works: Setting: 1.21", {
  expect_equal(PsProg_tte(HRgo = 0.8 ,n2 = 50 ,alpha = 0.05, beta = 0.1,
                          ec = 0.6, hr1 = 0.7, hr2 = 0.8, step1 = 1, step2 = 0.95,
                          strategy = 1, case = 21), 
               0.033617034)
})

test_that("PsProg_tte works: Setting: 1.22", {
  expect_equal(PsProg_tte(HRgo = 0.8 ,n2 = 50 ,alpha = 0.05, beta = 0.1,
                          ec = 0.6, hr1 = 0.7, hr2 = 0.8, step1 = 1, step2 = 0.95,
                          strategy = 1, case = 22), 
               0.025050005)
})

test_that("PsProg_tte works: Setting: 2.21", {
  expect_equal(PsProg_tte(HRgo = 0.8 ,n2 = 50 ,alpha = 0.05, beta = 0.1,
                          ec = 0.6, hr1 = 0.7, hr2 = 0.8, step1 = 1, step2 = 0.95,
                          strategy = 2, case = 21), 
               0.0179849417)
})

test_that("PsProg_tte works: Setting: 2.22", {
  expect_equal(PsProg_tte(HRgo = 0.8 ,n2 = 50 ,alpha = 0.05, beta = 0.1,
                          ec = 0.6, hr1 = 0.7, hr2 = 0.8, step1 = 1, step2 = 0.95,
                          strategy = 2, case = 22), 
               0.0144852983)
})

test_that("PsProg_tte works: Setting: 2.31", {
  expect_equal(PsProg_tte(HRgo = 0.8 ,n2 = 50 ,alpha = 0.05, beta = 0.1,
                          ec = 0.6, hr1 = 0.7, hr2 = 0.8, step1 = 1, step2 = 0.95,
                          strategy = 2, case = 31), 
               0.0204096092)
})

test_that("PsProg_tte works: Setting: 2.32", {
  expect_equal(PsProg_tte(HRgo = 0.8 ,n2 = 50 ,alpha = 0.05, beta = 0.1,
                          ec = 0.6, hr1 = 0.7, hr2 = 0.8, step1 = 1, step2 = 0.95,
                          strategy = 2, case = 32), 
               0.013853675)
})

