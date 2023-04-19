test_that("Type binary", {
  expect_type(optimal_return_doc(type = "binary"), 
              "character")
  })

test_that("multitrial_tte",{
  expect_type(optimal_return_doc(type = "binary",setting = "bias"),
              "character")
})

test_that("multitrial_tte",{
  expect_type(optimal_return_doc(type = "tte",setting = "multitrial"),
              "character")
})

test_that("multitrial_tte",{
  expect_type(optimal_return_doc(type = "tte",setting = "multiarm"),
               "character")
})

test_that("multiple_normal",{
  expect_type(optimal_return_doc(type = "normal",setting = "multiple"),
              "character")
})