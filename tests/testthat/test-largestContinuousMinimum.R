library(testthat)
context("largestContinuousMinimum")
library(PeakSegJoint)

test_that("(-Inf, Inf) when zero on both sides", {
  indices <- largestContinuousMinimum(c(0, 15, 0), c(1, 1, 1))
  expect_equal(indices$start, 1)
  expect_equal(indices$end, 3)
})
