library(testthat)
context("largestContinuousMinimum")
library(PeakSegJoint)

test_that("(start, end) when zero on both sides", {
  indices <- largestContinuousMinimum(
    c(0,   0, 1, 15, 2, 0),
    c(Inf, 1, 1,  1, 1, Inf))
  expect_equal(indices$start, 1)
  expect_equal(indices$end, 6)
})

test_that("(1, 2) when two zeros on left", {
  indices <- largestContinuousMinimum(
    c(0,   0, 1, 15, 2, 1),
    c(Inf, 1, 1,  1, 1, Inf))
  expect_equal(indices$start, 1)
  expect_equal(indices$end, 2)
})

test_that("(end, end) when one zero on right", {
  indices <- largestContinuousMinimum(
    c(1,   1, 1, 15, 2, 0),
    c(Inf, 1, 1,  1, 1, Inf))
  expect_equal(indices$start, 6)
  expect_equal(indices$end, 6)
})

test_that("(mid, mid) when one zero in middle", {
  indices <- largestContinuousMinimum(
    c(1,   1, 1, 0, 2, 2),
    c(Inf, 1, 1, 1, 1, Inf))
  expect_equal(indices$start, 4)
  expect_equal(indices$end, 4)
})

test_that("(minStart, minEnd) when two zeros in middle", {
  indices <- largestContinuousMinimum(
    c(1,   1, 0, 0, 2, 2),
    c(Inf, 1, 1, 1, 1, Inf))
  expect_equal(indices$start, 3)
  expect_equal(indices$end, 4)
})
