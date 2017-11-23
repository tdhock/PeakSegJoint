library(testthat)
context("mean mat")
library(PeakSegJoint)

data(H3K4me3.TDH.other.chunk8, envir=environment())
fit <- PeakSegJointHeuristic(H3K4me3.TDH.other.chunk8, 3)
test_that("model10 mean equal to model mean", {
  expect_equal(fit$models[[11]]$seg1_mean_vec, fit$mean_mat[,1])
  expect_equal(fit$models[[11]]$seg2_mean_vec, fit$mean_mat[,2])
  expect_equal(fit$models[[11]]$seg3_mean_vec, fit$mean_mat[,3])
})
