library(testthat)
library(PeakSegJoint)
context("Faster")
data(H3K36me3.TDH.other.chunk1, envir=environment())
some.counts <- subset(
  H3K36me3.TDH.other.chunk1$counts,
  43000000 < chromEnd &
  chromStart < 43200000)
some.counts$sample.group <- some.counts$cell.type

fit <- PeakSegJointFaster(some.counts, 2:7)
max.samples <- fit$sample.modelSelection$complexity[1]
sample.id.vec <- names(fit$sample.loss.diff.vec[1:max.samples])
mean.mat <- fit$mean_mat[sample.id.vec,]
is.feasible <- mean.mat[,1] < mean.mat[,2] & mean.mat[,2] > mean.mat[,3]
test_that("all selectable samples are feasible", {
  expect_true(all(is.feasible))
})
test_that("sample loss diff is sorted", {
  expect_identical(sort(fit$sample.loss.diff.vec), fit$sample.loss.diff.vec)
})

fit <- PeakSegJointFaster(some.counts, 2:7)
max.groups <- fit$group.modelSelection$complexity[1]
group.id.vec <- names(fit$group.loss.diff.vec[1:max.groups])
sample.id.vec <- unlist(fit$group.list[group.id.vec])
mean.mat <- fit$mean_mat[sample.id.vec,]
is.feasible <- mean.mat[,1] < mean.mat[,2] & mean.mat[,2] > mean.mat[,3]
test_that("all selectable samples are feasible for max groups", {
  expect_true(all(is.feasible))
})
