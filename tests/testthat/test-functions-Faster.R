library(testthat)
library(PeakSegJoint)
context("Faster")
data(H3K36me3.TDH.other.chunk1, envir=environment())
nogroup <- subset(
  H3K36me3.TDH.other.chunk1$counts,
  43000000 < chromEnd &
  chromStart < 43200000)
test_that("Faster needs group", {
  expect_error({
    fit <- PeakSegJointFaster(nogroup, 2:7)
  }, "no sample.group")
})

some.counts <- data.frame(
  nogroup,
  sample.group=nogroup$cell.type)
plist <- ProfileList(some.counts)
fit <- PeakSegJointFaster(plist, 2:7)
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
test_that("biggest selectable model has min.loss", {
  expect_equal(fit$sample.modelSelection$loss[1], fit$min.loss)
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

## what if there are no samples up?
pathological <- data.frame(
  sample.id="foo",
  sample.group="bar",
  chromStart=0:9,
  chromEnd=1:10,
  count=1:10)
fit <- PeakSegJointFaster(pathological, 2L)
test_that("no feasible samples is fine", {
  expect_identical(fit$sample.modelSelection$complexity, 0L)
  expect_identical(fit$group.modelSelection$complexity, 0L)
})

test_that("informative error", {
  expect_error({
    fit <- PeakSegJointFaster(pathological, 3L)
  }, "bin factor too large")
})
