library(testthat)
library(PeakSegJoint)
context("qsub-pipeline")

exampleData <-
  system.file("exampleData",
              mustWork=TRUE,
              package="PeakSegJoint")
labels.txt <- file.path(exampleData, "manually_annotated_region_labels.txt")
overlapping.txt <- file.path(exampleData, "overlapping_labels.txt")
other.txt <- file.path(exampleData, "other_labels.txt")
Step0 <-
  system.file("exec", "Step0-convert-labels.R",
              mustWork=TRUE,
              package="PeakSegJoint")
cmd.args <- paste(Step0, labels.txt, overlapping.txt)

test_that("overlapping chunks in different files is an error", {
  out.lines <- system2("Rscript", cmd.args, stderr=TRUE, stdout=TRUE)
  out.txt <- paste(out.lines, collapse="\n")
  expect_match(out.txt, "chunks in different label files should not overlap")
})

pred.RData <- file.path(exampleData, "PeakSegJoint.predictions.RData")
AllSteps <-
  system.file("exec", "00_AllSteps_qsub.R",
              mustWork=TRUE,
              package="PeakSegJoint")

unlink(pred.RData)
cmd <- paste("QSUB='echo SEQUENTIAL && bash' Rscript", AllSteps, labels.txt)
system(cmd)

test_that("pipeline trained on 4 samples predicts for 8 samples", {
  load(pred.RData)
  expect_equal(nrow(all.peaks.mat), 8)
  starts <- unique(all.peaks.df$chromStart)
  expect_equal(length(starts), ncol(all.peaks.mat))
})

test.errors.dir <- 
  system.file("exampleData", "PeakSegJoint-chunks", "figure-test-errors",
              mustWork=TRUE,
              package="PeakSegJoint")

test.errors.files <- dir(test.errors.dir)

test_that("pipeline trained on 3 chunks generates 3 test error plots", {
  expected.files <- paste0(1:3, ".png")
  expect_true(all(expected.files %in% test.errors.files))
})

test_that("pipeline trained on 3 chunks generates test error summary", {
  expect_true("figure-test-error-decreases.png" %in% test.errors.files)
})

test_that("pipeline trained on 1 file does 3 fold CV", {
  index.html <- file.path(test.errors.dir, "index.html")
  index.lines <- readLines(index.html)
  cv.line <- grep("cross-validation", index.lines, value=TRUE)
  expect_match(cv.line, "3 fold")
})

unlink(pred.RData)
cmd <-
  paste("QSUB='echo INTERACTIVE && bash' Rscript",
        AllSteps, labels.txt, other.txt)
system(cmd)

test_that("pipeline trained on 8 samples predicts for 8 samples", {
  load(pred.RData)
  expect_equal(nrow(all.peaks.mat), 8)
  starts <- unique(all.peaks.df$chromStart)
  expect_equal(length(starts), ncol(all.peaks.mat))
})

test.errors.files <- dir(test.errors.dir)

test_that("pipeline trained on 6 chunks generates 6 test error plots", {
  expected.files <- paste0(1:6, ".png")
  expect_true(all(expected.files %in% test.errors.files))
})

test_that("pipeline trained on 6 chunks generates test error summary", {
  expect_true("figure-test-error-decreases.png" %in% test.errors.files)
})

test_that("pipeline trained on 2 files does 2 fold CV", {
  index.html <- file.path(test.errors.dir, "index.html")
  index.lines <- readLines(index.html)
  cv.line <- grep("cross-validation", index.lines, value=TRUE)
  expect_match(cv.line, "2 fold")
})
