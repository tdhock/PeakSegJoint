library(testthat)
context("qsub-pipeline")

exampleData <-
  system.file("exampleData",
              mustWork=TRUE,
              package="PeakSegJoint")
labels.txt <- file.path(exampleData, "manually_annotated_region_labels.txt")
overlapping.txt <- file.path(exampleData, "overlapping_labels.txt")
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
unlink(pred.RData)

AllSteps <-
  system.file("exec", "00_AllSteps_qsub.R",
              mustWork=TRUE,
              package="PeakSegJoint")
cmd <- paste("QSUB='echo INTERACTIVE && bash' Rscript", AllSteps, labels.txt)
system(cmd)

test_that("pipeline trained on 4 samples predicts for 8 samples", {
  load(pred.RData)
  expect_equal(nrow(all.peaks.mat), 8)
  starts <- unique(all.peaks.df$chromStart)
  expect_equal(length(starts), ncol(all.peaks.mat))
})

test_that("pipeline generates test error plots", {
  test.errors.dir <- 
    system.file("exampleData", "PeakSegJoint-chunks", "figure-test-errors",
                mustWork=TRUE,
                package="PeakSegJoint")
  test.errors.files <- dir(test.errors.dir)
  expected.files <- paste0(1:3, ".png")
  expect_true(all(expected.files %in% test.errors.files))
})
