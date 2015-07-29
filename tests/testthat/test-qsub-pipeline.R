context("qsub-pipeline")

exampleData <-
  system.file("exampleData",
              mustWork=TRUE,
              package="PeakSegJoint")

pred.RData <- file.path(exampleData, "PeakSegJoint.predictions.RData")
unlink(pred.RData)

AllSteps <-
  system.file("exec", "00_AllSteps_qsub.R",
              mustWork=TRUE,
              package="PeakSegJoint")
labels.txt <- file.path(exampleData, "manually_annotated_region_labels.txt")
cmd <- paste("QSUB='echo INTERACTIVE && bash' Rscript", AllSteps, labels.txt)
system(cmd)

test_that("pipeline generates predicted peaks", {
  load(pred.RData)
  expect_equal(nrow(all.peaks.mat), 4)
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
