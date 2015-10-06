library(testthat)
library(PeakSegJoint)
context("qsub-pipeline")

AllSteps <-
  system.file("exec", "00_AllSteps_qsub.R",
              mustWork=TRUE,
              package="PeakSegJoint")
Rscript <- "QSUB='echo INTERACTIVE && bash' JOBS=2 Rscript"

orig.exampleData <-
  system.file("exampleData",
              mustWork=TRUE,
              package="PeakSegJoint")
exampleDir <- function(){
  tdir <- tempfile()
  dir.create(tdir)
  file.copy(orig.exampleData, tdir, recursive=TRUE)
  file.path(tdir, "exampleData")
}

test_that("overlapping chunks in different files is an error", {
  overlapping.chunks <- exampleDir()
  labels.txt <-
    file.path(overlapping.chunks, "manually_annotated_region_labels.txt")
  overlapping.txt <- file.path(overlapping.chunks, "overlapping_labels.txt")
  Step0 <-
    system.file("exec", "Step0-convert-labels.R",
                mustWork=TRUE,
                package="PeakSegJoint")
  cmd.args <- paste(Step0, labels.txt, overlapping.txt)
  out.lines <- system2("Rscript", cmd.args, stderr=TRUE, stdout=TRUE)
  out.txt <- paste(out.lines, collapse="\n")
  expect_match(out.txt, "chunks in different label files should not overlap")
})

test_that("pipeline trained on 4 samples predicts for 8 samples", {
  three.chunks <- exampleDir()
  labels.txt <-
    file.path(three.chunks, "manually_annotated_region_labels.txt")
  cmd <- paste(Rscript, AllSteps, labels.txt)
  system(cmd)
  ## There should be summary bed files for labels and peaks:
  labels.bed.gz <- file.path(three.chunks, "all_labels.bed.gz")
  bed.labels <- read.table(labels.bed.gz, skip=1)
  expect_equal(dim(bed.labels), c(16, 9))
  summary.bed.gz <- file.path(three.chunks, "PeakSegJoint.summary.bed.gz")
  bed.summary <- read.table(summary.bed.gz, skip=1)
  expect_equal(ncol(bed.summary), 5)
  ## There should be predictions for 8 samples:
  bed.gz.vec <- Sys.glob(file.path(three.chunks, "*", "*.bed.gz"))
  expect_equal(length(bed.gz.vec), 8)
  pred.RData <- file.path(three.chunks, "PeakSegJoint.predictions.RData")
  pred.objs <- load(pred.RData)
  expect_equal(nrow(all.peaks.mat), 8)
  starts <- unique(all.peaks.df$chromStart)
  expect_equal(length(starts), ncol(all.peaks.mat))
  test.errors.dir <- 
    file.path(three.chunks, "PeakSegJoint-chunks", "figure-test-errors")
  test.errors.files <- dir(test.errors.dir)
  ## pipeline trained on 3 chunks generates 3 test error plots:
  expected.files <- paste0(1:3, ".png")
  expect_true(all(expected.files %in% test.errors.files))
  ## pipeline trained on 3 chunks generates test error summary:
  expect_true("figure-test-error-decreases.png" %in% test.errors.files)
  ## pipeline trained on 1 file does 3 fold CV:
  index.html <- file.path(test.errors.dir, "index.html")
  index.lines <- readLines(index.html)
  cv.line <- grep("cross-validation", index.lines, value=TRUE)
  expect_match(cv.line, "3 fold")
  ## check for no peak filtering.
  expect_true("specific.error" %in% pred.objs)
  expect_null(specific.error)
})

test_that("pipeline trained on 8 samples predicts for 8 samples", {
  six.chunks <- exampleDir()
  labels.txt <- file.path(six.chunks, "manually_annotated_region_labels.txt")
  other.txt <- file.path(six.chunks, "other_labels.txt")
  cmd <- paste(Rscript, AllSteps, labels.txt, other.txt)
  system(cmd)
  ## There should be summary bed files for labels and peaks:
  labels.bed.gz <- file.path(six.chunks, "all_labels.bed.gz")
  bed.labels <- read.table(labels.bed.gz, skip=1)
  expect_equal(dim(bed.labels), c(28, 9))
  summary.bed.gz <- file.path(six.chunks, "PeakSegJoint.summary.bed.gz")
  bed.summary <- read.table(summary.bed.gz, skip=1)
  expect_equal(ncol(bed.summary), 5)
  ## There should be predictions for 8 samples:
  bed.gz.vec <- Sys.glob(file.path(six.chunks, "*", "*.bed.gz"))
  expect_equal(length(bed.gz.vec), 8)
  pred.RData <- file.path(six.chunks, "PeakSegJoint.predictions.RData")
  pred.objs <- load(pred.RData)
  expect_equal(nrow(all.peaks.mat), 8)
  starts <- unique(all.peaks.df$chromStart)
  expect_equal(length(starts), ncol(all.peaks.mat))
  test.errors.dir <- 
    file.path(six.chunks, "PeakSegJoint-chunks", "figure-test-errors")
  test.errors.files <- dir(test.errors.dir)
  ## pipeline trained on 6 chunks generates 6 test error plots:
  expected.files <- paste0(1:6, ".png")
  expect_true(all(expected.files %in% test.errors.files))
  ## pipeline trained on 6 chunks generates test error summary:
  expect_true("figure-test-error-decreases.png" %in% test.errors.files)
  ## pipeline trained on 2 files does 2 fold CV: 
  index.html <- file.path(test.errors.dir, "index.html")
  index.lines <- readLines(index.html)
  cv.line <- grep("cross-validation", index.lines, value=TRUE)
  expect_match(cv.line, "2 fold")
  ## check for no peak filtering.
  expect_true("specific.error" %in% pred.objs)
  expect_null(specific.error)
})

