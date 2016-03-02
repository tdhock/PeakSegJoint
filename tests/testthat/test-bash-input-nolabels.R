library(testthat)
library(PeakSegJoint)
context("qsub input")

AllSteps <-
  system.file("exec", "00_AllSteps_qsub.R",
              mustWork=TRUE,
              package="PeakSegJoint")
Rscript <- "QSUB='echo INTERACTIVE && bash' JOBS=2 Rscript"

## This needs to be a function so we can call on.exit(setwd())
download.zip <- function(){
  require(httr)
  zip.url <- "https://github.com/tdhock/input-test-data/archive/master.zip"
  owd <- setwd(tempdir())
  on.exit(setwd(owd))
  request <- GET(zip.url)
  stop_for_status(request)
  writeBin(content(request), "input-test-data-master.zip")
  unlink("input-test-data-master", recursive=TRUE)
  unzip("input-test-data-master.zip")
}

test_that("4 un-labeled input + 4 labeled H3K36me3", {
  download.zip()
  data.dir <- file.path(tempdir(), "input-test-data-master") 
  labels.txt <- file.path(data.dir, "kidney_bcell_labels.txt")
  cmd <- paste(Rscript, AllSteps, labels.txt)
  system(cmd)
  ## There should be summary bed files for labels and peaks:
  labels.bed.gz <- file.path(data.dir, "all_labels.bed.gz")
  bed.labels <- read.table(labels.bed.gz, skip=1)
  expect_equal(dim(bed.labels), c(25, 9))
  summary.bed.gz <- file.path(data.dir, "PeakSegJoint.summary.bed.gz")
  bed.summary <- read.table(summary.bed.gz, skip=1)
  expect_equal(ncol(bed.summary), 5)
  ## There should be predictions for 8 samples:
  pred.RData <- file.path(data.dir, "PeakSegJoint.predictions.RData")
  pred.objs <- load(pred.RData)
  expect_equal(nrow(all.peaks.mat), 8)
  starts <- unique(all.peaks.df$chromStart)
  expect_equal(length(starts), ncol(all.peaks.mat))
  test.errors.dir <- 
    file.path(data.dir, "PeakSegJoint-chunks", "figure-test-errors")
  test.errors.files <- dir(test.errors.dir)
  ## pipeline trained on 4 chunks generates 4 test error plots:
  expected.files <- paste0(1:4, ".png")
  expect_true(all(expected.files %in% test.errors.files))
  ## pipeline trained on 4 chunks generates test error summary:
  expect_true("figure-test-error-decreases.png" %in% test.errors.files)
  ## pipeline trained on 1 file does 4 fold CV:
  index.html <- file.path(test.errors.dir, "index.html")
  index.lines <- readLines(index.html)
  cv.line <- grep("cross-validation", index.lines, value=TRUE)
  expect_match(cv.line, "4 fold")
  ## check for no peak filtering when there are no labeled Inputs.
  expect_true("specific.error" %in% pred.objs)
  expect_null(specific.error$errors)
  ## Check for scatter viz.
  bars.vec <- Sys.glob(file.path(
    data.dir, "PeakSegJoint-predictions-viz", "*bars*.tsv"))
  expect_equal(length(bars.vec), 0)
  scatter.vec <- Sys.glob(file.path(
    data.dir, "PeakSegJoint-predictions-viz", "*scatter*.tsv"))
  expect_more_than(length(scatter.vec), 0)
})

