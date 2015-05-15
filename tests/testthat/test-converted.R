context("ConvertModelList")

test_that("segments for model with 0 peaks", {
  data(H3K4me3.TDH.other.chunk8)
  fit <- PeakSegJointHeuristic(H3K4me3.TDH.other.chunk8, 3)
  converted <- ConvertModelList(fit)
  expect_true(0 %in% converted$segments$peaks)
})
