context("trivial")

count <-
  as.integer(c(rep(c(0, 1), each=5),
               rep(c(10, 11), each=10),
               rep(c(0, 1), each=5)))

test_that("with 1 sample we refine peaks" , {
  chromEnd <- seq_along(count)
  profiles <-
    data.frame(chromStart=chromEnd-1L,
               chromEnd,
               sample.id="sample1",
               count)
  fit <- PeakSegJointHeuristic(profiles)
  converted <- ConvertModelList(fit)
  peak <- converted$peaks
  expect_equal(peak$chromStart, 10)
  expect_equal(peak$chromEnd, 30)
})

test_that("chromEnd <= chromStart is an error", {
  bad <- data.frame(chromStart=as.integer(c(0, 100)),
                    chromEnd=as.integer(c(100, 50)),
                    count=0L,
                    sample.id="foo")
  expect_error({
    PeakSegJointHeuristic(bad)
  }, "chromStart not less than chromEnd")
})

test_that("chromEnd[i-1] != chromStart[i] is an error", {
  bad <- data.frame(chromStart=as.integer(c(0, 100)),
                    chromEnd=as.integer(c(150, 200)),
                    count=0L,
                    sample.id="foo")
  expect_error({
    PeakSegJointHeuristic(bad)
  }, "chromStart[i] != chromEnd[i-1]", fixed=TRUE)
})

