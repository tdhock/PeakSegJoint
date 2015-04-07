context("multiSampleSeg")

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
  peak <- multiSampleSegHeuristic(profiles)
  expect_equal(peak$chromStart, 10)
  expect_equal(peak$chromEnd, 30)
})

test_that("chromEnd <= chromStart is an error", {
  bad <- data.frame(chromStart=as.integer(c(0, 100)),
                    chromEnd=as.integer(c(100, 50)),
                    count=0L,
                    sample.id="foo")
  expect_error({
    multiSampleSegHeuristic(bad)
  }, "chromStart not less than chromEnd")
})

test_that("chromEnd[i-1] != chromStart[i] is an error", {
  bad <- data.frame(chromStart=as.integer(c(0, 100)),
                    chromEnd=as.integer(c(150, 200)),
                    count=0L,
                    sample.id="foo")
  expect_error({
    multiSampleSegHeuristic(bad)
  }, "chromStart[i] != chromEnd[i-1]", fixed=TRUE)
})

data(H3K4me3.TDH.immune.chunk12.cluster4)
many <- H3K4me3.TDH.immune.chunk12.cluster4

test_that("optimal results for real data set", {
  ## These numbers obtained from running the optimal, slow DP algo.
  optimal <- data.frame(chromStart=27998215L, chromEnd=27999159L)
  ## optimal.seconds <- system.time({
  ##   optimal <- multiSampleSegOptimal(many)
  ## })[["elapsed"]]

  ## For some reason bin.factor=100 gives optimal results for these
  ## data.
  peak <- multiSampleSegHeuristic(many, 100L)

  expect_identical(peak, optimal)
})

test_that("reasonable results for fastest bin.factor=2", {
  peak <- multiSampleSegHeuristic(many, 2L)

  ## These numbers come from running the heuristic algorithm after
  ## extensive interactive testing, so I am pretty sure they are
  ## correct (even though they are not optimal).
  observed <- data.frame(chromStart=27998215L, chromEnd=27999172L)
  
  expect_identical(peak, observed)
})
