context("binSum")

test_that("binSum 1bp with constant 1 profile", {
  profile <- data.frame(chromStart=0L,
                        chromEnd=10000L,
                        count=1L)
  bins <- binSum(profile,
                 bin.chromStart=0L,
                 bin.size=1L,
                 n.bins=2000L)
  expect_identical(bins$chromStart, 0:1999)
  expect_identical(bins$chromEnd, 1:2000)
  expect_identical(bins$count, rep(1L, 2000))
  expect_equal(bins$mean, rep(1, 2000))
})

test_that("binSum 2bp with constant 1 profile", {
  profile <- data.frame(chromStart=0L,
                        chromEnd=10000L,
                        count=1L)
  bins <- binSum(profile,
                 bin.chromStart=0L,
                 bin.size=2L,
                 n.bins=2000L)
  expect_identical(bins$chromStart, as.integer(seq(0, by=2, l=2000)))
  expect_identical(bins$chromEnd, as.integer(seq(2, by=2, l=2000)))
  expect_identical(bins$count, rep(2L, 2000))
  expect_equal(bins$mean, rep(1, 2000))
})

test_that("binSum 3bp with non-constant profile", {
  ## bins of size 3bp.
  ## -1-   -3-   -5-
  ##    -2-   -4-
  ## 123456789012345 base index.
  ## --2---
  ##       --1-
  ##           --0-------
  ## Coverage profile.
  profile <- data.frame(chromStart=as.integer(c(0, 6, 10)),
                        chromEnd=as.integer(c(6, 10, 10000)),
                        count=as.integer(c(2, 1, 0)))
  bins <- binSum(profile,
                 bin.chromStart=0L,
                 bin.size=3L,
                 n.bins=2000L)
  expect_identical(bins$chromStart, as.integer(seq(0, by=3, l=2000)))
  expect_identical(bins$chromEnd, as.integer(seq(3, by=3, l=2000)))
  expected.total <- rep(0L, 2000)
  expected.total[1:4] <- as.integer(c(6, 6, 3, 1))
  expect_identical(bins$count, expected.total)
  expected.mean <- expected.total/3
  expect_equal(bins$mean, expected.mean)
})

test_that("binSum 3bp with non-constant profile + 1000", {
  ## bins of size 3bp.
  ## -1-   -3-   -5-
  ##    -2-   -4-
  ## 123456789012345 base index + 1000.
  ## --2---
  ##       --1-
  ##           --0-------
  ## Coverage profile.
  profile <- data.frame(chromStart=as.integer(c(0, 1000, 1006, 1010)),
                        chromEnd=as.integer(c(1000, 1006, 1010, 10000)),
                        count=as.integer(c(0, 2, 1, 0)))
  bins <- binSum(profile,
                 bin.chromStart=1000L,
                 bin.size=3L,
                 n.bins=2000L)
  expect_identical(bins$chromStart, as.integer(seq(1000, by=3, l=2000)))
  expect_identical(bins$chromEnd, as.integer(seq(1003, by=3, l=2000)))
  expected.total <- rep(0L, 2000)
  expected.total[1:4] <- as.integer(c(6, 6, 3, 1))
  expect_identical(bins$count, expected.total)
  expected.mean <- expected.total/3
  expect_equal(bins$mean, expected.mean)
})

test_that("binSum returns short data.frame if coverage ends", {
  profile <- data.frame(chromStart=as.integer(c(0, 100)),
                        chromEnd=as.integer(c(100, 200)),
                        count=as.integer(c(10, 5)))
  bins <- binSum(profile,
                 bin.chromStart=0L,
                 bin.size=10L,
                 n.bins=10L)
  expect_equal(bins$mean, rep(10, 10))

  bins <- binSum(profile,
                 bin.chromStart=0L,
                 bin.size=10L,
                 n.bins=20L)
  expect_equal(bins$mean, rep(c(10, 5), each=10))

  bins <- binSum(profile,
                 bin.chromStart=0L,
                 bin.size=10L,
                 n.bins=30L)
  expect_equal(bins$mean, rep(c(10, 5), each=10))

  bins <- binSum(profile,
                 bin.chromStart=0L,
                 bin.size=30L,
                 n.bins=10L)
  expect_equal(bins$mean,
               c(10, 10, 10, 200/30, 5, 5, 100/30))
  expect_identical(bins$chromStart,
                   as.integer(c(0, 30, 60, 90, 120, 150, 180)))
  expect_identical(bins$chromEnd,
                   as.integer(c(30, 60, 90, 120, 150, 180, 210)))
})

test_that("chromEnd <= chromStart is an error", {
  bad <- data.frame(chromStart=as.integer(c(0, 10)),
                    chromEnd=as.integer(c(10, 5)),
                    count=0L)
  expect_error({
    binSum(bad, 0L, 30L, 10L)
  }, "chromStart not less than chromEnd")
})

test_that("chromEnd[i-1] != chromStart[i] is an error", {
  bad <- data.frame(chromStart=as.integer(c(0, 10)),
                    chromEnd=as.integer(c(15, 20)),
                    count=0L)
  expect_error({
    binSum(bad, 0L, 30L, 10L)
  }, "chromStart[i] != chromEnd[i-1]", fixed=TRUE)
})

