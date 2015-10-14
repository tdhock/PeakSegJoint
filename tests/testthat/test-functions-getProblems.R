context("getProblems")

regions <-
  data.frame(chromStart=c(10, 50),
             chromEnd=c(110, 150))

test_that("getProblems ends after last region", {
  problems <- 
    getProblems("chrFake",
                min(regions$chromStart),
                max(regions$chromEnd),
                70)
  expect_true(150 < max(problems$problemEnd))
})

last.base.on.chrFake <- 150

test_that("getProblems ends at last base on chrom", {
  problems <- 
    getProblems("chrFake",
                min(regions$chromStart),
                max(regions$chromEnd),
                70,
                chrom.size=last.base.on.chrFake)
  expect_equal(max(problems$problemEnd), 150)
})
