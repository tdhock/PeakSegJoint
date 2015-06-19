library(testthat)
context("sparse data without counts equal to 0")

library(PeakSegJoint)

data(peak1.infeasible)
data.list <-
  list(with.zero=peak1.infeasible,
       without.zero=subset(peak1.infeasible, 0 < count))
test_that("PeakSegJointHeuristic loss same with or without zeros", {
  for(bp in 2:7){
    loss.list <- list()
    for(data.name in names(data.list)){
      L <- data.list[[data.name]]
      fit <- PeakSegJointHeuristic(L, bp)
      loss.list[[data.name]] <- sapply(fit$models, "[[", "loss")
    }
    with(loss.list, expect_equal(with.zero, without.zero))
  }
})

test_that("PeakSegJointSeveral loss same with or without zeros", {
  loss.list <- list()
  for(data.name in names(data.list)){
    L <- data.list[[data.name]]
    fit <- PeakSegJointSeveral(L)
    loss.vec <- sapply(fit$models, "[[", "loss")
    expect_true(all(is.finite(loss.vec)))
    loss.list[[data.name]] <- loss.vec
  }
  with(loss.list, expect_equal(with.zero, without.zero))
})

