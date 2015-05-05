context("fista")

data(H3K4me3.PGP.immune.4608)

test_that("coefficients satisfy subdifferential optimality", {
  train.problems <- H3K4me3.PGP.immune.4608

  positive.part <- function(x){
    ifelse(x<0, 0, x)
  }

  squared.hinge.deriv <- function(x){
    ifelse(x<1,2*(x-1),0)
  }  

  calc.grad <- function(x){
    grad.list <- list()
    for(problem.i in seq_along(train.problems)){
      problem <- train.problems[[problem.i]]
      raw.features <- problem$features[, fit$train.feature.names]
      mean.mat <- matrix(fit$mean.vec, nrow(raw.features), ncol(raw.features), byrow=TRUE)
      sd.mat <- matrix(fit$sd.vec, nrow(raw.features), ncol(raw.features), byrow=TRUE)
      norm.features <- (raw.features-mean.mat)/sd.mat
      intercept.features <- cbind(1, norm.features)
      pred.log.lambda <- sum(intercept.features %*% x)
      left.term <- squared.hinge.deriv(pred.log.lambda - problem$target[1])
      right.term <- squared.hinge.deriv(problem$target[2] - pred.log.lambda)
      full.grad <- colSums(intercept.features) * (left.term-right.term)
      grad.list[[problem.i]] <- full.grad
    }
    grad.mat <- do.call(rbind, grad.list)
    colMeans(grad.mat)
  }

  dist2subdiff.opt <- function(w,g){
    ifelse(w==0,positive.part(abs(g)-regularization),
           ifelse(w<0,abs(-regularization+g),abs(regularization+g)))
  }

  get.crit <- function(x){
    after.grad <- calc.grad(x)
    w.dist <- dist2subdiff.opt(x[-1], after.grad[-1])
    intercept.dist <- abs(after.grad[1])
    zero.at.optimum <- c(intercept.dist, w.dist)
    max(zero.at.optimum)
  }  

  thresh <- 1e-3
  fit <- IntervalRegressionProblems(train.problems, factor.regularization=1.2, threshold=thresh)
  computed.diff <- diff(log10(fit$regularization.vec))
  expected.diff <- log10(1.2)
  expect_true(all(abs(computed.diff-expected.diff) < 1e-10))
  for(regularization.i in seq_along(fit$regularization.vec)){
    regularization <- fit$regularization.vec[[regularization.i]]
    weight.vec <- fit$weight.mat[, regularization.i]
    crit <- get.crit(weight.vec)
    expect_that(crit, is_less_than(thresh))
  }
})
