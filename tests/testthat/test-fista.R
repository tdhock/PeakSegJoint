context("fista")

library(PeakSegJoint)
data(H3K4me3.PGP.immune.4608)
train.problems <- H3K4me3.PGP.immune.4608

positive.part <- function(x){
  ifelse(x<0, 0, x)
}
squared.hinge.deriv <- function(x){
  ifelse(x<1,2*(x-1),0)
}  
calc.grad <- function(x){
  linear.predictor <- as.numeric(features %*% x)
  left.term <- squared.hinge.deriv(linear.predictor-targets[,1])
  right.term <- squared.hinge.deriv(targets[,2]-linear.predictor)
  full.grad <- features * (left.term-right.term)
  colSums(full.grad)/nrow(full.grad)
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

test_that("coefficients satisfy subdifferential optimality", {
  fit <- IntervalRegressionProblems(train.problems, factor.regularization=1.2, threshold=1e-3)
  computed.diff <- diff(log10(fit$regularization.vec))
  expected.diff <- log10(1.2)
  expect_true(all(abs(computed.diff-expected.diff) < 1e-10))
  for(regularization.i in seq_along(fit$regularization.vec)){
    regularization <- fit$regularization.vec[[regularization.i]]
    weight.vec <- fit$weight.mat[, regularization.i]
    ## TODO: compute the optimality condition of the computed weight vector.
  }
})
