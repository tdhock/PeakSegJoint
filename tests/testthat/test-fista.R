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

  ## Compute the gradient of the smooth part of the loss
  ## function for a weight vector w \in \RR^d. This is
  ## \nabla L(w) = \sum_{i=1}^n \nabla l[y_i, f_w(X_i)]/n
  ## where i is the index of problems (n problems total),
  ## y_i is the target \in \RR^2
  ## X_i is the feature matrix \in \RR^{S_i \times d}
  ## S_i is the number of samples for problem i,
  ## d is the number of features (including intercept).
  calc.grad <- function(w){
    grad.list <- list()
    for(problem.i in seq_along(train.problems)){
      problem <- train.problems[[problem.i]]
      raw.features <- problem$features[, fit$train.feature.names]
      mean.mat <-
        matrix(fit$mean.vec, nrow(raw.features), ncol(raw.features), byrow=TRUE)
      sd.mat <-
        matrix(fit$sd.vec, nrow(raw.features), ncol(raw.features), byrow=TRUE)
      norm.features <- (raw.features-mean.mat)/sd.mat
      intercept.features <- cbind(1, norm.features)
      ## Predicted log penalty f_w(X) = 1_S^T X w for a weight vector
      ## w, a S-vector of ones 1_S, and a scaled feature matrix X (the
      ## first column is all ones).
      pred.log.lambda <- sum(intercept.features %*% w)
      left.term <- squared.hinge.deriv(pred.log.lambda - problem$target[1])
      right.term <- squared.hinge.deriv(problem$target[2] - pred.log.lambda)
      ## The gradient of one observation is \nabla l[y_i, f_w(X_i)] =
      ## X_i 1_S [ \phi'(f_w(X_i)-\underline y_i) - \phi'(\overline y_i-f_w(X_i))]
      ## where \phi' is the derivative of the squared hinge loss.
      one.grad <- colSums(intercept.features) * (left.term-right.term)
      grad.list[[problem.i]] <- one.grad
    }
    grad.mat <- do.call(rbind, grad.list)
    colMeans(grad.mat)
  }

  ## Compute the sub-differential optimality condition of weight
  ## vector w with gradient g. Each of these coefficients should be
  ## L1-regularized. Derivation of this expression in
  ## http://jmlr.org/proceedings/papers/v28/hocking13.html.
  dist2subdiff.opt <- function(w,g){
    ifelse(w==0,positive.part(abs(g)-regularization),
           ifelse(w<0,abs(-regularization+g),abs(regularization+g)))
  }

  ## Compute the optimality criterion for the weight vector w. Note
  ## that w[1] is the un-regularized intercept and w[-1] are the
  ## L1-regularized coefficients.
  get.crit <- function(w){
    after.grad <- calc.grad(w)
    w.dist <- dist2subdiff.opt(w[-1], after.grad[-1])
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
