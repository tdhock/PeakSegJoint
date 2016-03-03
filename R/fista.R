IntervalRegressionMatrixCV <- function
### Use cross-validation to estimate the optimal regularization, by
### picking the value that minimize the number of incorrectly
### predicted target intervals.
(feature.mat,
### Numeric feature matrix.
 target.mat,
### Numeric target matrix.
 n.folds=5,
### Number of cross-validation folds.
 fold.vec=sample(rep(1:n.folds, l=nrow(feature.mat)))
### Integer vector of fold id numbers.
 ){
  stopifnot(is.numeric(feature.mat))
  stopifnot(is.matrix(feature.mat))
  n.observations <- nrow(feature.mat)
  stopifnot(is.numeric(target.mat))
  stopifnot(is.matrix(target.mat))
  stopifnot(nrow(target.mat) == n.observations)
  stopifnot(ncol(target.mat) == 2)
  stopifnot(is.integer(fold.vec))
  stopifnot(length(fold.vec) == n.observations)
  best.reg.list <- list()
  for(validation.fold in unique(fold.vec)){
    ##print(validation.fold)
    is.validation <- fold.vec == validation.fold
    is.train <- !is.validation
    train.features <- feature.mat[is.train, ]
    train.targets <- target.mat[is.train, ]
    fit <- IntervalRegressionMatrixPath(
      train.features, train.targets, max.iterations=1e3, verbose=0)
    validation.features <- feature.mat[is.validation, ]
    pred.log.lambda <- fit$predict(validation.features)
    validation.targets <- target.mat[is.validation, ]
    too.small <- pred.log.lambda < validation.targets[, 1]
    too.big <- validation.targets[, 2] < pred.log.lambda
    is.error <- too.small | too.big
    error.vec <- colSums(is.error)
    best.reg.vec <- fit$regularization.vec[error.vec == min(error.vec)]
    best.reg.list[[paste(validation.fold)]] <- max(best.reg.vec) #simplest model
  }
  mean.reg <- mean(unlist(best.reg.list))
  IntervalRegressionMatrixPath(
    feature.mat, target.mat,
    initial.regularization=mean.reg,
    factor.regularization=NULL,
    verbose=0)
}

IntervalRegressionProblemsCV <- function
### Use cross-validation to estimate the optimal regularization, by
### picking the value that minimize the number of incorrectly
### predicted target intervals.
(problem.list,
### List of lists with features and targets.
 n.folds=5,
### Number of cross-validation folds.
 fold.vec=sample(rep(1:n.folds, l=length(problem.list)))
### Integer vector of fold id numbers.
 ){
  stopifnot(is.integer(fold.vec))
  stopifnot(length(fold.vec) == length(problem.list))
  best.reg.list <- list()
  for(validation.fold in unique(fold.vec)){
    ##print(validation.fold)
    is.validation <- fold.vec == validation.fold
    is.train <- !is.validation
    train.list <- problem.list[is.train]
    fit <- IntervalRegressionProblems(
      train.list, max.iterations=1e3, verbose=0)
    validation.list <- problem.list[is.validation]
    error.mat <- matrix(
      NA, length(validation.list), length(fit$regularization.vec))
    for(i in seq_along(validation.list)){
      info <- validation.list[[i]]
      pred.log.lambda <- fit$predict(info$features)
      too.small <- pred.log.lambda < info$target[1]
      too.big <- info$target[2] < pred.log.lambda
      is.error <- too.small | too.big
      error.mat[i, ] <- is.error
    }
    error.vec <- colSums(error.mat)
    best.reg.vec <- fit$regularization.vec[error.vec == min(error.vec)]
    best.reg.list[[paste(validation.fold)]] <- max(best.reg.vec) #simplest model
  }
  mean.reg <- mean(unlist(best.reg.list))
  IntervalRegressionProblems(
    problem.list,
    initial.regularization=mean.reg,
    factor.regularization=NULL,
    verbose=0)
}

IntervalRegressionMatrixPath <- function
### Interval regression using the squared hinge loss, for a path of
### regularization parameters.
(feature.mat,
### Numeric feature matrix.
 target.mat,
### Numeric target matrix.
 initial.regularization=0.001,
### Initial regularization parameter.
 factor.regularization=1.5,
### Increase regularization by this factor after finding an optimal
### solution. Or NULL to compute just one model
### (initial.regularization).
 verbose=1,
### Print messages if >= 1.
 ...
### Other parameters to pass to IntervalRegressionMatrix.
 ){
  stopifnot(is.numeric(feature.mat))
  stopifnot(is.matrix(feature.mat))
  n.observations <- nrow(feature.mat)
  stopifnot(is.numeric(target.mat))
  stopifnot(is.matrix(target.mat))
  stopifnot(nrow(target.mat) == n.observations)
  stopifnot(ncol(target.mat) == 2)
  stopifnot(is.numeric(initial.regularization))
  stopifnot(length(initial.regularization)==1)

  is.trivial.target <- apply(!is.finite(target.mat), 1, all)
  nontrivial.features <- feature.mat[!is.trivial.target, ]
  nontrivial.targets <- target.mat[!is.trivial.target, ]
  is.finite.feature <- apply(is.finite(nontrivial.features), 2, all)
  finite.features <- nontrivial.features[, is.finite.feature, drop=FALSE]
  all.mean.vec <- colMeans(finite.features)
  all.sd.vec <- apply(finite.features, 2, sd)
  is.invariant <- all.sd.vec == 0
  train.feature.i <- which(!is.invariant)
  train.feature.names <- colnames(finite.features)[train.feature.i]
  mean.vec <- all.mean.vec[train.feature.names]
  sd.vec <- all.sd.vec[train.feature.names]
  invariant.features <- finite.features[, train.feature.names, drop=FALSE]
  mean.mat <- matrix(
    mean.vec, nrow(invariant.features), ncol(invariant.features), byrow=TRUE)
  sd.mat <- matrix(
    sd.vec, nrow(invariant.features), ncol(invariant.features), byrow=TRUE)
  norm.features <- (invariant.features-mean.mat)/sd.mat
  intercept.features <- cbind("(Intercept)"=1, norm.features)
  apply(intercept.features, 2, mean)
  apply(intercept.features, 2, sd)

  regularization <- initial.regularization
  n.features <- ncol(intercept.features)
  param.vec <- rep(0, n.features)
  n.nonzero <- n.features

  param.vec.list <- list()
  regularization.vec.list <- list()
  while(n.nonzero > 1){
    param.vec <-
      IntervalRegressionMatrix(intercept.features, nontrivial.targets,
                               param.vec,
                               regularization,
                               verbose=verbose,
                               ...)
    n.zero <- sum(param.vec == 0)
    n.nonzero <- sum(param.vec != 0)
    l1.norm <- sum(abs(param.vec[-1]))
    if(verbose >= 1){
      cat(sprintf("regularization=%8.4f L1norm=%8.4f zeros=%d\n",
                  regularization, l1.norm, n.zero))
    }
    param.vec.list[[paste(regularization)]] <- param.vec
    regularization.vec.list[[paste(regularization)]] <- regularization
    if(is.null(factor.regularization)){
      n.nonzero <- 1 #stops while loop.
    }else{
      stopifnot(is.numeric(factor.regularization))
      stopifnot(length(factor.regularization)==1)
      regularization <- regularization * factor.regularization
    }
  }
  param.mat <- do.call(cbind, param.vec.list)
  if(verbose >= 1){
    cat(paste0("Done computing parameter matrix (",
               nrow(param.mat), " features x ",
               ncol(param.mat), " regularization parameters)\n"))
  }
  feature.not.used <- apply(param.mat[-1, ,drop=FALSE] == 0, 1, all)
  pred.feature.names <- train.feature.names[!feature.not.used]
  pred.param.mat <-
    param.mat[c("(Intercept)", pred.feature.names),,drop=FALSE]
  list(param.mat=param.mat,
       regularization.vec=do.call(c, regularization.vec.list),
       mean.vec=mean.vec,
       sd.vec=sd.vec,
       train.feature.names=train.feature.names,
       pred.feature.names=pred.feature.names,
       pred.param.mat=pred.param.mat,
       predict=function(mat){
         stopifnot(is.matrix(mat))
         stopifnot(is.numeric(mat))
         stopifnot(pred.feature.names %in% colnames(mat))
         raw.mat <- mat[, pred.feature.names, drop=FALSE]
         raw.mat[!is.finite(raw.mat)] <- 0 
         mean.mat <- matrix(
           mean.vec[pred.feature.names],
           nrow(raw.mat), ncol(raw.mat), byrow=TRUE)
         sd.mat <- matrix(
           sd.vec[pred.feature.names],
           nrow(raw.mat), ncol(raw.mat), byrow=TRUE)
         norm.mat <- (raw.mat-mean.mat)/sd.mat
         intercept.mat <- cbind("(Intercept)"=1, norm.mat)
         intercept.mat %*% pred.param.mat
       })
### List representing fit model. You can do
### fit$predict(feature.matrix) to get a predicted log penalty
### value. The mean.vec and sd.vec were used for scaling the training
### data matrices. The param.mat is the n.features * n.regularization
### numeric matrix of optimal coefficients.
}

IntervalRegressionProblems <- structure(function
### Compute a sequence of interval regression models for increasingly
### more L1 regularization, until we get to a regularization parameter
### that gives an optimal weight vector of zero. The problem is w* =
### argmin_w L(w) + regularization * ||w||_1 where L(w) is the mean
### squared hinge loss and ||w||_1 is the L1-norm of the non-intercept
### coefficients. We first scale the input features and then
### repeatedly call IntervalRegressionMatrix, using warm restarts.
(problem.list,
### List of problems with features (numeric matrix), target (numeric
### vector of length 2), and optionally weight (numeric length 1).
 initial.regularization=0.001,
### Initial regularization parameter.
 factor.regularization=1.5,
### Increase regularization by this factor after finding an optimal
### solution. Or NULL to compute just one model
### (initial.regularization).
 verbose=1,
### Print messages if >= 1.
 ...
### Other parameters to pass to IntervalRegressionMatrix.
 ){
  stopifnot(is.list(problem.list))
  stopifnot(is.numeric(initial.regularization))
  stopifnot(length(initial.regularization) == 1)
  stopifnot(initial.regularization > 0)
  if(!is.null(factor.regularization)){
    stopifnot(is.numeric(factor.regularization))
    stopifnot(length(factor.regularization) == 1)
    stopifnot(factor.regularization > 1)
  }
  stopifnot(length(problem.list) > 0)
  n.input.features <- ncol(problem.list[[1]]$features)

  featureSum.list <- list()
  targets.list <- list()
  for(problem.i in seq_along(problem.list)){
    problem <- problem.list[[problem.i]]
    stopifnot(is.matrix(problem$features))
    stopifnot(is.numeric(problem$features))
    stopifnot(ncol(problem$features) == n.input.features)
    stopifnot(is.numeric(problem$target))
    if(length(problem$target) != 2){
      print(problem)
      stop("target should be numeric vector of length 2")
    }
    problem <- problem.list[[problem.i]]
    featureSum.list[[problem.i]] <- colSums(problem$features)
    targets.list[[problem.i]] <- problem$target
  }
  all.targets <- do.call(rbind, targets.list)
  is.trivial.target <- apply(!is.finite(all.targets), 1, all)
  nontrivial.i.vec <- which(!is.trivial.target)

  all.featureSum <- do.call(rbind, featureSum.list)
  is.finite.feature <- apply(is.finite(all.featureSum), 2, all)

  features.list <- list()
  for(problem.i in nontrivial.i.vec){
    problem <- problem.list[[problem.i]]
    features.list[[paste(problem.i)]] <-
      problem$features[, is.finite.feature, drop=FALSE]
  }
  finite.features <- do.call(rbind, features.list)
  all.mean.vec <- colMeans(finite.features)
  all.sd.vec <- apply(finite.features, 2, sd)
  is.invariant <- all.sd.vec == 0
  train.feature.i <- which(!is.invariant)
  train.feature.names <- colnames(finite.features)[train.feature.i]
  mean.vec <- all.mean.vec[train.feature.names]
  sd.vec <- all.sd.vec[train.feature.names]

  norm.featureSum.list <- list()
  for(problem.i in nontrivial.i.vec){
    problem <- problem.list[[problem.i]]
    raw.features.mat <- problem$features[, train.feature.names, drop=FALSE]
    train.mean.mat <-
      matrix(mean.vec, nrow(raw.features.mat), ncol(raw.features.mat),
             byrow=TRUE)
    train.sd.mat <-
      matrix(sd.vec, nrow(raw.features.mat), ncol(raw.features.mat),
             byrow=TRUE)
    norm.mat <- (raw.features.mat-train.mean.mat)/train.sd.mat
    intercept.mat <- cbind("(Intercept)"=1, norm.mat)
    norm.featureSum.list[[paste(problem.i)]] <- colSums(intercept.mat)
  }
  norm.featureSum.mat <- do.call(rbind, norm.featureSum.list)
  targets.mat <- do.call(rbind, targets.list[nontrivial.i.vec])
  rownames(norm.featureSum.mat) <- rownames(targets.mat) <-
    names(problem.list)[nontrivial.i.vec]

  ## Do we need to take the feature sd again and divide by that for
  ## optimization purposes?
  apply(norm.featureSum.mat, 2, sd)

  regularization <- initial.regularization
  n.features <- ncol(norm.featureSum.mat)
  param.vec <- rep(0, n.features)
  n.nonzero <- n.features

  weight.vec <- sapply(problem.list[nontrivial.i.vec], function(p){
    if(is.numeric(p$weight)){
      p$weight[1]
    }else{
      1
    }
  })
  total.weight <- sum(weight.vec)
  norm.weight.vec <- weight.vec/total.weight * length(weight.vec)
  
  param.vec.list <- list()
  regularization.vec.list <- list()
  while(n.nonzero > 1){
    param.vec <-
      IntervalRegressionMatrix(norm.featureSum.mat, targets.mat,
                               param.vec,
                               regularization,
                               verbose=verbose,
                               weight.vec=norm.weight.vec,
                               ...)
    n.zero <- sum(param.vec == 0)
    n.nonzero <- sum(param.vec != 0)
    l1.norm <- sum(abs(param.vec[-1]))
    if(verbose >= 1){
      cat(sprintf("regularization=%8.4f L1norm=%8.4f zeros=%d\n",
                  regularization, l1.norm, n.zero))
    }
    param.vec.list[[paste(regularization)]] <- param.vec
    regularization.vec.list[[paste(regularization)]] <- regularization
    if(is.null(factor.regularization)){
      n.nonzero <- 1 #stops while loop.
    }else{
      regularization <- regularization * factor.regularization
    }
  }
  param.mat <- do.call(cbind, param.vec.list)
  if(verbose >= 1){
    cat(paste0("Done computing parameter matrix (",
               nrow(param.mat), " features x ",
               ncol(param.mat), " regularization parameters)\n"))
  }
  feature.not.used <- apply(param.mat[-1, ,drop=FALSE] == 0, 1, all)
  pred.feature.names <- train.feature.names[!feature.not.used]
  pred.param.mat <-
    param.mat[c("(Intercept)", pred.feature.names),,drop=FALSE]
  list(param.mat=param.mat,
       regularization.vec=do.call(c, regularization.vec.list),
       mean.vec=mean.vec,
       sd.vec=sd.vec,
       train.feature.names=train.feature.names,
       pred.feature.names=pred.feature.names,
       pred.param.mat=pred.param.mat,
       predict=function(mat){
         stopifnot(is.matrix(mat))
         stopifnot(is.numeric(mat))
         stopifnot(pred.feature.names %in% colnames(mat))
         raw.mat <- mat[, pred.feature.names, drop=FALSE]
         raw.mat[!is.finite(raw.mat)] <- 0 
         mean.mat <- matrix(
           mean.vec[pred.feature.names],
           nrow(raw.mat), ncol(raw.mat), byrow=TRUE)
         sd.mat <- matrix(
           sd.vec[pred.feature.names],
           nrow(raw.mat), ncol(raw.mat), byrow=TRUE)
         norm.mat <- (raw.mat-mean.mat)/sd.mat
         intercept.mat <- cbind("(Intercept)"=1, norm.mat)
         colSums(intercept.mat %*% pred.param.mat)
       })
### List representing fit model. You can do
### fit$predict(feature.matrix) to get a predicted log penalty
### value. The mean.vec and sd.vec were used for scaling the training
### data matrices. The param.mat is the n.features * n.regularization
### numeric matrix of optimal coefficients.
}, ex=function(){
  library(PeakSegJoint)
  data(H3K4me3.PGP.immune.4608)
  chrom.vec <- sub(":.*", "", names(H3K4me3.PGP.immune.4608))
  table(chrom.vec)
  train.chroms <- c("chr1", "chr9")
  sets <-
    list(train=chrom.vec %in% train.chroms,
         validation=! chrom.vec %in% train.chroms)
  train.problems <- H3K4me3.PGP.immune.4608[sets$train]
  fit <- IntervalRegressionProblems(train.problems)
  set.error.list <- list()
  for(set.name in names(sets)){
    in.set <- sets[[set.name]]
    problem.list <- H3K4me3.PGP.immune.4608[in.set]
    error.mat.list <- list()
    for(problem.name in names(problem.list)){
      problem <- problem.list[[problem.name]]
      pred.log.lambda <- fit$predict(problem$features)
      too.hi <- problem$target[2] < pred.log.lambda
      too.lo <- pred.log.lambda < problem$target[1]
      is.error <- too.hi | too.lo
      error.mat.list[[problem.name]] <- is.error
    }
    error.mat <- do.call(rbind, error.mat.list)
    percent.error <- colMeans(error.mat) * 100
    set.error.list[[set.name]] <-
      data.frame(set.name,
                 regularization=fit$regularization,
                 percent.error)
  }

  ## Coefficients of the best models.
  min.validation <- 
    subset(set.error.list$validation,
           percent.error==min(percent.error))
  best.models <- fit$param.mat[, rownames(min.validation)]
  best.nonzero <- best.models[apply(best.models!=0, 1, any), ]
  print(best.nonzero)

  ## Plot train/validation error curves.
  set.error <- do.call(rbind, set.error.list)
  library(ggplot2)
  ggplot()+
    geom_point(aes(-log10(regularization), percent.error),
               data=min.validation)+
    geom_line(aes(-log10(regularization), percent.error,
                  group=set.name, linetype=set.name),
              data=set.error)

  ## Fit model with the chosen regularization to the full
  ## train+validation set.
  chosen.regularization <- max(min.validation$regularization)
  full.fit <-
    IntervalRegressionProblems(H3K4me3.PGP.immune.4608,
                               initial.regularization=chosen.regularization,
                               factor.regularization=NULL)
  print(full.fit$param.mat)
  
})

IntervalRegressionMatrix <- function
### Solve the squared hinge loss interval regression problem for one
### regularization parameter: w* = argmin_w L(w) + regularization *
### ||w||_1 where L(w) is the average squared hinge loss with respect
### to the targets, and ||w||_1 is the L1-norm of the weight vector
### (excluding the first element, which is the un-regularized
### intercept or bias term).
(features,
### Scaled numeric feature matrix (problems x features). The first
### column/feature should be all ones and will not be regularized.
 targets,
### Numeric target matrix (problems x 2).
 initial.param.vec,
### initial guess for weight vector (features).
 regularization,
### Degree of L1-regularization.
 threshold=1e-2,
### When the stopping criterion gets below this threshold, the
### algorithm stops and declares the solution as optimal.
 max.iterations=1e5,
### Error if the algorithm has not found an optimal solution after
### this many iterations.
 weight.vec=NULL,
### A numeric vector of weights for each training example.
 Lipschitz=NULL,
### A numeric scalar or NULL, which means to compute Lipschitz as the
### mean of the squared L2-norms of the rows of the feature matrix.
 verbose=2
### Cat messages: for restarts and at the end if >= 1, and for every
### iteration if >= 2.
 ){
  stopifnot(is.matrix(features))
  stopifnot(is.numeric(features))
  n.features <- ncol(features)
  n.problems <- nrow(features)

  stopifnot(is.matrix(targets))
  stopifnot(nrow(targets) == n.problems)
  stopifnot(ncol(targets) == 2)

  if(is.null(weight.vec)){
    weight.vec <- rep(1, n.problems)
  }
  stopifnot(is.numeric(weight.vec))
  stopifnot(length(weight.vec) == n.problems)

  if(is.null(Lipschitz)){
    Lipschitz <- mean(rowSums(features * features) * weight.vec)
  }
  stopifnot(is.numeric(Lipschitz))
  stopifnot(length(Lipschitz) == 1)

  stopifnot(is.numeric(max.iterations))
  stopifnot(length(max.iterations) == 1)

  stopifnot(is.numeric(threshold))
  stopifnot(length(threshold) == 1)

  stopifnot(is.numeric(initial.param.vec))
  stopifnot(length(initial.param.vec) == n.features)

  ## Return 0 for a negative number and the same value otherwise.
  positive.part <- function(x){
    ifelse(x<0, 0, x)
  }
  squared.hinge <- function(x){
    ifelse(x<1,(x-1)^2,0)
  }
  squared.hinge.deriv <- function(x){
    ifelse(x<1,2*(x-1),0)
  }  
  calc.loss <- function(x){
    linear.predictor <- as.numeric(features %*% x)
    left.term <- squared.hinge(linear.predictor-targets[,1])
    right.term <- squared.hinge(targets[,2]-linear.predictor)
    both.terms <- left.term+right.term
    weighted.loss.vec <- both.terms * weight.vec
    mean(weighted.loss.vec)
  }
  calc.grad <- function(x){
    linear.predictor <- as.numeric(features %*% x)
    left.term <- squared.hinge.deriv(linear.predictor-targets[,1])
    right.term <- squared.hinge.deriv(targets[,2]-linear.predictor)
    full.grad <- features * (left.term-right.term) * weight.vec
    colSums(full.grad)/nrow(full.grad)
  }    
  calc.penalty <- function(x){
    regularization * sum(abs(x[-1]))
  }
  calc.cost <- function(x){
    calc.loss(x) + calc.penalty(x)
  }
  soft.threshold <- function(x,thresh){
    ifelse(abs(x) < thresh, 0, x-thresh*sign(x))
  }
  ## do not threshold the intercept.
  prox <- function(x,thresh){
    x[-1] <- soft.threshold(x[-1],thresh)
    x
  }
  ## p_L from the fista paper.
  pL <- function(x,L){
    grad <- calc.grad(x)
    prox(x - grad/L, regularization/L)
  }
  dist2subdiff.opt <- function(w,g){
    ifelse(w==0,positive.part(abs(g)-regularization),
           ifelse(w<0,abs(-regularization+g),abs(regularization+g)))
  }

  iterate.count <- 1
  stopping.crit <- threshold
  last.iterate <- this.iterate <- y <- initial.param.vec
  this.t <- 1
  while({
    ##browser(expr=is.na(stopping.crit))
    ##str(stopping.crit)
    stopping.crit >= threshold
  }){
    ## here we implement the FISTA method with constant step size, as
    ## described by in the Beck and Tebolle paper.
    last.iterate <- this.iterate
    this.iterate <- pL(y, Lipschitz)
    last.t <- this.t
    this.t <- (1+sqrt(1+4*last.t^2))/2
    y <- this.iterate + (last.t - 1)/this.t*(this.iterate-last.iterate)
    ## here we calculate the subgradient optimality condition, which
    ## requires 1 more gradient evaluation per iteration.
    after.grad <- calc.grad(this.iterate)
    w.dist <- dist2subdiff.opt(this.iterate[-1],after.grad[-1])
    zero.at.optimum <- c(abs(after.grad[1]),w.dist)
    stopping.crit <- max(zero.at.optimum)

    if(verbose >= 2){
      cost <- calc.cost(this.iterate)
      cat(sprintf("%10d cost %10f crit %10.7f\n",
                  iterate.count,
                  cost,
                  stopping.crit))
    }
    iterate.count <- iterate.count + 1
    if(iterate.count > max.iterations){
      Lipschitz <- Lipschitz * 1.5
      iterate.count <- 1
      if(verbose >= 1){
        cat(max.iterations, "iterations, increasing Lipschitz.",
            "crit =", stopping.crit, "\n")
      }
    }
    if(any(!is.finite(this.iterate))){
      if(verbose >= 1){
        cat("infinite parameter, restarting with bigger Lipschitz.\n")
      }
      iterate.count <- 1
      stopping.crit <- threshold
      last.iterate <- this.iterate <- y <- initial.param.vec
      this.t <- 1
      Lipschitz <- Lipschitz * 1.5
    }
  }
  if(verbose >= 1){
    cat("solution with crit =", stopping.crit, "\n")
  }
  this.iterate
### Numeric vector of scaled weights w of the affine function f_w(X) =
### X %*% w for a scaled feature matrix X with the first row entirely
### ones.
}
