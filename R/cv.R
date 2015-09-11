tv.curves <- function
### Fit regularized interval regression models and compute train and
### validation error curves.
(chunk.problem.list,
### chunk.problem.list[[chunk.name]][[problem.name]] is a list with
### peaks, segments, loss, modelSelection, features, target.
 regions.list,
### regions.list[[chunk.name]] is a data.frame of annotated region
### labels for one chunk.
 n.folds=4,
### Number of cross-validation folds to use to estimate the
### regularization parameter.
 seed=1
### Random seed to use to choose cross-validation fold ID numbers for
### each chunk.
 ){
  ## First we need to filter only the labeled problems, since
  ## IntervalRegressionProblems stops with an error if its inputs do
  ## not have targets.
  labeled.chunk.problem.list <- list()
  for(chunk.name in names(chunk.problem.list)){
    one.chunk <- chunk.problem.list[[chunk.name]]
    for(problem.name in names(one.chunk)){
      prob <- one.chunk[[problem.name]]
      if(is.numeric(prob$target)){
        labeled.chunk.problem.list[[chunk.name]][[problem.name]] <- prob
      }
    }
  }
  train.validation <- names(chunk.problem.list)
  stopifnot(is.character(train.validation))
  stopifnot(2 <= length(train.validation))
  if(length(train.validation) < n.folds){
    n.folds <- length(train.validation)
  }
  ## If the chunk numbers are not sorted, then two permutations of the
  ## same train.validation vector will end up with two different
  ## models.
  train.validation <- sort(train.validation)
  set.seed(seed)
  fold.id <- sample(rep(1:n.folds, l=length(train.validation)))

  one.fold <- function(validation.fold){
    is.validation <- fold.id == validation.fold
    sets <- list(validation=train.validation[is.validation],
                 train=train.validation[!is.validation])
    t.chunk.data <- do.call(c, labeled.chunk.problem.list[sets$train])
    fit <-
      IntervalRegressionProblems(t.chunk.data,
                                 initial.regularization=0.005,
                                 factor.regularization=1.1,
                                 verbose=0)
    error.by.tv <- list()
    for(tv in names(sets)){
      chunk.name.vec <- sets[[tv]]
      tv.chunk.data <- chunk.problem.list[chunk.name.vec]
      tv.regions <- regions.list[chunk.name.vec]
      error.results <- error.metrics(tv.chunk.data, tv.regions, fit)
      error.df <- error.results$metrics
      error.by.tv[[tv]] <- data.frame(tv, validation.fold, error.df)
    }#tv
    do.call(rbind, error.by.tv)
  }#validation.fold

  error.by.fold <- mclapply.or.stop(1:n.folds, one.fold)
  do.call(rbind, error.by.fold)
### Data.frame of train/validation error curves.
}

best.on.validation <- function
### For each validation fold, pick the most regularized (least
### complex) model that has minimum error.
(err.curves
### Result data.frame of train/validation error curves from tv.curves.
 ){
  v.err <-
    subset(err.curves,
           metric.name=="incorrect.regions" &
             tv=="validation")
  v.list <- split(v.err, v.err$validation.fold)
  picked.list <- list()
  for(validation.fold in names(v.list)){
    v <- v.list[[validation.fold]]
    min.err.params <- subset(v, metric.value == min(metric.value))
    least.complex <- which.max(min.err.params$regularization)
    picked.list[[validation.fold]] <- min.err.params[least.complex, ]
  }
  do.call(rbind, picked.list)
### Data.frame with one row for each validation fold.
}

error.metrics <- function
### Use a model to predict peaks, and evaluate several error metrics:
### incorrect.targets, incorrect.regions, false.positives,
### false.negatives.
(problems.by.chunk,
### List[[chunk.name]][[problem.name]], a list which must have
### features, modelSelection, peaks (and if it also has target, it
### will be used to compute number of incorrect targets).
 regions.by.chunk,
### List[[chunk.name]], a data.frame of annotated region labels for
### each chunk on which the error metrics will be evaluated.
 fit
### Model fit list from IntervalRegressionProblems.
 ){
  chunk.name.vec <- names(regions.by.chunk)
  error.vec.list <- list()
  regions.list <- list()
  outside.target.list <- list()
  result.list <- list()
  for(chunk.name in chunk.name.vec){
    chunk.regions <- regions.by.chunk[[chunk.name]]
    regions.list[[chunk.name]] <- nrow(chunk.regions)
    chunk.problems <- problems.by.chunk[[chunk.name]]
    if(is.null(chunk.problems)){
      stop("no problem data for chunk ", chunk.name)
    }
    peaks.by.regularization <- list()
    for(problem.name in names(chunk.problems)){
      problem <- chunk.problems[[problem.name]]
      if(is.data.frame(problem$peaks) && 0 < nrow(problem$peaks)){
        log.lambda.vec <- fit$predict(problem$features)
        if(is.numeric(problem$target)){
          too.hi <- problem$target[2] < log.lambda.vec
          too.lo <- log.lambda.vec < problem$target[1]
          outside.target.list[[problem.name]] <- too.hi | too.lo
        }
        for(regularization.i in seq_along(log.lambda.vec)){
          log.lambda <- log.lambda.vec[[regularization.i]]
          selected <- 
            subset(problem$modelSelection,
                   min.log.lambda < log.lambda &
                     log.lambda < max.log.lambda)
          stopifnot(nrow(selected) == 1)
          reg.str <- paste(fit$regularization.vec[[regularization.i]])
          peaks.by.regularization[[reg.str]][[problem.name]] <-
            subset(problem$peaks, peaks == selected$peaks)
        }#log.lambda
      }#if(is.data.frame(problem$peaks)
    }#problem.name
    metric.vec.list <- list()
    for(metric.name in c("fp", "fn", "possible.fp", "possible.tp")){
      metric.vec.list[[metric.name]] <-
        rep(NA, length(fit$regularization.vec))
    }
    for(regularization.i in seq_along(peaks.by.regularization)){
      chunk.peaks <-
        do.call(rbind, peaks.by.regularization[[regularization.i]])
      chunk.error <- PeakErrorSamples(chunk.peaks, chunk.regions)
      if(regularization.i == 1){ # first one for test error.
        result.list$peaks[[chunk.name]] <- chunk.peaks
        result.list$error.regions[[chunk.name]] <- chunk.error
      }
      for(metric.name in names(metric.vec.list)){
        metric.vec.list[[metric.name]][[regularization.i]] <-
          sum(chunk.error[, metric.name])
      }
    }
    error.vec.list[[chunk.name]] <- metric.vec.list
  }#chunk.name
  fp.mat <- do.call(cbind, lapply(error.vec.list, "[[", "fp"))
  fn.mat <- do.call(cbind, lapply(error.vec.list, "[[", "fn"))
  fn.possible.mat <- do.call(cbind, lapply(error.vec.list, "[[", "possible.tp"))
  fp.possible.mat <- do.call(cbind, lapply(error.vec.list, "[[", "possible.fp"))
  false.positives <- rowSums(fp.mat)
  false.negatives <- rowSums(fn.mat)
  false.positives.possible <- rowSums(fp.possible.mat)
  false.negatives.possible <- rowSums(fn.possible.mat)
  regions.vec <- do.call(c, regions.list)
  outside.target.mat <- do.call(rbind, outside.target.list)
  incorrect.regions <- false.positives + false.negatives
  incorrect.regions.possible <- sum(regions.vec)
  incorrect.targets <- colSums(outside.target.mat)
  incorrect.targets.possible <- nrow(outside.target.mat)
  metrics <- function(...){
    df.list <- list()
    for(metric.name in c(...)){
      possible <- get(paste0(metric.name, ".possible"))
      df.list[[metric.name]] <- 
        data.frame(metric.name,
                   metric.value=get(metric.name),
                   possible,
                   regularization=fit$regularization.vec,
                   row.names=NULL)
    }
    do.call(rbind, df.list)
  }
  result.list$metrics <-
    metrics("incorrect.targets", "incorrect.regions",
            "false.positives", "false.negatives")
  result.list
### list of metrics, error.regions, peaks.
}
