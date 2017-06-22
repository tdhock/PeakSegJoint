library(PeakSegJoint)

model.path <-
  system.file("exampleData",
              "PeakSegJoint-chunks",
              "trained.model.RData",
              package="PeakSegJoint")
model.path <- "~/exampleData/PeakSegJoint-chunks/trained.model.RData"
model.path <- "~/projects/H3K27ac_TDH/PeakSegJoint-chunks/trained.model.RData"
argv <- c(model.path, "1", "1")

argv <- commandArgs(trailingOnly=TRUE)

print(argv)

ppn <- PPN.cores()
if(!is.na(ppn))options(mc.cores=ppn/2)

if(length(argv) != 3){
  stop("usage: Step3e.R path/to/PeakSegJoint-chunks/trained.model.RData chunk.order.seed test.fold")
}

trained.model.RData <- normalizePath(argv[1], mustWork=TRUE)
chunk.order.seed <- as.integer(argv[2])
test.fold <- as.integer(argv[3])
chunks.dir <- dirname(trained.model.RData)
data.dir <- dirname(chunks.dir)

model.objs <- load(trained.model.RData)

## Divide chunks into train+validation/test.
map.RData <- file.path(data.dir, "chunk.file.map.RData")
load(map.RData)
chunks.by.file <- split(chunk.file.map, chunk.file.map$labels.file)
all.chunk.names <- names(problems.by.chunk)
names(all.chunk.names) <- basename(dirname(all.chunk.names))
cv.chunk.names <- all.chunk.names[paste(chunk.file.map$chunk.id)]
is.bad <- is.na(cv.chunk.names)
if(any(is.bad)){
  print(chunk.file.map[is.bad,])
  stop("some chunks have not been computed")
}

## For each test fold, hold it out and train a sequence of models with
## increasingly more data, and compute test error of each.
test.metrics.curve.list <- list()

is.test <- outer.fold.id == test.fold
set.seed(chunk.order.seed)
sets <-
  list(test=cv.chunk.names[is.test],
       train.validation=sample(cv.chunk.names[!is.test]))
test.regions <- regions.by.chunk[sets$test]
test.problems <- problems.by.chunk[sets$test]
test.features <- feature.mat[chunk.vec %in% sets$test,]

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
### Model fit list from IntervalRegressionCV.
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
        fmat <- rbind(colSums(problem$features))
        log.lambda.vec <- fit$predict(fmat)
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

if(length(sets$train.validation) < 2){
  print(sets$train.validation)
  stop("need at least 2 train chunks, please add more labels")
}
## to estimate if we have enough labeled chunks, we fit a sequence of
## models from 2 to N chunks.
train.chunks.vec <- 2:length(sets$train.validation)
for(train.chunks in train.chunks.vec){
  print(system.time({
    cat("estimating test error:",
        train.chunks, "/", length(sets$train.validation), "chunks.\n")
    some.train.validation <- sets$train.validation[1:train.chunks]
    row.is.train <- chunk.vec %in% some.train.validation
    fit <- penaltyLearning::IntervalRegressionCV(
      feature.mat[row.is.train,], target.mat[row.is.train,],
      min.observations=sum(row.is.train))
    test.results <- error.metrics(test.problems, test.regions, fit)
    test.results$metrics$test.fold <- test.fold
    test.results$metrics$chunk.order.seed <- chunk.order.seed
    test.results$metrics$train.chunks <- train.chunks
    stopifnot(test.results$metrics["incorrect.regions", "possible"] ==
                sum(sapply(test.regions, nrow)))
    test.metrics.curve.list[[paste(train.chunks)]] <- test.results$metrics
  }))
}
test.metrics.curve <- do.call(rbind, test.metrics.curve.list)
rownames(test.metrics.curve) <- NULL

test.out.dir <- file.path(chunks.dir, "figure-test-errors")
seed.RData <-
  file.path(test.out.dir,
            paste0("seed", chunk.order.seed, "fold", test.fold, ".RData"))

save(test.results,
     test.metrics.curve,
     file=seed.RData)

## ggplot()+
##   ggtitle(paste("test error for one",
##                 "random ordering of the labeled train chunks"))+
##   geom_text(aes(train.chunks, metric.value/possible*100,
##                 label=sprintf("%.1f%%", metric.value/possible*100)),
##             data=test.results$metrics,
##             vjust=-1)+
##   ylab("percent incorrect (test error)")+
##   scale_x_continuous("number of labeled chunks in train set",
##                      breaks=function(lim.vec){
##                        ##print(lim.vec)
##                        ceiling(lim.vec[1]):floor(lim.vec[2])
##                      })+
##   geom_line(aes(train.chunks, metric.value/possible*100,
##                 group=chunk.order.seed),
##             data=test.metrics.curve)+
##   geom_point(aes(train.chunks, metric.value/possible*100,
##                  group=chunk.order.seed),
##              data=test.metrics.curve)+
##   theme_bw()+
##   theme(panel.margin=grid::unit(0, "cm"))+
##   facet_grid(metric.name ~ test.fold, scales="free")

