library(data.table)
library(ggplot2)
library(xtable)
library(PeakSegJoint)
library(parallel)

argv <-
  system.file(file.path("exampleData",
                        "PeakSegJoint-chunks"),
              package="PeakSegDP")

argv <- "PeakSegJoint-chunks/H3K36me3_TDH_immune"
argv <- "~/exampleData/PeakSegJoint-chunks"

argv <- commandArgs(trailingOnly=TRUE)

print(argv)

if(length(argv) != 1){
  stop("usage: Step3.R path/to/PeakSegJoint-chunks")
}

chunks.dir <- normalizePath(argv[1])

problems.RData.vec <- Sys.glob(file.path(chunks.dir, "*", "problems.RData"))

err.mat.list <- list()
chunk.best.list <- list()
bpp.list <- list()
for(chunk.i in seq_along(problems.RData.vec)){
  problems.RData <- problems.RData.vec[[chunk.i]]
  chunk.dir <- dirname(problems.RData)
  chunk.id <- basename(chunk.dir)
  index.html <-
    file.path("..", chunk.id, "figure-train-errors", "index.html")
  first.png <-
    file.path("..", chunk.id, "figure-train-errors.png")
  objs <- load(problems.RData)
  err.vec <- res.error$errors
  min.df <- subset(res.error, errors==min(errors))
  best.errors <- min(res.error$errors)
  bpp.list[[chunk.id]] <- min.df$bases.per.problem
  href <- sprintf('<a href="%s">
  <img src="%s" alt="chunk%s" />
</a>', index.html, first.png, chunk.id)
  chunk.best.list[[chunk.id]] <-
    data.frame(selected.model=href,
               best.errors,
               regions=res.error$regions[1])
  names(err.vec) <- res.error$bases.per.problem
  err.mat.list[[chunk.i]] <- err.vec
}
chunk.best <- do.call(rbind, chunk.best.list)
err.mat <- do.call(rbind, err.mat.list)
err.vec <- colSums(err.mat)
train.errors <-
  data.frame(bases.per.problem=as.integer(names(err.vec)),
             errors=err.vec,
             regions=sum(chunk.best$regions))
chosen.i <- pick.best.index(err.vec)
res.str <- names(err.vec)[chosen.i]
bases.per.problem <- as.integer(res.str)
train.errors.picked <- train.errors[chosen.i, ]

chunk.best$selected.errors <- err.mat[, res.str]
chunk.best$best.bases.per.problem <- NA
for(chunk.i in seq_along(bpp.list)){
  bpp <- bpp.list[[chunk.i]]
  bpp[bpp == res.str] <-
    paste0("<b>", bpp[bpp==res.str], "</b>")
  chunk.best$best.bases.per.problem[chunk.i] <- paste0(bpp, collapse="<br />")
}

chunk.ordered <-
  chunk.best[order(chunk.best$selected.errors, decreasing=TRUE),]

resCurve <- 
ggplot()+
  geom_line(aes(bases.per.problem, errors),
            data=train.errors)+
  scale_x_log10()+
  ylab("minimum incorrect labels (train error)")+
  ggtitle(paste(res.str, "bases/problem"))+
  geom_point(aes(bases.per.problem, errors),
             data=train.errors.picked,
             pch=1)

png.name <-
  file.path(chunks.dir, "figure-train-errors", "figure-train-errors.png")
png.dir <- dirname(png.name)
dir.create(png.dir, showWarnings=FALSE, recursive=TRUE)
png(png.name, width=400, h=300, units="px")
print(resCurve)
dev.off()

res.xt <- xtable(train.errors)
res.html <-
  print(res.xt, type="html",
        include.rownames=FALSE)

xt <- xtable(chunk.ordered)
html.table <-
  print(xt, type="html",
        include.rownames=FALSE,
        sanitize.text.function=identity)
html.out <-
  paste("<h1>Train error totals per problem size</h1>",
        "<table><tr>",
        '<td> <img src="figure-train-errors.png" /> </td>',
        "<td>", res.html, "</td>",
        "</tr></table>",
        "<h1>Train error details for each chunk of labels</h1>",
        html.table)

out.file <- file.path(chunks.dir, "figure-train-errors", "index.html")

cat(html.out, file=out.file)

problems.by.chunk <- list()
regions.by.chunk <- list()
for(chunk.i in seq_along(problems.RData.vec)){
  problems.RData <- problems.RData.vec[[chunk.i]]
  objs <- load(problems.RData)
  if(! "step2.data.list" %in% objs){
    stop("step.data.list not found in ", problems.RData)
  }
  res.data <- step2.data.list[[res.str]]
  for(problem.name in names(res.data$regions)){
    target <- step2.error.list[[paste(res.str, problem.name)]]$problem$target
    if(is.numeric(target)){
      n.finite <- sum(is.finite(target))
      if(n.finite > 0){
        mlist <- step2.model.list[[problem.name]]
        mlist$target <- target
        problems.by.chunk[[problems.RData]][[problem.name]] <- mlist
      }
    }
  }
  chunk.dir <- dirname(problems.RData)
  regions.RData <- file.path(chunk.dir, "regions.RData")
  load(regions.RData)
  regions.by.chunk[[problems.RData]] <- regions
}

train.problem.counts <- sapply(problems.by.chunk, length)
print(train.problem.counts)

stopifnot(train.problem.counts > 0)

set.seed(1)

my.mclapply <- function(...){
  result.list <- mclapply(...)
  is.error <- sapply(result.list, inherits, "try-error")
  if(any(is.error)){
    print(result.list[is.error])
    stop("errors in mclapply")
  }
  result.list
}

tv.curves <- function(train.validation, n.folds=4){
  stopifnot(is.character(train.validation))
  stopifnot(2 <= length(train.validation))
  if(length(train.validation) < n.folds){
    n.folds <- length(train.validation)
  }
  fold.id <- sample(rep(1:n.folds, l=length(train.validation)))

  one.fold <- function(validation.fold){
    is.validation <- fold.id == validation.fold
    sets <- list(validation=train.validation[is.validation],
                 train=train.validation[!is.validation])
    train.list <- do.call(c, problems.by.chunk[sets$train])
    fit <-
      IntervalRegressionProblems(train.list,
                                 initial.regularization=0.005,
                                 factor.regularization=1.1,
                                 verbose=0)
    error.by.tv <- list()
    for(tv in names(sets)){
      chunk.name.vec <- sets[[tv]]
      error.results <- error.metrics(chunk.name.vec, fit)
      error.df <- error.results$metrics
      error.by.tv[[tv]] <- data.frame(tv, validation.fold, error.df)
    }#tv
    do.call(rbind, error.by.tv)
  }#validation.fold

  error.by.fold <- my.mclapply(1:n.folds, one.fold)
  do.call(rbind, error.by.fold)
}

error.metrics <- function(chunk.name.vec, fit){
  error.vec.list <- list()
  regions.list <- list()
  outside.target.list <- list()
  result.list <- list()
  for(chunk.name in chunk.name.vec){
    chunk.regions <- regions.by.chunk[[chunk.name]]
    regions.list[[chunk.name]] <- nrow(chunk.regions)
    chunk.problems <- problems.by.chunk[[chunk.name]]
    peaks.by.regularization <- list()
    for(problem.name in names(chunk.problems)){
      problem <- chunk.problems[[problem.name]]
      log.lambda.vec <- fit$predict(problem$features)
      too.hi <- problem$target[2] < log.lambda.vec
      too.lo <- log.lambda.vec < problem$target[1]
      outside.target.list[[problem.name]] <- too.hi | too.lo
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
  fp.mat <- sapply(error.vec.list, "[[", "fp")
  fn.mat <- sapply(error.vec.list, "[[", "fn")
  fn.possible.mat <- sapply(error.vec.list, "[[", "possible.tp")
  fp.possible.mat <- sapply(error.vec.list, "[[", "possible.fp")
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
                   regularization=fit$regularization.vec)
    }
    do.call(rbind, df.list)
  }
  result.list$metrics <-
    metrics("incorrect.targets", "incorrect.regions",
            "false.positives", "false.negatives")
  result.list
}

estimate.regularization <- function(train.validation){
  full.curves <- tv.curves(train.validation)

  v.err <-
    subset(full.curves,
           metric.name=="incorrect.regions" &
             tv=="validation")
  v.list <- split(v.err, v.err$validation.fold)
  picked.list <- list()
  for(validation.fold in names(v.list)){
    v <- v.list[[validation.fold]]
    picked.i <- pick.best.index(v$metric.value)
    picked.list[[validation.fold]] <- v[picked.i, ]
  }
  picked.error <- do.call(rbind, picked.list)

  tvPlot <- 
    ggplot()+
      geom_point(aes(-log10(regularization), metric.value, color=tv),
                 pch=1,
                 data=picked.error)+
      geom_line(aes(-log10(regularization), metric.value, color=tv),
                data=full.curves)+
      theme_bw()+
      theme(panel.margin=grid::unit(0, "cm"))+
      facet_grid(metric.name ~ validation.fold, scales="free")
  print(tvPlot)

  mean(picked.error$regularization)
}

all.chunk.names <- names(problems.by.chunk)
mean.reg <- estimate.regularization(all.chunk.names)
train.list <- do.call(c, problems.by.chunk)
fit <-
  IntervalRegressionProblems(train.list,
                             initial.regularization=mean.reg,
                             factor.regularization=1.1,
                             verbose=0)

outer.folds <- 4
if(length(all.chunk.names) < outer.folds){
  outer.folds <- length(all.chunk.names)
}
outer.fold.id <- sample(rep(1:outer.folds, l=length(all.chunk.names)))
test.error.list <- list()
test.peaks.list <- list()
test.regions.list <- list()
for(test.fold in 1:outer.folds){
  is.test <- outer.fold.id == test.fold
  sets <- list(train.validation=all.chunk.names[!is.test],
               test=all.chunk.names[is.test])
  mean.reg <- estimate.regularization(sets$train.validation)
  tv.list <- do.call(c, problems.by.chunk[sets$train.validation])
  fit <-
    IntervalRegressionProblems(tv.list,
                               initial.regularization=mean.reg,
                               factor.regularization=10000,
                               verbose=0)
  test.results <- error.metrics(sets$test, fit)
  test.regions.list[names(test.results$error.regions)] <-
    test.results$error.regions
  test.peaks.list[names(test.results$peaks)] <- test.results$peaks
  test.metrics <-
    subset(test.results$metrics, regularization == regularization[1])
  rownames(test.metrics) <- NULL
  test.error.list[[paste("test fold", test.fold)]] <-
    data.frame(test.fold, test.metrics)
}

incorrect <- rowSums(sapply(test.error.list, "[[", "metric.value"))
possible <- rowSums(sapply(test.error.list, "[[", "possible"))
percent.incorrect <- incorrect / possible * 100
test.error.summary <- 
  data.frame(metric.name=test.error.list[[1]]$metric.name,
             incorrect, possible, percent.incorrect)

print("TODO: plot predictions and test errors for each chunk")

coverage.RData.vec <- Sys.glob(file.path(chunks.dir, "*", "*", "*.RData"))
coverage.RData <- coverage.RData.vec[1]
cobjs <- load(coverage.RData)

## These calculations try to split the genome into manageable
## subsets/jobs, each of which can fit several samples' coverage data
## into RAM at once.
rows.per.megabyte <- 53000 
bases.per.row <- 20
max.samples <- 100
max.megabytes <- 4000
max.megabytes.per.sample <- max.megabytes/max.samples
max.rows.per.sample <- max.megabytes.per.sample * rows.per.megabyte
max.bases.per.sample <- max.rows.per.sample * bases.per.row
hg19.bases <- 3137161264
bases.per.job <- as.integer(max.bases.per.sample)
estimated.jobs <- hg19.bases/bases.per.job
problems.per.job <- bases.per.job/bases.per.problem
problems.by.chrom <- list()
jobs.by.chrom <- list()
for(chrom.i in 1:nrow(chrom.ranges)){
  chrom.range <- chrom.ranges[chrom.i, ]
  cat(sprintf("%4d / %4d %s\n", chrom.i, nrow(chrom.ranges), chrom.range$chrom))
  chrom.bases <- with(chrom.range, max.chromEnd - min.chromStart)
  chrom.problems <- with(chrom.range, {
    getProblems(chrom, min.chromStart, max.chromEnd, bases.per.problem)
  })
  jobStart <- with(chrom.range, {
    seq(min.chromStart, max.chromEnd, by=bases.per.job)
  })
  jobEnd <- jobStart + bases.per.job + bases.per.problem
  jobEnd[chrom.range$max.chromEnd < jobEnd] <- chrom.range$max.chromEnd
  job.name <- sprintf("%s:%d-%d", chrom.range$chrom, jobStart, jobEnd)
  chrom.jobs <- data.table(job.name, jobStart, jobEnd)
  chrom.problems$job.name <- NA
  for(job.i in 1:nrow(chrom.jobs)){
    job <- chrom.jobs[job.i, ]
    to.assign <-
      job$jobStart <= chrom.problems$problemStart &
        chrom.problems$problemEnd <= job$jobEnd
    chrom.problems$job.name[to.assign] <- job$job.name
  }
  stopifnot(!is.na(chrom.problems$job.name))
  problems.by.chrom[[paste(chrom.range$chrom)]] <- chrom.problems
  jobs.by.chrom[[paste(chrom.range$chrom)]] <- chrom.jobs
}
test.problems <- do.call(rbind, problems.by.chrom)
test.jobs <- do.call(rbind, jobs.by.chrom)
problem.job.tab <- table(test.problems$job.name)
cat(length(problem.job.tab), "jobs with",
    nrow(test.problems), "problems.",
    "Problems per job:\n")
print(quantile(problem.job.tab))
tab.sorted <- sort(problem.job.tab)
print(head(tab.sorted))
print(tail(tab.sorted))

trained.model.RData <- file.path(chunks.dir, "trained.model.RData")
save(train.errors, train.errors.picked,
     test.problems, test.jobs,
     file=trained.model.RData)

problems.by.job <- split(test.problems, test.problems$job.name)
Step4 <-
  system.file(file.path("exec", "Step4-test-segmentation.R"),
              mustWork=TRUE,
              package="PeakSegJoint")
R.bin <- R.home("bin")
Rscript <- file.path(R.bin, "Rscript")
for(job.name in names(problems.by.job)){
  job.problems <- problems.by.job[[job.name]]
  cmd <- paste("qsub", Rscript, Step4, trained.model.RData, job.name)
  print(cmd)
}
