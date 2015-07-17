library(ggplot2)
library(xtable)
library(PeakSegJoint)
options(xtable.print.results=FALSE)

argv <-
  system.file(file.path("exampleData",
                        "PeakSegJoint-chunks"),
              package="PeakSegDP")

argv <- "PeakSegJoint-chunks/H3K36me3_TDH_immune"
argv <- "~/exampleData/PeakSegJoint-chunks"

argv <- commandArgs(trailingOnly=TRUE)

print(argv)

PPN.cores()

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
labeled.problems.by.chunk <- list()
regions.by.chunk <- list()
for(chunk.i in seq_along(problems.RData.vec)){
  problems.RData <- problems.RData.vec[[chunk.i]]
  objs <- load(problems.RData)
  if(! "step2.data.list" %in% objs){
    stop("step.data.list not found in ", problems.RData)
  }
  res.data <- step2.data.list[[res.str]]
  for(problem.i in 1:nrow(res.data$problems)){
    prob.info <- res.data$problems[problem.i, ]
    problem.name <- paste(prob.info$problem.name)
    target <- step2.error.list[[paste(res.str, problem.name)]]$problem$target
    mlist <- step2.model.list[[problem.name]]
    stopifnot(prob.info$problemStart < mlist$peaks$chromStart)
    stopifnot(mlist$peaks$chromEnd < prob.info$problemEnd)
    mlist$target <- target
    problems.by.chunk[[problems.RData]][[problem.name]] <- mlist
    if(is.numeric(target)){
      labeled.problems.by.chunk[[problems.RData]][[problem.name]] <- mlist
    }
  }
  chunk.dir <- dirname(problems.RData)
  regions.RData <- file.path(chunk.dir, "regions.RData")
  load(regions.RData)
  regions.by.chunk[[problems.RData]] <- regions
}

train.problem.counts <- sapply(labeled.problems.by.chunk, length)
print(train.problem.counts)
stopifnot(train.problem.counts > 0)

set.seed(1)

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
    train.list <- do.call(c, labeled.problems.by.chunk[sets$train])
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

  error.by.fold <- mclapply.or.stop(1:n.folds, one.fold)
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

## Divide chunks into train+validation/test, compute test error.
all.chunk.names <- names(problems.by.chunk)
outer.folds <- 4
if(length(all.chunk.names) < outer.folds){
  outer.folds <- length(all.chunk.names)
}
outer.fold.id <- sample(rep(1:outer.folds, l=length(all.chunk.names)))
test.error.list <- list()
test.peaks.list <- list()
test.regions.list <- list()
for(test.fold in 1:outer.folds){
  cat(sprintf("estimating test error for fold %4d / %4d\n",
              test.fold, outer.folds))
  is.test <- outer.fold.id == test.fold
  sets <- list(train.validation=all.chunk.names[!is.test],
               test=all.chunk.names[is.test])
  mean.reg <- estimate.regularization(sets$train.validation)
  tv.list <- do.call(c, labeled.problems.by.chunk[sets$train.validation])
  tv.fit <-
    IntervalRegressionProblems(tv.list,
                               initial.regularization=mean.reg,
                               factor.regularization=10000,
                               verbose=0)
  test.results <- error.metrics(sets$test, tv.fit)
  test.regions.list[names(test.results$error.regions)] <-
    test.results$error.regions
  test.peaks.list[names(test.results$peaks)] <- test.results$peaks
  test.metrics <-
    subset(test.results$metrics, regularization == regularization[1])
  rownames(test.metrics) <- NULL
  test.error.list[[paste("test fold", test.fold)]] <-
    data.frame(test.fold, test.metrics)
}

incorrect <- as.integer(rowSums(sapply(test.error.list, "[[", "metric.value")))
possible <- as.integer(rowSums(sapply(test.error.list, "[[", "possible")))
percent.incorrect <- incorrect / possible * 100
test.error.summary <- 
  data.frame(row.names=test.error.list[[1]]$metric.name,
             incorrect, possible, percent.incorrect)

test.out.dir <- file.path(chunks.dir, "figure-test-errors")
dir.create(test.out.dir, showWarnings=FALSE)
test.index.html <- file.path(test.out.dir, "index.html")
test.row.list <- list()
for(problems.RData in names(test.regions.list)){
  chunk.dir <- dirname(problems.RData)
  chunk.id <- basename(chunk.dir)
  test.error.figure <-
    sprintf('<img src="%s.png" alt="chunk%s" />', chunk.id, chunk.id)
  error.regions <- test.regions.list[[problems.RData]]
  chunk.peaks <- test.peaks.list[[problems.RData]]
  out.RData <- file.path(test.out.dir, paste0(chunk.id, ".RData"))
  save(error.regions, chunk.peaks, file=out.RData)
  test.row.list[[chunk.id]] <- with(error.regions, {
    data.frame(test.error.figure,
               errors=sprintf("<pre>%4d / %4d</pre>",
                 sum(fp+fn), length(fp)),
               fp=sprintf("<pre>%4d / %4d</pre>",
                 sum(fp), sum(possible.fp)),
               fn=sprintf("<pre>%4d / %4d</pre>",
                 sum(fn), sum(possible.tp)))
  })
}
test.row.df <- do.call(rbind, test.row.list)
test.xt <- xtable(test.row.df)
test.html.table <-
  print(test.xt, type="html",
        include.rownames=FALSE,
        sanitize.text.function=identity)
summary.xt <- xtable(test.error.summary)
summary.html <- print(summary.xt, type="html",
                      include.rownames=TRUE)
test.html.out <-
  paste("<h1>Test error summary</h1>",
        summary.html,
        "<p>Targets counts examples (genomic regions to segment)",
        "in the interval regression problem.</p>",
        "<p>Regions, false postives, and false negatives",
        "count labels (peakStart, peakEnd, peaks, noPeaks).</p>",
        "<p>Test error was estimated using",
        outer.folds, "fold cross-validation.",
        "</p>",
        "<h1>Test error details for each chunk of labels</h1>",
        test.html.table)
cat(test.html.out, file=test.index.html)

stopifnot(test.error.summary["incorrect.regions", "possible"] ==
            sum(sapply(regions.by.chunk, nrow)))

## Fit model to all training data.
mean.reg <- estimate.regularization(all.chunk.names)
train.list <- do.call(c, labeled.problems.by.chunk)
full.fit <-
  IntervalRegressionProblems(train.list,
                             initial.regularization=mean.reg,
                             factor.regularization=1.1,
                             verbose=0)

data.dir <- dirname(chunks.dir)
bigwig.file.vec <- Sys.glob(file.path(data.dir, "*", "*.bigwig"))
bigwig.file <- bigwig.file.vec[1]
chrom.ranges <- bigWigInfo(bigwig.file)

hg19.bases <- 3137161264
max.jobs <- 100 ## soft limit.
min.bases.per.job <- 1e7
total.bases <- sum(chrom.ranges$chromEnd)
bases.per.job <- total.bases / max.jobs
if(bases.per.job < min.bases.per.job){
  bases.per.job <- min.bases.per.job
}
bases.per.job <- as.integer(bases.per.job)

problems.by.chrom <- list()
jobs.by.chrom <- list()
for(chrom.i in 1:nrow(chrom.ranges)){
  chrom.range <- chrom.ranges[chrom.i, ]
  cat(sprintf("%4d / %4d %s\n", chrom.i, nrow(chrom.ranges), chrom.range$chrom))
  chrom.bases <- with(chrom.range, chromEnd - chromStart)
  chrom.problems <- with(chrom.range, {
    getProblems(chrom, chromStart, chromEnd, bases.per.problem)
  })
  jobStart <- with(chrom.range, {
    seq(chromStart, chromEnd, by=bases.per.job)
  })
  jobEnd <- jobStart + bases.per.job + bases.per.problem
  jobEnd[chrom.range$chromEnd < jobEnd] <- chrom.range$chromEnd
  job.name <- sprintf("%s:%d-%d", chrom.range$chrom, jobStart, jobEnd)
  chrom.jobs <- data.table(chrom=chrom.range$chrom, job.name, jobStart, jobEnd)
  ## chrom.problems$job.name <- NA
  ## for(job.i in 1:nrow(chrom.jobs)){
  ##   job <- chrom.jobs[job.i, ]
  ##   to.assign <-
  ##     job$jobStart <= chrom.problems$problemStart &
  ##       chrom.problems$problemEnd <= job$jobEnd
  ##   chrom.problems$job.name[to.assign] <- job$job.name
  ## }
  ## stopifnot(!is.na(chrom.problems$job.name))
  problems.by.chrom[[paste(chrom.range$chrom)]] <- chrom.problems
  jobs.by.chrom[[paste(chrom.range$chrom)]] <- chrom.jobs
}
test.problems <- do.call(rbind, problems.by.chrom)
test.jobs <- do.call(rbind, jobs.by.chrom)
problem.job.tab <- table(test.problems$job.name)
## cat(length(problem.job.tab), "jobs with",
##     nrow(test.problems), "problems.",
##     "Problems per job:\n")
## print(quantile(problem.job.tab))
## tab.sorted <- sort(problem.job.tab)
## print(head(tab.sorted))
## print(tail(tab.sorted))

trained.model.RData <- file.path(chunks.dir, "trained.model.RData")
save(train.errors, train.errors.picked,
     full.fit,
     test.problems,
     ##test.jobs,
     file=trained.model.RData)

