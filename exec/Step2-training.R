library(data.table)
library(ggplot2)
library(xtable)
library(PeakSegJoint)

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

estimate.regularization <- function(train.validation, n.folds=4){
  stopifnot(is.character(train.validation))
  stopifnot(2 <= length(train.validation))
  n.folds <- if(length(train.validation) < n.folds) length(train.validation)
  fold.id <- sample(rep(1:n.folds, l=length(train.validation)))
  picked.by.fold <- list()
  for(validation.fold in 1:n.folds){
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
      error.vec.list <- list()
      regions.list <- list()
      outside.target.list <- list()
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
              subset(problem$peaks, peaks == selected)
          }#log.lambda
        }#problem.name
        error.vec <- rep(NA, length(fit$regularization.vec))
        for(regularization.i in seq_along(peaks.by.regularization)){
          chunk.peaks <-
            do.call(rbind, peaks.by.regularization[[regularization.i]])
          chunk.error <- PeakErrorSamples(chunk.peaks, chunk.regions)
          error.vec[[regularization.i]] <- with(chunk.error, sum(fp+fn))
        }
        error.vec.list[[chunk.name]] <- error.vec
      }#chunk.name
      error.mat <- do.call(rbind, error.vec.list)
      regions.vec <- do.call(c, regions.list)
      outside.target.mat <- do.call(rbind, outside.target.list)
      regions <- sum(regions.vec)
      incorrect.regions <- colSums(error.mat)
      targets <- nrow(outside.target.mat)
      incorrect.targets <- colSums(outside.target.mat)
      tvdf <- function(metric.name){
        data.frame(tv,
                   metric.name,
                   metric.value=get(metric.name),
                   regularization=fit$regularization.vec)
      }
      error.by.tv[[tv]] <-
        rbind(tvdf("incorrect.targets"), tvdf("incorrect.regions"))
    }#tv
    tv.error <- do.call(rbind, error.by.tv)
    v.err <- subset(error.by.tv$validation, metric.name=="incorrect.regions")
    picked.i <- pick.best.index(v.err$metric.value)
    picked.error <- v.err[picked.i, ]
    tvPlot <- 
      ggplot()+
        geom_point(aes(-log10(regularization), metric.value, color=tv),
                   pch=1,
                   data=picked.error)+
        geom_line(aes(-log10(regularization), metric.value, color=tv),
                  data=tv.error)+
        theme_bw()+
        theme(panel.margin=grid::unit(0, "cm"))+
        facet_grid(metric.name ~ .)
    print(tvPlot)
    print(tv.error)
    picked.by.fold[[validation.fold]] <- picked.error
  }#validation.fold
  picked <- do.call(rbind, picked.by.fold)
  mean(picked$regularization)
}  

estimate.regularization(names(problems.by.chunk))

print("TODO: train it!")

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
