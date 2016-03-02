library(ggplot2)
library(xtable)
library(PeakSegJoint)
options(xtable.print.results=FALSE)

ex.dir <-
  system.file(file.path("exampleData",
                          "PeakSegJoint-chunks"),
                package="PeakSegJoint")
ex.dir <- "PeakSegJoint-chunks/H3K36me3_TDH_immune"
ex.dir <- "~/projects/H3K27ac_TDH/PeakSegJoint-chunks"
ex.dir <- "~/exampleData/PeakSegJoint-chunks"
argv <- file.path(ex.dir, c(1:3), "models.by.problem.RData")

argv <- commandArgs(trailingOnly=TRUE)

print(argv)

ppn <- PPN.cores()
if(!is.na(ppn))options(mc.cores=ppn/2)

if(length(argv) == 0){
  stop("usage: Step2.R chunk1/models.by.problem.RData [...]")
}

models.RData.vec <- normalizePath(argv, mustWork=TRUE)

models.by.res <- list()
for(models.RData in models.RData.vec){
  objs <- load(models.RData)
  chunk.id <- basename(dirname(models.RData))
  for(problem.i in 1:nrow(sample.problems.dt)){
    problem <- sample.problems.dt[problem.i,]
    res.str <- paste(problem$bases.per.problem)
    models.by.res[[res.str]][[paste(chunk.id, problem$problem.name)]] <-
      models.by.problem[[problem.i]]
  }
}

set.seed(1)
for(res.str in names(models.by.res)){
  models.by.problem <- models.by.res[[res.str]]
  feature.mat.list <- list()
  target.mat.list <- list()
  for(chunk.problem in names(models.by.problem)){
    models.by.sample <- models.by.problem[[chunk.problem]]
    for(id.group in names(models.by.sample)){
      model <- models.by.sample[[id.group]]
      if(any(is.finite(model$target))){
        target.mat.list[[paste(chunk.problem, id.group)]] <- model$target
        feature.mat.list[[paste(chunk.problem, id.group)]] <- model$features
      }
    }
  }
  feature.mat <- do.call(rbind, feature.mat.list)
  target.mat <- do.call(rbind, target.mat.list)
  fit <- IntervalRegressionMatrixCV(feature.mat, target.mat)
  for(chunk.problem in names(models.by.problem)){
    models.by.sample <- models.by.problem[[chunk.problem]]
    for(id.group in names(models.by.sample)){
      model <- models.by.sample[[id.group]]
      pred.log.lambda <- as.numeric(fit$predict(rbind(model$features)))
      selected <- subset(
        model$modelSelection,
        min.log.lambda < pred.log.lambda &
          pred.log.lambda < max.log.lambda)
      stopifnot(nrow(selected)==1)
      peaks.str <- paste(selected$peaks)
      pred.peaks <- model$fit$peaks[[peaks.str]]
      stop("TODO store metadata, plot pred peaks")
    }
  }
}

if(all(sapply(sample.results.list, is.null))){
  print(sample.problems.dt)
  stop("no computable models for any uniform size segmentation problems")
}
sample.results <- do.call(rbind, sample.results.list)

## Plot step1 problems with detected peaks.
## ggplot()+
##   geom_segment(aes(problemStart/1e3, problem.name,
##                    xend=problemEnd/1e3, yend=problem.name),
##                data=sample.results)+
##   geom_segment(aes(chromStart/1e3, problem.name,
##                    xend=chromEnd/1e3, yend=problem.name),
##                size=2,
##                color="deepskyblue",
##                data=sample.results)+
##   theme_bw()+
##   theme(panel.margin=grid::unit(0, "cm"))+
##   facet_grid(bases.per.problem ~ ., scales="free", space="free")

step1.by.res <- split(step1.results, step1.results$bases.per.problem)
Step1Step2 <- function(res.str){
  bases.per.problem <- as.integer(res.str)
  step1.peaks <- step1.by.res[[res.str]]
  step2.problems <- clusterProblems(step1.peaks)
  
  problems.dt <- data.table(bases.per.problem, step2.problems)
  problems.dt[, problemStart1 := problemStart + 1L]
  setkey(problems.dt, problemStart1, problemEnd)
  setkey(regions, chromStart1, chromEnd)
  over.regions <- foverlaps(regions, problems.dt, nomatch=0L)
  over.regions[,
               `:=`(overlapStart=ifelse(problemStart < chromStart,
                      chromStart, problemStart),
                    overlapEnd=ifelse(problemEnd < chromEnd,
                      problemEnd, chromEnd))]
  over.regions[,
               overlapBases := overlapEnd-overlapStart]
  wrong.direction <- with(over.regions, {
    (annotation=="peakEnd" & chromStart < problemStart) |
      (annotation=="peakStart" & problemEnd < chromEnd)
  })
  wrong.regions <-
    unique(over.regions[wrong.direction, .(problem.name, chromStart)])
  setkey(over.regions, problem.name, chromStart)
  over.regions[wrong.regions, overlapBases := 0]
  ## ggplot()+
  ##   geom_tallrect(aes(xmin=chromStart/1e3, xmax=chromEnd/1e3),
  ##                 data=over.regions,
  ##                 fill="grey")+
  ##   geom_text(aes(problemStart/1e3, factor(problem.i),
  ##                 label=overlapBases),
  ##             hjust=1,
  ##             data=over.regions)+
  ##   theme_bw()+
  ##   theme(panel.margin=grid::unit(0, "cm"))+
  ##   facet_grid(sample.id ~ ., scales="free")+
  ##   geom_segment(aes(problemStart/1e3, factor(problem.i),
  ##                    xend=problemEnd/1e3, yend=factor(problem.i)),
  ##                data=over.regions)
  region.i.problems <-
    over.regions[,
                 .(problem.name=problem.name[which.max(overlapBases)]),
                 by=region.i]
  stopifnot(nrow(region.i.problems) <= nrow(regions))
  setkey(regions, region.i)
  setkey(region.i.problems, region.i)
  assigned.regions <- regions[region.i.problems,]
  stopifnot(nrow(assigned.regions) == nrow(region.i.problems))
  regions.by.problem <-
    split(assigned.regions, assigned.regions$problem.name, drop=TRUE)

  list(problems=problems.dt,
       regions=regions.by.problem)
}
step2.data.list <- mclapply.or.stop(names(step1.by.res), Step1Step2)
names(step2.data.list) <- names(step1.by.res)

step2.problems.list <- list()
problems.with.regions.list <- list()
for(res.str in names(step2.data.list)){
  res.data <- step2.data.list[[res.str]]
  problems.dt <- res.data$problems[, .(problem.name, problemStart, problemEnd)]
  problems.by.name <- split(problems.dt, problems.dt$problem.name)
  step2.problems.list[names(problems.by.name)] <- problems.by.name
  if(length(res.data$regions)){
    problems.with.regions.list[[res.str]] <- 
      data.table(bases.per.problem=as.integer(res.str),
                 problem.name=names(res.data$regions))
  }
}
step2.problems <- do.call(rbind, step2.problems.list)
problems.with.regions <- do.call(rbind, problems.with.regions.list)

blank <- unique(counts[, .(sample.id)])
SegmentStep2 <- function(row.i){
  ##print(row.i)
  problem <- step2.problems[row.i, ]
  problem[, problemStart1 := problemStart + 1L]
  setkey(problem, problemStart1, problemEnd)
  problem.counts <-
    foverlaps(counts, problem, nomatch=0L, type="any")
  problem.counts[, .(chromStart=min(chromStart),
                     chromEnd=max(chromEnd)),
                 by=sample.id]
  problem.counts[, `:=`(
    chromStart=ifelse(chromStart < problemStart, problemStart, chromStart),
    chromEnd=ifelse(problemEnd < chromEnd, problemEnd, chromEnd)
    )]
  coverage.limits <- 
    problem.counts[, .(chromStart=min(chromStart),
                       chromEnd=max(chromEnd)),
                   by=sample.id]
  with(coverage.limits, {
    stopifnot(problem$problemStart <= chromStart)
    stopifnot(chromEnd <= problem$problemEnd)
  })
  with(problem.counts, stopifnot(chromStart < chromEnd))

  ## counts.by.sample <- table(problem.counts$sample.id)
  ## is.segment <- counts.by.sample == 1
  ## is.step <- !is.segment
  ## step.ids <- names(counts.by.sample)[is.step]
  ## segment.ids <- names(counts.by.sample)[is.segment]
  ## setkey(problem.counts, sample.id)
  ## segment.dt <- problem.counts[segment.ids]
  ## step.dt <- problem.counts[step.ids]
  ## ggplot()+
  ##   theme_bw()+
  ##   theme(panel.margin=grid::unit(0, "cm"))+
  ##   facet_grid(sample.id ~ ., labeller=function(var, val){
  ##     sub("McGill0", "", sub(" ", "\n", val))
  ##   }, scales="free")+
  ##   geom_segment(aes((problemStart+0.5)/1e3, 0,
  ##                    xend=problemEnd/1e3, yend=0),
  ##                data=problem,
  ##                color="black",
  ##                size=1)+
  ##   geom_segment(aes((chromStart+0.5)/1e3, count,
  ##                    xend=chromEnd/1e3, yend=count),
  ##                data=problem.counts,
  ##                color="grey50")+
  ##   ## geom_segment(aes((chromStart+0.5)/1e3, count,
  ##   ##                  xend=chromEnd/1e3, yend=count),
  ##   ##              data=segment.dt,
  ##   ##              color="grey50")+
  ##   ## geom_step(aes((chromStart+0.5)/1e3, count),
  ##   ##           data=step.dt,
  ##   ##           color="grey50")+
  ##   geom_blank(data=blank)

  profile.list <- ProfileList(problem.counts)
  fit <- tryCatch({
    PeakSegJointSeveral(profile.list)
  }, error=function(e){
    NULL
  })
  if(is.null(fit)){
    return(NULL)
  }
  info <- ConvertModelList(fit)
  info$features <- featureMatrix(profile.list)
  info
}
step2.model.list <-
  mclapply.or.stop(seq_along(step2.problems$problem.name), SegmentStep2)
## It is OK to index the model list on problem name (even though the
## same problem could occur in several resolutions), since anyways the
## model should not change between resolutions.
names(step2.model.list) <- step2.problems$problem.name
stopifnot(table(names(step2.model.list)) == 1)

ProblemError <- function(row.i){
  prob.meta <- problems.with.regions[row.i, ]
  one.res <- step2.data.list[[paste(prob.meta$bases.per.problem)]]
  problem.name <- paste(prob.meta$problem.name)
  problem.regions <- one.res$regions[[problem.name]]
  converted <- step2.model.list[[problem.name]]
  if(is.null(converted)){
    pred.peaks <- Peaks()
    prob.err.list <- PeakSegJointError(list(peaks=NULL), problem.regions)
  }else{
    prob.err.list <- PeakSegJointError(converted, problem.regions)
    best.models <-
      subset(prob.err.list$modelSelection, errors==min(errors))
    peaks.num <- min(best.models$peaks)
    pred.peaks <- if(peaks.num > 0){
      subset(converted$peaks, peaks == peaks.num)
    }
  }
  list(problem=prob.err.list,
       peaks=pred.peaks)
}
step2.error.list <-
  mclapply.or.stop(seq_along(problems.with.regions$problem.name), ProblemError)
names(step2.error.list) <- with(problems.with.regions, {
  paste(bases.per.problem, problem.name)
})

ResError <- function(res.str){
  res.data <- step2.data.list[[res.str]]
  error.row.name.vec <- with(problems.with.regions.list[[res.str]], {
    paste(bases.per.problem, problem.name)
  })
  pred.peaks.list <- list()
  for(error.row.name in error.row.name.vec){
    error.info <- step2.error.list[[error.row.name]]
    pred.peaks.list[[error.row.name]] <- error.info$peaks
  }
  pred.peaks <- if(length(pred.peaks.list) == 0){
    Peaks()
  }else{
    do.call(rbind, pred.peaks.list)
  }
  res.regions <- do.call(rbind, res.data$regions)
  error.regions <- PeakErrorSamples(pred.peaks, res.regions)
  with(error.regions, {
    data.table(bases.per.problem=as.integer(res.str),
               fp=sum(fp),
               fn=sum(fn),
               errors=sum(fp+fn),
               regions=length(fp))
  })
}
res.error.list <- mclapply.or.stop(names(problems.with.regions.list), ResError)
res.error <- do.call(rbind, res.error.list)
print(res.error)

## Check to make sure each peak occurs within its problem.
for(prob.i in 1:nrow(step2.problems)){
  prob <- step2.problems[prob.i, ]
  problem.name <- paste(prob$problem.name)
  model <- step2.model.list[[problem.name]]
  stopifnot(prob$problemStart < model$peaks$chromStart)
  stopifnot(model$peaks$chromEnd < prob$problemEnd)
}

## Save results for this chunk/resolution.
problems.RData <- file.path(chunk.dir, "problems.RData")
save(step1.results,
     ##step2.problems,
     ##regions.by.problem,
     step2.data.list,
     step2.model.list,
     step2.error.list,
     res.error,
     file=problems.RData)



err.mat.list <- list()
chunk.best.list <- list()
bpp.list <- list()
for(chunk.i in seq_along(problems.RData.vec)){
  problems.RData <- problems.RData.vec[[chunk.i]]
  message(sprintf("%4d / %4d chunks %s", chunk.i, length(problems.RData.vec),
                  problems.RData))
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
             errors=as.integer(err.vec),
             regions=sum(chunk.best$regions))
chosen.i <- max(which(err.vec == min(err.vec)))
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
print(train.errors)
cat("Picked the following problem size:\n")
print(train.errors.picked)

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
        html.table,
        "<h1>Train/validation error curves for selecting regularization</h1>",
        '<p><img src="figure-regularization.png" /></td>')

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
  if(! "step2.error.list" %in% objs){
    stop("step.error.list not found in ", problems.RData)
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

## Fit model to all training data. need (unlabeled) problems.by.chunk
## here since the predictions must be made on all problems in the
## validation set (not just the subset of problems for which we have
## assigned some regions and computed a target).
full.curves <- tv.curves(problems.by.chunk, regions.by.chunk)
picked.error <- best.on.validation(full.curves)

tvPlot <- 
  ggplot()+
    geom_point(aes(-log10(regularization), metric.value, color=tv),
               pch=1,
               data=picked.error)+
    geom_line(aes(-log10(regularization), metric.value, color=tv),
              data=full.curves)+
    theme_bw()+
    theme(panel.margin=grid::unit(0, "cm"))+
    facet_grid(metric.name ~ validation.fold, labeller=function(var, val){
      if(var=="validation.fold"){
        paste("fold", val)
      }else{
        paste(val)
      }
    }, scales="free")
reg.png <-
  file.path(chunks.dir, "figure-train-errors", "figure-regularization.png")
png(reg.png, width=600, h=400, units="px")
print(tvPlot)
dev.off()

mean.reg <- mean(picked.error$regularization)

## labeled.problems.by.chunk needed here since
## IntervalRegressionProblems does not accept problems without
## targets.
train.list <- do.call(c, labeled.problems.by.chunk)
full.fit <-
  IntervalRegressionProblems(train.list,
                             initial.regularization=mean.reg,
                             factor.regularization=NULL,
                             verbose=0)

## get chrom size info from a bigwig, so we can generate a list of
## segmentation problems and divide them into jobs.
bigwig.file.vec <- Sys.glob(file.path(data.dir, "*", "*.bigwig"))
bigwig.file <- bigwig.file.vec[1]
chrom.ranges <- bigWigInfo(bigwig.file)
ranges.by.chrom <- split(chrom.ranges, chrom.ranges$chrom)
problems.by.chrom <- list()
for(chrom.i in seq_along(ranges.by.chrom)){
  chrom <- names(ranges.by.chrom)[[chrom.i]]
  message(sprintf("%4d / %4d chroms %s", chrom.i, length(ranges.by.chrom),
                  chrom))
  chrom.range <- ranges.by.chrom[[chrom]]
  all.chrom.problems <- with(chrom.range, {
    getProblems(chrom, chromStart, chromEnd, bases.per.problem,
                chrom.size=chromEnd)
  })
  ## Look at a bigwig file to see where the first chromStart and last
  ## chromEnd are, and then only process the problems which have some
  ## data.
  cmd <-
    sprintf("bigWigToBedGraph %s /dev/stdout -chrom=%s",
            bigwig.file, chrom)
  head.cmd <- paste(cmd, "| head -1")
  tail.cmd <- paste(cmd, "| tail -1")
  head.dt <- fread.or.null(head.cmd)
  tail.dt <- fread.or.null(tail.cmd)
  two <- rbind(head.dt, tail.dt)
  setnames(two, c("chrom", "chromStart", "chromEnd", "count"))
  first.chromStart <- two$chromStart[1]
  last.chromEnd <- two$chromEnd[2]
  problems.by.chrom[[chrom]] <- 
    all.chrom.problems[problemStart < last.chromEnd &
                         first.chromStart < problemEnd, ]
}
overlapping.problems <- do.call(rbind, problems.by.chrom)
overlapping.problems$job <- sort(rep(1:numJobs, l=nrow(overlapping.problems)))
table(overlapping.problems$job)
problems.by.job <- split(overlapping.problems, overlapping.problems$job)

trained.model.RData <- file.path(chunks.dir, "trained.model.RData")
save(train.errors, train.errors.picked,
     full.fit,
     ## for estimating test error later:
     problems.by.chunk,
     regions.by.chunk,
     ## for parallelizing prediction on jobs:
     problems.by.job,
     file=trained.model.RData)

