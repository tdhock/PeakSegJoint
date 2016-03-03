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
regions.by.chunk <- list()
chunks.list <- list()
for(models.RData in models.RData.vec){
  objs <- load(models.RData)
  chunk.dir <- dirname(models.RData)
  chunk.id <- basename(chunk.dir)
  regions.RData <- file.path(chunk.dir, "regions.RData")
  objs <- load(regions.RData)
  regions.by.chunk[[chunk.id]] <- regions
  chunks.list[[chunk.id]] <- chunk
  for(problem.i in 1:nrow(sample.problems.dt)){
    problem <- sample.problems.dt[problem.i,]
    res.str <- paste(problem$bases.per.problem)
    models.by.res[[res.str]][[paste(chunk.id, problem$problem.name)]] <-
      models.by.problem[[problem.i]]
  }
}

chunks.dir <- dirname(dirname(models.RData))
data.dir <- dirname(chunks.dir)
bigwig.file.vec <- Sys.glob(file.path(data.dir, "*", "*.bigwig"))

SeparatePeaks <- function(res.str, seed=1){
  set.seed(1)
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
  peaks.by.chunk.list <- list()
  bins.by.chunk.list <- list()
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
      chunk.str <- paste(model$meta$chunk.id)
      bins.by.chunk.list[[chunk.str]][[paste(chunk.problem, id.group)]] <-
        data.table(model$meta, model$bins)
      if(0 < selected$peaks){
        peaks.str <- paste(selected$peaks)
        pred.peaks <- model$fit$peaks[[peaks.str]]
        peak.cols <- c(
          "sample.id", "sample.group",
          "chromStart", "chromEnd")
        peaks.by.chunk.list[[chunk.str]][[paste(chunk.problem, id.group)]] <- 
          data.table(model$meta, pred.peaks)[, peak.cols, with=FALSE]
      }#if
    }#for(id.group
  }#for(chunk.problem
  problems.by.chunk <- list()
  for(chunk.str in names(peaks.by.chunk.list)){
    peaks.list <- peaks.by.chunk.list[[chunk.str]]
    regions <- regions.by.chunk[[chunk.str]]
    prop.noPeaks <- regions[, list(
      prop=mean(annotation=="noPeaks")
      ), by=.(chromStart, chromEnd)]
    prop.noPeaks[, bases := chromEnd - chromStart]
    peaks.list$joint <- prop.noPeaks[prop==1, {
      data.table(sample.id="joint", sample.group="joint",
                 chromStart=as.integer(chromStart+bases/3),
                 chromEnd=as.integer(chromEnd-bases/3))
    }]
    peaks <- do.call(rbind, peaks.list)
    clustered <- clusterPeaks(peaks)
    joint.problems <- clustered2problems(clustered)
    setkey(joint.problems, problemStart, problemEnd)
    setkey(regions, chromStart, chromEnd)
    over.dt <- foverlaps(joint.problems, regions, nomatch=0L)
    joint.problems[, overlaps.regions := ifelse(
      cluster %in% over.dt$cluster, "some", "none")]
    problems.by.chunk[[chunk.str]] <- joint.problems
    if(FALSE){
      bins <- do.call(rbind, bins.list)
      bins.list <- bins.by.chunk.list[[chunk.str]]
      problems <- unique(bins[, .(problemStart, problemEnd)])
      problems[, problem.i := 1:.N]
      ggplot()+
        geom_tallrect(aes(xmin=chromStart/1e3, xmax=chromEnd/1e3,
                         fill=annotation),
                      alpha=0.5,
                      color="grey",
                      data=regions)+
        scale_fill_manual(values=ann.colors)+
        geom_point(aes((chromStart+chromEnd)/2e3, mean),
                   color="grey50",
                   shape=1,
                   data=bins)+
        theme_bw()+
        theme(panel.margin=grid::unit(0, "lines"))+
        facet_grid(sample.group + sample.id ~ ., scales="free")+
        geom_segment(aes(problemStart/1e3, problem.i,
                         xend=problemEnd/1e3, yend=problem.i),
                     data=data.table(
                       sample.id="problems",
                       sample.group="problems",
                       problems))+
        geom_segment(aes(problemStart/1e3, cluster,
                         xend=problemEnd/1e3, yend=cluster),
                     data=data.table(
                       sample.id="joint",
                       sample.group="joint",
                       joint.problems))+
        geom_segment(aes(chromStart/1e3, 0,
                         xend=chromEnd/1e3, yend=0),
                     color="deepskyblue",
                     data=peaks)+
        geom_point(aes(chromStart/1e3, 0),
                   color="deepskyblue",
                   data=peaks)
    }#if(FALSE
  }#for(chunk.str
  list(problems=problems.by.chunk,
       fit=fit)
}
separate.by.res <- mclapply.or.stop(names(models.by.res), SeparatePeaks)
names(separate.by.res) <- names(models.by.res)

problems.by.chunk <- list()
for(res.str in names(separate.by.res)){
  problem.list <- separate.by.res[[res.str]]$problems
  for(chunk.str in names(problem.list)){
    problems.by.chunk[[chunk.str]][[res.str]] <-
      data.table(sample.id=res.str, sample.group=res.str,
                 problem.list[[chunk.str]])
  }
}

joint.by.problem <- list()
joint.by.chunk <- list()
labeled.by.chunk <- list()
for(chunk.str in names(problems.by.chunk)){
  joint.problems.by.res <- problems.by.chunk[[chunk.str]]
  regions <- regions.by.chunk[[chunk.str]]
  chunk.problems <- do.call(rbind, joint.problems.by.res)
  labeled.problems <- chunk.problems[overlaps.regions=="some",]
  zoom <- labeled.problems[, list(
    zoomStart=min(problemStart),
    zoomEnd=max(problemEnd))]
  counts.by.sample <- list()
  for(bigwig.file in bigwig.file.vec){
    sample.counts <-
      readBigWig(bigwig.file, chunk$chrom, zoom$zoomStart, zoom$zoomEnd)
    sample.id <- sub("[.]bigwig$", "", basename(bigwig.file))
    sample.group <- basename(dirname(bigwig.file))
    counts.by.sample[[paste(sample.id, sample.group)]] <-
      data.table(sample.id, sample.group, sample.counts)
  }
  counts <- do.call(rbind, counts.by.sample)
  setkey(counts, chromStart, chromEnd)
  if(FALSE){
    ggplot()+
      geom_segment(aes(problemStart/1e3, cluster,
                       color=overlaps.regions,
                       xend=problemEnd/1e3, yend=cluster),
                   data=chunk.problems)+
      geom_tallrect(aes(xmin=chromStart/1e3, xmax=chromEnd/1e3,
                        fill=annotation),
                    alpha=0.5,
                    color="grey",
                    data=regions)+
      scale_fill_manual(values=ann.colors)+
      geom_rect(aes(xmin=chromStart/1e3, xmax=chromEnd/1e3,
                    ymin=0, ymax=count),
                color="grey50",
                data=counts)+
      theme_bw()+
      theme(panel.margin=grid::unit(0, "lines"))+
      facet_grid(sample.group + sample.id ~ ., scales="free")
  }
  labeled.problems[, problem.name := sprintf(
    "%s:%d-%d", chunk$chrom, problemStart, problemEnd)]
  uniq.problems <- unique(
    labeled.problems[, .(problem.name, problemStart, problemEnd)])
  SegmentJoint <- function(row.i){
    ##print(row.i)
    problem <- uniq.problems[row.i, ]
    problem[, problemStart1 := problemStart + 1L]
    setkey(problem, problemStart1, problemEnd)
    problem.regions <- foverlaps(problem, regions, nomatch=0L)
    problem.counts <-
      foverlaps(problem, counts, nomatch=0L, type="any")
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
    seg.counts <- problem.counts[chromStart < chromEnd,]
    if(FALSE){
      ggplot()+
        theme_bw()+
        theme(panel.margin=grid::unit(0, "cm"))+
        facet_grid(sample.id ~ ., scales="free")+
        geom_segment(aes((problemStart+0.5)/1e3, 0,
                         xend=problemEnd/1e3, yend=0),
                     data=problem,
                     color="black",
                     size=1)+
        geom_rect(aes(xmin=chromStart/1e3, ymin=0,
                      xmax=chromEnd/1e3, ymax=count),
                  data=seg.counts,
                  color="grey50")
    }
    profile.list <- ProfileList(seg.counts)
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
    info$error <- PeakSegJointError(info, problem.regions)
    info
  }
  labeled.model.list <-
    mclapply.or.stop(uniq.problems[, 1:.N], SegmentJoint)
  ## It is OK to index the model list on problem name (even though the
  ## same problem could occur in several resolutions), since anyways the
  ## model should not change between resolutions.
  names(labeled.model.list) <- uniq.problems$problem.name
  joint.by.problem[uniq.problems$problem.name] <- labeled.model.list
  stopifnot(table(names(labeled.model.list)) == 1)
  joint.by.chunk[[chunk.str]] <- labeled.model.list
  labeled.by.chunk[[chunk.str]] <- data.table(chunk.str, labeled.problems)
}
all.labeled.problems <- do.call(rbind, labeled.by.chunk)
setkey(all.labeled.problems, sample.id)

ResError <- function(res.str){
  res.problems <- all.labeled.problems[res.str]
  train.by.problem <- list()
  for(problem.name in res.problems$problem.name){
    info <- joint.by.problem[[problem.name]]
    train.by.problem[[problem.name]] <- list(
      features=info$features,
      target=info$error$target)
  }
  joint.fit <- IntervalRegressionProblemsCV(train.by.problem)
  peaks.by.chunk <- list()
  for(problem.i in 1:nrow(res.problems)){
    problem <- res.problems[problem.i,]
    info <- joint.by.problem[[problem$problem.name]]
    pred.log.lambda <- joint.fit$predict(info$features)
    pred.row <- subset(
      info$error$modelSelection,
      min.log.lambda < pred.log.lambda &
        pred.log.lambda < max.log.lambda)
    stopifnot(nrow(pred.row)==1)
    peaks.by.chunk[[problem$chunk.str]][[problem$problem.name]] <- 
      subset(info$peaks, peaks==pred.row$peaks)
  }
  error.by.chunk <- list()
  for(chunk.str in names(peaks.by.chunk)){
    peaks.by.problem <- peaks.by.chunk[[chunk.str]]
    chunk.peaks <- do.call(rbind, peaks.by.problem)
    chunk.regions <- regions.by.chunk[[chunk.str]]
    error.by.chunk[[chunk.str]] <- PeakErrorSamples(chunk.peaks, chunk.regions)
  }
  error <- do.call(rbind, error.by.chunk)
  with(error, {
    data.table(res.str,
               bases.per.problem=as.integer(res.str),
               fp=sum(fp),
               fn=sum(fn),
               errors=sum(fp+fn),
               regions=length(fp))
  })
}
res.vec <- unique(all.labeled.problems$sample.id)
res.error.list <- mclapply.or.stop(res.vec, ResError)
res.error <- do.call(rbind, res.error.list)
print(res.error)
stop("check peaks")

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

