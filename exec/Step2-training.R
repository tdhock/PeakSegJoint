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
argv <- file.path(ex.dir, c(1:3))

argv <- commandArgs(trailingOnly=TRUE)

print(argv)

ppn <- PPN.cores()
if(!is.na(ppn))options(mc.cores=ppn/2)

if(length(argv) == 0){
  stop("usage: Step2.R PeakSegJoint-chunks/1 [...]")
}

ann.colors <-
  c(noPeaks="#f6f4bf",
    peakStart="#ffafaf",
    peakEnd="#ff4c4c",
    peaks="#a445ee")

chunk.dir.vec <- normalizePath(argv, mustWork=TRUE)
chunks.dir <- dirname(chunk.dir.vec[1])
data.dir <- dirname(chunks.dir)
bigwig.file.vec <- Sys.glob(file.path(data.dir, "*", "*.bigwig"))

models.by.res <- list()
regions.by.chunk <- list()
chunks.list <- list()
for(chunk.dir in chunk.dir.vec){
  chunk.id <- basename(chunk.dir)
  bins.per.problem <- 500L
  chunk.id <- basename(chunk.dir)
  regions.RData <- file.path(chunk.dir, "regions.RData")
  objs <- load(regions.RData)
  regions$region.i <- 1:nrow(regions)
  chrom <- paste(regions$chrom[1])
  regions[, chromStart1 := chromStart + 1L]
  regions[, id.group := paste(sample.id, sample.group)]

  counts.by.sample <- list()
  for(bigwig.file in bigwig.file.vec){
    sample.counts <-
      readBigWig(bigwig.file, chunk$chrom, chunk$chunkStart, chunk$chunkEnd)
    sample.id <- sub("[.]bigwig$", "", basename(bigwig.file))
    sample.group <- basename(dirname(bigwig.file))
    counts.by.sample[[paste(sample.id, sample.group)]] <-
      data.table(sample.id, sample.group, sample.counts)
  }
  counts <- do.call(rbind, counts.by.sample)
  counts[, chromStart1 := chromStart + 1L]
  setkey(counts, chromStart1, chromEnd)

  sample.problems.dt <- do.call(rbind, problems.by.res)

  ## ggplot()+
  ##   ggtitle(regions.RData)+
  ##   geom_tallrect(aes(xmin=chromStart/1e3, xmax=chromEnd/1e3, fill=annotation),
  ##                 color="grey",
  ##                 alpha=0.5,
  ##                 data=regions)+
  ##   geom_segment(aes(problemStart/1e3, seq_along(bases.per.problem),
  ##                    xend=problemEnd/1e3, yend=seq_along(bases.per.problem)),
  ##                data=data.table(sample.problems.dt,
  ##                                sample.group="problems", sample.id="problems"))+
  ##   scale_fill_manual(values=ann.colors)+
  ##   geom_rect(aes(xmin=chromStart/1e3, xmax=chromEnd/1e3,
  ##                 ymin=0, ymax=count),
  ##             fill="grey50",
  ##             data=counts)+
  ##   theme_bw()+
  ##   theme(panel.margin=grid::unit(0, "cm"))+
  ##   facet_grid(sample.group + sample.id ~ ., scales="free")

  SampleProblems <- function(problem.i){
    cat(sprintf(
      "%4d / %4d problems %s\n",
      problem.i, nrow(sample.problems.dt),
      chunk.dir))
    problem <- sample.problems.dt[problem.i, ]
    bases.per.bin <- as.integer(problem$bases.per.problem/bins.per.problem)
    problem.name <- paste(problem$problem.name)
    problem.regions <-
      regions[! (chromEnd < problem$problemStart |
                   problem$problemEnd < chromStart), ]
    setkey(problem.regions, id.group)
    
    ## ggplot()+
    ##   ggtitle(problem.name)+
    ##   geom_tallrect(aes(xmin=chromStart/1e3, xmax=chromEnd/1e3, fill=annotation),
    ##                 color="grey",
    ##                 alpha=0.5,
    ##                 data=regions)+
    ##   geom_segment(aes(problemStart/1e3, 0,
    ##                    xend=problemEnd/1e3, yend=0),
    ##                data=problem)+
    ##   scale_fill_manual(values=ann.colors)+
    ##   geom_rect(aes(xmin=chromStart/1e3, xmax=chromEnd/1e3,
    ##                 ymin=0, ymax=count),
    ##             fill="grey50",
    ##             data=counts)+
    ##   theme_bw()+
    ##   theme(panel.margin=grid::unit(0, "cm"))+
    ##   facet_grid(sample.group + sample.id ~ ., scales="free")

    models.by.sample <- list()
    for(id.group in names(counts.by.sample)){
      sample.counts <- counts.by.sample[[id.group]]
      problem.counts <-
        sample.counts[! (chromEnd < problem$problemStart |
                           problem$problemEnd < chromStart), ]
      start <- as.integer(problem$problemStart)
      model <- segmentBins(
        problem.counts, start, bases.per.bin, bins.per.problem)
      model$meta <- data.table(
        problem.counts[1, .(sample.id, sample.group)],
        problem,
        chunk.id)
      sample.regions <- problem.regions[id.group]
      if(!is.na(sample.regions$sample.id[1])){
        
        ## ggplot()+
        ##   ggtitle(problem.name)+
        ##   geom_tallrect(aes(xmin=chromStart/1e3, xmax=chromEnd/1e3,
        ##                     fill=annotation),
        ##                 color="grey",
        ##                 alpha=0.5,
        ##                 data=sample.regions)+
        ##   geom_segment(aes(problemStart/1e3, -1,
        ##                    xend=problemEnd/1e3, yend=-1),
        ##                size=2,
        ##                data=problem)+
        ##   scale_fill_manual(values=ann.colors)+
        ##   geom_rect(aes(xmin=chromStart/1e3, xmax=chromEnd/1e3,
        ##                 ymin=0, ymax=count),
        ##             fill="grey50",
        ##             data=problem.counts)+
        ##   geom_line(aes((chromStart+chromEnd)/2e3, mean),
        ##             data=model$bins)+
        ##   theme_bw()+
        ##   theme(panel.margin=grid::unit(0, "cm"))+
        ##   facet_grid(sample.group + sample.id ~ ., scales="free")
        
        model$fit$error$incorrect.regions <- NA
        for(peaks.str in names(model$fit$peaks)){
          peaks <- model$fit$peaks[[peaks.str]]
          error.regions <- PeakErrorChrom(peaks, sample.regions)
          model$fit$error[peaks.str, "incorrect.regions"] <-
            with(error.regions, sum(fp+fn))
        }
        model$modelSelection$incorrect.regions <-
          model$fit$error[paste(model$modelSelection$peaks), "incorrect.regions"]
        target.indices <- with(model$modelSelection, largestContinuousMinimum(
          incorrect.regions, max.log.lambda-min.log.lambda))
        model$target <- with(model$modelSelection, c(
          min.log.lambda[target.indices$start],
          max.log.lambda[target.indices$end]))
      }
      models.by.sample[[id.group]] <- model
    }#id.group
    models.by.sample
  }#problem.i
  models.by.problem <-
    mclapply.or.stop(seq_along(sample.problems.dt$problem.name), SampleProblems)

  regions.by.chunk[[chunk.id]] <- regions
  chunks.list[[chunk.id]] <- chunk
  for(problem.i in 1:nrow(sample.problems.dt)){
    problem <- sample.problems.dt[problem.i,]
    res.str <- paste(problem$bases.per.problem)
    models.by.res[[res.str]][[paste(chunk.id, problem$problem.name)]] <-
      models.by.problem[[problem.i]]
  }
}

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
    chunk <- chunks.list[[chunk.str]]
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
    over.dt[, problem.name := sprintf(
      "%s:%d-%d", chunk$chrom, problemStart, problemEnd)]
    over.dt[, region.name := sprintf(
      "%s:%d-%d", chunk$chrom, chromStart, chromEnd)]
    ## Edit assignment of labels to problems -- there should be only
    ## one segmentation problem with a peak for each peakStart/peakEnd
    ## label.
    problems.by.chromStart <- split(over.dt, over.dt$chromStart)
    assigned.by.chromStart <- list()
    for(start.str in names(problems.by.chromStart)){
      region.problems <- problems.by.chromStart[[start.str]]
      multiple.problems <- any(1 < table(region.problems$sample.id))
      peakStart <- any(region.problems$annotation=="peakStart")
      peakEnd <- any(region.problems$annotation=="peakEnd")
      if(multiple.problems){
        if(peakStart){
          region.problems[problemEnd < chromEnd, annotation := "noPeaks"]
        }
        if(peakEnd){
          region.problems[chromStart < problemStart, annotation := "noPeaks"]
        }
      }
      assigned.by.chromStart[[start.str]] <- region.problems
    }
    problems.by.chunk[[chunk.str]] <- do.call(rbind, assigned.by.chromStart)
    if(FALSE){
      bins.list <- bins.by.chunk.list[[chunk.str]]
      bins <- do.call(rbind, bins.list)
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
      data.table(res.str,
                 problem.list[[chunk.str]])
  }
}

joint.by.res <- list()
labeled.by.chunk <- list()
for(chunk.str in names(problems.by.chunk)){
  joint.problems.by.res <- problems.by.chunk[[chunk.str]]
  labeled.problems <- do.call(rbind, joint.problems.by.res)
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
    show.problems <- data.table(labeled.problems)
    show.problems[, sample.id := res.str]
    show.problems[, sample.group := res.str]
    ggplot()+
      geom_segment(aes(problemStart/1e3, cluster,
                       xend=problemEnd/1e3, yend=cluster),
                   data=show.problems)+
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
  uniq.problems <- unique(
    labeled.problems[, .(res.str, problem.name, problemStart, problemEnd)])
  setkey(labeled.problems, res.str, problem.name)
  SegmentJoint <- function(row.i){
    ##print(row.i)
    problem <- uniq.problems[row.i, ]
    problem[, problemStart1 := problemStart + 1L]
    setkey(problem, res.str, problem.name)
    problem.regions <- labeled.problems[problem]
    setkey(problem, problemStart1, problemEnd)
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
    ##lapply(uniq.problems[, 1:.N], SegmentJoint)
    mclapply.or.stop(uniq.problems[, 1:.N], SegmentJoint)
  for(problem.i in 1:nrow(uniq.problems)){
    problem <- uniq.problems[problem.i, ]
    joint.by.res[[problem$res.str]][[problem$problem.name]] <- 
      labeled.model.list[[problem.i]]
  }
  labeled.by.chunk[[chunk.str]] <- data.table(chunk.str, labeled.problems)
}
all.labeled.problems <- do.call(rbind, labeled.by.chunk)
setkey(all.labeled.problems, res.str)

ResError <- function(res.str){
  set.seed(1)
  res.problems <- all.labeled.problems[res.str]
  joint.by.problem <- joint.by.res[[res.str]]
  train.by.problem <- list()
  oracle.by.chunk <- list()
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
    min.err <- subset(info$error$modelSelection, errors==min(errors))
    simplest <- subset(min.err, peaks==min(peaks))
    oracle.by.chunk[[problem$chunk.str]][[problem$problem.name]] <-
      subset(info$peaks, peaks==simplest$peaks)
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
  pred.peaks.by.chunk <- list()
  for(chunk.str in names(peaks.by.chunk)){
    peaks.by.problem <- oracle.by.chunk[[chunk.str]]
    peaks.by.problem <- peaks.by.chunk[[chunk.str]]
    chunk.peaks <- do.call(rbind, peaks.by.problem)
    chunk.regions <- regions.by.chunk[[chunk.str]]
    error.regions <- PeakErrorSamples(chunk.peaks, chunk.regions)
    if(FALSE){# plot errors.
      chunk <- chunks.list[[chunk.str]]
      labeled.problems <- labeled.by.chunk[[chunk.str]]
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
      ggplot()+
        geom_tallrect(aes(xmin=chromStart/1e3, xmax=chromEnd/1e3,
                          fill=annotation),
                      alpha=0.5,
                      color="grey",
                      data=chunk.regions)+
        scale_y_continuous("aligned read coverage",
                           breaks=function(limits){
                             floor(limits[2])
                           })+
        scale_linetype_manual("error type",
                              limits=c("correct", 
                                "false negative",
                                "false positive"
                                       ),
                              values=c(correct=0,
                                "false negative"=3,
                                "false positive"=1))+
        scale_x_continuous(paste("position on", chunk$chrom,
                                 "(kilo bases = kb)"))+
        geom_tallrect(aes(xmin=chromStart/1e3, xmax=chromEnd/1e3,
                          linetype=status),
                      fill=NA,
                      size=1,
                      color="black",
                      data=error.regions)+
        scale_fill_manual(values=ann.colors)+
        geom_rect(aes(xmin=chromStart/1e3, xmax=chromEnd/1e3,
                      ymin=0, ymax=count),
                  color="grey50",
                  data=counts)+
        geom_segment(aes(chromStart/1e3, 0,
                         xend=chromEnd/1e3, yend=0),
                     data=chunk.peaks,
                     color="deepskyblue",
                     size=2)+
        theme_bw()+
        theme(panel.margin=grid::unit(0, "lines"))+
        facet_grid(sample.group + sample.id ~ ., scales="free")
    }
    bases.per.problem <- as.integer(res.str)
    pred.peaks.by.chunk[[chunk.str]] <- data.table(
      bases.per.problem, chunk.str, chunk.peaks)
    error.by.chunk[[chunk.str]] <- data.table(
      bases.per.problem, chunk.str, error.regions)
  }#for(chunk.str
  error <- do.call(rbind, error.by.chunk)
  list(error=error,
       peaks=do.call(rbind, pred.peaks.by.chunk),
       fit=joint.fit)
}
res.vec <- unique(all.labeled.problems$res.str)
res.result.list <- mclapply.or.stop(res.vec, ResError)
names(res.result.list) <- res.vec
res.error.list <- lapply(res.result.list, "[[", "error")
res.error.regions <- do.call(rbind, res.error.list)
res.peaks.list <- lapply(res.result.list, "[[", "peaks")
res.peaks <- do.call(rbind, res.peaks.list)
res.error <- res.error.regions[, data.table(
  fp=sum(fp),
  fn=sum(fn),
  errors=sum(fp+fn),
  regions=length(fp)), by=bases.per.problem]
print(train.errors <- res.error[order(bases.per.problem),])
min.errors <- res.error[errors==min(errors),]
train.errors.picked <- min.errors[which.max(bases.per.problem),]
cat("Picked the following problem size:\n")
print(train.errors.picked)
res.chunk.error <-
  res.error.regions[bases.per.problem==train.errors.picked$bases.per.problem,
                    data.table(
                      fp=sum(fp),
                      fn=sum(fn),
                      errors=sum(fp+fn),
                      regions=length(fp)
                      ), by=chunk.str]
(chunks.ordered <- res.chunk.error[order(-errors), ])
joint.fit <-
  res.result.list[[paste(train.errors.picked$bases.per.problem)]]$fit
separate.fit <-
  separate.by.res[[paste(train.errors.picked$bases.per.problem)]]$fit
resCurve <- 
  ggplot()+
    geom_line(aes(bases.per.problem, errors),
              data=train.errors)+
    scale_x_log10()+
    ylab("incorrect labels (train error)")+
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
xt <- xtable(chunks.ordered)
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

## get chrom size info from a bigwig, so we can generate a list of
## segmentation problems and divide them into jobs.
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
    getProblems(chrom, chromStart, chromEnd,
                train.errors.picked$bases.per.problem,
                overlap.count=2L,
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
load(file.path(data.dir, "jobs.RData"))
overlapping.problems$job <- sort(rep(1:n.jobs, l=nrow(overlapping.problems)))
table(overlapping.problems$job)
problems.by.job <- split(overlapping.problems, overlapping.problems$job)

trained.model.RData <- file.path(chunks.dir, "trained.model.RData")
save(train.errors, train.errors.picked,
     separate.fit, joint.fit,
     res.peaks,
     res.error.regions,
     ## for parallelizing prediction on jobs:
     problems.by.job,
     file=trained.model.RData)

