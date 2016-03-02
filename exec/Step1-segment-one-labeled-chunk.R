library(PeakSegJoint)
library(PeakError)

## Compute PeakSegJoint segmentations for one labeled chunk in the
## train data set.

argv <-
  system.file(file.path("exampleData",
                        "PeakSegJoint-chunks",
                        "1"),
              package="PeakSegJoint")

argv <- "~/exampleData/PeakSegJoint-chunks/3"
argv <- "~/projects/PeakSegJoint-paper/PeakSegJoint-chunks/H3K36me3_AM_immune/21"
argv <- "~/projects/PeakSegJoint-paper/PeakSegJoint-chunks/H3K4me3_PGP_immune/2"

argv <- commandArgs(trailingOnly=TRUE)

ppn <- PPN.cores()
if(!is.na(ppn))options(mc.cores=ppn/2)

print(argv)

if(length(argv) != 1){
  stop("usage: Step1.R path/to/PeakSegJoint-chunks/012354")
}

chunk.dir <- argv[1]
bins.per.problem <- 500L
chunk.id <- basename(chunk.dir)
chunks.dir <- dirname(chunk.dir)
data.dir <- dirname(chunks.dir)
bigwig.file.vec <- Sys.glob(file.path(data.dir, "*", "*.bigwig"))

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

ann.colors <-
  c(noPeaks="#f6f4bf",
    peakStart="#ffafaf",
    peakEnd="#ff4c4c",
    peaks="#a445ee")
step1.problems.dt <- do.call(rbind, problems.by.res)

ggplot()+
  ggtitle(regions.RData)+
  geom_tallrect(aes(xmin=chromStart/1e3, xmax=chromEnd/1e3, fill=annotation),
                color="grey",
                alpha=0.5,
                data=regions)+
  geom_segment(aes(problemStart/1e3, seq_along(bases.per.problem),
                   xend=problemEnd/1e3, yend=seq_along(bases.per.problem)),
               data=data.table(step1.problems.dt,
                               sample.group="problems", sample.id="problems"))+
  scale_fill_manual(values=ann.colors)+
  geom_rect(aes(xmin=chromStart/1e3, xmax=chromEnd/1e3,
                ymin=0, ymax=count),
            fill="grey50",
            data=counts)+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "cm"))+
  facet_grid(sample.group + sample.id ~ ., scales="free")

Step1Problem <- function(problem.i){
  problem <- step1.problems.dt[problem.i, ]
  bases.per.bin <- as.integer(problem$bases.per.problem/bins.per.problem)
  problem.name <- paste(problem$problem.name)
  problem.regions <-
    regions[! (chromEnd < problem$problemStart |
                 problem$problemEnd < chromStart), ]
  setkey(problem.regions, id.group)
  
  ggplot()+
    ggtitle(problem.name)+
    geom_tallrect(aes(xmin=chromStart/1e3, xmax=chromEnd/1e3, fill=annotation),
                  color="grey",
                  alpha=0.5,
                  data=regions)+
    geom_segment(aes(problemStart/1e3, 0,
                     xend=problemEnd/1e3, yend=0),
                 data=problem)+
    scale_fill_manual(values=ann.colors)+
    geom_rect(aes(xmin=chromStart/1e3, xmax=chromEnd/1e3,
                  ymin=0, ymax=count),
              fill="grey50",
              data=counts)+
    theme_bw()+
    theme(panel.margin=grid::unit(0, "cm"))+
    facet_grid(sample.group + sample.id ~ ., scales="free")

  for(id.group in names(counts.by.sample)){
    sample.regions <- problem.regions[id.group]
    if(nrow(sample.regions)){
      sample.counts <- counts.by.sample[[id.group]]
      problem.counts <-
        sample.counts[! (chromEnd < problem$problemStart |
                           problem$problemEnd < chromStart), ]
      start <- as.integer(problem$problemStart)
      model <- segmentBins(
        problem.counts, start, bases.per.bin, bins.per.problem)

      ggplot()+
        ggtitle(problem.name)+
        geom_tallrect(aes(xmin=chromStart/1e3, xmax=chromEnd/1e3,
                          fill=annotation),
                      color="grey",
                      alpha=0.5,
                      data=sample.regions)+
        geom_segment(aes(problemStart/1e3, -1,
                         xend=problemEnd/1e3, yend=-1),
                     size=2,
                     data=problem)+
        scale_fill_manual(values=ann.colors)+
        geom_rect(aes(xmin=chromStart/1e3, xmax=chromEnd/1e3,
                      ymin=0, ymax=count),
                  fill="grey50",
                  data=problem.counts)+
        geom_line(aes((chromStart+chromEnd)/2e3, mean),
                  data=model$bins)+
        theme_bw()+
        theme(panel.margin=grid::unit(0, "cm"))+
        facet_grid(sample.group + sample.id ~ ., scales="free")
      
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
    }#if(nrow(sample.regions))
  }#id.group
  stop("TODO return something")
}#problem.i
step1.results.list <-
  mclapply.or.stop(seq_along(step1.problems.dt$problem.name), Step1Problem)
if(all(sapply(step1.results.list, is.null))){
  print(step1.problems.dt)
  stop("no computable models for any uniform size segmentation problems")
}
step1.results <- do.call(rbind, step1.results.list)

## Plot step1 problems with detected peaks.
## ggplot()+
##   geom_segment(aes(problemStart/1e3, problem.name,
##                    xend=problemEnd/1e3, yend=problem.name),
##                data=step1.results)+
##   geom_segment(aes(chromStart/1e3, problem.name,
##                    xend=chromEnd/1e3, yend=problem.name),
##                size=2,
##                color="deepskyblue",
##                data=step1.results)+
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
