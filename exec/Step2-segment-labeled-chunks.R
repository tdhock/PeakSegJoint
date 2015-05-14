require(data.table)
require(ggplot2)
require(xtable)
require(PeakSegJoint)
require(PeakError)

## Compute PeakSegJoint segmentations for one labeled chunk in the
## train data set.

argv <-
  system.file(file.path("exampleData",
                        "PeakSegJoint-chunks",
                        "3"),
              package="PeakSegDP")

argv <- "~/exampleData/PeakSegJoint-chunks/1"
argv <- "PeakSegJoint-chunks/H3K36me3_AM_immune/21"

argv <- commandArgs(trailingOnly=TRUE)

print(argv)

if(length(argv) != 1){
  stop("usage: Step2.R path/to/PeakSegJoint-chunks/012354")
}

chunk.dir <- argv[1]
chunk.id <- basename(chunk.dir)

regions.RData <- file.path(chunk.dir, "regions.RData")

objs <- load(regions.RData)
regions$region.i <- 1:nrow(regions)
chrom <- paste(regions$chrom[1])

ann.colors <-
  c(noPeaks="#f6f4bf",
    peakStart="#ffafaf",
    peakEnd="#ff4c4c",
    peaks="#a445ee")

counts.RData.vec <- Sys.glob(file.path(chunk.dir, "*", "*.RData"))

counts.by.sample <- list()
for(counts.RData.path in counts.RData.vec){
  objs <- load(counts.RData.path)
  sample.id <- sub(".RData$", "", basename(counts.RData.path))
  counts.by.sample[[sample.id]] <- data.table(sample.id, counts)
}
counts <- do.call(rbind, counts.by.sample)
setkey(counts, chromStart, chromEnd)

step2.error.list <- list()
step2.by.res <- list()
for(res.str in names(problems.by.res)){
  bases.per.problem <- as.integer(res.str)
  res.problems <- problems.by.res[[res.str]]

  step1.peaks.list <- list()
  for(problem.i in 1:nrow(res.problems)){
    problem <- res.problems[problem.i, ]
    problem.name <- paste(problem$problem.name)
    problem.counts <-
      counts[! (chromEnd < problem$problemStart |
                  problem$problemEnd < chromStart), ]
    profile.list <- ProfileList(problem.counts)
    problem.peaks <- tryCatch({
      fit <- PeakSegJointHeuristic(profile.list)
      ConvertModelList(fit)$peaks
    }, error=function(e){
      NULL
    })
    if(is.numeric(problem.peaks$peaks)){
      peaks.by.peaks <- split(problem.peaks, problem.peaks$peaks)
      peaks.num <- as.numeric(names(peaks.by.peaks))
      peaks.df <- peaks.by.peaks[[which.max(peaks.num)]]
      step1.peaks.list[[problem.name]] <-
        data.table(problem, peaks.df[1,])
    }
  }#problem.name
  if(length(step1.peaks.list) == 0){
    pred.peaks <- Peaks()
    step2.by.problem <- list()
  }else{
    step1.peaks <- do.call(rbind, step1.peaks.list)
    step2.overlap.list <- list()
    clustered.peaks <- clusterPeaks(step1.peaks)
    peaks.by.cluster <- split(clustered.peaks, clustered.peaks$cluster)
    pred.by.cluster <- list()
    for(cluster.name in names(peaks.by.cluster)){
      cluster <- peaks.by.cluster[[cluster.name]]
      merged.peak <- with(cluster, {
        data.frame(chromStart=min(chromStart),
                   chromEnd=max(chromEnd))
      })
      pred.by.cluster[[cluster.name]] <-
        data.frame(merged.peak,
                   sample.id=unique(cluster$sample.id))
      cluster.chromStart <- min(cluster$chromStart)
      cluster.chromEnd <- max(cluster$chromEnd)
      cluster.mid <-
        as.integer((cluster.chromEnd + cluster.chromStart)/2)
      half.bases <- as.integer(bases.per.problem/2)
      cluster.num <- as.numeric(cluster.name)
      before.name <- paste(cluster.num-1)
      chromEnd.before <- if(before.name %in% names(peaks.by.cluster)){
        max(peaks.by.cluster[[before.name]]$chromEnd)
      }else{
        0
      }
      after.name <- paste(cluster.num+1)
      chromStart.after <- if(after.name %in% names(peaks.by.cluster)){
        min(peaks.by.cluster[[after.name]]$chromStart)
      }else{
        Inf
      }
      problemStart <- as.integer(cluster.chromStart - half.bases)
      if(problemStart < chromEnd.before){
        problemStart <-
          as.integer((chromEnd.before+cluster.chromStart)/2)
      }
      problemEnd <- as.integer(cluster.chromEnd + half.bases)
      if(chromStart.after < problemEnd){
        problemEnd <- as.integer((chromStart.after+cluster.chromEnd)/2)
      }
      stopifnot(problemStart <= cluster.chromStart)
      stopifnot(cluster.chromEnd <= problemEnd)
      problem.i <- as.numeric(cluster.name)+1
      step2.overlap.list[[problem.i]] <-
        data.frame(problem.i,
                   problem.name=sprintf("%s:%d-%d",
                     chrom, problemStart, problemEnd),
                   problemStart, problemEnd,
                   peakStart=merged.peak$chromStart,
                   peakEnd=merged.peak$chromEnd)
    }
    
    step2.overlap <- do.call(rbind, step2.overlap.list)
    step2.problems <- with(step2.overlap, {
      prev.problemEnd <- problemEnd[-length(problemEnd)]
      next.problemStart <- problemStart[-1]
      overlaps.next <- which(next.problemStart <= prev.problemEnd)
      mid <- as.integer((prev.problemEnd+next.problemStart)/2)
      problemEnd[overlaps.next] <- mid[overlaps.next]
      problemStart[overlaps.next+1] <- mid[overlaps.next]+1L
      data.frame(problem.i=seq_along(problemStart),
                 problem.name=sprintf("%s:%d-%d",
                   chrom, problemStart, problemEnd),
                 problemStart, problemEnd)
    })
    stopifnot(with(step2.problems, {
      problemEnd[-length(problemEnd)] < problemStart[-1]
    }))
    
    ## ggplot()+
    ##   geom_segment(aes(problemStart/1e3, problem.i,
    ##                    color=what,
    ##                    xend=problemEnd/1e3, yend=problem.i),
    ##                data=data.frame(step2.overlap, what="overlap"))+
    ##   geom_segment(aes(problemStart/1e3, problem.i,
    ##                    color=what,
    ##                    xend=problemEnd/1e3, yend=problem.i),
    ##                data=data.frame(step2.problems, what="corrected"))

    problems.dt <- data.table(step2.problems)
    setkey(problems.dt, problemStart, problemEnd)
    setkey(regions, chromStart, chromEnd)
    over.regions <- foverlaps(regions, problems.dt, nomatch=0L)
    over.regions[,
                 `:=`(overlapStart=ifelse(problemStart < chromStart,
                        chromStart, problemStart),
                      overlapEnd=ifelse(problemEnd < chromEnd,
                        problemEnd, chromEnd))]
    over.regions[, overlapBases := overlapEnd-overlapStart]
    region.i.problems <-
      over.regions[,
                   .(problem.name=problem.name[which.max(overlapBases)]),
                   by=region.i]
    stopifnot(nrow(region.i.problems) == nrow(regions))
    setkey(regions, region.i)
    setkey(region.i.problems, region.i)
    assigned.regions <- regions[region.i.problems,]
    stopifnot(nrow(assigned.regions) == nrow(regions))
    regions.by.problem <-
      split(assigned.regions, assigned.regions$problem.name, drop=TRUE)
    setkey(problems.dt, problem.name)
    peaks.by.problem <- list()
    step2.by.problem <- list()
    step2.peak.list <- list()
    saved.problems.list <- list()
    saved.i <- 1
    for(problem.name in names(regions.by.problem)){
      problem.regions <- regions.by.problem[[problem.name]]
      problem <- problems.dt[problem.name]
      problem$saved.i <- saved.i
      saved.i <- saved.i+1
      saved.problems.list[[problem.name]] <-
        data.frame(problem, sample.id="step 2")
      setkey(problem, problemStart, problemEnd)
      problem.i <- problem$problem.i
      problem.counts <-
        foverlaps(counts, problem, nomatch=0L, type="within")

      ## ggplot()+
      ##   theme_bw()+
      ##   theme(panel.margin=grid::unit(0, "cm"))+
      ##   facet_grid(sample.id ~ ., labeller=function(var, val){
      ##     sub("McGill0", "", sub(" ", "\n", val))
      ##   }, scales="free")+
      ##   geom_step(aes(chromStart/1e3, count),
      ##             data=problem.counts,
      ##             color="grey50")
      
      tryCatch({
        profile.list <- ProfileList(problem.counts)
        fit <- PeakSegJointHeuristic(profile.list)
        converted <- ConvertModelList(fit)
        prob.err.list <- PeakSegJointError(converted, problem.regions)
        step2.by.problem[[problem.name]] <-
          list(converted=converted,
               error=prob.err.list,
               features=featureMatrix(profile.list))
        best.models <-
          subset(prob.err.list$modelSelection, errors==min(errors))
        peaks.num <- min(best.models$peaks)
        if(peaks.num > 0){
          show.peaks <- subset(converted$peaks, peaks == peaks.num)
          peaks.by.problem[[problem.name]] <- show.peaks
          peak.row <- show.peaks[1,]
          peak.row$sample.id <- "step 2"
          step2.peak.list[[problem.name]] <-
            data.frame(problem, peak.row)
        }
      }, error=function(e){
        paste(e) #ignore.
      })
    }#problem.name
    saved.problems <- do.call(rbind, saved.problems.list)
    pred.peaks <- do.call(rbind, peaks.by.problem)
    step2.peaks <- do.call(rbind, step2.peak.list)
    ## ggplot()+
    ##   geom_point(aes(chromStart/1e3, sample.id),
    ##              data=pred.peaks,
    ##              pch=1)+
    ##   geom_segment(aes(chromStart/1e3, sample.id,
    ##                    xend=chromEnd/1e3, yend=sample.id),
    ##                data=pred.peaks)
  }#if any peaks
  step2.by.res[[res.str]] <- step2.by.problem
  if(is.null(pred.peaks))pred.peaks <- Peaks()
  error.regions <- PeakErrorSamples(pred.peaks, regions)

  limits <- regions[, .(chromStart=min(chromStart),
                        chromEnd=max(chromEnd))]
  limits[, bases := chromEnd - chromStart ]
  limits[, expand := as.integer(bases/10) ]
  limits[, min := chromStart - expand]
  limits[, max := chromEnd + expand]
  lim.vec <- with(limits, c(min, max))/1e3
  some.counts <- counts[limits$min < chromEnd &
                          chromStart < limits$max,]

  label.df <-
    data.frame(problemEnd=max(step1.peaks$problemEnd),
               problem.i=1,
               label=paste(res.str, "bases/problem"),
               sample.id="step 1")

  step1.peaks$sample.id <- "step 1"
  step1.peaks$problem.i <- 1:nrow(step1.peaks)

  problemPlot <- 
    ggplot()+
      ggtitle(chunk.dir)+
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
      scale_x_continuous("position on chromosome (kilo bases = kb)")+
      coord_cartesian(xlim=lim.vec)+
      geom_tallrect(aes(xmin=chromStart/1e3, xmax=chromEnd/1e3,
                        fill=annotation),
                    alpha=0.5,
                    color="grey",
                    data=regions)+
      scale_fill_manual(values=ann.colors)+
      geom_text(aes(problemEnd/1e3, problem.i, label=label),
                data=label.df,
                size=3,
                vjust=0,
                hjust=1)+
      theme_bw()+
      theme(panel.margin=grid::unit(0, "cm"))+
      facet_grid(sample.id ~ ., labeller=function(var, val){
        sub("McGill0", "", sub(" ", "\n", val))
      }, scales="free")+
      geom_step(aes(chromStart/1e3, count),
                data=some.counts,
                color="grey50")+
      geom_tallrect(aes(xmin=chromStart/1e3, xmax=chromEnd/1e3,
                        linetype=status),
                    color="black",
                    fill=NA,
                    size=1,
                    data=error.regions)+
      geom_segment(aes(problemStart/1e3, problem.i,
                       xend=problemEnd/1e3, yend=problem.i),
                   data=step1.peaks)+
      geom_segment(aes(problemStart/1e3, saved.i,
                       xend=problemEnd/1e3, yend=saved.i),
                   data=saved.problems)
  if(is.data.frame(pred.peaks) && nrow(pred.peaks) > 0){
    problemPlot <- problemPlot+
      geom_segment(aes(chromStart/1e3, 0,
                       xend=chromEnd/1e3, yend=0),
                   size=2,
                   color="deepskyblue",
                   data=pred.peaks)
  }
  if(is.data.frame(step2.peaks) && nrow(step2.peaks) > 0){
    problemPlot <- problemPlot+
      geom_segment(aes(chromStart/1e3, saved.i,
                       xend=chromEnd/1e3, yend=saved.i),
                   size=2,
                   color="deepskyblue",
                   data=step2.peaks)
  }
  if(is.data.frame(step1.peaks) && nrow(step1.peaks) > 0){
    problemPlot <- problemPlot+
      geom_segment(aes(chromStart/1e3, problem.i,
                       xend=chromEnd/1e3, yend=problem.i),
                   size=2,
                   color="deepskyblue",
                   data=step1.peaks)
  }
  png.base <- 
    file.path(chunk.dir, "figure-train-errors", res.str)
  png.name <- paste0(png.base, ".png")
  
  png.dir <- dirname(png.name)
  dir.create(png.dir, showWarnings=FALSE, recursive=TRUE)

  png(png.name, width=14, h=10, res=100, units="in")
  print(problemPlot)
  dev.off()

  thumb <- paste0(png.base, "-thumb.png")
  cmd <- sprintf("convert %s -resize 230 %s", png.name, thumb)
  system(cmd)
  
  step2.error.list[[res.str]] <- 
    with(error.regions, {
      data.frame(chunk.id,
                 bases.per.problem,
                 fp=sum(fp),
                 fn=sum(fn),
                 errors=sum(fp+fn),
                 regions=length(fp))
    })
}#bases.per.problem
step2.error <- do.call(rbind, step2.error.list)
print(step2.error)
print(table(regions$annotation))

xcols <- c("bases.per.problem", "fp", "fn", "errors", "regions")
bpp <- step2.error$bases.per.problem
htable <-
  data.frame(thumb=sprintf('<a href="%s">
  <img src="%s" />
</a>', paste0(bpp, ".png"),
               paste0(bpp, "-thumb.png")),
             step2.error[, xcols])
xt <- xtable(htable)
print(xt, type="html", include.rownames=FALSE, sanitize.text.function=identity,
      file=file.path(chunk.dir, "figure-train-errors", "index.html"))

## Save results for this chunk/resolution.
problems.RData <- file.path(chunk.dir, "problems.RData")
save(step2.error, step2.by.res,
     file=problems.RData)
