if(!require(data.table))
  install.packages("data.table")
if(!require(PeakSegJoint))
  devtools::install_github("tdhock/PeakSegJoint")
if(!require(PeakError))
  devtools::install_github("tdhock/PeakError")

## Compute PeakSegJoint segmentations for one labeled chunk in the
## train data set.

argv <-
  system.file(file.path("exampleData",
                        "PeakSegJoint-chunks",
                        "3"),
              package="PeakSegDP")

argv <- "~/exampleData/PeakSegJoint-chunks/1"

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

  peaks.by.problem <- list()
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
      peaks.by.problem[[problem.name]] <- peaks.df[1,]
    }
  }#problem.name
  if(length(peaks.by.problem) == 0){
    pred.peaks <- Peaks()
    step2.by.problem <- list()
  }else{
    step1.peaks <- do.call(rbind, peaks.by.problem)
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
    saved.problem.list <- list()
    for(problem.name in names(regions.by.problem)){
      problem.regions <- regions.by.problem[[problem.name]]
      problem <- problems.dt[problem.name]
      setkey(problem, problemStart, problemEnd)
      problem.i <- problem$problem.i
      problem.counts <-
        foverlaps(counts, problem, nomatch=0L, type="within")
      
      tryCatch({
        profile.list <- ProfileList(problem.counts)
        fit <- PeakSegJointHeuristic(profile.list)
        converted <- ConvertModelList(fit)
        prob.err.list <- PeakSegJointError(converted, problem.regions)
        step2.by.problem[[problem.name]] <-
          list(converted=converted,
               error=prob.err.list,
               features=featureMatrix(profile.list))
        saved.problem.list[[problem.name]] <-
          data.frame(problem, sample.id="step 2")
        best.models <-
          subset(prob.err.list$modelSelection, errors==min(errors))
        peaks.num <- min(best.models$peaks)
        if(peaks.num > 0){
          show.peaks <- subset(converted$peaks, peaks == peaks.num)
          peaks.by.problem[[problem.name]] <- show.peaks
          peak.row <- show.peaks[1,]
          peak.row$sample.id <- "step 2"
          step2.peak.list[[problem.name]] <-
            data.frame(problem.i, peak.row)
        }
      }, error=function(e){
        paste(e) #ignore.
      })
    }#problem.name
    pred.peaks <- do.call(rbind, peaks.by.problem)
    step2.peaks <- do.call(rbind, step2.peak.list)
    saved.problems <- do.call(rbind, saved.problem.list)
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

## Save results for this chunk/resolution.
problems.RData <- file.path(chunk.dir, "problems.RData")
save(step2.error, step2.by.res,
     file=problems.RData)
