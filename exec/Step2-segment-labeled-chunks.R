require(data.table)
require(PeakSegJoint)
require(PeakError)
require(parallel)

## Compute PeakSegJoint segmentations for one labeled chunk in the
## train data set.

argv <-
  system.file(file.path("exampleData",
                        "PeakSegJoint-chunks",
                        "3"),
              package="PeakSegDP")

argv <- "~/exampleData/PeakSegJoint-chunks/1"
argv <- "~/projects/PeakSegJoint-paper/PeakSegJoint-chunks/H3K36me3_AM_immune/21"
argv <- "~/projects/PeakSegJoint-paper/PeakSegJoint-chunks/H3K4me3_PGP_immune/2"

argv <- commandArgs(trailingOnly=TRUE)

print(argv)

if(length(argv) != 1){
  stop("usage: Step2.R path/to/PeakSegJoint-chunks/012354")
}

options(mc.cores=detectCores())
ppn <- as.integer(Sys.getenv("PBS_NUM_PPN"))
if(!is.finite(ppn))ppn <- 2
options(mc.cores=ppn)
print(options("mc.cores"))

chunk.dir <- argv[1]
chunk.id <- basename(chunk.dir)

regions.RData <- file.path(chunk.dir, "regions.RData")

objs <- load(regions.RData)
regions$region.i <- 1:nrow(regions)
chrom <- paste(regions$chrom[1])

counts.RData.vec <- Sys.glob(file.path(chunk.dir, "*", "*.RData"))

counts.by.sample <- list()
for(counts.RData.path in counts.RData.vec){
  objs <- load(counts.RData.path)
  sample.id <- sub(".RData$", "", basename(counts.RData.path))
  counts.by.sample[[sample.id]] <- data.table(sample.id, counts)
}
counts <- do.call(rbind, counts.by.sample)
setkey(counts, chromStart, chromEnd)

step1.problems.dt <- do.call(rbind, problems.by.res)
Step1Problem <- function(problem.i){
  problem <- step1.problems.dt[problem.i, ]
  problem.name <- paste(problem$problem.name)
  problem.counts <-
    counts[! (chromEnd < problem$problemStart |
                problem$problemEnd < chromStart), ]
  profile.list <- ProfileList(problem.counts)
  fit <- PeakSegJointSeveral(profile.list)
  models <- fit$models[-1]
  loss.vec <- sapply(models, "[[", "loss")
  is.feasible <- is.finite(loss.vec)
  peak.mat <- sapply(models, "[[", "peak_start_end")
  rownames(peak.mat) <- c("chromStart", "chromEnd")  
  feasible.mat <- peak.mat[, is.feasible]
  if(ncol(feasible.mat) > 0){
    peak.vec <- feasible.mat[, ncol(feasible.mat)]
    data.table(problem, t(peak.vec))
  }
}
step1.results.list <- mclapply(1:nrow(step1.problems.dt), Step1Problem)
step1.results <- do.call(rbind, step1.results.list)
setkey(step1.results, chromStart, chromEnd)

step1.by.res <- split(step1.results, step1.results$bases.per.problem)
Step1Step2 <- function(res.str){
  bases.per.problem <- as.integer(res.str)
  step1.peaks <- step1.by.res[[res.str]]
  step2.overlap.list <- list()
  clustered.peaks <- clusterPeaks(step1.peaks)
  peaks.by.cluster <- split(clustered.peaks, clustered.peaks$cluster)
  for(cluster.name in names(peaks.by.cluster)){
    cluster <- peaks.by.cluster[[cluster.name]]
    merged.peak <- with(cluster, {
      data.frame(chromStart=min(chromStart),
                 chromEnd=max(chromEnd))
    })
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
    problemStart <- as.integer(cluster.mid - half.bases)
    if(cluster.chromStart < problemStart){
      problemStart <- as.integer(cluster.chromStart - half.bases) #old
      problemStart <- as.integer(cluster.chromStart - half.bases/2) 
    }
    if(problemStart < chromEnd.before){
      problemStart <-
        as.integer((chromEnd.before+cluster.chromStart)/2)
    }
    problemEnd <- as.integer(cluster.mid + half.bases)
    if(problemEnd < cluster.chromEnd){
      problemEnd <- as.integer(cluster.chromEnd + half.bases) #old
      problemEnd <- as.integer(cluster.chromEnd + half.bases/2)
    }
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

  problems.dt <- data.table(bases.per.problem, step2.problems)
  setkey(problems.dt, problemStart, problemEnd)
  setkey(regions, chromStart, chromEnd)
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
step2.data.list <- mclapply(names(step1.by.res), Step1Step2)
regions.by.problem <-
  do.call(c, lapply(step2.data.list, "[[", "regions"))
step2.problems <-
  do.call(rbind, lapply(step2.data.list, "[[", "problems"))
## Careful to set names after, so that the previous lapply calls get
## correct names.
names(step2.data.list) <- names(step1.by.res)

SegmentStep2 <- function(row.i){
  problem <- step2.problems[row.i, ]
  setkey(problem, problemStart, problemEnd)
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
  loss <- info$loss
  loss$cummin <- cummin(loss$loss)
  cummin.reduced <- c(TRUE, diff(loss$cummin) < 0)
  some.loss <- loss[cummin.reduced, ]
  info$modelSelection <- with(some.loss, {
    exactModelSelection(loss, peaks, peaks)
  })
  info
}
step2.model.list <- mclapply(1:nrow(step2.problems), SegmentStep2)
names(step2.model.list) <- step2.problems$problem.name

setkey(step2.problems, problem.name)
ProblemError <- function(problem.name){
  problem <- step2.problems[problem.name]
  problem.regions <- regions.by.problem[[problem.name]]
  converted <- step2.model.list[[problem.name]]
  pred.peaks <- if(!is.null(converted)){
    prob.err.list <- PeakSegJointError(converted, problem.regions)
    best.models <-
      subset(prob.err.list$modelSelection, errors==min(errors))
    peaks.num <- min(best.models$peaks)
    if(peaks.num > 0){
      subset(converted$peaks, peaks == peaks.num)
    }
  }
  list(problem=prob.err.list,
       peaks=pred.peaks)
}
step2.error.list <- mclapply(names(regions.by.problem), ProblemError)
names(step2.error.list) <- names(regions.by.problem)

ResError <- function(res.str){
  res.data <- step2.data.list[[res.str]]
  problems.with.regions <-
    res.data$problems[problem.name %in% names(step2.error.list), ]
  pred.peaks.list <- list()
  for(problem.name in problems.with.regions$problem.name){
    pred.peaks.list[[problem.name]] <- step2.error.list[[problem.name]]$peaks
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
res.error.list <- mclapply(names(step2.data.list), ResError)
res.error <- do.call(rbind, res.error.list)
print(res.error)

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
