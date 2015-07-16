getProblems <- function
### Tile problems over an entire chromosome.
(chrom,
### chrom for problem names.
 min.chromStart,
### First base of problems.
 max.chromEnd,
### Last base of problems.
 bases.per.problem,
### each problem has this many bases.
 overlap.count=3
### how many problems overlap at any given genome position?
 ){
  stopifnot(is.character(chrom))
  stopifnot(length(chrom) == 1)
  stopifnot(is.numeric(min.chromStart))
  stopifnot(length(min.chromStart) == 1)
  stopifnot(is.numeric(max.chromEnd))
  stopifnot(length(max.chromEnd) == 1)
  stopifnot(is.numeric(bases.per.problem))
  stopifnot(length(bases.per.problem) == 1)
  stopifnot(is.numeric(overlap.count))
  stopifnot(length(overlap.count) == 1)
  
  problemSeq <- seq(0, max.chromEnd, by=bases.per.problem)
  
  ## mix small and big problems for step 1?
  bigSize <- bases.per.problem*1.5
  bigSeq <- seq(0, max.chromEnd, by=bigSize)
  bigStart <- as.integer(sort(c(bigSeq, bigSeq + bases.per.problem)))
  bigEnd <- bigStart + bigSize
  smallEnd <- problemSeq + bases.per.problem
  problemStart <- c(problemSeq, bigStart)
  problemEnd <- c(smallEnd, bigEnd)
  
  ## overlapping problems of the same size?
  overlap.list <- list()
  for(overlap.i in 1:overlap.count){
    offset <- bases.per.problem*(overlap.i-1)/overlap.count
    overlap.list[[overlap.i]] <- problemSeq + offset
  }
  problemStart <-
    as.integer(sort(do.call(c, overlap.list)))
  problemEnd <- problemStart + bases.per.problem

  problemEnd[max.chromEnd < problemEnd] <- max.chromEnd
  problemStart[problemStart < 0] <- 0
  
  problem.name <- sprintf("%s:%d-%d", chrom, problemStart, problemEnd)
  problems <-
    data.table(chrom, problem.name, 
               bases.per.problem, problemStart, problemEnd)
  problems[min.chromStart < problemEnd &
             problemStart < max.chromEnd, ]
}

clusterProblems <- function
### Compute a set of non-overlapping segmentation problems.
(step1.peaks
### Maximum number of peaks from overlapping segmentation problems.
 ){
  step2.overlap.list <- list()
  bases.per.problem <- step1.peaks$bases.per.problem[[1]]
  stopifnot(step1.peaks$bases.per.problem == bases.per.problem)
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
  step2.problems
### data.table of non-overlapping segmentation problems.
}

