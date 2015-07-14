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
  ## TDH 12 June 2015 do not make the problems shorter.
  
  ##problemEnd[max.chromEnd < problemEnd] <- max.chromEnd
  problemStart[problemStart < 0] <- 0
  
  problem.name <- sprintf("%s:%d-%d", chrom, problemStart, problemEnd)
  problems <-
    data.table(problem.name, 
               bases.per.problem, problemStart, problemEnd)
  problems[min.chromStart < problemEnd &
             problemStart < max.chromEnd, ]
}
