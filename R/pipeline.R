problem.joint <- function
### Fit a joint model.
(jointProblem.dir
### path/to/jointProblem
){
  segmentations.RData <- file.path(jointProblem.dir, "segmentations.RData")
  if(file.exists(segmentations.RData)){
    load(segmentations.RData)
  }else{
    problem.bed <- file.path(jointProblem.dir, "problem.bed")
    problem <- fread(problem.bed)
    setnames(problem, c("chrom",  "problemStart", "problemEnd", "problem.name"))
    problem[, problemStart1 := problemStart + 1L]
    setkey(problem, chrom, problemStart1, problemEnd)
    jointProblems <- dirname(jointProblem.dir)
    data.dir <- dirname(jointProblems)
    samples.dir <- file.path(data.dir, "samples")
    coverage.bedGraph.vec <- Sys.glob(file.path(
      samples.dir, "*", "*", "problems",
      problem$problem.name, "coverage.bedGraph"))
    coverage.list <- list()
    for(coverage.i in seq_along(coverage.bedGraph.vec)){
      coverage.bedGraph <- coverage.bedGraph.vec[[coverage.i]]
      cat(sprintf(
        "%4d / %4d %s\n",
        coverage.i, length(coverage.bedGraph.vec), coverage.bedGraph))
      sample.coverage <- fread(coverage.bedGraph)
      setnames(sample.coverage, c("chrom", "chromStart", "chromEnd", "count"))
      sample.coverage[, chromStart1 := chromStart + 1L]
      setkey(sample.coverage, chrom, chromStart1, chromEnd)
      problem.coverage <- foverlaps(sample.coverage, problem, nomatch=0L)
      problem.dir <- dirname(coverage.bedGraph)
      problems.dir <- dirname(problem.dir)
      sample.dir <- dirname(problems.dir)
      sample.id <- basename(sample.dir)
      group.dir <- dirname(sample.dir)
      sample.group <- basename(group.dir)
      coverage.list[[coverage.bedGraph]] <- data.table(
        sample.id, sample.group, problem.coverage)
    }
    coverage <- do.call(rbind, coverage.list)
    setkey(coverage, sample.id, chrom, chromStart, chromEnd)
    fit <- PeakSegJointSeveral(coverage)
    segmentations <- ConvertModelList(fit)
    save(segmentations, file=segmentations.RData)
  }
  segmentations
### Model from ConvertModelList.
}
