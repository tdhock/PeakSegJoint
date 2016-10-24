problem.joint <- function
### Fit a joint model.
(jointProblem.dir
### path/to/jointProblem
){
  segmentations.RData <- file.path(jointProblem.dir, "segmentations.RData")
  problem.bed <- file.path(jointProblem.dir, "problem.bed")
  problem <- fread(problem.bed)
  setnames(problem, c("chrom",  "problemStart", "problemEnd", "problem.name"))
  problem[, problemStart1 := problemStart + 1L]
  setkey(problem, chrom, problemStart1, problemEnd)
  jointProblems <- dirname(jointProblem.dir)
  prob.dir <- dirname(jointProblems)
  probs.dir <- dirname(prob.dir)
  data.dir <- dirname(probs.dir)
  samples.dir <- file.path(data.dir, "samples")
  coverage.bedGraph.vec <- Sys.glob(file.path(
    samples.dir, "*", "*", "problems",
    problem$problem.name, "coverage.bedGraph"))
  cat("Found", length(coverage.bedGraph.vec), "samples to jointly segment.\n")
  coverage.list <- list()
  for(coverage.i in seq_along(coverage.bedGraph.vec)){
    coverage.bedGraph <- coverage.bedGraph.vec[[coverage.i]]
    ## cat(sprintf(
    ##   "%4d / %4d %s\n",
    ##   coverage.i, length(coverage.bedGraph.vec), coverage.bedGraph))
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
    sample.path <- paste0(sample.group, "/", sample.id)
    coverage.list[[sample.path]] <- data.table(
      sample.id, sample.group, problem.coverage)
  }
  coverage <- do.call(rbind, coverage.list)
  setkey(coverage, sample.id, chrom, chromStart, chromEnd)
  profile.list <- ProfileList(coverage)
  fit <- PeakSegJointSeveral(coverage)
  segmentations <- ConvertModelList(fit)
  segmentations$features <- featureMatrix(profile.list)
  cat("Writing segmentation and features to", segmentations.RData, "\n")
  save(segmentations, file=segmentations.RData)
  segmentations$coverage <- coverage
  segmentations
### Model from ConvertModelList.
}

problem.joint.predict <- function
### Compute peak predictions for a joint problem.
(joint.model.RData,
### path/to/joint.model.RData file (which contains an object called joint.model)
 jointProblem.dir
### problem to predict peaks.
){
  load(joint.model.RData)
  converted <- problem.joint(jointProblem.dir)
  log.penalty <- joint.model$predict(converted$features)
  stopifnot(length(log.penalty)==1)
  selected <- subset(
    converted$modelSelection,
    min.log.lambda < log.penalty & log.penalty < max.log.lambda)
  loss.tsv <- file.path(jointProblem.dir, "loss.tsv")
  pred.dt <- if(selected$peaks == 0){
    unlink(loss.tsv)
    data.table()
  }else{
    selected.loss <- converted$loss[paste(selected$peaks), "loss"]
    flat.loss <- converted$loss["0", "loss"]
    loss.dt <- data.table(
      loss.diff=flat.loss-selected.loss)
    write.table(
      loss.dt,
      loss.tsv,
      quote=FALSE,
      sep="\t",
      col.names=FALSE,
      row.names=FALSE)
    pred.df <- subset(converted$peaks, peaks==selected$peaks)
    chrom <- paste(converted$coverage$chrom[1])
    with(pred.df, data.table(
      chrom,
      chromStart,
      chromEnd,
      name=paste0(sample.group, "/", sample.id),
      mean))
  }
  peaks.bed <- file.path(jointProblem.dir, "peaks.bed")
  cat("Writing ",
      nrow(pred.dt), " peaks to ",
      peaks.bed,
      "\n", sep="")
  write.table(
    pred.dt, peaks.bed,
    quote=FALSE,
    sep="\t",
    col.names=FALSE,
    row.names=FALSE)
### Nothing.
}

problem.joint.target <- function
### Compute target interval for a joint problem.
(jointProblem.dir
### Joint problem directory.
){
  converted <- problem.joint(jointProblem.dir)
  labels.bed <- file.path(jointProblem.dir, "labels.tsv")
  labels <- fread(labels.bed)
  setnames(labels, c(
    "chrom", "chromStart", "chromEnd", "annotation",
    "sample.id", "sample.group"))
  fit.error <- PeakSegJointError(converted, labels)
  if(FALSE){
    show.peaks <- 8
    show.peaks.df <- subset(converted$peaks, peaks==show.peaks)
    show.errors <- fit.error$error.regions[[paste(show.peaks)]]
    ann.colors <-
      c(noPeaks="#f6f4bf",
        peakStart="#ffafaf",
        peakEnd="#ff4c4c",
        peaks="#a445ee")
    ggplot()+
      theme_bw()+
      theme(panel.margin=grid::unit(0, "lines"))+
      facet_grid(sample.group + sample.id ~ ., scales="free")+
      scale_fill_manual(values=ann.colors)+
      geom_tallrect(aes(
        xmin=chromStart/1e3,
        xmax=chromEnd/1e3,
        fill=annotation),
        color="grey",
        alpha=0.5,
        data=labels)+
      geom_tallrect(aes(
        xmin=chromStart/1e3,
        xmax=chromEnd/1e3,
        linetype=status),
        color="black",
        size=1,
        fill=NA,
        data=show.errors)+
      scale_linetype_manual(
        "error type",
        limits=c("correct", 
                 "false negative",
                 "false positive"),
        values=c(correct=0,
                 "false negative"=3,
                 "false positive"=1))+
      geom_step(aes(chromStart/1e3, count),
                data=converted$coverage,
                color="grey50")+
      geom_segment(aes(chromStart/1e3, 0,
                       xend=chromEnd/1e3, yend=0),
                   data=show.peaks.df,
                   color="deepskyblue",
                   size=2)
    ## geom_segment(aes(chromStart/1e3, mean,
    ##                  xend=chromEnd/1e3, yend=mean),
    ##              data=
    ##              color="green")
  }
  cat("Train error:\n")
  data.table(fit.error$modelSelection)[, .(
    min.log.lambda, max.log.lambda, peaks, errors)]
  target.tsv <- file.path(jointProblem.dir, "target.tsv")
  cat(
    "Writing target interval (",
    paste(fit.error$target, collapse=", "),
    ") to ", 
    target.tsv,
    "\n", sep="")
  write(fit.error$target, target.tsv, sep="\t")
### Nothing.
}
