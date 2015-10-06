library(PeakSegJoint)

getenv.or <- function(env.var, default){
  env.value <- Sys.getenv(env.var)
  if(env.value == ""){
    default
  }else{
    env.value
  }
}

## Some arbitrary parameters that affect how long (and how much
## embarrassing paralellelism) the computation will take.
n.jobs <- as.integer(getenv.or("JOBS", 1000))
qsub <- getenv.or("QSUB", "qsub")
n.chunk.order.seeds <- 2 # for estimating test error.

## Make and run qsub scripts for all steps of the PeakSegJoint pipeline.
R.bin <- R.home("bin")
Rscript <- file.path(R.bin, "Rscript")
argv <- # interactive default for debugging.
  system.file("exampleData", "manually_annotated_region_labels.txt",
              package="PeakSegJoint")
argv <- # interactive default for debugging.
  system.file("exampleData",
              c("other_labels.txt", "manually_annotated_region_labels.txt"),
              package="PeakSegJoint")
argv <- commandArgs(trailingOnly=TRUE)
if(length(argv) == 0){
  stop("usage: AllSteps_qsub.R path/to/labels.txt
where there are path/to/*/*.bigwig files")
}
##NEEDS TO BE ABSOLUTE!
labels.file.vec <- normalizePath(argv, mustWork=TRUE)
data.dir <- dirname(labels.file.vec[1])

## Step0 generates some files for every chunk, which we need to
## examine to make the step2 commands, so we run it interactively.
Step0 <-
  system.file(file.path("exec", "Step0-convert-labels.R"),
              mustWork=TRUE,
              package="PeakSegJoint")
labels.files <- paste(labels.file.vec, collapse=" ")
cmd <- paste(Rscript, Step0, labels.files)
status <- system(cmd)
if(status != 0){
  stop("error in Step0, most likely problem with labels")
}
load(file.path(data.dir, "chunk.file.map.RData"))

Step1 <-
  system.file(file.path("exec", "Step1-segment-one-labeled-chunk.R"),
              mustWork=TRUE,
              package="PeakSegJoint")
regions.RData.vec <-
  Sys.glob(file.path(data.dir, "PeakSegJoint-chunks", "*", "regions.RData"))
chunk.dir.vec <- dirname(regions.RData.vec)
problems.RData.vec <- file.path(chunk.dir.vec, "problems.RData")

Step2 <- 
  system.file(file.path("exec", "Step2-training.R"),
              mustWork=TRUE,
              package="PeakSegJoint")
chunks.dir <- file.path(data.dir, "PeakSegJoint-chunks")
trained.model.RData <- file.path(chunks.dir, "trained.model.RData")

## Step3e and Step4v make test error estimates and visualizations
## using the labeled data.
Step3e <-
  system.file(file.path("exec", "Step3e-estimate-test-error.R"),
              mustWork=TRUE,
              package="PeakSegJoint")
test.error.params <-
  expand.grid(outer.fold=1:outer.folds,
              chunk.order.seed=1:n.chunk.order.seeds)
Step4e <-
  system.file(file.path("exec", "Step4e-plot-test-error.R"),
              mustWork=TRUE,
              package="PeakSegJoint")
oJob.dir <- file.path(data.dir, "PeakSegJoint-overlapping")
dir.create(oJob.dir, showWarnings=FALSE)
Step4v <-
  system.file(file.path("exec", "Step4v-viz-one-labeled-chunk.R"),
              mustWork=TRUE,
              package="PeakSegJoint")

job.vec <- 1:n.jobs
Step3 <-
  system.file(file.path("exec", "Step3-overlapping-problems.R"),
              mustWork=TRUE,
              package="PeakSegJoint")
Step4 <-
  system.file(file.path("exec", "Step4-combine-overlapping.R"),
              mustWork=TRUE,
              package="PeakSegJoint")
Step5 <-
  system.file(file.path("exec", "Step5-final-problems.R"),
              mustWork=TRUE,
              package="PeakSegJoint")
pred.dir <- file.path(data.dir, "PeakSegJoint-predictions")
dir.create(pred.dir, showWarnings=FALSE)
combined.problems.RData <- file.path(data.dir, "combined.problems.RData")
Step6 <-
  system.file(file.path("exec", "Step6-write-bed-files.R"),
              mustWork=TRUE,
              package="PeakSegJoint")

step <- function(step.name, walltime, jobs, depends=NULL){
  list(step.name=step.name,
       walltime=walltime,
       jobs=jobs,
       depends=depends)
}

job <- function(command.line, name, produces){
  prefix <- sub("[.][^.]*$", "", produces)
  data.table(command.line, name, prefix)
}

cmd.list <-
  list(step("Step1", "02:00:00",
            job(paste(Rscript, Step1, chunk.dir.vec),
                paste0("chunk", basename(chunk.dir.vec)),
                problems.RData.vec)),
       step("Step2", "02:00:00",
            job(paste(Rscript, Step2, chunks.dir, n.jobs),
                "training",
                trained.model.RData),
            "Step1"),
       step("Step3e", "20:00:00",
            with(test.error.params, job(
              paste(Rscript, Step3e, trained.model.RData,
                    chunk.order.seed, outer.fold),
              paste0("seed", chunk.order.seed, "fold", outer.fold),
              file.path(chunks.dir,
                        "figure-test-errors",
                        paste0("seed", chunk.order.seed,
                               "fold", outer.fold,
                               ".RData")))),
            "Step2"),
       step("Step4e", "02:00:00",
            job(paste(Rscript, Step4e, chunks.dir),
                "vizTest",
                file.path(chunks.dir, "figure-test-errors",
                          "test.metrics.curves.RData")),
            "Step3e"),
       step("Step4v", "02:00:00",
            job(paste(Rscript, Step4v, chunk.dir.vec),
                paste0("chunk", basename(chunk.dir.vec), "viz"),
                file.path(chunk.dir.vec, "figure-train-errors.png")),
            "Step3e"),
       step("Step3", "20:00:00",
            job(paste(Rscript, Step3, trained.model.RData, job.vec),
                paste0("oJob", job.vec),
                file.path(oJob.dir, paste0(job.vec, ".RData"))),
            "Step2"),
       step("Step4", "02:00:00",
            job(paste(Rscript, Step4, oJob.dir, n.jobs),
                "combine",
                file.path(data.dir, "combined.problems.RData")),
            "Step3"),
       step("Step5", "20:00:00",
            job(paste(Rscript, Step5, combined.problems.RData, job.vec),
                paste0("finalJob", job.vec),
                file.path(pred.dir, paste0(job.vec, ".RData"))),
            "Step4"),
       step("Step6", "02:00:00",
            job(paste(Rscript, Step6, pred.dir),
                "bed",
                file.path(data.dir, "PeakSegJoint.predictions.RData")),
            "Step5"))

depend.list <- list()
for(step.list in cmd.list){
  depend.txt <- if(length(step.list$depends)==1 && qsub == "qsub"){
    depend.vec <- depend.list[[step.list$depends]]
    pid.txt <- paste(depend.vec, collapse=":")
    paste0("-W depend=afterok:", pid.txt)
  }else{
    ""
  }
  for(cmd.i in 1:nrow(step.list$jobs)){
    cmd.row <- step.list$jobs[cmd.i, ]
    script.txt <-
      paste0("#!/bin/bash
#PBS -l nodes=1:ppn=4
#PBS -l walltime=", step.list$walltime, "
#PBS -A bws-221-ae
#PBS -m ae
#PBS -M tdhock5@gmail.com
#PBS -o ", cmd.row$prefix, ".out
#PBS -e ", cmd.row$prefix, ".err
#PBS -V                                        
#PBS -N ", cmd.row$name, "\n", cmd.row$command.line, "\n")
    script.file <- paste0(cmd.row$prefix, ".sh")
    script.dir <- dirname(script.file)
    dir.create(script.dir, showWarnings=FALSE, recursive=TRUE)
    cat(script.txt, file=script.file)
    qsub.cmd <- paste(qsub, depend.txt, script.file)
    qsub.out <- system(qsub.cmd, intern=TRUE)
    status.code <- attr(qsub.out, "status")
    ## status.code is NULL if qsub.cmd exited with status 0.
    if(length(status.code) == 1 && status.code != 0){
      stop(qsub.cmd, " exited with status ", status.code)
    }
    qsub.id <- sub("[.].*", "", qsub.out[1])
    cat(step.list$step.name, " ",
        cmd.row$name, " ",
        "submitted as job ",
        qsub.id, "\n", sep="")
    depend.list[[step.list$step.name]][cmd.row$name] <- qsub.id
  }
}

