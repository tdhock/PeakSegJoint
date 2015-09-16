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
n.jobs <- as.integer(getenv.or("JOBS", 400))
short.time <- "02:00:00"
long.time <- "100:00:00"
qsub <- getenv.or("QSUB", "qsub")

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

## Starting with Step1, we add vectors of commands that will wait for
## each other.
cmd.list <- list()

data.dir <- normalizePath(data.dir, mustWork=TRUE)#NEEDS TO BE ABSOLUTE!

Step1 <-
  system.file(file.path("exec", "Step1-segment-one-labeled-chunk.R"),
              mustWork=TRUE,
              package="PeakSegJoint")
regions.RData.vec <-
  Sys.glob(file.path(data.dir, "PeakSegJoint-chunks", "*", "regions.RData"))
chunk.dir.vec <- dirname(regions.RData.vec)
cmd.list$Step1 <-
  structure(paste(Rscript, Step1, chunk.dir.vec),
            names=paste0("chunk", basename(chunk.dir.vec)))

Step2 <- 
  system.file(file.path("exec", "Step2-training.R"),
              mustWork=TRUE,
              package="PeakSegJoint")
chunks.dir <- file.path(data.dir, "PeakSegJoint-chunks")
cmd.list$Step2 <-
  c(training=paste(Rscript, Step2, chunks.dir, n.jobs))

Step3e <-
  system.file(file.path("exec", "Step3e-estimate-test-error.R"),
              mustWork=TRUE,
              package="PeakSegJoint")
Step3 <-
  system.file(file.path("exec", "Step3-overlapping-problems.R"),
              mustWork=TRUE,
              package="PeakSegJoint")
trained.model.RData <- file.path(chunks.dir, "trained.model.RData")
## We want Step3-overlapping-problems and Step3e-estimate-test-error
## to start immediately after Step2-training, which is what the next
## line does:
job.vec <- 1:n.jobs
cmd.list$Step3 <-
  c(structure(paste(Rscript, Step3e, trained.model.RData),
              names="testJob"),
    structure(paste(Rscript, Step3, trained.model.RData, job.vec),
              names=paste0("oJob", job.vec)))

Step4 <-
  system.file(file.path("exec", "Step4-combine-overlapping.R"),
              mustWork=TRUE,
              package="PeakSegJoint")
oJob.dir <- file.path(data.dir, "PeakSegJoint-overlapping")
dir.create(oJob.dir, showWarnings=FALSE)
cmd.list$Step4 <-
 structure(paste(Rscript, Step4, oJob.dir, n.jobs),
           names="combine")

Step5 <-
  system.file(file.path("exec", "Step5-final-problems.R"),
              mustWork=TRUE,
              package="PeakSegJoint")
Step5v <-
  system.file(file.path("exec", "Step5v-viz-one-labeled-chunk.R"),
              mustWork=TRUE,
              package="PeakSegJoint")
pred.dir <- file.path(data.dir, "PeakSegJoint-predictions")
dir.create(pred.dir, showWarnings=FALSE)
combined.problems.RData <- file.path(data.dir, "combined.problems.RData")
cmd.list$Step5 <-
  c(structure(paste(Rscript, Step5, combined.problems.RData, job.vec),
              names=paste0("finalJob", job.vec)),
    structure(paste(Rscript, Step5v, chunk.dir.vec),
              names=paste0("chunk", basename(chunk.dir.vec), "viz")))

Step6 <-
  system.file(file.path("exec", "Step6-write-bed-files.R"),
              mustWork=TRUE,
              package="PeakSegJoint")
cmd.list$Step6 <-
  c(structure(paste(Rscript, Step6, pred.dir),
              names="bed"))

depend.list <- list()
for(step.name in names(cmd.list)){
  depend.txt <- if(length(depend.list)==0){
    ""
  }else{
    depend.vec <- do.call(c, depend.list)
    pid.txt <- paste(depend.vec, collapse=":")
    paste0("\n#PBS -W depend=afterok:", pid.txt)
  }
  depend.list <- list()
  cmd.vec <- cmd.list[[step.name]]
  for(cmd.name in names(cmd.vec)){
    cmd <- cmd.vec[[cmd.name]]
    walltime <- if(grepl("Job", cmd.name)){
      long.time
    }else{
      short.time
    }
    last.args <- sub(".*[.]R ", "", cmd)
    last.file <- sub(" ", "-", last.args)
    prefix <- paste0(last.file, "-", step.name)
    script.txt <-
      paste0("#!/bin/bash
#PBS -l nodes=1:ppn=4
#PBS -l walltime=", walltime, "
#PBS -A bws-221-ae", depend.txt, "
#PBS -m ae
#PBS -M tdhock5@gmail.com
#PBS -o ", prefix, ".out
#PBS -e ", prefix, ".err
#PBS -V                                        
#PBS -N ", cmd.name, "\n", cmd, "\n")
    script.file <- paste0(prefix, ".sh")
    cat(script.txt, file=script.file)
    qsub.cmd <- paste(qsub, script.file)
    qsub.out <- system(qsub.cmd, intern=TRUE)
    qsub.id <- sub("[.].*", "", qsub.out[1])
    cat(step.name, " ",
        cmd.name, " ",
        "submitted as job ",
        qsub.id, "\n", sep="")
    is.viz <- grepl("viz", cmd.name)
    if(!is.viz){
      depend.list[[cmd.name]] <- qsub.id
    }
  }
}

