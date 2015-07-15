## Make and run qsub scripts for all steps of the PeakSegJoint pipeline.
R.bin <- R.home("bin")
Rscript <- file.path(R.bin, "Rscript")
labels.txt.file <- # interactive default for debugging.
  system.file(file.path("exampleData", "manually_annotated_region_labels.txt"),
              package="PeakSegDP")
argv <- commandArgs(trailingOnly=TRUE)
if(length(argv) != 1){
  stop("usage: AllSteps_qsub.R path/to/labels.txt
where there are path/to/*/*.bigwig files")
}
labels.txt.file <- normalizePath(argv[1], mustWork=TRUE)
data.dir <- dirname(labels.txt.file)

## Step0 generates some files for every chunk, which we need to
## examine to make the step2 commands, so we run it interactively.
Step0 <-
  system.file(file.path("exec", "Step0-convert-labels.R"),
              mustWork=TRUE,
              package="PeakSegJoint")
cmd <- paste(Rscript, Step0, labels.txt.file)
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
  c(training=paste(Rscript, Step2, chunks.dir))

Step3v <-
  system.file(file.path("exec", "Step3v-viz-train-errors.R"),
              mustWork=TRUE,
              package="PeakSegJoint")

cmd.list$Step3 <-
  structure(paste(Rscript, Step3v, chunk.dir.vec),
            names=paste0("chunk", basename(chunk.dir.vec), "viz"))

qsub <- "echo 1 && bash"
qsub <- "qsub"

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
    last.file <- sub(".* ", "", cmd)
    last.base <- basename(last.file)
    last.dir <- dirname(last.file)
    prefix.only <- sub("[.].*?$", "", last.base)
    prefix.nostep <- file.path(last.dir, prefix.only)
    prefix <- paste0(prefix.nostep, "-", step.name)
    script.txt <-
      paste0("#!/bin/bash
#PBS -l nodes=1:ppn=5
#PBS -l walltime=02:00:00
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
    qsub.id <- sub("[.].*", "", qsub.out)
    cat(step.name, " ",
        cmd.name, " ",
        "submitted as job ",
        qsub.id, "\n", sep="")
    if(!grepl("viz", cmd.name)){
      depend.list[[cmd.name]] <- qsub.id
    }
  }
}

problems.by.job <- split(test.problems, test.problems$job.name)
Step4 <-
  system.file(file.path("exec", "Step4-test-segmentation.R"),
              mustWork=TRUE,
              package="PeakSegJoint")
R.bin <- R.home("bin")
Rscript <- file.path(R.bin, "Rscript")
for(job.name in names(problems.by.job)){
  job.problems <- problems.by.job[[job.name]]
  cmd <- paste("qsub", Rscript, Step4, trained.model.RData, job.name)
  print(cmd)
}
