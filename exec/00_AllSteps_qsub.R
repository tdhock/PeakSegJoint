library(PeakSegJoint)

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
  c(training=paste(Rscript, Step2, chunks.dir))

Step3e <-
  system.file(file.path("exec", "Step3e-estimate-test-error.R"),
              mustWork=TRUE,
              package="PeakSegJoint")

bigwig.file.vec <- Sys.glob(file.path(data.dir, "*", "*.bigwig"))
chrom.ranges <- bigWigInfo(bigwig.file.vec[1])
Step3 <-
  system.file(file.path("exec", "Step3-segment-one-chrom.R"),
              mustWork=TRUE,
              package="PeakSegJoint")
trained.model.RData <- file.path(chunks.dir, "trained.model.RData")

## We want Step3-predict and Step3e-estimate-test-error to start
## immediate after Step2-training, which is what the next line does:
cmd.list$Step3 <-
  c(structure(paste(Rscript, Step3e, trained.model.RData),
              names="testError"),
    structure(paste(Rscript, Step3, trained.model.RData, chrom.ranges$chrom),
              names=paste0(chrom.ranges$chrom, "predict")))

Step4 <-
  system.file(file.path("exec", "Step4-write-bed-files.R"),
              mustWork=TRUE,
              package="PeakSegJoint")
Step4v <-
  system.file(file.path("exec", "Step4v-viz-one-labeled-chunk.R"),
              mustWork=TRUE,
              package="PeakSegJoint")
## TODO: start Step4v immediately after Step3e finishes!
pred.dir <- file.path(data.dir, "PeakSegJoint-predictions")
cmd.list$Step4 <-
  c(structure(paste(Rscript, Step4, pred.dir),
              names="bed"),
    structure(paste(Rscript, Step4v, chunk.dir.vec),
              names=paste0("chunk", basename(chunk.dir.vec), "viz")))

qsub <- Sys.getenv("QSUB")
if(qsub == ""){
  qsub <- "qsub"
}

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
    is.viz <- grepl("viz", cmd.name)
    cmd <- cmd.vec[[cmd.name]]
    is.prediction <- grepl(Step3, cmd)
    walltime <- if(step.name == "Step3"){
      "10:00:00"
    }else{
      "01:00:00"
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
    if(!is.viz){
      depend.list[[cmd.name]] <- qsub.id
    }
  }
}

