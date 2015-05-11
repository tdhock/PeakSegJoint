## Make and run qsub scripts for all steps of the PeakSegJoint pipeline.
R.bin <- R.home("bin")
Rscript <- file.path(R.bin, "Rscript")
labels.txt.file <- # interactive default for debugging.
  system.file(file.path("exampleData", "manually_annotated_region_labels.txt"),
              package="PeakSegDP")
argv <- commandArgs(trailingOnly=TRUE)
if(length(argv) != 1){
  stop("usage: AllSteps_qsub.R path/to/labels.txt
where there are path/to/*/*.bedGraph files")
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
system(cmd)

## Starting with Step1, we add vectors of commands that will wait for
## each other.
cmd.list <- list()

data.dir <- normalizePath(data.dir, mustWork=TRUE)#NEEDS TO BE ABSOLUTE!
Step1 <-
  system.file(file.path("exec", "Step1-chunk-coverage.R"),
              mustWork=TRUE,
              package="PeakSegJoint")
bedGraph.path.vec <- Sys.glob(file.path(data.dir, "*", "*.bedGraph"))
Step1.list <- list()
for(bedGraph.i in seq_along(bedGraph.path.vec)){
  bedGraph.path <- bedGraph.path.vec[[bedGraph.i]]
  Step1.list[[basename(bedGraph.path)]] <-
    paste0(Rscript, " ", Step1, " ", bedGraph.path)
}
cmd.list$Step1 <- do.call(c, Step1.list)

Step2 <-
  system.file(file.path("exec", "Step2-segment-labeled-chunks.R"),
              mustWork=TRUE,
              package="PeakSegJoint")
regions.RData.vec <-
  Sys.glob(file.path(data.dir, "PeakSegJoint-chunks", "*", "regions.RData"))
Step2.list <- list()
for(chunk.i in seq_along(regions.RData.vec)){
  regions.RData <- regions.RData.vec[[chunk.i]]
  chunk.dir <- dirname(regions.RData)
  chunk.str <- paste0("chunk", basename(chunk.dir))
  Step2.list[[chunk.str]] <-
    paste0(Rscript, " ", Step2, " ", chunk.dir)
}
cmd.list$Step2 <- do.call(c, Step2.list)

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
    prefix <- file.path(last.dir, prefix.only)
    script.txt <-
      paste0("#!/bin/bash
#PBS -l nodes=1:ppn=4
#PBS -l walltime=12:00:00                      
#PBS -A bws-221-ae", depend.txt, "
#PBS -o ", prefix, ".out
#PBS -e ", prefix, ".err
#PBS -V                                        
#PBS -N ", cmd.name, "\n", cmd, "\n")
    script.file <- paste0(prefix, ".sh")
    cat(script.txt, file=script.file)
    qsub.cmd <- paste(qsub, script.file)
    qsub.out <- system(qsub.cmd, intern=TRUE)
    qsub.id <- sub("[.].*", "", qsub.out)
    cat("submitted job ", qsub.id, "\n", sep="")
    depend.list[[cmd.name]] <- qsub.id
  }
}
