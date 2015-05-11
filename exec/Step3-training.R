if(!require(data.table))
  install.packages("data.table")
if(!require(PeakSegJoint))
  install.packages("PeakSegJoint")

argv <-
  system.file(file.path("exampleData",
                        "PeakSegJoint-chunks"),
              package="PeakSegDP")

argv <- commandArgs(trailingOnly=TRUE)

print(argv)

if(length(argv) != 1){
  stop("usage: Step2.R path/to/PeakSegJoint-chunks")
}

chunks.dir <- normalizePath(argv[1])

problems.RData.vec <- Sys.glob(file.path(chunks.dir, "*", "problems.RData"))

err.mat.list <- list()
for(chunk.i in seq_along(problems.RData.vec)){
  problems.RData <- problems.RData.vec[[chunk.i]]
  objs <- load(problems.RData)
  err.vec <- step2.error$errors
  names(err.vec) <- rownames(step2.error)
  err.mat.list[[chunk.i]] <- err.vec
}
err.mat <- do.call(rbind, err.mat.list)
err.vec <- colSums(err.mat)
res.str <- names(err.vec)[which.min(err.vec)]

stop("train supervised model")
