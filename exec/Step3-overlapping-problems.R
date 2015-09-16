library(PeakSegJoint)

argv <-
  c("~/exampleData/PeakSegJoint-chunks/trained.model.RData",
    "1")

argv <-
  c("~/genomepipelines/H3K4me3_TDH_immune/PeakSegJoint-chunks/trained.model.RData",
    "1")

argv <- commandArgs(trailingOnly=TRUE)

print(argv)

ppn <- PPN.cores()
if(!is.na(ppn))options(mc.cores=ppn/2)

if(length(argv) != 2){
  stop("usage: Step3.R PeakSegJoint-chunks/trained.model.RData jobNum")
}

trained.model.RData <- normalizePath(argv[1])
jobNum <- argv[2]

objs <- load(trained.model.RData)

job.problems <- problems.by.job[[jobNum]]

bases.per.problem <- train.errors.picked$bases.per.problem

chunks.dir <- dirname(trained.model.RData)
data.dir <- dirname(chunks.dir)
bigwig.file.vec <- Sys.glob(file.path(data.dir, "*", "*.bigwig"))
sample.group.vec <- basename(dirname(bigwig.file.vec))
sample.id.vec <- sub("[.]bigwig$", "", basename(bigwig.file.vec))
names(bigwig.file.vec) <- paste0(sample.group.vec, "/", sample.id.vec)

OverlappingProblem <- function(problem.i){
  message(sprintf("%10d / %10d overlapping problems",
                  problem.i, nrow(job.problems)))
  problem <- job.problems[problem.i, ]
  profile.list <- readBigWigSamples(problem, bigwig.file.vec)
  peak.only <- peak.or.null(profile.list)
  if(!is.null(peak.only)){
    data.table(problem, peak.only)
  }
}

overlapping.peaks.list <-
  mclapply.or.stop(1:nrow(job.problems), OverlappingProblem)
overlapping.peaks <- do.call(rbind, overlapping.peaks.list)

job.RData <-
  file.path(data.dir,
            "PeakSegJoint-overlapping",
            paste0(jobNum, ".RData"))

save(overlapping.peaks, file=job.RData)
