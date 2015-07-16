library(PeakSegJoint)

argv <-
  c("~/exampleData/PeakSegJoint-chunks/trained.model.RData",
    "chr11")

argv <- commandArgs(trailingOnly=TRUE)

print(argv)

if(length(argv) != 2){
  stop("usage: Step3.R PeakSegJoint-chunks/trained.model.RData chrom")
}

trained.model.RData <- argv[1]
chrom <- argv[2]

objs <- load(trained.model.RData)

chunks.dir <- dirname(trained.model.RData)
data.dir <- dirname(chunks.dir)
bigwig.file.vec <- Sys.glob(file.path(data.dir, "*", "*.bigwig"))

setkey(test.problems, chrom)
chrom.problems <- test.problems[chrom]

readBigWigSamples <- function(problem){
  counts.by.sample <- list()
  for(bigwig.file in bigwig.file.vec){
    sample.counts <- 
      readBigWig(bigwig.file, chrom,
                 problem$problemStart, problem$problemEnd)
    sample.id <- sub("[.]bigwig$", "", basename(bigwig.file))
    counts.by.sample[[sample.id]] <- with(sample.counts, {
      data.frame(chromStart, chromEnd, count)
    })
  }
  counts.by.sample
}

Step1Problem <- function(problem.i){
  cat(sprintf("%10d / %10d problems\n", problem.i, nrow(chrom.problems)))
  problem <- chrom.problems[problem.i, ]
  profile.list <- readBigWigSamples(problem)
  peak.only <- peak.or.null(profile.list)
  if(!is.null(peak.only)){
    data.table(problem, peak.only)
  }
}

step1.results.list <-
  mclapply.or.stop(seq_along(chrom.problems$problem.name), Step1Problem)
