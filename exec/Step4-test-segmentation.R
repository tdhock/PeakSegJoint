library(data.table)

argv <-
  c("~/exampleData/PeakSegJoint-chunks/trained.model.RData",
    "jobName")

argv <- commandArgs(trailingOnly=TRUE)

print(argv)

if(length(argv) != 2){
  stop("usage: Step4.R PeakSegJoint-chunks/trained.model.RData jobName")
}

trained.model.RData <- argv[1]
job.name <- argv[2]

objs <- load(trained.model.RData)
