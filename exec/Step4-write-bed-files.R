library(PeakSegJoint)

argv <-
  c("~/exampleData/PeakSegJoint-predictions/trained.model.RData",
    "chr11")

argv <- commandArgs(trailingOnly=TRUE)

print(argv)
