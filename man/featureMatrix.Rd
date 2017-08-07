\name{featureMatrix}
\alias{featureMatrix}
\title{featureMatrix}
\description{Compute the feature matrix for this joint segmentation problem.}
\usage{featureMatrix(profile.list)}
\arguments{
  \item{profile.list}{
}
}

\value{Numeric feature matrix (samples x features).}

\author{Toby Dylan Hocking}




\examples{
library(PeakSegJoint)
data(H3K36me3.TDH.other.chunk1, envir=environment())
lims <- c(43000000, 43200000) # left
some.counts <-
  subset(H3K36me3.TDH.other.chunk1$counts,
         lims[1] < chromEnd & chromStart < lims[2])
profile.list <- ProfileList(some.counts)
featureMatrix(profile.list)
}
