\name{PeakSegJointHeuristicStep2}
\alias{PeakSegJointHeuristicStep2}
\title{PeakSegJointHeuristicStep2}
\description{Run the first and second steps of the PeakSegJoint fast heuristic
optimization algorithm. Step2 the SearchNearPeak subroutine
described in the JointZoom Algorithm of arXiv:1506.01286, and it
is guaranteed to return feasible segmentations (seg1 < seg2 >
seg3). NB: this function is only for testing the C code against
the R implementation (search tests/testthat/*.R files for
Step2). For real data see \code{\link{PeakSegJointSeveral}}.}
\usage{PeakSegJointHeuristicStep2(profiles, 
    bin.factor = 2L)}
\arguments{
  \item{profiles}{List of data.frames with columns chromStart, chromEnd, count, or
single data.frame with additional column sample.id.}
  \item{bin.factor}{Size of bin pyramid. Bigger values result in slower computation.}
}

\value{List of model fit results, which can be passed to \code{\link{ConvertModelList}}
for easier interpretation.}

\author{Toby Dylan Hocking}





