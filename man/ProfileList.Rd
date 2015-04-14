\name{ProfileList}
\alias{ProfileList}
\title{ProfileList}
\description{Convert a data.frame or list of profiles to a list that can be
passed to .Call("PeakSegJointHeuristic...").}
\usage{ProfileList(profiles)}
\arguments{
  \item{profiles}{List of data.frames with columns chromStart, chromEnd, count, or
single data.frame with additional column sample.id.}
}

\value{Named list of data.frames with columns chromStart, chromEnd,
count.}

\author{Toby Dylan Hocking}





