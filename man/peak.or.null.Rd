\name{peak.or.null}
\alias{peak.or.null}
\title{peak or null}
\description{Run PeakSegJointSeveral in segmentation problems that overlap
across the genome. Will later use the results of this function in
clusterProblems to determine a set of non-overlapping segmentation
problems.}
\usage{peak.or.null(count.data)}
\arguments{
  \item{count.data}{data.frame of counts to segment.}
}

\value{data.table(chromStart, chromEnd) = the peak of the most complex
model, or NULL if there are no feasible models.}

\author{Toby Dylan Hocking}





