\name{error.metrics}
\alias{error.metrics}
\title{error metrics}
\description{Use a model to predict peaks, and evaluate several error metrics:
incorrect.targets, incorrect.regions, false.positives,
false.negatives.}
\usage{error.metrics(problems.by.chunk, regions.by.chunk, 
    fit)}
\arguments{
  \item{problems.by.chunk}{List[[chunk.name]][[problem.name]], a list which must have
features, modelSelection, peaks (and if it also has target, it
will be used to compute number of incorrect targets).}
  \item{regions.by.chunk}{List[[chunk.name]], a data.frame of annotated region labels for
each chunk on which the error metrics will be evaluated.}
  \item{fit}{Model fit list from IntervalRegressionProblems.}
}

\value{list of metrics, error.regions, peaks.}

\author{Toby Dylan Hocking}





