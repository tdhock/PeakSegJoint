\name{tv.curves}
\alias{tv.curves}
\title{tv curves}
\description{Fit regularized interval regression models and compute train and
validation error curves.}
\usage{tv.curves(chunk.problem.list, regions.list, n.folds = 4, 
    seed = 1)}
\arguments{
  \item{chunk.problem.list}{chunk.problem.list[[chunk.name]][[problem.name]] is a list with
peaks, segments, loss, modelSelection, features, target.}
  \item{regions.list}{regions.list[[chunk.name]] is a data.frame of annotated region
labels for one chunk.}
  \item{n.folds}{Number of cross-validation folds to use to estimate the
regularization parameter.}
  \item{seed}{Random seed to use to choose cross-validation fold ID numbers for
each chunk.}
}

\value{Data.frame of train/validation error curves.}

\author{Toby Dylan Hocking}





