\name{best.on.validation}
\alias{best.on.validation}
\title{best on validation}
\description{For each validation fold, pick the most regularized (least
complex) model that has minimum error.}
\usage{best.on.validation(err.curves)}
\arguments{
  \item{err.curves}{Result data.frame of train/validation error curves from tv.curves.}
}

\value{Data.frame with one row for each validation fold.}

\author{Toby Dylan Hocking}





