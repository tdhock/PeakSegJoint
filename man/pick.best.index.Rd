\name{pick.best.index}
\alias{pick.best.index}
\title{pick best index}
\description{Find a minimum near the middle.}
\usage{pick.best.index(err)}
\arguments{
  \item{err}{Vector of errors to minimize.}
}

\value{Integer index of the minimal error.}

\author{Toby Dylan Hocking}




\examples{
stopifnot(pick.best.index(rep(0,100))==50)

err <- rep(1,100)
err[5] <- 0
stopifnot(pick.best.index(err)==5)

## should pick the middle
err <- rep(1,100)
err[40:60] <- 0
stopifnot(pick.best.index(err)==50)

## should pick the biggest
err <- rep(1,100)
err[1:60] <- 0
stopifnot(pick.best.index(err)==60)

## should pick the smallest
err <- rep(1,100)
err[50:100] <- 0
stopifnot(pick.best.index(err)==50)
}
