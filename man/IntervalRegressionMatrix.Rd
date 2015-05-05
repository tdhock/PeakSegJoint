\name{IntervalRegressionMatrix}
\alias{IntervalRegressionMatrix}
\title{IntervalRegressionMatrix}
\description{Solve the squared hinge loss interval regression problem for one
regularization parameter: w* = argmin_w L(w) + regularization *
||w||_1 where L(w) is the average squared hinge loss with respect
to the targets, and ||w||_1 is the L1-norm of the weight vector
(excluding the first element, which is the un-regularized
intercept or bias term).}
\usage{IntervalRegressionMatrix(features, targets, initial.weight.vec, 
    regularization, threshold = 0.01, max.iterations = 10000, 
    Lipschitz = NULL, verbose = 2)}
\arguments{
  \item{features}{Scaled numeric feature matrix (problems x features). The first
column/feature should be all ones and will not be regularized.}
  \item{targets}{Numeric target matrix (problems x 2).}
  \item{initial.weight.vec}{initial guess for weight vector (features).}
  \item{regularization}{Degree of L1-regularization.}
  \item{threshold}{When the stopping criterion gets below this threshold, the
algorithm stops and declares the solution as optimal.}
  \item{max.iterations}{Error if the algorithm has not found an optimal solution after
this many iterations.}
  \item{Lipschitz}{A numeric scalar or NULL, which means to compute Lipschitz as the
mean of the squared L2-norms of the rows of the feature matrix.}
  \item{verbose}{Print messages if >= 2.}
}

\value{Numeric vector of scaled weights w of the affine function f_w(X) =
X \%*\% w for a scaled feature matrix X with the first row entirely
ones.}

\author{Toby Dylan Hocking}





