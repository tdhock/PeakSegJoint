\name{IntervalRegressionProblems}
\alias{IntervalRegressionProblems}
\title{IntervalRegressionProblems}
\description{Compute a sequence of interval regression models for increasingly
more L1 regularization, until we get to a regularization parameter
that gives an optimal weight vector of zero. The problem is w* =
argmin_w L(w) + regularization * ||w||_1 where L(w) is the mean
squared hinge loss and ||w||_1 is the L1-norm of the non-intercept
coefficients. We first scale the input features and then
repeatedly call IntervalRegressionMatrix, using warm restarts.}
\usage{IntervalRegressionProblems(problem.list, initial.regularization = 0.001, 
    factor.regularization = 1.5, verbose = 1, ...)}
\arguments{
  \item{problem.list}{List of problems with features (numeric matrix), target (numeric
vector of length 2), and optionally weight (numeric length 1).}
  \item{initial.regularization}{Initial regularization parameter.}
  \item{factor.regularization}{Increase regularization by this factor after finding an optimal
solution. Or NULL to compute just one model
(initial.regularization).}
  \item{verbose}{Print messages if >= 1.}
  \item{\dots}{Other parameters to pass to IntervalRegressionMatrix.}
}

\value{List representing fit model. You can do
fit$predict(feature.matrix) to get a predicted log penalty
value. The mean.vec and sd.vec were used for scaling the training
data matrices. The param.mat is the n.features * n.regularization
numeric matrix of optimal coefficients.}

\author{Toby Dylan Hocking}




\examples{
library(PeakSegJoint)
data(H3K4me3.PGP.immune.4608)
chrom.vec <- sub(":.*", "", names(H3K4me3.PGP.immune.4608))
table(chrom.vec)
train.chroms <- c("chr1", "chr9")
sets <-
  list(train=chrom.vec \%in\% train.chroms,
       validation=! chrom.vec \%in\% train.chroms)
train.problems <- H3K4me3.PGP.immune.4608[sets$train]
fit <- IntervalRegressionProblems(train.problems)
set.error.list <- list()
for(set.name in names(sets)){
  in.set <- sets[[set.name]]
  problem.list <- H3K4me3.PGP.immune.4608[in.set]
  error.mat.list <- list()
  for(problem.name in names(problem.list)){
    problem <- problem.list[[problem.name]]
    pred.log.lambda <- fit$predict(problem$features)
    too.hi <- problem$target[2] < pred.log.lambda
    too.lo <- pred.log.lambda < problem$target[1]
    is.error <- too.hi | too.lo
    error.mat.list[[problem.name]] <- is.error
  }
  error.mat <- do.call(rbind, error.mat.list)
  percent.error <- colMeans(error.mat) * 100
  set.error.list[[set.name]] <-
    data.frame(set.name,
               regularization=fit$regularization,
               percent.error)
}

## Coefficients of the best models.
min.validation <- 
  subset(set.error.list$validation,
         percent.error==min(percent.error))
best.models <- fit$param.mat[, rownames(min.validation)]
best.nonzero <- best.models[apply(best.models!=0, 1, any), ]
print(best.nonzero)

## Plot train/validation error curves.
set.error <- do.call(rbind, set.error.list)
library(ggplot2)
ggplot()+
  geom_point(aes(-log10(regularization), percent.error),
             data=min.validation)+
  geom_line(aes(-log10(regularization), percent.error,
                group=set.name, linetype=set.name),
            data=set.error)

## Fit model with the chosen regularization to the full
## train+validation set.
chosen.regularization <- max(min.validation$regularization)
full.fit <-
  IntervalRegressionProblems(H3K4me3.PGP.immune.4608,
                             initial.regularization=chosen.regularization,
                             factor.regularization=NULL)
print(full.fit$param.mat)

}
