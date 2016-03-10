cDPA <- structure(function
### A constrained dynamic programming algorithm (cDPA) can be used to
### compute the best segmentation with respect to the Poisson
### likelihood, subject to a constraint on the number of segments, and
### the changes which must alternate: up, down, up, down, ...
(count,
### Integer vector of count data to segment.
 weight=rep(1, length(count)),
### Data weights (normally this is the number of base pairs).
 maxSegments
### Maximum number of segments to consider.
 ){
  stopifnot(is.numeric(count))
  stopifnot(is.numeric(weight))
  stopifnot(length(count) == length(weight))
  stopifnot(length(maxSegments) == 1)
  nData <- length(count)
  A <- .C("cDPA",
          count=as.integer(count),
          weight=as.integer(weight),
          nData=as.integer(nData),
          maxSegments=as.integer(maxSegments),
          loss=double(maxSegments*nData),
          ends=integer(maxSegments*nData),
          mean=double(maxSegments*nData),
          PACKAGE="PeakSegJoint")
  A$loss <- matrix(A$loss, nrow=maxSegments, byrow=TRUE)
  A$ends <- matrix(A$ends, nrow=maxSegments, byrow=TRUE)
  A$mean <- matrix(A$mean, nrow=maxSegments, byrow=TRUE)
  A
}, ex=function(){
  fit <- cDPA(c(0, 10, 11, 1), maxSegments=3)
  stopifnot(fit$ends[3,4] == 3)
  stopifnot(fit$ends[2,3] == 1)
})

### Extract endpoint matrix from cDPA result.
getPath <- function(A){
  n <- ncol(A$ends)
  end.mat <- matrix(NA, nrow=nrow(A$ends), ncol=nrow(A$ends))
  end.mat[1, 1] <- 0
  for(i in 2: nrow(A$ends)){
    end.mat[i, i-1] <- A$ends[i, n]
    if(end.mat[i, i-1] > 0){
      for(k in 1:(i-1)){
        end.mat[i, i-1-k] <- A$ends[i-k, end.mat[i, i-k]]
      }
    }
  }
  diag(end.mat) <- ncol(A$ends)
  return(end.mat)
}

