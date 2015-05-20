pick.best.index <- structure(function
### Find a minimum near the middle.
(err
### Vector of errors to minimize.
 ){
  stopifnot(is.numeric(err))
  nparam <- length(err)
  candidates <- which(err==min(err))
  if(length(err)==1)return(candidates)
  st <- abs(median(candidates)-candidates)
  middle <- candidates[which.min(st)]
  if(all(diff(err)==0))return(middle)
  if(nparam %in% candidates && 1 %in% candidates){
    cat("Warning: strange error profile, picking something near the center\n")
    print(as.numeric(err))
    d <- diff(candidates)>1
    if(any(d)){
      which(d)[1]
    }else{
      middle
    }
  }else if(1 %in% candidates){
    max(candidates)
  }else if(nparam %in% candidates){
    min(candidates)
  }else {
    middle
  }
### Integer index of the minimal error.
},ex=function(){
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
})
