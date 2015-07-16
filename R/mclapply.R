### mclapply with error checking.
mclapply.or.stop <- function(...){
  result.list <- mclapply(...)
  is.error <- sapply(result.list, inherits, "try-error")
  if(any(is.error)){
    print(result.list[is.error])
    stop("errors in mclapply")
  }
  result.list
}

