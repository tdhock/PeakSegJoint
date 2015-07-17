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

### Set mc.cores option from the PBS_NUM_PPN environment variable.
PPN.cores <- function(variable="PBS_NUM_PPN"){
  stopifnot(is.character(variable))
  stopifnot(length(variable) == 1)
  ppn.txt <- Sys.getenv(variable)
  ppn <- as.integer(ppn.txt)
  if(is.finite(ppn)){
    options(mc.cores=ppn)
  }
  ppn
}

