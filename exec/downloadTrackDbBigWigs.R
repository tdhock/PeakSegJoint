argv <- "~/projects/blueprint/labels/H3K27ac_TDH_MonoMacroMyeloid/trackDb.txt"

argv <- commandArgs(trailingOnly=TRUE)

trackDb.txt <- argv[1]

trackDb.lines <- readLines(trackDb.txt)

trackDb.str <- paste(trackDb.lines, collapse="\n")

## Parse the first occurance of pattern from each of several strings
## using (named) capturing regular expressions, returning a matrix
## (with column names).
str_match_perl <- function(string,pattern){
  stopifnot(is.character(string))
  stopifnot(is.character(pattern))
  stopifnot(length(pattern)==1)
  parsed <- regexpr(pattern,string,perl=TRUE)
  captured.text <- substr(string,parsed,parsed+attr(parsed,"match.length")-1)
  captured.text[captured.text==""] <- NA
  captured.groups <- if(is.null(attr(parsed, "capture.start"))){
    NULL
  }else{
    do.call(rbind,lapply(seq_along(string),function(i){
      st <- attr(parsed,"capture.start")[i,]
      if(is.na(parsed[i]) || parsed[i]==-1)return(rep(NA,length(st)))
      substring(string[i],st,st+attr(parsed,"capture.length")[i,]-1)
    }))
  }
  result <- cbind(captured.text,captured.groups)
  colnames(result) <- c("",attr(parsed,"capture.names"))
  result
}

## Parse several occurances of pattern from each of several strings
## using (named) capturing regular expressions, returning a list of
## matrices (with column names).
str_match_all_perl <- function(string,pattern){
  stopifnot(is.character(string))
  stopifnot(is.character(pattern))
  stopifnot(length(pattern)==1)
  parsed <- gregexpr(pattern,string,perl=TRUE)
  lapply(seq_along(parsed),function(i){
    r <- parsed[[i]]
    full <- substring(string[i],r,r+attr(r,"match.length")-1)
    starts <- attr(r,"capture.start")
    if(is.null(starts)){
      m <- cbind(full)
      colnames(m) <- ""
      m
    }else{
      if(r[1]==-1)return(matrix(nrow=0,ncol=1+ncol(starts)))
      names <- attr(r,"capture.names")
      lengths <- attr(r,"capture.length")
      subs <- substring(string[i],starts,starts+lengths-1)
      m <- matrix(c(full,subs),ncol=length(names)+1)
      colnames(m) <- c("",names)
      m
    }
  })
}

track.pattern <-
  paste0("    ",
         "(?<variable>.*?)",
         " ",
         "(?<value>.*?)",
         "\n")

track.vec <- strsplit(trackDb.str, "\n\n")[[1]]
subtrack.vec <- track.vec[-1]

match.list <- str_match_all_perl(subtrack.vec, track.pattern)

subGroup.pattern <-
  paste0("(?<variable>[^ ]+)",
         "=",
         "(?<value>[^ ]+)")

data.dir <- dirname(trackDb.txt)

for(match.mat in match.list){
  rownames(match.mat) <- match.mat[, "variable"]
  value.vec <- match.mat[, "value"]
  if(value.vec[["type"]] == "bigWig"){
    match <-
      str_match_all_perl(value.vec[["subGroups"]],
                         subGroup.pattern)[[1]]
    rownames(match) <- match[, "variable"]
    cell.type <- match["sampleType", "value"]
    u <- value.vec[["bigDataUrl"]]
    bigwig.base <- sub("[.]bw$", ".bigwig", basename(u))
    bigwig.dir <- file.path(data.dir, cell.type)
    dir.create(bigwig.dir, showWarnings=FALSE)
    bigwig.path <- file.path(bigwig.dir, bigwig.base)
    if(!file.exists(bigwig.path)){
      download.file(u, bigwig.path)
    }
  }
}
