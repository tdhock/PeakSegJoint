library(PeakSegJoint)

argv <- # interactive default for debugging.
  system.file("exampleData", "manually_annotated_region_labels.txt",
              package="PeakSegJoint")
argv <- # interactive default for debugging.
  system.file("exampleData",
              c("other_labels.txt", "manually_annotated_region_labels.txt"),
              package="PeakSegJoint")

argv <- "~/exampleData/manually_annotated_region_labels.txt"

argv <- "/gs/project/mugqic/epigenome/pipelines/atac_seq/v_1/EMC_Temporal_Change/hg19/toby_peak_calling/toby_chunks"

argv <- commandArgs(trailingOnly=TRUE)

print(argv)

if(length(argv) == 0){
  stop("usage: convert-labels.R path/to/labels.txt
where there are path/to/sample_group/*.bigwig files")
}

g.pos.pattern <-
  paste0("(?<chrom>chr.+?)",
         ":",
         "(?<chromStart>[0-9 ,]+)",
         "-",
         "(?<chromEnd>[0-9 ,]+)",
         " ",
         "(?<annotation>[a-zA-Z]+)",
         "(?<sample_groups>.*)")

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
  captured.groups <- do.call(rbind,lapply(seq_along(string),function(i){
    st <- attr(parsed,"capture.start")[i,]
    if(is.na(parsed[i]) || parsed[i]==-1)return(rep(NA,length(st)))
    substring(string[i],st,st+attr(parsed,"capture.length")[i,]-1)
  }))
  result <- cbind(captured.text,captured.groups)
  colnames(result) <- c("",attr(parsed,"capture.names"))
  result
}

valid.labels <- c("peakStart", "peakEnd", "peaks", "noPeaks")
regions.by.file <- list()
regions.by.chunk.file <- list()
chunk.limits.list <- list()
for(labels.file in argv){
  ## there should be bigwig files in subdirectories under the same
  ## directory as labels.file.
  bigwig.files <- Sys.glob(file.path(dirname(labels.file), "*", "*.bigwig"))
  sample.id <- sub("[.]bigwig$", "", basename(bigwig.files))
  sample.group <- basename(dirname(bigwig.files))
  sample.df <- data.frame(sample.id, sample.group)
  samples.by.group <- split(sample.df, sample.df$sample.group)

  cat("Reading ", labels.file, "\n", sep="")
  labels.lines <- readLines(labels.file)
  is.blank <- labels.lines == ""
  chunk.id <- cumsum(is.blank)+1L
  label.df <- data.frame(chunk.id, line=labels.lines)[!is.blank, ]
  cat(length(unique(label.df$chunk.id)), " chunks, ",
      nrow(label.df), " label lines\n", sep="")

  ## Error checking.
  raw.vec <- paste(label.df$line)
  line.vec <- gsub(",", "", raw.vec)
  match.mat <- str_match_perl(line.vec, g.pos.pattern)
  stopifnot(!is.na(match.mat[,1]))
  not.recognized <-
    ! match.mat[, "annotation"] %in% valid.labels
  if(any(not.recognized)){
    print(raw.vec[not.recognized])
    print(match.mat[not.recognized, , drop=FALSE])
    stop("unrecognized annotation")
  }
  match.df <-
    data.frame(chrom=match.mat[, "chrom"],
               chromStart=as.integer(match.mat[, "chromStart"]),
               chromEnd=as.integer(match.mat[, "chromEnd"]),
               annotation=match.mat[, "annotation"],
               sample.groups=match.mat[, "sample_groups"],
               chunk.id=label.df$chunk.id,
               stringsAsFactors=FALSE)
  match.by.chrom <- split(match.df, match.df$chrom)
  for(chrom in names(match.by.chrom)){
    chrom.df <- match.by.chrom[[chrom]]
    sorted <- chrom.df[with(chrom.df, order(chromStart, chromEnd)), ]
    same.as.next <- diff(sorted$chromStart) <= 0
    if(any(same.as.next)){
      bad.i <- which(same.as.next)
      print(sorted[c(bad.i, bad.i+1), ])
      stop("chromStart not increasing")
    }
    if(any(with(sorted, chromStart >= chromEnd))){
      print(sorted)
      stop("chromStart >= chromEnd")
    }
    overlaps.next <-
      with(sorted, chromStart[-1] < chromEnd[-length(chromEnd)])
    if(any(overlaps.next)){
      print(data.frame(sorted, overlaps.next=c(overlaps.next, FALSE)))
      stop("overlapping regions")
    }
  }

  ## determine total set of sample groups with positive=Peak
  ## annotations.
  stripped <- gsub(" *$", "", gsub("^ *", "", match.df$sample.groups))
  is.noPeaks <- stripped == ""
  sample.group.list <- strsplit(stripped, split=" ")
  names(sample.group.list) <- rownames(match.df)
  sample.group.vec <- unique(unlist(sample.group.list))
  cat("sample groups with peak annotations: ",
      paste(sample.group.vec, collapse=", "),
      "\n",
      sep="")

  match.by.chunk <- split(match.df, match.df$chunk.id)
  for(chunk.id in names(match.by.chunk)){
    ## Check that all regions are on the same chrom.
    chunk.df <- match.by.chunk[[chunk.id]]
    chunkChrom <- paste(chunk.df$chrom[1])
    if(any(chunk.df$chrom != chunkChrom)){
      print(chunk.df)
      stop("each chunk must span only 1 chrom")
    }
    regions.list <- list()
    for(ann.i in 1:nrow(chunk.df)){
      chunk.row <- chunk.df[ann.i, ]
      groups.up.vec <- sample.group.list[[rownames(chunk.row)]]
      is.observed <- sample.group.vec %in% groups.up.vec
      observed <- sample.group.vec[is.observed]
      not.observed <- sample.group.vec[!is.observed]
      to.assign <- list()
      ann <- chunk.row$annotation
      to.assign[observed] <- ann
      to.assign[not.observed] <- "noPeaks"
      for(sample.group in names(to.assign)){
        relevant.samples <- samples.by.group[[sample.group]]
        if(length(relevant.samples) == 0){
          glob.str <- file.path(sample.group, "*.bigwig")
          stop("no ", glob.str, " files (but labels are present)")
        }
        annotation <- to.assign[[sample.group]]
        regions.list[[paste(ann.i, sample.group)]] <- 
          data.table(sample.id=paste(relevant.samples$sample.id),
                     sample.group,
                     chrom=chunk.row$chrom,
                     chromStart=chunk.row$chromStart,
                     chromEnd=chunk.row$chromEnd,
                     annotation)
      }
    }
    one.chunk <- do.call(rbind, regions.list)
    setkey(one.chunk, chromStart, chromEnd)
    file.and.chunk <- paste0(basename(labels.file), "-chunk", chunk.id)
    regions.by.chunk.file[[file.and.chunk]] <- one.chunk
    regions.by.file[[labels.file]][[chunk.id]] <- one.chunk
    chunk.limits.list[[file.and.chunk]] <- with(one.chunk, {
      data.frame(file.and.chunk,
                 chrom=chrom[1],
                 chromStart=min(chromStart),
                 chromEnd=max(chromEnd))
    })
  }
}
chunk.limits <- do.call(rbind, chunk.limits.list)
rownames(chunk.limits) <- NULL
limits.by.chrom <- split(chunk.limits, chunk.limits$chrom)
for(chrom in names(limits.by.chrom)){
  chrom.limits <- limits.by.chrom[[chrom]]
  ## Find overlapping chunks, and join them:
  clustered <- clusterPeaks(chrom.limits)
  limits.by.cluster <- split(clustered, clustered$cluster)
  chunks.per.cluster <- sapply(limits.by.cluster, nrow)
  not.ok <- 1 < chunks.per.cluster
  if(any(not.ok)){
    print(limits.by.cluster[not.ok])
    stop("chunks in different label files should not overlap")
  }
}

## Use positive region sizes to limit problem size grid search.
all.regions <- do.call(rbind, regions.by.chunk.file)
positive.df <- subset(all.regions, annotation %in% c("peakStart", "peakEnd"))
positive.bases <- unique(with(positive.df, chromEnd-chromStart))
cat("peak region size in bases:\n")
print(positive.q <- quantile(positive.bases))
min.bases.per.problem <- positive.q[["25%"]]

## Define possible problem sizes.
bases.per.problem.all <- as.integer(4.5 * 2^seq(3, 20, by=0.5))
bases.per.problem.vec <-
  bases.per.problem.all[min.bases.per.problem < bases.per.problem.all &
                          bases.per.problem.all < min.bases.per.problem * 100]

## Write chunks to separate RData files.
chunk.id <- 1
chunk.file.list <- list()
for(labels.file in names(regions.by.file)){
  regions.by.chunk <- regions.by.file[[labels.file]]
  for(chunk.str in names(regions.by.chunk)){
    regions <- regions.by.chunk[[chunk.str]]

    max.chromEnd <- max(regions$chromEnd)
    min.chromStart <- min(regions$chromStart)

    problems.by.res <- list()
    for(bases.per.problem in bases.per.problem.vec){
      res.str <- paste(bases.per.problem)

      chunk.problems <-
        getProblems(chrom, min.chromStart, max.chromEnd, bases.per.problem)

      problems.by.res[[res.str]] <- chunk.problems
    }

    problems <- do.call(rbind, problems.by.res)
    ## ggplot()+
    ##   geom_tallrect(aes(xmin=chromStart/1e3, xmax=chromEnd/1e3),
    ##                 data=data.frame(chromStart=min.chromStart,
    ##                   chromEnd=max.chromEnd),
    ##                 color="black",
    ##                 size=2,
    ##                 fill="grey")+
    ##   geom_tallrect(aes(xmin=chromStart/1e3, xmax=chromEnd/1e3,
    ##                     fill=annotation),
    ##                 data=regions)+
    ##   scale_fill_manual(values=ann.colors)+
    ##   geom_segment(aes(problemStart/1e3, problem.name,
    ##                    xend=problemEnd/1e3, yend=problem.name),
    ##                data=problems)

    chunk <- 
      data.table(chrom=regions$chrom[1],
                 chunkStart=min(problems$problemStart),
                 chunkEnd=max(problems$problemEnd))
    setkey(chunk, chrom, chunkStart, chunkEnd)
    
    RData.file <-
      file.path(dirname(labels.file),
                "PeakSegJoint-chunks",
                chunk.id,
                "regions.RData")
    RData.dir <- dirname(RData.file)
    dir.create(RData.dir, showWarnings=FALSE, recursive=TRUE)
    cat("Writing ", nrow(regions), " labels to ", RData.file, "\n", sep="")
    save(regions, chunk, problems.by.res, file=RData.file)

    chunk.file.list[[paste(chunk.id)]] <-
      data.frame(labels.file=basename(labels.file),
                 chunk.id)

    chunk.id <- chunk.id+1
  }
}
chunk.file.map <- do.call(rbind, chunk.file.list)
rownames(chunk.file.map) <- NULL

## Divide chunks into train+validation/test.
chunks.by.file <- split(chunk.file.map, chunk.file.map$labels.file)
all.chunk.names <- chunk.file.map$chunk.id
if(length(chunks.by.file) == 1){
  outer.folds <- 4
  if(length(all.chunk.names) < outer.folds){
    outer.folds <- length(all.chunk.names)
  }
  set.seed(1)
  outer.fold.id <- sample(rep(1:outer.folds, l=length(all.chunk.names)))
  names(outer.fold.id) <- all.chunk.names
  fold.msg <- "randomly selected folds"
}else{
  outer.folds <- length(chunks.by.file)
  outer.fold.id <- rep(NA, l=length(all.chunk.names))
  names(outer.fold.id) <- all.chunk.names
  for(file.i in seq_along(chunks.by.file)){
    file.chunks <- chunks.by.file[[file.i]]
    outer.fold.id[paste(file.chunks$chunk.id)] <- file.i
  }
  fold.msg <- "one fold for each labels file"
}
test.error.msg <-
  paste0(outer.folds, " fold cross-validation (", fold.msg, ").")

RData.file <- file.path(dirname(labels.file), "chunk.file.map.RData")
save(chunk.file.map,
     outer.folds,
     outer.fold.id,
     fold.msg,
     test.error.msg,
     file=RData.file)
