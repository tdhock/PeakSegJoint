if(!require(data.table))
  install.packages("data.table")

argv <-
  system.file(file.path("exampleData", "manually_annotated_region_labels.txt"),
              package="PeakSegDP")

argv <- "/gs/project/mugqic/epigenome/pipelines/atac_seq/v_1/EMC_Temporal_Change/hg19/toby_peak_calling/toby_chunks"

argv <- commandArgs(trailingOnly=TRUE)

print(argv)

if(length(argv) != 1){
  stop("usage: convert-labels.R path/to/labels.txt
where there are path/to/celltype/*.bedGraph files")
}

labels.file <- argv[1]

g.pos.pattern <-
  paste0("(?<chrom>chr.+?)",
         ":",
         "(?<chromStart>[0-9 ,]+)",
         "-",
         "(?<chromEnd>[0-9 ,]+)",
         " ",
         "(?<annotation>[a-zA-Z]+)",
         "(?<types>.*)")

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

## there should be bedGraph files in subdirectories under the same
## directory as labels.file.
bedGraph.files <- Sys.glob(file.path(dirname(labels.file), "*", "*.bedGraph"))
sample.id <- sub("[.]bedGraph$", "", basename(bedGraph.files))
cell.type <- basename(dirname(bedGraph.files))
sample.df <- data.frame(sample.id, cell.type)
samples.by.type <- split(sample.df, sample.df$cell.type)

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
  ! match.mat[, "annotation"] %in% c("peakStart", "peakEnd", "peaks", "noPeaks")
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
             types=match.mat[, "types"],
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

## determine total set of cell types with positive=Peak annotations.
stripped <- gsub(" *$", "", gsub("^ *", "", match.df$types))
is.noPeaks <- stripped == ""
type.list <- strsplit(stripped, split=" ")
names(type.list) <- rownames(match.df)
cell.types <- unique(unlist(type.list))
cat("cell types with peak annotations: ",
    paste(cell.types, collapse=", "),
    "\n",
    sep="")

## Use positive region sizes to limit problem size grid search.
positive.df <- match.df[!is.noPeaks, ]
positive.bases <- with(positive.df, chromEnd-chromStart)
cat("peak region size in bases:\n")
print(positive.q <- quantile(positive.bases))
min.bases.per.problem <- positive.q[["25%"]]

match.by.chunk <- split(match.df, match.df$chunk.id)
regions.by.chunk <- list()
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
    type.vec <- type.list[[rownames(chunk.row)]]
    is.observed <- cell.types %in% type.vec
    observed <- cell.types[is.observed]
    not.observed <- cell.types[!is.observed]
    to.assign <- list()
    ann <- chunk.row$annotation
    to.assign[observed] <- ann
    to.assign[not.observed] <- "noPeaks"
    for(cell.type in names(to.assign)){
      relevant.samples <- samples.by.type[[cell.type]]
      annotation <- to.assign[[cell.type]]
      regions.list[[paste(ann.i, cell.type)]] <- 
        data.table(sample.id=paste(relevant.samples$sample.id),
                   cell.type,
                   chrom=chunk.row$chrom,
                   chromStart=chunk.row$chromStart,
                   chromEnd=chunk.row$chromEnd,
                   annotation)
    }
  }
  one.chunk <- do.call(rbind, regions.list)
  setkey(one.chunk, chromStart, chromEnd)
  regions.by.chunk[[paste(chunk.id)]] <- one.chunk
}

## Define possible problem sizes.
bases.per.problem.all <- as.integer(4.5 * 2^seq(3, 20, by=0.5))
bases.per.problem.vec <-
  bases.per.problem.all[min.bases.per.problem < bases.per.problem.all &
                          bases.per.problem.all < min.bases.per.problem * 100]

## Write chunks to separate RData files.
for(chunk.str in names(regions.by.chunk)){
  regions <- regions.by.chunk[[chunk.str]]

  max.chromEnd <- max(regions$chromEnd)
  min.chromStart <- min(regions$chromStart)

  problems.by.res <- list()
  for(bases.per.problem in bases.per.problem.vec){
    res.str <- paste(bases.per.problem)
    problemSeq <- seq(0, max.chromEnd, by=bases.per.problem)
    problemStart <-
      as.integer(sort(c(problemSeq,
                        problemSeq+bases.per.problem/2,
                        problemSeq+bases.per.problem*3/4,
                        problemSeq+bases.per.problem*1/4)))
    problemEnd <- problemStart+bases.per.problem
    is.overlap <- min.chromStart < problemEnd &
      problemStart < max.chromEnd
    problem.name <- sprintf("%s:%d-%d", chrom, problemStart, problemEnd)
    res.problems <- 
      data.table(problem.name, 
                 bases.per.problem, problemStart, problemEnd)[is.overlap,]
    setkey(res.problems, problemStart, problemEnd)
    problems.by.res[[res.str]] <- res.problems
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
              chunk.str,
              "regions.RData")
  RData.dir <- dirname(RData.file)
  dir.create(RData.dir, showWarnings=FALSE, recursive=TRUE)
  cat("Writing ", nrow(regions), " labels to ", RData.file, "\n", sep="")
  save(regions, chunk, problems.by.res, file=RData.file)
}
