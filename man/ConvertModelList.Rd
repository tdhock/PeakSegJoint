\name{ConvertModelList}
\alias{ConvertModelList}
\title{ConvertModelList}
\description{Convert a model list from the non-repetitive format that we get
from the C code to the repetitive format that is more useful for
plotting.}
\usage{ConvertModelList(model.list)}
\arguments{
  \item{model.list}{List from PeakSegJointHeuristic(...) or PeakSegJointSeveral(...).}
}

\value{List of data.frames: segments has 1 row for each segment mean,
sample, and model size (peaks, sample.id, sample.group,
chromStart, chromEnd, mean); peaks is the same kind of data.frame
as segments, but with only the second/peak segments; loss has one
row for each model size; modelSelection has one row for each model
size that can be selected, see exactModelSelection.}

\author{Toby Dylan Hocking}





