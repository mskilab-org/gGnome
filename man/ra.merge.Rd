% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{ra.merge}
\alias{ra.merge}
\title{Merges rearrangements represented by \code{GRangesList} objects}
\usage{
ra.merge(..., pad = 0, ignore.strand = FALSE)
}
\arguments{
\item{...}{GRangesList representing rearrangements to be merged}

\item{pad}{non-negative integer specifying padding}

\item{ignore.strand}{whether to ignore strand (implies all strand information will be ignored, use at your own risk)}

\item{ind}{logical flag (default FALSE) specifying whether the "seen.by" fields should contain indices of inputs (rather than logical flags) and NA if the given junction is missing}
}
\value{
\code{GRangesList} of merged junctions with meta data fields specifying which of the inputs each outputted junction was "seen.by"
}
\description{
Determines overlaps between two or more piles of rearrangement junctions (as named or numbered arguments) +/- padding
and will merge those that overlap into single junctions in the output, and then keep track for each output junction which
of the input junctions it was "seen in" using logical flag  meta data fields prefixed by "seen.by." and then the argument name
(or "seen.by.ra" and the argument number)
}
\examples{

# generate some junctions
gr1 <- GRanges(1, IRanges(1:10, width = 1), strand = rep(c('+', '-'), 5))
gr2 <- GRanges(1, IRanges(4 + 1:10, width = 1), strand = rep(c('+', '-'), 5))
ra1 = split(gr1, rep(1:5, each = 2))
ra2 = split(gr2, rep(1:5, each = 2))

ram = ra.merge(ra1, ra2)
values(ram) # shows the metadata with TRUE / FALSE flags

ram2 = ra.merge(ra1, ra2, pad = 5) # more inexact matching results in more merging
values(ram2)

ram3 = ra.merge(ra1, ra2) #indices instead of flags
values(ram3)
}
