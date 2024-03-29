% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/eventCallers.R
\name{microhomology}
\alias{microhomology}
\title{get microhomology}
\usage{
microhomology(
  gg,
  hg,
  prefix_only = FALSE,
  pad = c(5, 10, 50, 100),
  ignore_missing = FALSE
)
}
\arguments{
\item{gg}{gGraph or Junctions}

\item{hg}{DNAStringSet or path to reference fasta}

\item{prefix_only}{(logical) default FALSE. if TRUE, considers only the longest 
common prefix. if FALSE, considers the longest local alignment.}

\item{pad}{(numeric) default NA (use the default window lengths of 5, 10, 50, 
and 100). otherwise, an integer specifying window length.}

\item{ignore_missing}{(logical) ignore junctions where at least one breakend 
is not found on the reference, and return -1 for microhomology. default FALSE, 
which will cause an error.}
}
\value{
gGraph with $pyrgo marking on nodes and edges labeling unique "events"
}
\description{
Computes microhomology at 5bp, 10bp, 50bp, and 100bp windows around ALT junctions 
of input gGraph (or Junction object) gg and adds these as an edge 
annotation to the appropriate edges.
}
\details{
The default behavior is to compute the maximum microhomology using local alignment 
across the entire window. However, the longest common prefix within each window 
can be specified by setting the argument prefix_only to TRUE.

Care should be taken that the sequence names of junctions are consistent with 
those provided in the reference. There will be an error if the sequence names 
of the junction are not a subset of those of the reference, if ignore_missing 
is FALSE (default). If ignore_missing is TRUE, then those junctions with missing 
seqnames will be assigned score -1.

Requires Biostrings.
}
