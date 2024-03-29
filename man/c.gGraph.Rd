% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gGnome.R
\name{c.gGraph}
\alias{c.gGraph}
\title{c.gGraph}
\usage{
\method{c}{gGraph}(...)
}
\arguments{
\item{...}{set of gGraph arguments (names will be incorporated as graph metadata $parent.graph)}
}
\value{
A new gGraph object that is the union of the nodes and edges in the input gGraphs
}
\description{
Concatenates gGraphs without doing any merging or aggregating
}
\author{
Marcin Imielinski
}
