[![Build Status](https://travis-ci.org/mskilab/gGnome.svg?branch=master)](https://travis-ci.org/mskilab/gGnome)
[![codecov.io](https://img.shields.io/codecov/c/github/mskilab/gGnome.svg)](https://codecov.io/github/mskilab/gGnome?branch=master)

# gGnome

The **gGnome** package provides a flexible, queriable `R` interface to graphs
and walks of reference genomic intervals.  **gGnome** is written in the `R6` object
oriented standard and built around a powerful `GenomicRanges`, `data.table`, and
`igraph` backend, and thus supports agile interaction with graphs consisting of
hundreds of thousands of nodes and edges.  

## Install

1. Install devtools

```{r}
install.packages('devtools')
install.packages('testthat')
```
2. Install gGnome and dependent packages

```{r}
## set this for R-3.5.2 
Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS = TRUE)

devtools::install_github('mskilab/gGnome)
```

Documentation 
------------

[gGnome Tutorial](http://mskilab.com/gGnome/tutorial.html)

[![alttext](https://github.com/mskilab/gGnome/raw/master/docs/gGnome.png) ](http://mskilab.com/gGnome/tutorial.html)

<!---
[gGnome Developer Reference](docs/reference.md)
-->

<div id="attributions"/>

Attributions
------------
> Marcin Imielinski - Assistant Professor, Weill Cornell Medicine
> Core Member, New York Genome Center.

> Xiaotong Yao - Graduate Research Assistant, Weill Cornell Medicine, New York
> Genome Center.

> Joseph DeRose - Undergraduate Research Assistant, New York Genome Center.

> Rick Mortensen - Undergraduate Research Assistant, New York Genome Center,
> Memorial Sloan-Kettering Cancer Center

Funding sources
------------

<img
src="https://static1.squarespace.com/static/562537a8e4b0bbf0e0b819f1/5ad81984575d1f7d69517350/5ad819f02b6a28750f79597c/1524111879079/DDCF.jpeg?format=1500w"
height="150" class ="center"> <img
src="https://static1.squarespace.com/static/562537a8e4b0bbf0e0b819f1/5ad81984575d1f7d69517350/5ad819b8aa4a996c2d584594/1524111841815/BWF.png?format=500w"
height="150" class ="center">




```
