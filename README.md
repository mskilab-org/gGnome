[![Build Status](https://travis-ci.org/mskilab/gGnome.svg?branch=master)](https://travis-ci.org/mskilab/gGnome)
[![codecov.io](https://img.shields.io/codecov/c/github/mskilab/gGnome.svg)](https://codecov.io/github/mskilab/gGnome?branch=master)

# gGnome

The **gGnome** package provides a flexible, queriable `R` interface to graphs
and walks of reference genomic intervals.  **gGnome** is written in the `R6` object
oriented standard and built around a powerful `GenomicRanges`, `data.table`, and
`igraph` backend, and thus supports agile interaction with graphs consisting of
hundreds of thousands of nodes and edges.  

## Install

1. Install dependent packages and latest Bioconductor (if you haven't already)

```{r}
source("https://bioconductor.org/biocLite.R")
biocLite("GenomicRanges")
```
install.packages('devtools')
install.packages('testthat')


2. Install dependent mskilab R packages

```{r}
devtools::install_github('mskilab/gUtils')
devtools::install_github('mskilab/gTrack')
```

4. Install gGnome

```{r}
devtools::install_github('mskilab/gGnome)
```

Documentation 
------------

[gGnome Tutorial](http://mskilab.com/gGnome/tutorial.html)

[gGnome Developer Reference](docs/reference.md)

<div id="attributions"/>

![alttext](https://github.com/mskilab/gGnome/raw/master/docs/gGnome.png) 

Attributions
------------
> Marcin Imielinski - Assistant Professor, Weill-Cornell Medical College. 
> Core Member, New York Genome Center.

> Xiaotong Yao - Graduate Research Assistant, Weill Cornell Medicine, New York
> Genome Center.

> Joseph DeRose - Undergraduate Research Assistant, New York Genome Center.

> Rick Mortensen - Undergraduate Research Assistant, New York Genome Center,
> Memorial Sloan-Kettering Cancer Center





```
