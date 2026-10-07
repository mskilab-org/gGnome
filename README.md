[![R-CMD-check](https://github.com/zining01/gGnome/actions/workflows/check-release.yaml/badge.svg)](https://github.com/zining01/gGnome/actions/workflows/check-release.yaml)
[![codecov.io](https://img.shields.io/codecov/c/github/mskilab/gGnome.svg)](https://codecov.io/github/mskilab/gGnome?branch=master)

# gGnome

```
            █████████                                             
           ███░░░░░███                                            
  ███████ ███     ░░░ ████████    ██████  █████████████    ██████ 
 ███░░███░███        ░░███░░███  ███░░███░░███░░███░░███  ███░░███
░███ ░███░███    █████░███ ░███ ░███ ░███ ░███ ░███ ░███ ░███████ 
░███ ░███░░███  ░░███ ░███ ░███ ░███ ░███ ░███ ░███ ░███ ░███░░░  
░░███████ ░░█████████ ████ █████░░██████  █████░███ █████░░██████ 
 ░░░░░███  ░░░░░░░░░ ░░░░ ░░░░░  ░░░░░░  ░░░░░ ░░░ ░░░░░  ░░░░░░  
 ███ ░███                                                         
░░██████                                                          
 ░░░░░░                                                           
```

**R API for genome graphs**

The **gGnome** package provides a flexible, queriable `R` interface to graphs
and walks of reference genomic intervals.  **gGnome** is written in the `R6` object
oriented standard and built around a powerful `GenomicRanges`, `data.table`, and
`igraph` backend, and thus supports agile interaction with graphs consisting of
hundreds of thousands of nodes and edges.  

## VCF compatibility (0.1.1)

`read.juncs()` supports Esvee paired-breakend ALT records with biological
`SVTYPE=DEL`, `DUP`, `INV`, `TRA`, or `BND`. A valid paired ALT takes precedence
over symbolic type: source POS coordinates and reciprocal mate pairing produce
one adjacency per pair, without symbolic offsets or extra inversion boundaries.
Original SVTYPE, FILTER, and sample metadata are retained; no FILTER selection
policy is applied. Older GRIDSS PARID pairing remains supported.

Symbolic records still require valid END values (TRA also requires CHR2/CT).
Missing or invalid symbolic endpoints are warned about and skipped independently
of valid paired records. INS and SGL remain excluded, with skipped type counts
reported. Empty/header-only and all-unsupported VCFs return an empty GRangesList.

The fix requires reinstalling gGnome in the runtime that reads the VCFs. A host
installation does not update a JaBbA container; rebuild under a new image identity
and verify `packageVersion("gGnome")` and `find.package("gGnome")` inside it.
Local parser regressions do not establish JaBbA optimization success.

## Install

1. Install R-3.6 or up

2. Install devtools

```{r}
install.packages('devtools')
install.packages('testthat')
```
3. Install gGnome and dependent packages

```{r}
## allows dependencies that throw warnings to install
Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS = TRUE)


devtools::install_github('mskilab/gGnome')
```


Documentation 
------------

[gGnome Tutorial](http://mskilab.com/gGnome/tutorial.html)

[alttext](https://github.com/mskilab/gGnome/articles/gGnome.png)

[PkgDown Site](https://tanubrata.github.io/gGnome/index.html)

<!---
[gGnome Developer Reference](docs/reference.md)
-->

<div id="attributions"/>

Attributions
------------
> Marcin Imieliński - Associate Professor, Weill Cornell Medicine
> Core Member, New York Genome Center.

> Xiaotong Yao - (as) Graduate Research Assistant, Weill Cornell Medicine, New York
> Genome Center.

> Zi-Ning Choo - Graduate Research Assistant, Weill Cornell Medicine, New York
> Genome Center.

> Alon Shaiber - Genomics Data Scientist, Weill Cornell Medicine, New York
> Genome Center.

> Joseph DeRose - (as) Undergraduate Research Assistant, New York Genome Center.

> Rick Mortensen - (as) Undergraduate Research Assistant, New York Genome Center,
> Memorial Sloan-Kettering Cancer Center

> Max Chao - Bioinformatic Analyst, Weill Cornell Medicine, New York Genome
> Center

Funding sources
------------

<img
src="https://static1.squarespace.com/static/562537a8e4b0bbf0e0b819f1/5ad81984575d1f7d69517350/5ad819f02b6a28750f79597c/1524111879079/DDCF.jpeg?format=1500w"
height="150" class ="center"> <img
src="https://static1.squarespace.com/static/562537a8e4b0bbf0e0b819f1/5ad81984575d1f7d69517350/5ad819b8aa4a996c2d584594/1524111841815/BWF.png?format=500w"
height="150" class ="center">




```
