[![Build Status](https://travis-ci.org/mskilab/gGnome.svg?branch=master)](https://travis-ci.org/mskilab/gGnome)
[![codecov.io](https://img.shields.io/codecov/c/github/mskilab/gGnome.svg)](https://codecov.io/github/mskilab/gGnome?branch=master)

# gGnome
Reference-based graph representation of structurally-altered genome based on GenomicRanges.

## Install
```
library(devtools)
devtools::install_github("mskilab/gGnome")
```

## Quick start
```
library(gGnome)
```

* The object `gGraph` is instantaited from an actual JaBbA output inferred from HCC1143 cell line whole genome sequencing as a gGraph. We see now in the graph that there are 350 aberrant junctions (i.e. a somatic adjacency that is not present in reference or germline), 310 loose ends (i.e. false negative junction calls), and 1000 reference connections which help connect the adjacencies consistent with the reference genome.


```R
jab = readRDS(system.file("extdata", "jabba.simple.rds", package = "gGnome"))
win = readRDS(system.file("extdata", "win_17.21.rds", package = "gGnome"))

g1 = gGraph$new(jabba=jab)
g1
```

    only use 'jabba' or 'weaver' field, not both



    A gGraph object.
    Based on reference genome: GENOME
    
    Total non-loose segmentation:1025
    
    Junction counts:
    type
     aberrant     loose reference 
          350       310      1000 



```R
* Gray rectangle: DNA segment
* red links: aberrant junction
* blue dashed: loose ends
* gray link: reference adjacency.
```

* Take a quick look around, say chr 17 through 22

```R
plot(g1$td, c(as.character(17:22)), xaxis.chronly=T, labels.suppress=T, gap=1e7, xaxis.cex.tick=0.5)
```

![vis4](../master/inst/extdata/images/output_13_0.png)


* We can also extract the subgraph containing the 1Mbp neighborhood around a TRA between chr17 and chr21.


```R
juncs = g1$junctions
win = unlist(juncs[grl.in(juncs, GRanges(c("17", "21")), only=TRUE)])
g2 = g1$hood(win, d=1e6)
plot(g2$td, g2$window(1e5), xaxis.chronly=T, labels.suppress=T, gap=1e5, xaxis.cex.tick=0.5)
```

![vis6](../master/inst/extdata/images/output_17_0.png)
