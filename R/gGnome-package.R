#' gGnome: An R API to manipulate genome graphs
#' 
#' @importFrom parallel mclapply
#' @importFrom reshape2 melt
#' @importFrom VariantAnnotation readVcf info geno
#' @importFrom RCurl url.exists
#' @importMethodsFrom S4Vectors do.call Rle
#' @importMethodsFrom BiocGenerics width
#' @importFrom gUtils %$% %Q% %*%
#' @importFrom GenomeInfoDb seqlengths seqlevels seqnames isCircular Seqinfo seqlengths<- seqlevels<- seqnames<-
#' @importFrom stats setNames
#' @import methods
#' @import R6
#' @import data.table
#' @import Matrix
#' @import jsonlite
#' @import GenomicRanges
#' @import igraph
#' @import gUtils
#' @import gTrack
#' @import fishHook
#' @useDynLib gGnome
"_PACKAGE"