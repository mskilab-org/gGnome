#' gGnome: An R API to manipulate genome graphs
#' 
#' @importFrom parallel mclapply
#' @importFrom reshape2 melt
#' @importFrom VariantAnnotation readVcf info geno
#' @importFrom RCurl url.exists
#' @importFrom S4Vectors cbind.DataFrame
#' @importMethodsFrom S4Vectors do.call elementNROWS Rle mcols mcols<- values values<- elementMetadata elementMetadata<- from to
#' @importMethodsFrom BiocGenerics width cbind
#' @importFrom gUtils %$% %Q% %*% gr_construct_by gr_deconstruct_by MULTIDIM DIM NCOL2 append_by_field_from_seqnames do.assign rleseq BY.SEP1 BY.SEP2
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