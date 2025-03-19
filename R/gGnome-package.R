#' gGnome: An R API to manipulate genome graphs
#' 
#' @importFrom parallel mclapply
#' @importFrom reshape2 melt
#' @importFrom VariantAnnotation readVcf info geno
#' @importFrom RCurl url.exists
#' @importFrom S4Vectors cbind.DataFrame
#' @importMethodsFrom S4Vectors do.call elementNROWS Rle mcols mcols<- values values<- elementMetadata elementMetadata<- from to
#' @importMethodsFrom BiocGenerics width cbind union
#' @importFrom GenomeInfoDb seqlengths seqlevels seqnames isCircular Seqinfo seqinfo seqlengths<- seqlevels<- seqnames<-
#' @importFrom stats setNames
#' @importFrom jsonlite read_json
#' @importFrom grDevices col2rgb rgb
#' @importFrom graphics title
#' @importFrom stats aggregate as.dist chisq.test cutree dist hclust median predict quantile runif var
#' @importFrom utils fix read.delim str
#' @importFrom MatrixGenerics rowRanges
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