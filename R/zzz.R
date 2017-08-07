## preset environmental variables
GENOME = readRDS(system.file("extdata", "hg19.broad.BSgenome.rds", package="gGnome"))
regularChr = c(as.character(1:22), "X", "Y") ## 24 regular chrs
Sys.setenv(DEFAULT_BSGENOME=system.file("extdata", "human_g1k_v37.chrom.sizes", package="gGnome"))
