## preset environmental variables
GENOME = readRDS(system.file("extdata", "hg19.broad.BSgenome.rds", package="gGnome"))
regularChr = c(as.character(1:22), "X", "Y") ## 24 regular chrs
