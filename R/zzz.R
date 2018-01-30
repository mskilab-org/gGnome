.onAttach <- function(libname, pkgname) {
  packageStartupMessage("gGnome: genomic structural variations as graph")
}

.onLoad <- function(libname, pkgname) {
  op <- options()
  op.gGnome <- list(
      gGnome.verbose = FALSE
      ## devtools.path = "~/R-dev",
      ## devtools.install.args = "",
      ## devtools.name = "Your name goes here",
      ## devtools.desc.author = '"First Last <first.last@example.com> [aut, cre]"',
      ## devtools.desc.license = "What license is it under?",
      ## devtools.desc.suggests = NULL,
      ## devtools.desc = list()
  )
  toset <- !(names(op.gGnome) %in% names(op))
  if(any(toset)) options(op.gGnome[toset])

  invisible()
}

## preset environmental variables
## GENOME = readRDS(system.file("extdata", "hg19.broad.BSgenome.rds", package="gGnome"))
## regularChr = c(as.character(1:22), "X", "Y") ## 24 regular chrs
Sys.setenv(DEFAULT_BSGENOME=system.file("extdata", "human_g1k_v37.chrom.sizes", package="gGnome"))
strmap = setNames(c("+", "-"), c("-", "+"))
