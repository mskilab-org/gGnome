.onAttach <- function(libname, pkgname) {
  packageStartupMessage("gGnome: genomic structural variations as graph")
}

.onLoad <- function(libname, pkgname) {
    Sys.setenv(_R_CHECK_FORCE_SUGGESTS_=FALSE)
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
## just DEFAULT values that user could set
Sys.setenv(DEFAULT_BSGENOME=system.file("extdata", "human_g1k_v37.chrom.sizes", package="gGnome"))
Sys.setenv(DEFAULT_REGULAR_CHR=system.file("extdata", "human_g1k_v37.regular.chrom.sizes", package="gGnome"))
