## .onAttach <- function(libname, pkgname) {
##   packageStartupMessage("loaded gGnome")
## }

## .onLoad <- function(libname, pkgname) {
##   op <- options()
##   op.gGnome <- list(
##       gGnome.verbose = TRUE,
##       gGnome.debug = FALSE
##       ## devtools.path = "~/R-dev",
##       ## devtools.install.args = "",
##       ## devtools.name = "Your name goes here",
##       ## devtools.desc.author = '"First Last <first.last@example.com> [aut, cre]"',
##       ## devtools.desc.license = "What license is it under?",
##       ## devtools.desc.suggests = NULL,
##       ## devtools.desc = list()
##   )
##   toset <- !(names(op.gGnome) %in% names(op))
##   if(any(toset)) options(op.gGnome[toset])

##   ## preset environmental variables
##   ## just DEFAULT values that user could set
##   Sys.setenv(DEFAULT_BSGENOME=system.file("extdata", "human_g1k_v37.chrom.sizes", package="gGnome"))
##   Sys.setenv(DEFAULT_REGULAR_CHR=system.file("extdata", "human_g1k_v37.regular.chrom.sizes", package="gGnome"))
##   Sys.setenv(DEFAULT_GENE_ANNOTATION="ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_27/GRCh37_mapping/gencode.v27lift37.basic.annotation.gff3.gz")
##   Sys.setenv(DEFAULT_GGNOMEJS="~/git/gGnome.js")
##   invisible()
## }
