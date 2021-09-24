library(testthat)
library(gUtils)

setDTthreads(1)

test_that('cougar2gg', {
    # download cougar files from cougar repo

    system('wget https://github.com/compbio-UofT/CouGaR-viz/raw/master/TCGA.samples.all_Contigs.tar.gz')
    system('tar -zxf TCGA.samples.all_Contigs.tar.gz TCGA-LN-A4A3/')

    chrom.sizes = system.file("extdata", "human_g1k_v37.regular.chrom.sizes", package="gGnome")
    Sys.setenv(DEFAULT_GENOME = chrom.sizes)
    cougar_list = cougar2gg('TCGA-LN-A4A3')
    cougar_gg = gG(breaks = cougar_list$breaks, junctions = cougar_list$juncs)
    expect_true(inherits(cougar_gg, 'gGraph'))

    expect_error(cougar2gg('not a dir'))
    expect_error(cougar2gg('./')) # not a cougar output directory
}
