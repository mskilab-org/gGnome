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
})

test_that('rck2gg', {
    # download cougar files from cougar repo

    rck.tar = system.file('extdata', "rck.hcc1954.tar.gz", package = "gGnome")
    system(paste0('tar -xzvf ', rck.tar))
    rck.dir = 'rck.hcc1954'

    haploid.res = rck2gg(rck.dir, haploid = TRUE)
    haploid.gg = gG(nodes = haploid.res$nodes, edges = haploid.res$edges)
    expect_true(length(haploid.res$edges) == 10)
    expect_true(inherits(haploid.gg, 'gGraph'))

    haploid.res = rck2gg(rck.dir, haploid = FALSE)
    haploid.gg.not.haploid = gG(nodes = haploid.res$nodes, edges = haploid.res$edges)
    expect_true(length(haploid.res$edges) == 14)
    expect_true(inherits(haploid.gg.not.haploid, 'gGraph'))

    expect_error(rck2gg('not a dir'))
    expect_error(rck2gg('./')) # not a cougar output directory

    # make bad rck files
    tmp = paste0(tempdir(), '/bad.rck')
    dir.create(tmp)
    fwrite(data.table(V1 = 'bad', V2 = 'not better'), paste0(tmp, '/rck.scnt.tsv'))
    fwrite(data.table(V1 = 'bad', V2 = 'not better'), paste0(tmp, '/rck.acnt.tsv'))
    expect_error(rck2gg(tmp))

    # good scnt file, bad acnt file
    good_scnt = paste0(rck.dir, '/rck.scnt.tsv')
    system(paste0('cp ', good_scnt, ' ', tmp))
    expect_error(rck2gg(tmp))
})
