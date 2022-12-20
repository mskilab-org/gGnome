library(testthat)
library(gUtils)

setDTthreads(1)

test_that(desc = "test read.juncs for svaba vcf",
          code = {
              fn = system.file("extdata", "svaba.subset.vcf", package = "gGnome")
              svaba.all = gGnome:::read.juncs(fn)
              svaba.std = gGnome:::read.juncs(fn, standard.only = TRUE)
              svaba.chr = gGnome:::read.juncs(fn, standard.only = TRUE, chr.convert = TRUE)
              ## make sure correct number of junctions is read in
              expect_true(length(svaba.all) == 20)
              expect_true(length(svaba.std) == 20)
              expect_true(length(svaba.chr) == 20)
              ## make sure seqlengths only includes standard chr
              expect_true(all(names(seqlengths(svaba.std)) %in% paste0("chr", c(as.character(1:22), "X", "Y"))))
              ## make sure chr prefex is removed
              expect_true(!any(grepl("chr", names(seqlengths(svaba.chr)))))
          })

test_that(desc = "test ONT VCF",
          code = {
              fn = system.file("extdata", "ont.subset.vcf", package = "gGnome")
              ont.all = gGnome:::read.juncs(fn)
              ont.std = gGnome:::read.juncs(fn, standard.only = TRUE)
              expect_true(length(ont.all) == 25)
              expect_true(length(ont.std) < 25)
              expect_true(all(names(seqlengths(ont.std)) %in% c(as.character(1:22), "X", "Y")))
              bp1 = grl.pivot(ont.std)[[1]]
              bp2 = grl.pivot(ont.std)[[2]]
              expect_true(all(as.character(seqnames(bp1)) %in% c(as.character(1:22), "X", "Y")))
              expect_true(all(as.character(seqnames(bp2)) %in% c(as.character(1:22), "X", "Y")))
          })

test_that(desc = "test reading Junction object",
          code = {
              fn = system.file("extdata", "svaba.subset.junction.rds", package = "gGnome")
              svaba.all = gGnome:::read.juncs(fn)
              expect_true(length(svaba.all) == 20)
          })

test_that(desc = "test reading bedpe",
          code = {
              fn = system.file("extdata", "ont.subset.bedpe", package = "gGnome")
              ont.bpe = gGnome:::read.juncs(fn)
              ont.std = gGnome:::read.juncs(fn, standard.only = TRUE)
              expect_true(length(ont.bpe) == 25)
              expect_true(length(ont.std) < 25)
          })


test_that('jab2gg', {
    # test allelic annotation of jabba
    jabba = readRDS(system.file('extdata/hcc1954', 'jabba.rds', package="gGnome"))
    ccols = c("cn", "edges.in", "edges.out", "tile.id")
    # random values for cn.low and cn.high
    asegstats.low = gr2dt(jabba$segstats[, ccols])
    asegstats.low[, parent := .I]
    asegstats.high = asegstats.low
    asegstats.low[!is.na(cn), cn := floor(runif(1, 0, pmax(0, cn/2 - 1)))]
    asegstats.high$cn = jabba$segstats$cn - asegstats.low$cn
    asegstats.low[, type := 'low']
    asegstats.high[, type := 'high']
    # notice that this fake asegstats does not have all the fields that a real asegstats will have (for example "phased" is missing)
    jabba$asegstats = c(dt2gr(asegstats.low), dt2gr(asegstats.high))
    gg = gG(jabba = jabba)
    expect_true(gg$nodes$dt[!is.na(cn), all(cn == cn.low + cn.high)])
})

test_that('jab2gg with purity', {
    # test allelic annotation of jabba
    jabba = gG(jabba = system.file("extdata", "hcc1954.rigma.sg.with.pp.rds", package = "gGnome"))
    expect_true(!is.null(jabba$meta$purity))
    expect_true(!is.null(jabba$meta$ploidy))
    expect_true(!is.na(jabba$meta$purity))
    expect_true(!is.na(jabba$meta$ploidy))
})

test_that('cougar1gg', {
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


test_that('haplograph', {
    setDTthreads(1)
    bionano.fn = system.file('extdata', 'H838_rare_variant_pipeline_EXP_REFINEFINAL1_bionano.xmap.gz', package = 'gGnome')
    # use utils.R read_xmap to read the xmap to GRanges and generate a gWalk
    bgw = bionano.fn %>% read_xmap(merge = FALSE) %>% gW(grl = .)
    # generate haplograph from the walks
    bgg = gG(walks = bgw)
    expect_true(inherits(bgg, 'gGraph'))
})

test_that('jab2gg', {
    setDTthreads(1)
    # almost everything is tested elsewhere. just testing some edge cases here
    expect_error(jab2gg(list()))
    expect_error(jab2gg('no.such.file.rds'))
    expect_error(jab2gg(123))
    empty_gg = jab2gg(gG()) # this should work
})

test_that('read.juncs', {
    setDTthreads(1)
    # almost everything is tested elsewhere. just testing some edge cases here
    expect_error(read.juncs(NA))
    expect_error(jab2gg('no.such.file.rds'))
    expect_error(jab2gg(123))
    empty_gg = jab2gg(gG()) # this should work
})

test_that('karyotype', {
    setDTthreads(1)
    # almost everything is tested elsewhere. just testing some edge cases here
    expect_true(inherits(karyotype(), 'gGraph'))
    expect_error(karyotype(karyo = 'not null')) # since this is not implemented yet we expect error
})

test_that('alignments2gg', {
    setDTthreads(1)
    # download lightweight (~40Mb) hg19 BAM file
    system('wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeUwRepliSeq/wgEncodeUwRepliSeqK562G1AlnRep1.bam')
    system('wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeUwRepliSeq/wgEncodeUwRepliSeqK562G1AlnRep1.bam.bai')
    reads = bamUtils::read.bam('wgEncodeUwRepliSeqK562G1AlnRep1.bam', intervals = 'chr1')
    expect_true(inherits(gG(alignments = reads), 'gGraph'))
    # cleanup
    system('rm wgEncodeUwRepliSeqK562G1AlnRep1.bam*')
})
