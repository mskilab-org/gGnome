library(testthat)
library(gUtils)

setDTthreads(1)

test_that('read_xmap', {
    # one test for read_xmap is already included in the haplograph tests
    # here we only test compliment functions not tested there
    bionano.fn = system.file('extdata', 'H838_rare_variant_pipeline_EXP_REFINEFINAL1_bionano.xmap.gz', package = 'gGnome')

    # test win and grl = FALSE #
    win = '1:627999-1213255'
    gr = bionano.fn %>% read_xmap(merge = FALSE, verbose = TRUE, win = win, grl = FALSE)
    expect_true(gr.string(gr.stripstrand(gr %>% reduce)) == win)

    # test seqlevels #
    seqlevels.dict = c(1:22, 'X', 'Y')
    names(seqlevels.dict) = 1:24
    gr1 = bionano.fn %>% read_xmap(merge = FALSE, verbose = TRUE, seqlevels = seqlevels.dict)
    expect_true(all(sort(seqlevels(gr1)) == sort(as.character(seqlevels.dict))))

}
