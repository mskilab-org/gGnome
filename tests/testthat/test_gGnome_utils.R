library(testthat)
library(gUtils)

setDTthreads(1)

test_that('read_xmap, read_cmap', {
    # one test for read_xmap is already included in the haplograph tests
    # here we only test compliment functions not tested there
    system('wget https://mskilab.s3.amazonaws.com/gGnome/H838-bionano.tar.gz && tar -xzvf H838-bionano.tar.gz')
    xmap.fn = './H838-bionano/EXP_REFINEFINAL1.xmap'

    # test win and grl = FALSE #
    win = '1:627999-1213255'
    gr = xmap.fn %>% read_xmap(merge = FALSE, verbose = TRUE, win = win, grl = FALSE)
    expect_true(gr.string(gr.stripstrand(gr %>% reduce)) == win)

    # test seqlevels #
    seqlevels.dict = c(1:22, 'X', 'Y')
    names(seqlevels.dict) = 1:24
    gr1 = xmap.fn %>% read_xmap(merge = FALSE, verbose = TRUE, seqlevels = seqlevels.dict)
    expect_true(all(sort(seqlevels(gr1)) == sort(as.character(seqlevels.dict))))

    expect_warning(read_cmap('/dev/null'))
    expect_warning(read_cmap('/dev/null', gr = FALSE))

    cmap.fn = './H838-bionano/EXP_REFINEFINAL1.cmap'
    cmap.gr = read_cmap(cmap.fn, gr = TRUE)
    expect_true(inherits(cmap.gr, 'GRanges'))

    # test read_xmap with merge and lift
    gr = xmap.fn %>% read_xmap(merge = TRUE, lift = TRUE, verbose = TRUE, grl = FALSE)

    expect_warning(read_xmap('/dev/null'))
    expect_warning(read_xmap('/dev/null', grl = FALSE))

})
