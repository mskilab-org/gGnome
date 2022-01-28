library(gUtils)
library(gTrack)
library(testthat)

## read input
jab = system.file('extdata', 'jabba.simple.rds', package="gGnome")
jab.gg = gG(jabba = jab)

## use just a subgraph for faster walks
sg = jab.gg$copy$subgraph(seed = GRanges(seqnames = "1",
                                         ranges = IRanges(start = 1,
                                                          end = seqlengths(jab.gg$gr)["1"])
                                         ))

test_that("testing peel function", {

    ## test that circular walks are returned if embed.loops is FALSE
    wks = peel(sg, field = NULL, embed.loops = FALSE, verbose = FALSE)
    expect_true(any(wks$dt$circular == TRUE, na.rm = TRUE))

    ## test that loops are embedded if embed.loops is TRUE
    wks.embed = peel(sg, field = NULL, embed.loops = TRUE, verbose = FALSE)
    expect_true(all(wks.embed$dt$circular == FALSE, na.rm = TRUE))

    ## test that walk copy numbers without embed loops sum to original node copy numbers
    wks$set(walk.cn = wks$dt$cn)
    gr = stack(wks$grl) %>% gr.sum(field = "walk.cn")
    gr = gr %*% sg$nodes$gr[, "cn"]
    expect_true(all(gr$cn == gr$walks.cn))

    ## same with embedded loops
    wks.embed$set(walk.cn = wks.embed$dt$cn)
    gr = stack(wks.embed$grl) %>% gr.sum(field = "walk.cn")
    gr = gr %*% sg$nodes$gr[, "cn"]
    expect_true(all(gr$cn == gr$walks.cn))
    
})

