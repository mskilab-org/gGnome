library(testthat)
library(gUtils)

setDTthreads(1)

###########################
## test binstats and unphased balance
##
###########################

message("Reading test subgraph")
sg = readRDS(system.file("extdata", "hcc1954.rigma.sg.rds", package = "gGnome"))

message("Reading test coverage")
sg.cov = readRDS(system.file("extdata", "hcc1954.rigma.sg.cov.rds", package = "gGnome"))

message("Reading expected results")
sg.gr = readRDS(system.file("extdata", "hcc1954.rigma.sg.cn.rds", package = "gGnome"))
sg.gr$expected = sg.gr$cn

message("Testing binstats")
pl = weighted.mean(sg$nodes$dt$cn, sg$nodes$dt$width, na.rm = TRUE)
binstats.sg = binstats(sg, bins = sg.cov, field = "ratio", purity = 1, ploidy = pl, lp = TRUE)

res = binstats.sg$nodes$gr[, "cn"] %$% sg.gr[, "expected"]
expect_equal(round(res$cn), res$expected, tolerance = 1e-2)

qp.binstats.sg = binstats(sg, bins = sg.cov, field = "ratio", purity = 1, ploidy = pl, lp = FALSE) ## this should only change node weights and not cnmle

res = qp.binstats.sg$nodes$gr[, "cn"] %$% sg.gr[, "expected"]
expect_equal(round(res$cn), res$expected, tolerance = 1e-2)

message("Testing unphased LP balance")
lp.bal.sg = balance(binstats.sg, lambda = 10, epgap = 1e-6, tilim = 60, lp = TRUE, verbose = 2)

res = lp.bal.sg$nodes$gr[, "cn"] %$% sg.gr[, "expected"]
expect_equal(res$cn, res$expected, tolerance = 1e-2)

message("Testing unphased QP balance")
qp.bal.sg = balance(binstats.sg, lambda = 10, epgap = 1e-6, tilim = 60, lp = FALSE, verbose = 2)

res = qp.bal.sg$nodes$gr[, "cn"] %$% sg.gr[, "expected"]
expect_equal(res$cn, res$expected, tolerance = 1e-2)

####################################
## test binstats and phased balance
##
####################################

message("Reading subgraph hets")
sg.hets = readRDS(system.file("extdata", "hcc1954.rigma.sg.hets.rds", package = "gGnome"))

message("Reading phased graph expected CN")
phased.gr = readRDS(system.file("extdata", "hcc1954.rigma.phased.cn.rds", package = "gGnome"))
phased.gr$expected = phased.gr$cn

message("Checking phased binstats")
phased.binstats.sg = phased.binstats(sg, bins = sg.hets, count.field = "count", allele.field = "allele", min.bins = 3, min.var = 1e-3, purity = 0.94, ploidy = pl)

res = gr.findoverlaps(phased.binstats.sg$nodes$gr %Q% (!is.na(weight)),
                      phased.gr,
                      qcol = c("cn", "allele"),
                      scol = c("expected", "allele"),
                      by = "allele")

expect_equal(round(res$cn), res$expected)

message("Checking phased balance with marginals")

mg = sg$nodes$gr[, "cn"]
mg$fix = 1
mg$weight = 1

phased.bal.sg = balance(phased.binstats.sg, lambda = 10, epgap = 1e-6, tilim = 60, lp = TRUE, verbose = 2, phased = TRUE, marginal = mg, ism = FALSE)

res = gr.findoverlaps(phased.bal.sg$nodes$gr %Q% (!is.na(weight)),
                      phased.gr,
                      qcol = c("cn", "allele"),
                      scol = c("expected", "allele"),
                      by = "allele")

expect_equal(res$cn, res$expected)

message("Checking phased balance without marginals")
phased.bal.sg.nomarginal = balance(phased.binstats.sg, lambda = 10, epgap = 1e-6, tilim = 60, lp = TRUE, verbose = 2, phased = TRUE, ism = FALSE)

res = gr.findoverlaps(phased.bal.sg.nomarginal$nodes$gr %Q% (!is.na(weight)),
                      phased.gr,
                      qcol = c("cn", "allele"),
                      scol = c("expected", "allele"),
                      by = "allele")

expect_equal(res$cn, res$expected)
