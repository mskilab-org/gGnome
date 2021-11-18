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

pl = weighted.mean(sg$nodes$dt$cn, sg$nodes$dt$width, na.rm = TRUE)
binstats.sg = binstats(sg, bins = sg.cov,
                       field = "ratio", purity = 1,
                       ploidy = pl, lp = TRUE)

## this should only change node weights and not cnmle
qp.binstats.sg = binstats(sg, bins = sg.cov,
                          field = "ratio", purity = 1,
                          ploidy = pl, lp = FALSE) 


test_that(desc = "Testing binstats",
          code = {
              res = binstats.sg$nodes$gr[, "cn"] %$% sg.gr[, "expected"]
              expect_equal(round(res$cn), res$expected, tolerance = 1e-2)
              res = qp.binstats.sg$nodes$gr[, "cn"] %$% sg.gr[, "expected"]
              expect_equal(round(res$cn), res$expected, tolerance = 1e-2)
          })

## run LP balance
lp.bal.sg = balance(binstats.sg,
                    lambda = 10, epgap = 1e-6,
                    tilim = 60, lp = TRUE,
                    verbose = 2)

## check LP balanced with node ub/lb
bounded.sg = binstats.sg$copy
bounded.sg$nodes$mark(ub = 1, lb = 1)

bounded.bal.sg = balance(bounded.sg,
                    lambda = 10, epgap = 1e-6,
                    tilim = 60, lp = TRUE,
                    verbose = 2)


test_that(desc = "Testing unphased LP balance with node upper and lower bounds",
          code = {
              res = bounded.bal.sg$nodes$gr[, "cn"] %$% sg.gr[, "expected"]
              expect_true(all(res$cn == 1))
          })

test_that(desc = "Testing unphased LP balance",
          code = {
              res = lp.bal.sg$nodes$gr[, "cn"] %$% sg.gr[, "expected"]
              expect_equal(res$cn, res$expected, tolerance = 1e-2)
          })


message("Testing unphased QP balance")
qp.bal.sg = balance(binstats.sg, lambda = 10, epgap = 1e-6, tilim = 60, lp = FALSE, verbose = 2)

res = qp.bal.sg$nodes$gr[, "cn"] %$% sg.gr[, "expected"]
expect_equal(res$cn, res$expected, tolerance = 1e-2)

####################################
## test implementation of edge reward
##
####################################

ereward.sg = binstats.sg$copy

ereward.sg$edges[type == "ALT"]$mark(reward = 10)

test_that(desc = "Testing edge reward implementation with ISM",
          code = {
              expect_message(balance(ereward.sg,
                                     lambda = 10,
                                     epgap = 1e-6,
                                     tilim = 60,
                                     lp = TRUE,
                                     ism = TRUE,
                                     verbose = 2),
                             regexp = "reward")
          })

test_that(desc = "Testing edge reward implementation without ISM",
          code = {
              expect_message(balance(ereward.sg,
                                     lambda = 10,
                                     epgap = 1e-6,
                                     tilim = 60,
                                     lp = TRUE,
                                     ism = TRUE,
                                     verbose = 2),
                             regexp = "reward")
          })
          
                             
## make a pair of reciprocal edges and make sure they get added
tiny.grl = GRangesList(
    GRanges(seqnames = c('1', '8'),
            ranges = IRanges(start = c(1e6, 1e6), width = 1),
            strand = c('+', '+')),
    GRanges(seqnames = c('1', '8'),
            ranges = IRanges(start = c(1e6, 1e6), width = 1),
            strand = c('-', '-')))


## use a positive reward
reciprocal.jj = jJ(tiny.grl)
reciprocal.jj$set(reward = 100)
reciprocal.gg = gG(junctions = reciprocal.jj)
reciprocal.gg$nodes$mark(cn = 2)

test_that(desc = "Testing edge reward with positive reward",
          code = {
              res = balance(reciprocal.gg,
                            lambda = 10,
                            epgap = 1e-8,
                            tilim = 60,
                            lp = TRUE,
                            ism = TRUE,
                            verbose = 2)
              altedges = res$edges$dt[type == "ALT", cn]
              expect_true(all(altedges > 0))
          })


reciprocal.jj = jJ(tiny.grl)
reciprocal.jj$set(reward = -100)
reciprocal.gg = gG(junctions = reciprocal.jj)
reciprocal.gg$nodes$mark(cn = 2)

test_that(desc = "Testing edge reward implementation with negative reward",
          code = {
              res = balance(reciprocal.gg,
                            lambda = 10,
                            epgap = 1e-8,
                            tilim = 60,
                            lp = TRUE,
                            ism = TRUE,
                            verbose = 2)
              altedges = res$edges$dt[type == "ALT", cn]
              expect_true(all(altedges == 0))
          })

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

# don't provide purity ploidy
phased.binstats.sg = phased.binstats(sg, bins = sg.hets, count.field = "count", allele.field = "allele", min.bins = 3, min.var = 1e-3, purity = NA, ploidy = NA, verbose = TRUE)

expect_error(phased.binstats(sg, bins = 'not a GRanges', count.field = "count", allele.field = "allele", min.bins = 3, min.var = 1e-3, purity = 0.94, ploidy = pl, verbose = TRUE))

expect_error(phased.binstats(sg, bins = sg.hets, count.field = "NOT A VALID FIELD", allele.field = "allele", min.bins = 3, min.var = 1e-3, purity = 0.94, ploidy = pl, verbose = TRUE))

sg.hets.no.major = sg.hets %Q% (allele != 'major')
expect_error(phased.binstats(sg, bins = sg.hets.no.major , count.field = "count", allele.field = "allele", min.bins = 3, min.var = 1e-3, purity = 0.94, ploidy = pl, verbose = TRUE))

expect_error(phased.binstats(phase.blocks = 'NOT A GRanges', sg, bins = sg.hets, count.field = "count", allele.field = "allele", min.bins = 3, min.var = 1e-3, purity = 0.94, ploidy = pl, verbose = TRUE))

expect_warning(phased.binstats(sg, edge.phase.dt = 'Not a data.table', bins = sg.hets, count.field = "count", allele.field = "allele", min.bins = 3, min.var = 1e-3, purity = 0.94, ploidy = pl, verbose = TRUE))

data.table.with.wrong.columns = data.table()
expect_warning(phased.binstats(sg, edge.phase.dt = data.table.with.wrong.columns, bins = sg.hets, count.field = "count", allele.field = "allele", min.bins = 3, min.var = 1e-3, purity = 0.94, ploidy = pl, verbose = TRUE))

phased.binstats.sg = phased.binstats(sg, bins = sg.hets, count.field = "count", allele.field = "allele", min.bins = 3, min.var = 1e-3, purity = 0.94, ploidy = pl, verbose = TRUE)

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
