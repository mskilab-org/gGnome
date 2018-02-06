library(gGnome)
library(testthat)
library(gUtils)

## HCC1143 real data
jab = readRDS(system.file('extdata', 'jabba.simple.rds', package="gGnome"))  ## HCC1143
segments = jab$segs
junctions = jab$junc

## small example, nested tDUP
test_segs = readRDS(system.file('extdata', 'testing.segs.rds', package="gGnome"))
test_es = readRDS(system.file('extdata', 'testing.es.rds', package="gGnome"))

##-------------------------------------------------------##
test_that('gGraph constructor, initalize', {
    expect_error(gGraph$new(), NA)  ## test it works
    foo = gGraph$new(segs=test_segs, es=test_es)
    expect_equal(dim(foo$edges)[1], 12)
    expect_equal(dim(foo$edges)[2], 16)
    expect_equal(max((foo$edges)$cn), 3)
    expect_equal(max((foo$edges)$fromStart), 18593415)
    expect_equal(max((foo$edges)$fromEnd), 18793414)
})

##-------------------------------------------------------##
test_that('gGraph, nullGGraph', {
    foo = gGraph$new(segs=test_segs, es=test_es)
    expect_equal(dim(foo$edges)[1], 12)
    expect_equal(dim(foo$edges)[2], 16)
    expect_equal(max((foo$edges)$cn), 3)
    expect_equal(max((foo$edges)$fromStart), 18593415)
    expect_equal(max((foo$edges)$fromEnd), 18793414)
    expect_equal(dim((foo$nullGGraph())$edges)[1], 0)
    expect_equal(dim((foo$nullGGraph())$edges)[2], 3)
})


## ── 1. Failure: gGraph, dipGraph (@test_gGnome_ops.R#40)  ───────────────────────
## `gGraph$new()$dipGraph()` threw an error.
## Message: unable to find an inherited method for function 'seqinfo' for signature '"NULL"'
## Class:   simpleError/error/condition

##
## ── 2. Error: gGraph, dipGraph (@test_gGnome_ops.R#41)  ─────────────────────────
## unable to find an inherited method for function 'seqinfo' for signature '"NULL"'

##-------------------------------------------------------##
## dipGraph = function(genome = NULL, chr=FALSE, regular=TRUE)
## test_that('gGraph, dipGraph', {
##     expect_error(gGraph$new()$dipGraph(), NA)
##    expect_equal(nrow(gGraph$new()$dipGraph()$edges), 0)
##    expect_equal(length(gGraph$new()$dipGraph()$segstats), 50)
##})


## ── 3. Error: karyograph (@test_gGnome_ops.R#49)  ───────────────���───────────────
## unable to find an inherited method for function 'seqinfo' for signature '"NULL"'

##-------------------------------------------------------##
## test_that('karyograph', {
##     ## init with only tile
##     kag.tile = gGraph$new(tile = segments)
##     expect_equal(length(kag.tile$segstats), 2222)
##     expect_equal(kag.tile$edges[, sum(type=="reference")/2], 1086)
## 
##    ## init with only junc
##    kag.junc = gGraph$new(junctions = junctions)
## 
##    ## init with both
##    kag = gGraph$new(tile = segments, junctions = junctions)
## })

## ── 1. Error: gread (@test_gGnome_ops.R#80)  ────────────────────────────────────
## Input is either empty or fully whitespace after the skip or autostart. Run again with verbose=TRUE.
## 1: expect_true(inherits(gread(prego), "gGraph")) at testthat/test_gGnome_ops.R:80

##-------------------------------------------------------##
## test_that('gread', {
##     jab = system.file('extdata', 'jabba.simple.rds', package='gGnome')
##   prego = system.file('extdata', 'intervalFile.results', package='gGnome')
##   weaver = system.file('extdata', 'weaver', package='gGnome')
##   expect_error(gread('no_file_here'))
##   expect_true(inherits(gread(jab), "bGraph"))
##   expect_true(inherits(gread(prego), "gGraph"))
##   expect_true(inherits(gread(weaver), "gGraph"))
##})

##-------------------------------------------------------##
test_that('gtf2json', {
    expect_error(gread('no_file_here'))
    expect_equal(gtf2json(system.file('extdata', 'test.gtf', package='gGnome')), "./gtf.json")
    system(paste('rm', "./gtf.json"))
})

##-------------------------------------------------------##
test_that('get.ploidy', {
    expect_error(get.ploidy(GRangesList()))
    expect_true(round(get.ploidy(segments), 3)-3.817<0.2)
})

##-------------------------------------------------------##
test_that('setxor', {
    A = c(1, 2, 3)
    B = c(1, 4, 5)
    expect_equal(setxor(A, B), c(2, 3, 4, 5))
})

##-------------------------------------------------------##
test_that('etype', {
    ## default
    expect_equal(dim(etype(test_segs, test_es))[1], 12)
    expect_equal(dim(etype(test_segs, test_es))[2], 16)
    expect_equal(unique(as.integer(etype(test_segs, test_es)$toChr)), 5)
    expect_equal(any(etype(test_segs, test_es)$fromLoose), FALSE)
    expect_equal(any(etype(test_segs, test_es)$toLoose), FALSE)
    expect_error(etype(GRangesList(), GRangesList()))  ## Error in etype(GRangesList(), GRangesList()) : Error:segs must be GRanges
    expect_error(etype(GRanges(), GRangesList()))      ## Error in etype(GRanges(), GRangesList()) : Error:es must be data.frame
    expect_error(etype(GRanges(), data.table()))       ## Error: 'from' & 'to' must be in es!
})

##-------------------------------------------------------##
test_that('seg.fill', {
    expect_equal(length(seg.fill(GRanges())), 0)
    expect_equal(length(seg.fill(segments)), 2346)

    ## check 'verbose'
    expect_equal(length(seg.fill(segments, verbose=TRUE)),
                 2346)  ## doesn't hit verbose statements
        expect_equal(length(seg.fill(segments %Q% (strand=="+"), verbose=TRUE)),
                 2346)  ## doesn't hit verbose statements
})
##-------------------------------------------------------##
test_that('hydrogenBonds', {
    expect_equal(dim(hydrogenBonds(segments))[1], length(segments))
    expect_equal(dim(hydrogenBonds(segments))[2], 3)
    expect_equal(unique(hydrogenBonds(segments)$type), 'hydrogen')
})

##-------------------------------------------------------##
