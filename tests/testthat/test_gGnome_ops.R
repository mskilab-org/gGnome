context('testing gGnome')

library(gGnome)
library(testthat)
library(gUtils)

##-------------------------------------------------------##
test_that('constructors and essential functions', {
    ## small example, nested tDUP
    message("Toy segments: ", system.file('extdata', 'testing.segs.rds', package="gGnome"))
    test_segs = readRDS(system.file('extdata', 'testing.segs.rds', package="gGnome"))
    message("Toy edges: ", system.file('extdata', 'testing.es.rds', package="gGnome"))
    test_es = readRDS(system.file('extdata', 'testing.es.rds', package="gGnome"))
    ## default
    expect_equal(dim(etype(test_segs, test_es))[1], 12)
    expect_equal(dim(etype(test_segs, test_es))[2], 16)
    expect_equal(unique(as.integer(etype(test_segs, test_es)$toChr)), 5)
    expect_equal(any(etype(test_segs, test_es)$fromLoose), FALSE)
    expect_equal(any(etype(test_segs, test_es)$toLoose), FALSE)
    expect_error(etype(GRangesList(), GRangesList()))  ## Error in etype(GRangesList(), GRangesList()) : Error:segs must be GRanges
    expect_error(etype(GRanges(), GRangesList()))      ## Error in etype(GRanges(), GRangesList()) : Error:es must be data.frame
    expect_error(etype(GRanges(), data.table()))       ## Error: 'from' & 'to' must be in es!
    expect_error(gGraph$new(), NA)  ## test it works
    foo = gGraph$new(segs=test_segs, es=test_es)
    expect_equal(dim(foo$edges)[1], 12)
    expect_equal(dim(foo$edges)[2], 16)
    expect_equal(max((foo$edges)$cn), 3)
    expect_equal(max((foo$edges)$fromStart), 18593415)
    expect_equal(max((foo$edges)$fromEnd), 18793414)
    foo = gGraph$new(segs=test_segs, es=test_es)
    expect_equal(dim(foo$edges)[1], 12)
    expect_equal(dim(foo$edges)[2], 16)
    expect_equal(max((foo$edges)$cn), 3)
    expect_equal(max((foo$edges)$fromStart), 18593415)
    expect_equal(max((foo$edges)$fromEnd), 18793414)
    expect_equal(dim((foo$nullGGraph())$edges)[1], 0)
    expect_equal(dim((foo$nullGGraph())$edges)[2], 3)
})

##-------------------------------------------------------##
test_that('gGraph, dipGraph', {
    expect_error(gGraph$new()$dipGraph(), NA)
    expect_equal(nrow(gGraph$new()$dipGraph()$edges), 0)
    expect_equal(length(gGraph$new()$dipGraph()$segstats), length(gUtils::hg_seqlengths())*2)
})

##-------------------------------------------------------##
test_that('karyograph', {
    jab = system.file('extdata', 'jabba.simple.rds', package="gGnome")
    message("JaBbA result: ", jab)
    segments = readRDS(jab)$segstats
    junctions = readRDS(jab)$junctions
    ## init with only tile
    expect_true(inherits(kag.tile <<- gGraph$new(tile = segments), "gGraph"))
    ## expect_equal(length(kag.tile$segstats), 2220)
    ## expect_equal(kag.tile$edges[, sum(type=="reference")/2], 1086)
    ## ## init with only junc
    expect_true(inherits(kag.junc <<- gGraph$new(junc = junctions), "gGraph"))
})

##-------------------------------------------------------##
test_that('gread', {
    jab = system.file('extdata', 'jabba.simple.rds', package="gGnome")
    message("JaBbA result: ", jab)
    prego = system.file('extdata', 'intervalFile.results', package='gGnome')
    message("PREGO results: ", prego)
    weaver = system.file('extdata', 'weaver', package='gGnome')
    message("Weaver results: ", weaver)
    expect_error(gread('no_file_here'))
    expect_true(inherits(jab <<- gread(jab), "bGraph"))
    expect_true(inherits(prego <<- gread(prego), "gGraph"))
    expect_true(inherits(wv <<-gread(weaver), "gGraph"))
})

##-------------------------------------------------------##
test_that('gtf2json', {
    expect_error(gread('no_file_here'))
    expect_equal(gtf2json(system.file('extdata', 'test.gtf', package='gGnome')), "./gtf.json")
    system(paste('rm', "./gtf.json"))
})

##-------------------------------------------------------##
test_that('setxor', {
    A = c(1, 2, 3)
    B = c(1, 4, 5)
    expect_equal(setxor(A, B), c(2, 3, 4, 5))
})

##-------------------------------------------------------##
test_that('special ranges functions for skew-symmetric graph', {
    jab = system.file('extdata', 'jabba.simple.rds', package="gGnome")
    message("JaBbA result: ", jab)
    segments = readRDS(jab)$segstats
    junctions = readRDS(jab)$junctions
    expect_equal(length(seg.fill(GRanges())), 0)
    expect_equal(length(seg.fill(segments)), 2346)
    ## check 'verbose'
    expect_equal(length(seg.fill(segments, verbose=TRUE)), 2346)
    expect_equal(length(seg.fill(segments %Q% (strand=="+"), verbose=TRUE)), 2346)
    expect_equal(dim(hydrogenBonds(segments))[1], length(segments))
    expect_equal(dim(hydrogenBonds(segments))[2], 3)
    expect_equal(unique(hydrogenBonds(segments)$type), 'hydrogen')
})

##-------------------------------------------------------##
test_that('gWalks', {
    jab = system.file('extdata', 'jabba.simple.rds', package="gGnome")
    message("JaBbA result: ", jab)
    segments = readRDS(jab)$segstats
    junctions = readRDS(jab)$junctions
    grl = system.file("extdata", "gw.grl.rds", package="gGnome")
    message("Walks for testing:", grl)
    grl = readRDS(grl)
    expect_equal(length(gw <<- as(grl, "gWalks")), length(grl))
    expect_error(bg <<- as(gw, "bGraph"), NA)
    expect_equal(length(bg$junctions), sum(values(junctions)$cn>0))
    expect_true(inherits(gw.simp <<- gw$simplify(mod=FALSE), "gWalks"))
    expect_error(bg.simp <<- as(gw.simp, "bGraph"), NA)
    expect_error(bg.dc <<- bg.simp$decouple(), NA)
    expect_equal(length(bg.dc$junctions), length())
})
