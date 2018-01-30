library(gGnome)

library(testthat)
library(gUtils)


jab = readRDS('jabba.simple.rds')  ## HCC1143

segments = jab$segs

junctions = jab$junc

test_segs = readRDS('testing.segs.rds')

test_es = readRDS('testing.es.rds')

## gencode.v19.annotation.nochr.head1000.gtf




### begin by testing JaBbA functions:
### 'proximity', 'karyograph', 'karyoMIP', 'karyoMIP.to.path', and 'jabba.walk


## proximity
test_that('proximity', {

    gr.grl1 = grl.unlist(grl1)
    foo = proximity(query=gr.grl1, subject=example_genes, ra=junctions)
    expect_true(is(foo, 'list'))
    expect_equal(length(names(foo)), 9)
    expect_equal(dim(foo$sum)[1], 581)
    expect_equal(dim(foo$sum)[2], 9)
    expect_equal(max((foo$sum)$j), 18187)
    expect_equal(as.integer(max((foo$sum)$query.nm)), 82)
    expect_equal(as.integer(max((foo$sum)$subject.nm)), 9996)
    expect_equal(round(max((foo$sum)$rel)), 1)
    expect_equal(max((foo$sum)$ra), 998280)
    expect_equal((foo$sum)$wt[581], 134851)
    expect_equal(length((foo$sum)$paths[[581]]), 11)
    expect_equal((foo$sum)$ab.edges[[581]], 16)
    ## check 'jab'
    foobar = proximity(query=gr.grl1, subject=example_genes, jab=jab)
    expect_true(is(foobar, 'list'))
    expect_equal(length(names(foobar)), 9)
    expect_equal(dim(foobar$sum)[1], 498)
    expect_equal(dim(foobar$sum)[2], 9)
    expect_equal(max((foobar$sum)$j), 18187)
    expect_equal(as.integer(max((foobar$sum)$query.nm)), 82)
    expect_equal(as.integer(max((foobar$sum)$subject.nm)), 9996)
    expect_equal(round(max((foobar$sum)$rel)), 1)
    expect_equal(max((foobar$sum)$ra), 998280)
    expect_equal((foobar$sum)$wt[498], 134851)
    expect_equal(length((foobar$sum)$paths[[498]]), 11)
    expect_equal((foo$sum)$ab.edges[[498]], 257)
    ## check 'verbose'
    ## check 'mc.cores'
    foo1 = proximity(query=gr.grl1, subject=example_genes, ra=junctions, verbose=TRUE, mc.cores=2)
    expect_true(is(foo1, 'list'))
    expect_equal(length(names(foo1)), 9)
    expect_equal(dim(foo1$sum)[1], 581)
    expect_equal(dim(foo1$sum)[2], 9)

})



## karyograph
test_that('karyograph', {

    foo = karyograph(junctions)
    expect_equal(length(names(foo)), 6)
    expect_true(is(foo$tile, 'GRanges'))
    expect_equal(length(foo$tile), 1332)
    expect_equal(max((foo$tile)$tile.id), 1247)
    expect_equal(length(foo$adj), 1774224)
    expect_equal(length(foo$G), 10)
    expect_equal(length(foo$ab.adj), 1774224)
    expect_equal(length(foo$ab.edges), 1770)
    expect_equal(length(foo$junctions), 295)
    ## check 'tile'
    ## check 'label.edges'
    foobar = karyograph(junctions, label.edges=TRUE)
    expect_equal(length(names(foobar)), 6)
    expect_true(is(foobar$tile, 'GRanges'))

})




## karyoMIP
## test_that('karyoMIP', {
##
##
##
## })


## karyoMIP.to.path
## test_that('karyoMIP', {
##
##
##
## })





## jabba.walk
test_that('jabba.walk', {



})






### gGraph, initialize
## initialize = function(tile=NULL, junctions=NULL, cn = FALSE,
##                                        jabba=NULL, weaver=NULL, prego=NULL,
##                                       segs=NULL, es=NULL, ploidy=NULL, purity=NULL,
##                                        regular=TRUE, rescue.balance=FALSE){
##
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

##-------------------------------------------------------##
## dipGraph = function(genome = NULL, chr=FALSE, regular=TRUE)
test_that('gGraph, dipGraph', {
    expect_error(gGraph$new()$dipGraph(), NA)
    expect_equal(nrow(gGraph$new()$dipGraph()$edges), 0)
    expect_true(all(table(seqnames(gGraph$new()$dipGraph()$segstats))==2))
})

##-------------------------------------------------------##
test_that('gread', {
    jab = system.file('extdata', 'jabba.simple.rds', package='gGnome')
    prego = system.file('extdata', 'intervalFile.results', package='gGnome')
    weaver = system.file('extdata', 'weaver', package='gGnome')
    expect_error(gread('no_file_here'))
    expect_true(inherits(gread(jab), "bGraph"))
    expect_true(inherits(gread(prego), "gGraph"))
    expect_true(inherits(gread(weaver), "gGraph"))
})

##-------------------------------------------------------##
test_that('gtf2json', {
    expect_error(gread('no_file_here'))
    expect_equal(gtf2json(system.file('extdata', 'test.gtf', package='gGnome')), "./gtf.json")
    system(paste('rm', "./gtf.json"))
})

##-------------------------------------------------------##
test_that('getPloidy', {

    expect_error(getPloidy(GRangesList()))
    expect_equal(round(getPloidy(segments), 3), 3.817)
})

##-------------------------------------------------------##
test_that('grl.duplicated', {

    expect_error(grl.duplicated(GRangesList()))   ### I don't think this should give an error
    expect_equal(round(getPloidy(segments), 3), 3.817)
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
## currently doesn't export
seg.fill = function(segs, verbose=FALSE){
    if (length(segs)==0){
        return(segs)
    }
    segs = segs[!duplicated(segs)]

    ## collapse +/- strand
    ss = unique(gr.stripstrand(segs))
    idss = match(gr.stripstrand(segs), ss)
    if (all(table(idss)==2)){
        ## if the variety matches, check the copy numbers
        ## TODO: back to this later
    }
    else {
        if (verbose){
            warning("Some segments do not have matching strand.")
        }
        fill = tofill = segs[which(idss %in% which(table(idss)==1))]
        strmap = setNames(c("+", "-"), c("-", "+"))
        strand(fill) = strmap[as.character(strand(tofill))]
        mcols(fill) = mcols(tofill)
        segs = c(segs, fill)
    }
    return(segs)
}

test_that('seg.fill', {

    expect_equal(length(seg.fill(GRanges())), 0)
    expect_equal(length(seg.fill(segments)), 2296)
    ## check 'verbose'
    expect_equal(length(seg.fill(segments, verbose=TRUE)), 2296)  ## doesn't hit verbose statements
    expect_equal(max(seg.fill(segments, verbose=FALSE)$start.ix), 77862)
    expect_equal(max(seg.fill(segments, verbose=FALSE)$end.ix), 77862)
})
##-------------------------------------------------------##

### currently doesn't export
hydrogenBonds = function(segs){
    ## collapse +/- strand
    ss = unique(gr.stripstrand(segs))
    idss = match(gr.stripstrand(segs), ss)
    if (!all(table(idss)==2)){
        stop("Error: Malformed object. Suggest creation again.")
    }
    tmpDt = data.table(ssid = seq_along(ss))
    tmpDt[, ":="(n1 = which(idss==ssid)[1],
                 n2 = which(idss==ssid)[2]), by=ssid]
    hydrogenBs = tmpDt[, .(from = n1, to = n2,
                           type="hydrogen")]
    return(hydrogenBs)
}

test_that('hydrogenBonds', {

    expect_equal(dim(hydrogenBonds(segments))[1], 1148)
    expect_equal(dim(hydrogenBonds(segments))[2], 3)
    expect_equal(unique(hydrogenBonds(segments)$type), 'hydrogen')
    expect_equal(max(hydrogenBonds(segments)$from), 2006)
    expect_equal(max(hydrogenBonds(segments)$to), 2296)
    ## check 'if (!all(table(idss)==2)){stop("Error: Malformed object. Suggest creation again.")}'

})
##-------------------------------------------------------##









### Functions:
#

## junctions?

## Class:
## gGraph
## -- initialize
## -- nullGraph
## -- dipGraph
## -- addJuncs
## -- addSegs
## -- karyograph
## -- simplify
## -- decouple
## -- add
## -- jabba2gGraph
## -- weaver2gGraph
## -- prego2gGraph
## -- print
## -- plot
## -- layout
## -- summary
## -- length
## -- gGraph2gTrack
## -- json
## -- html
## -- gGraph2json
## -- qw
## -- hydrogenBonds
## -- components
## -- subgraph
## -- fillin
## -- trim
## -- getSeqInfo
## -- makeAbEdges
## -- getAdj
## -- hood
## -- dist
## -- e2j
## -- jGraph
## -- isJunctionBalanced
## -- getLooseEnds
## -- walk

## private methods


## CLASS:
## bGraph
## --
## --
## --
## --
## --
## --
## --
## --


### FUNCTIONS
## -- gread
## -- ul (not exported)
## -- get.constrained.shortest.path (not exported)
## -- gtf2json
## -- getPloidy
## -- vaggregate
## -- mmatch
## -- alpha
## -- rev.comp
## -- etype
## -- setxor
## -- write.tab
## -- dedup
## -- isInteger
## -- hydrogenBonds
## --
## --
## --
## --
## --
## --
## --
## --
## --





