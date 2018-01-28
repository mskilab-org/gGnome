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



test_that('gGraph, nullGGraph', {

    foo = gGraph$new(segs=test_segs, es=test_es)
    expect_equal(dim(foo$edges)[1], 12)
    expect_equal(dim(foo$edges)[2], 16)
    expect_equal(max((foo$edges)$cn), 3)
    expect_equal(max((foo$edges)$fromStart), 18593415)
    expect_equal(max((foo$edges)$fromEnd), 18793414)
    expect_equal(dim((foo$nullGGraph())$edges)[1], 0)
    expect_equal(dim((foo$nullGGraph())$edges)[2], 3)
    ## check 'regular'
    ## foo$nullGGraph(regular=FALSE)  

})



## dipGraph = function(genome = NULL, chr=FALSE, regular=TRUE)
test_that('gGraph, dipGraph', {

    foo = gGraph$new(segs=test_segs, es=test_es)
    expect_error(foo$dipGraph())  ## Error in match(x, table, nomatch = 0L) : object 'regularChr' not found


})









test_that('gread', {

    expect_error(gread('no_file_here'))
    ## gread('jabba.simple.rds')
    ## 
})


##
## gtf2json = function(gtf=NULL, gtf.rds=NULL, gtf.gr.rds=NULL, filename="./gtf.json",
##                    genes=NULL, grep=NULL, grepe=NULL, chrom.sizes=NULL, include.chr=NULL,
##                    gene.collapse=TRUE, verbose = TRUE)
##                    
test_that('gtf2json', {

    expect_error(gread('no_file_here'))
    ## gread('jabba.simple.rds')
    ## 
})




test_that('isInteger', {

    expect_equal(isInterger('hey'), FALSE)
    expect_equal(isInterger(2), TRUE)
    expect_equal(isInterger(2.2), FALSE)
})





test_that('getPloidy', {

    expect_error(getPloidy(GRangesList()))
    expect_equal(round(getPloidy(segments), 3), 3.817)
})




test_that('grl.duplicated', {

    expect_error(grl.duplicated(GRangesList()))   ### I don't think this should give an error
    expect_equal(round(getPloidy(segments), 3), 3.817)
})



test_that('setxor', {

    A = c(1, 2, 3)
    B = c(1, 4, 5)
    expect_equal(setxor(A, B), c(2, 3, 4, 5))
})



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



gencode2json = function(gencode=NULL, file="."){
    ## ASSUMPTION: gencode is a GR, presumably read from skidb function
    if (is.null(gencode)){
        require(skidb)
        ## ALERT: if you don't give me anything, I'm only including known genes
        gencode = skidb::read_gencode()
    }
}

test_that('gencode2json', {

    A = c(1, 2, 3)
    B = c(1, 4, 5)
    expect_equal(setxor(A, B), c(2, 3, 4, 5))
})




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





