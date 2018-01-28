library(gGnome)

library(testthat)
library(gUtils)


jab = readRDS('jabba.simple.rds')  ## HCC1143

segments = jab$segs












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





