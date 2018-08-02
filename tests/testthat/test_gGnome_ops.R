library(gGnome)
library(testthat)
library(gUtils)

context('testing gGnome')

message("Toy segments: ", system.file('extdata', 'testing.segs.rds', package="gGnome"))
test_segs = readRDS(system.file('extdata', 'testing.segs.rds', package="gGnome"))

message("Toy edges: ", system.file('extdata', 'testing.es.rds', package="gGnome"))
test_es = readRDS(system.file('extdata', 'testing.es.rds', package="gGnome"))

jab = system.file('extdata', 'jabba.simple.rds', package="gGnome")
message("JaBbA result: ", jab)

jab.gw.grl = system.file('extdata', 'gw.grl.rds', package = "gGnome")
message("Example JaBbA genome walks: ", jab.gw.grl)

prego = system.file('extdata', 'intervalFile.results', package='gGnome')
message("PREGO results: ", prego)

weaver = system.file('extdata', 'weaver', package='gGnome')
message("Weaver results: ", weaver)

## 

gr = GRanges(1, IRanges(c(3,7,13), c(5,9,16)), strand=c('+','-','-'), seqinfo=Seqinfo("1", 25), name=c("A","B","C"))
gr2 = GRanges(1, IRanges(c(1,9), c(6,14)), strand=c('+','-'), seqinfo=Seqinfo("1", 25), field=c(1,2))
dt = data.table(seqnames=1, start=c(2,5,10), end=c(3,8,15))


## FIXME: REQUIREMENTS FOR THE BELOW TESTS
## TEST FOR GGRAPH DEFAULT CONSTRUCTOR
## TEST FOR QUERYLOOKUP


test_that('Constructors', {
    expect_is(gGraph$new(jabba = jab), "gGraph")
    expect_is(gGraph$new(weaver = weaver), "gGraph")
    expect_is(gGraph$new(prego = prego), "gGraph")
})




test_that('gNode Class Constructor/length, gGraph length/active $nodes', {
     nodes1 = c(GRanges("1",IRanges(1,100),"*"), GRanges("1",IRanges(101,200),"*"),
                GRanges("1",IRanges(201,300),"*"), GRanges("1",IRanges(301,400),"*"),
                GRanges("1",IRanges(401,500),"*"))
     edges = data.table(n1 = c(3,2,4,1,3), n2 = c(3,4,2,5,4), n1.side = c(1,1,0,0,1), n2.side = c(0,0,0,1,0))

     gg = gGraph$new(nodes = nodes1, edges = edges)

     ## Testing some errors here
     expect_error(gNode$new(0, gg))
     expect_error(gNode$new(100, gg))
     expect_error(gNode$new(-30, gg))
     expect_error(gNode$new(4, gGraph$new()))
     expect_error(gNode$new(c(1,2,-10), gg))

     ## Testing with no snode.id
     gn = gNode$new(graph = gg)

     expect_equal(length(gn), 0)
     expect_identical(gn$graph, gg)
     expect_equal(length(gn$sid), 0)
     
     ## Testing building using active fields of gg
     gn = gg$nodes
     
     expect_is(gn, "gNode")
     expect_equal(length(gn), length(gg))
     expect_equal(gn$gr, gg$gr %Q% (strand == "+"))
     expect_equal(gn$id, (gg$gr %Q% (strand == "+"))$snode.id)
     expect_equal(gn$id, gn$sid)
})


test_that('gNode subsetting', {
    nodes1 = c(GRanges("1",IRanges(1,100),"*"), GRanges("1",IRanges(101,200),"*"),
                GRanges("1",IRanges(201,300),"*"), GRanges("1",IRanges(301,400),"*"),
                GRanges("1",IRanges(401,500),"*"))
    edges = data.table(n1 = c(3,2,4,1,3), n2 = c(3,4,2,5,4), n1.side = c(1,1,0,0,1), n2.side = c(0,0,0,1,0))

    gg = gGraph$new(nodes = nodes1, edges = edges)
    gn = gg$nodes

    ## Testing positive indicies
    expect_equal(gn[1:5]$gr, gn$gr)
    expect_equal(length(gn[1:5]), 5)
    expect_identical(gn[1:5]$graph, gg)
    
    expect_equal(gn[c(2,4)]$gr, gn$gr[c(2,4)])
    expect_equal(length(gn[c(2,4)]), 2)
    expect_equal(gn[c(2,4)]$id, c(2,4))
    
    expect_equal(gn[1]$gr, gn$gr[1])
    expect_equal(length(gn[1]), 1)
    expect_equal(gn[1]$graph, gg)

    expect_error(gn[6])
    expect_error(gn[-10])
    expect_error(gn[0])
    
    ## Indexing via snode.id name
    expect_equal(gn["4"]$gr, gg$gr[4])
    expect_equal(gn["-1"]$gr, gg$gr[6])
    expect_equal(gn[c("1","-3")]$gr, gg$gr[c(1,8)])

    expect_error(gn["10"])
    expect_error(gn["-6"])
    
    ## Some complex indexing
    expect_equal(gn["1"][-1], gn[-1])
    expect_equal(gn[2:4]["-2"][-1], gn[2])
    expect_equal(gn["4"][-1][-1], gn["4"])
    expect_equal(gn[c(-1,-4)]$gr, gg$gr[c(6,9)])
    
    ## Indexing via queries
    expect_equal(gn[loose == FALSE], gn)
    expect_equal(gn[snode.id > 3], gn[4:5])
    expect_equal(gn[start > 150 & end < 450 & loose == FALSE], gn[3:4])
})



test_that('Junction', {
    juncs = readRDS(jab)$junctions
    badjuncs = GRangesList(GRanges("1", IRanges(1,100), "+"))
    
    ## Some errors
    expect_error(Junction$new(GRanges("1", IRanges(1,100), "+")))
    expect_error(Junction$new(badjuncs))

    ## Empty Junctions
    jj = Junction$new(GRangesList())
    expect_equal(length(jj), 0)
    
    ## Build
    jj = Junction$new(juncs)

    expect_equal(length(jj), 500)
    expect_equal(unlist(jj$grl), unlist(juncs))
    expect_equal(jj$dt, as.data.table(values(juncs)))

    ## c() / +
    jj2 = c(jj, jj)
    expect_equal(length(jj2), length(jj)*2)
    expect_equal(jj2$grl, c(jj$grl, jj$grl))

    jj2 = jj + 100
    expect_true(all(all(width(jj2$grl) == 100)))
    
    ## $breakpoints
    bps = jj$breakpoints
    bps = c(GenomicRanges::shift(bps %Q% (strand == "+"), 1), bps %Q% (strand == "-"))
    expect_equal(length(findOverlaps(unlist(juncs), bps)), length(juncs)*2)
    
    bps = jj2$breakpoints
    bps = c(GenomicRanges::shift(bps %Q% (strand == "+"), 1), bps %Q% (strand == "-"))
    expect_equal(length(findOverlaps(unlist(juncs), bps)), length(juncs)*2)

    ## $gGraph
    gg = jj$graph
    gg1 = gGraph$new(juncs = juncs)

    expect_is(gg, "gGraph")
    expect_equal(length(gg), length(gg1))
})

test_that('gGraph, empty constructor/length', {
    ## Test with unspecified genome - seqinfo(default values)
    gg = gGraph$new()
    expect_true(is(gg, 'gGraph'))
    
    ## Check the active bindings
    expect_equal(length(seqinfo(gg)), 0)
    
    ## Test with specified genome - valid - just check seqinfo (set.seqinfo w/valid genome)
    ##gg = gGraph$new(genome = hg_seqlengths())
    ##expect_equal(length(gg$seqinfo), 25)
})

## TESTING TRIM
test_that('gGraph, trim', {
    ## 1) trim within single node
    ## 2) trim across multiple nodes
    ## 3) Trim is not continuous

    ## Build a gGraph
    gr = c(GRanges("1", IRanges(1001,2000), "*"), GRanges("1", IRanges(2001,3000), "*"),
           GRanges("1", IRanges(3001,4000), "*"), GRanges("1", IRanges(4001,5000), "*"))
    gr = gr.fix(gr, hg_seqlengths())
    
    es = data.table(n1 = c(1,2,3,4,4), n2 = c(2,3,4,1,3), n1.side = c(1,1,1,0,1), n2.side = c(0,0,0,1,0))
    
    graph = gGraph$new(nodes = gr, edges = es)

    ## CASE 1
    gr1 = GRanges("1", IRanges(1200,1500), "+")
    gr1 = gr.fix(gr1, hg_seqlengths())

    tmp = graph$trim(gr1)

    expect_equal(streduce(tmp$nodes$gr), streduce(gr1))
    expect_equal(granges(tmp$nodes$gr), granges(gr1))
    expect_equal(length(gg$edges), 0)
    expect_equal(length(tmp), 1)
    
    ## CASE 2
    gr2 = GRanges("1", IRanges(2200,4500), "+")
    gr2 = gr.fix(gr2, hg_seqlengths())
    edges = data.table(n1 = c(1,2), n2 = c(2,3), n1.side = c(1,1), n2.side = c(0,0))

    tmp2 = graph$trim(gr2)

    expect_equal(streduce(tmp2$nodes), streduce(gr2))

    es = tmp2$edges[order(n1,n2,n1.side,n2.side),][, c("type") := NULL]
    expect_equal(es, edges[order(n1,n2,n1.side,n2.side),])
    expect_equal(length(tmp2), 3)
   
    ## Case 3
    ranges(gr2) = IRanges(1800,2200)
    gr2 = c(gr2, GRanges("1", IRanges(2800,4500), "+"))
    edges = data.table(n1 = c(2,4,5,6), n2 = c(3,5,6,2), n1.side = c(1,1,1,0), n2.side = c(0,0,0,1))
    
    graph$trim(c(gr1,gr2), mod=T)
    
    expect_equal(streduce(graph$nodes), streduce(c(gr1,gr2)))

    es = graph$edges[order(n1,n2,n1.side,n2.side),][, c("type") := NULL]
    expect_equal(es, edges[order(n1,n2,n1.side,n2.side),])
    expect_equal(length(graph), 6)
    
})


## ## FIXME: Definitely could add more tests to catch errors and things in constructor but it seems to be working
## test_that('gGraph, gNodes and Edges Constructor/active bindings/loosegNodes', {
##     ## Make sure it add the right nodes and edges
##     ## 1) edges = NULL
##     ##     a) looseterm = FALSE
##     ##     b) looseterm = TRUE
##     ## 2) edges != NULL
##     ##     a) looseterm = FALSE
##     ##     b) looseterm = TRUE
##     ## - check edges, nodes, seqinfo, loosegNodes

##     ## CASE 1a
##     nodes1 = c(GRanges("1",IRanges(1,100),"*"), GRanges("1",IRanges(101,200),"*"),
##                GRanges("1",IRanges(201,300),"*"), GRanges("1",IRanges(301,400),"*"),
##                GRanges("1",IRanges(401,500),"*"))
    
##     ## Set the seqinfo to make sure it carries over
##     nodes1 = gUtils::gr.fix(nodes1, hg_seqlengths())
##     seq = hg_seqlengths()
##     names(seq) = NULL

##     ## Check the null case
##     gg = gGraph$new(nodes = GRanges(), looseterm = FALSE)

##     expect_equal(gg$nodes, GRanges())
##     expect_equal(dim(gg$edges)[1], 0)
##     expect_error(gg$graph, NA)
##     expect_equal(length(seqinfo(gg)), 0)
    
##     ## Check with looseterm = F
##     gg = gGraph$new(nodes = nodes1, looseterm = FALSE)

##     expect_equal(sort(granges(gg$nodes)), sort(granges(nodes1)))
##     expect_equal(dim(gg$edges)[1], 0)
##     expect_equal(length(gg$loosegNodes()), 0)
    
##     ## CASE 1b
##     gg = gGraph$new(nodes = nodes1, looseterm = TRUE)

##     loosenodes = c(GRanges("1",IRanges(1,1),"*"), GRanges("1",IRanges(100,100),"*"),
##                    GRanges("1",IRanges(101,101),"*"), GRanges("1",IRanges(200,200),"*"),
##                    GRanges("1",IRanges(201,201),"*"), GRanges("1",IRanges(300,300),"*"),
##                    GRanges("1",IRanges(301,301),"*"), GRanges("1",IRanges(400,400),"*"),
##                    GRanges("1",IRanges(401,401),"*"), GRanges("1",IRanges(500,500),"*"))
##     loosenodes = gUtils::gr.fix(loosenodes, hg_seqlengths())
    
##     expect_equal(sort(granges(nodes1)), sort(granges(gg$nodes)))
##     expect_equal(dim(gg$edges)[1], 0)
##     expect_equal(length(gg$loosegNodes()), 10)
##     expect_equal(sort(granges(gg$loosegNodes())), sort(granges(loosenodes)))
    
##     ## CASE 2a
##     edges = data.table(n1 = c(3,2,4,1,3), n2 = c(3,4,2,5,4), n1.side = c(1,1,0,0,1), n2.side = c(0,0,0,1,0))

##     gg = gGraph$new(nodes = nodes1, edges = edges, looseterm = FALSE)

##     expect_true(is(gg, 'gGraph'))
    
##     ## Check that the correct things are stored in our gGraph
##     expect_equal(sort(granges(nodes1)), sort(granges(gg$nodes)))

##     es = gg$edges[order(n1,n2,n1.side,n2.side),][, c("type") := NULL]
##     expect_equal(es, edges[order(n1,n2,n1.side,n2.side),])
##     expect_equal(length(gg), 5)
##     expect_equal(length(gg$loosegNodes()), 0)
##     expect_equal(seqinfo(gg)@seqlengths, seq)
    
##     ## CASE 2b
##     gg = gGraph$new(nodes = nodes1, edges = edges, looseterm = TRUE)

##     loosenodes = c(GRanges("1", IRanges(100,100), "*"),
##                    GRanges("1", IRanges(400,400), "*"), GRanges("1", IRanges(401,401), "*"))
##     loosenodes = gUtils::gr.fix(loosenodes, hg_seqlengths())
    
##     ## Check that the correct things are stored in our gGraph
##     expect_equal(sort(granges(nodes1)), sort(granges(gg$nodes)))

##     es = gg$edges[order(n1,n2,n1.side,n2.side),][, c("type") := NULL]
    
##     expect_equal(es, edges[order(n1,n2,n1.side,n2.side),])
##     expect_equal(length(gg), 5)
##     expect_equal(sort(granges(gg$loosegNodes())), sort(granges(loosenodes)))

##     gg1 = gGraph$new(nodes = gg$nodes, edges = gg$edges, looseterm = TRUE)

##     expect_equal(sort(granges(gg1$nodes)), sort(granges(gg$nodes)))
##     expect_equal(gg1$edges[order(n1,n2,n1.side,n2.side),], gg$edges[order(n1,n2,n1.side,n2.side),])
##     expect_equal(length(gg), length(gg1))
##     expect_equal(sort(granges(gg$loosegNodes())), sort(granges(loosenodes)))

## })


## ## Just checks that the jabba input is correct processed not that the result is correct
## test_that('gGraph, jabba input', {
##     expect_error(gGraph$new(jabba = "badinput"))
##     badinput = list(1,2,3)
##     names(badinput) = c("segstats", "ab.edges", "purity")
##     expect_error(gGraph$new(jabba = badinput))

##     gg = gGraph$new(jabba = jab)
##     expect_true(is(gg, "gGraph"))
    
##     gg = gGraph$new(jabba = jab, regular=T)
##     gg1 = gGraph$new()$simpleGraph()

##     expect_equal(length(gr.findoverlaps(gg$nodes, gg1$nodes)), length(gg$nodes))
## })


## ## FIXME: test the circular thing idk how to do that
## test_that('gGraph, simpleGraph', {
##     ggnew = gGraph$new()
    
##     ## simpleGraph = function(genome = NULL, chr=FALSE, include.junk=FALSE)
##     ## Testing simpleGraph(default values)
##     gg = ggnew$simpleGraph()
##     expect_true(is(gg, 'gGraph'))
##     expect_equal(length(gg$nodes), 25)
##     expect_equal(dim(gg$edges)[1], 0)
##     expect_equal(length(gg$junctions), 0)
##     expect_error(gg$graph, NA) ## check it works
##     ##expect_equal(length(gg$adj), 2500) -- FIXME: Don't know what this does but it doesn't work
##     ## expect_equal(gg$parts, NULL) -- FIXME: also broken, all don't know what it does
##     expect_equal(length(gg$seqinfo), 25)
##     expect_equal(length(gg$loosegNodes()), 50)

##     ## Don't know what the fuck is going on here
##     expect_true(is(gg$td, 'gTrack'))
##     expect_equal((gg$td)$ygap, 2)
##     expect_match((gg$td)$name, 'CN')
##     expect_equal(length(gg$win), 25)

##     ## Testing simpleGraph with genome = Something
##     gg = ggnew$simpleGraph(genome = hg_seqlengths())
##     expect_equal(length(gg$nodes), 25)
##     expect_equal(length(gg$loosegNodes()), 50)
## })





## ## Test the addSegs/makeSegs functionality
## test_that('gGraph, addSegs/karyograph no juncs', {

##     expect_equal(gGraph$new()$nodes, gGraph$new(tile=NULL)$nodes)
##     expect_equal(gGraph$new()$edges, gGraph$new(tile=NULL)$edges)
    
##     ## Check again at addSegs level
##     gg = gGraph$new()
##     expect_error(gg$addSegs(tile=NULL))
##     expect_error(gg$addSegs(tile=5))
    
##     gr = GRanges("fakeChr", IRanges(1,10000), "*")
##     expect_error(gGraph$new(tile=gr))

##     gr = c(GRanges("1", IRanges(1001,2000), "*"), GRanges("1", IRanges(2001,3000), "*"),
##            GRanges("2", IRanges(3001,4000), "*"), GRanges("3", IRanges(5000,6000), "*"),
##            GRanges("3", IRanges(5500,6500), "*"))

##     expect_error(gg$addSegs(tile=gr))

##     ## Check building a gGraph from a tile, no cn
##     gg = gGraph$new(tile = gr)
##     gg1 = gGraph$new()$simpleGraph()
##     nodes = gr.breaks(gr,gg1$nodes)    
    
##     expect_equal(length(gg), 34)
##     expect_equal(length(gg$loosegNodes()),50) 
##     expect_equal(sort(granges(gg$nodes)), sort(granges(nodes)))

##     pairs = table(gg$nodes$pair.id)
##     expect_true(all(pairs == 1))

##     edges = data.table(n1 = c(1,2,3,5,6,8,9,10,11),
##                        n2 = c(2,3,4,6,7,9,10,11,12),
##                        n1.side = c(1,1,1,1,1,1,1,1,1),
##                        n2.side = c(0,0,0,0,0,0,0,0,0))

##     expect_equal(nrow(gg$edges), 9)
##     es = gg$edges[order(n1,n2,n1.side,n2.side),][, c("type") := NULL]
##     expect_equal(es, edges[order(n1,n2,n1.side,n2.side),][, c("type") := NULL])

##     ## Try again with cn
    
## })


## test_that('gGraph, mergeGraphs', {
##     ## CASES
##     ## 1) Self is null, adding null
##     ## 2) Self is null, adding non-null ---- Could be caught with a clause a the beginning to save work
##     ## 3) Self is non-null, adding null
##     ## 4) Self is non-null, adding non-null
##     ##     a) Neither have edges
##     ##     b) Either has edges but not both
##     ##     b) Elaborate both have edges

##     gr = GRanges("1",IRanges(1,1000),"*")
    
##     graph = gGraph$new(nodes = gr)

##     ## CASE 1
##     nullgg = gGraph$new()
##     tmp = nullgg$mergeGraphs(nullgg, decouple=F)

##     expect_equal(nullgg, tmp)

##     ## CASE 2
##     tmp = nullgg$mergeGraphs(graph, decouple=F)

##     expect_equal(graph, tmp)
    
##     ## CASE 3
##     tmp = graph$mergeGraphs(nullgg, decouple=F)

##     expect_equal(graph, tmp)

##     ## CASE 4a
##     gr1 = GRanges("1", IRanges(4001,5000), "*")

##     graph1 = gGraph$new(nodes = gr1)
##     tmp = graph$mergeGraphs(graph1, decouple=F)
    
##     ## Check nodes and edges
##     expect_equal(sort(granges(c(graph$nodes,graph1$nodes))), sort(granges(tmp$nodes)))
##     expect_equal(dim(tmp$edges)[1], 0)
##     expect_equal(length(tmp), 2)

##     ## CASE 4b
##     nodes = c(GRanges("1",IRanges(1,100),"*"), GRanges("1",IRanges(101,200),"*"),
##               GRanges("1",IRanges(201,300),"*"), GRanges("1",IRanges(301,400),"*"),
##               GRanges("1",IRanges(401,500),"*"))

##     nodes1 = c(GRanges("1",IRanges(150,250),"*"), GRanges("1",IRanges(401,500),"*"),
##                GRanges("1",IRanges(601,1000),"*"))
    
##     edges = data.table(n1 = c(3,2,4,1,5), n2 = c(3,4,2,5,1), n1.side = c(1,1,0,0,1), n2.side = c(0,0,0,1,0))
##     edges1 = data.table(n1 = c(2,1,3,2), n2 = c(3,2,1,2), n1.side = c(1,1,0,0), n2.side = c(0,0,1,1))

##     ## graph has edges and graph1 does not
##     graph = gGraph$new(nodes = nodes, edges = edges)
##     graph1 = gGraph$new(nodes = nodes1)
##     tmp = graph$mergeGraphs(graph1, decouple=F)
    
##     expect_equal(sort(granges(tmp$nodes)), sort(granges(c(graph$nodes, graph1$nodes))))
    
##     es = graph$edges[order(n1,n2,n1.side,n2.side),][, c("type") := NULL]
##     expect_equal(es, tmp$edges[order(n1,n2,n1.side,n2.side),][, c("type") := NULL])
##     expect_equal(length(tmp), 8)    
    
##     ## graph does not have edges and graph1 does
##     graph = gGraph$new(nodes = nodes)
##     graph1 = gGraph$new(nodes = nodes1, edges = edges1)
##     tmp = graph$mergeGraphs(graph1, decouple=F)

##     expect_equal(sort(granges(tmp$nodes)), sort(granges(c(graph$nodes, graph1$nodes))))
    
##     es = graph1$edges[order(n1,n2,n1.side,n2.side),][, c("type") := NULL]
##     expect_equal(es, tmp$edges[order(n1,n2,n1.side,n2.side),][, c("type") := NULL])
##     expect_equal(length(tmp), 8)
    
##     ## CASE 4c
##     graph = gGraph$new(nodes = nodes, edges = edges)
##     graph1 = gGraph$new(nodes = nodes1, edges = edges1)
##     tmp = graph$clone()
##     tmp$mergeGraphs(graph1, mod=T)

##     expect_equal(sort(granges(tmp$nodes)), sort(granges(c(graph$nodes, graph1$nodes))))
    
##     es = rbind(graph$edges, graph1$edges[, ":="(n1 = n1 + length(graph), n2 = n2 + length(graph))])
##     es = es[order(n1,n2,n1.side,n2.side),][, c("type") := NULL]
##     expect_equal(es, tmp$edges[order(n1,n2,n1.side,n2.side),][, c("type") := NULL])
##     expect_equal(length(tmp), 8)

    
##     ## TESTS FOR SUBGRAPHS
##     ## Check that all the subgraphs are null
##     ##for(i in 1:length(tmp$subgraphs)) {
##     ##    expect_equal(tmp$subgraphs[[i]]$length(), 0)
##     ##}
    
##     ## Case 4b - We know the main level works, just check subgraphs
##     ##subgraphs = list(gGraph$new(tile = gr1[1])$trim(gr1[1]),
##     ##                 gGraph$new(tile = gr1[2])$trim(gr1[2]))

##     ##graph1$resetSubgraphs(subgraphs)
##     ##tmp = graph$mergeGraphs(graph1, decouple=F)
    
##     ## FIXME: insert subgraph test
    
##     ## Case 4c
##     ##subgraphs = list(gGraph$new(tile = gr[1])$trim(gr[1]),
##     ##                 gGraph$new(tile = gr[2])$trim(gr[2]))

##     ##graph$resetSubgraphs(subgraphs)
##     ##tmp = graph$mergeGraphs(graph1, decouple=F)

##     ## FIXME: insert subgraph test
    
## })


## test_that('gGraph, mergeOverlaps', {
    
## })


## test_that('gGraph, simplify', {

##     nodes = c(GRanges("1", IRanges(1001,2000), "*"), GRanges("1", IRanges(2001,3000), "*"),
##               GRanges("1", IRanges(3001,4000), "*"), GRanges("1", IRanges(4001,5000), "*"),
##               GRanges("1", IRanges(5001,6000), "*"), GRanges("1", IRanges(6001,7000), "*"),
##               GRanges("1", IRanges(7001,8000), "*"), GRanges("1", IRanges(8001,9000), "*"),
##               GRanges("1", IRanges(9001,10000), "*"))
    
##     edges = data.table(n1 = c(1,2,3,5,6,7,8,1,4,6,3),
##                        n2 = c(2,3,4,6,7,8,9,4,8,6,2),
##                        n1.side = c(1,1,1,1,1,1,1,1,1,0,1),
##                        n2.side = c(0,0,0,0,0,0,0,0,0,1,0))

##     graph = gGraph$new(nodes = nodes, edges = edges, looseterm=T)
##     graphSimple = graph$simplify(mod=F)

##     ## Check to make sure the simplifed version correctly merged everythign
##     nodes = c(GRanges("1", IRanges(1001,2000), "*"), GRanges("1", IRanges(2001,4000), "*"),
##               GRanges("1", IRanges(4001,5000), "*"),
##               GRanges("1", IRanges(5001,6000), "*"), GRanges("1", IRanges(6001,7000), "*"),
##               GRanges("1", IRanges(7001,8000), "*"), GRanges("1", IRanges(8001,10000), "*"))
##     nodes = gr.fix(nodes, hg_seqlengths())
    
##     edges = data.table(n1 = c(1,2,2,3,1,4,5,6,5),
##                        n2 = c(2,3,2,7,3,5,6,7,5),
##                        n1.side = c(1,1,1,1,1,1,1,1,0),
##                        n2.side = c(0,0,0,0,0,0,0,0,1))
    
##     expect_equal(length(graphSimple), 7)
##     expect_equal(nrow(graphSimple$edges), 9)
##     expect_equal(length(graphSimple$loosegNode s()), 3)
##     expect_equal(sort(granges(nodes)), sort(granges(graphSimple$nodes)))

##     expect_equal(graphSimple$edges[order(n1,n2,n1.side,n2.side), ][, c("type") := NULL],
##                  edges[order(n1,n2,n1.side,n2.side), ][, c("type") := NULL])
## })


## test_that('gGraph, window', {
##     ## Cases
##     ## 1) There is not a cn
##     ## 2) different pad values

##     ## CASE 1: pad = 0
##     nodes = c(GRanges("1", IRanges(400,600), "*"), GRanges("1", IRanges(500,800), "*"), GRanges("1", IRanges(1000,1400), "*"))
##     graph = gGraph$new(nodes = nodes)

##     result = c(GRanges("1", IRanges(400,800), "*"), GRanges("1", IRanges(1000,1400), "*"))
##     result = gr.fix(result, hg_seqlengths())
    
##     expect_equal(graph$window(), graph$win)
##     expect_equal(graph$win, result)
    
##     ## Case 2: pad != 0

##     result = GRanges("1", IRanges(200,1600), "*")
##     result = gr.fix(result, hg_seqlengths())
    
##     expect_equal(graph$window(200), result)
## })


## test_that('gGnome, gg2js/json', {

##     ## Make sure it throws an error when the graph is empty
##     gg = gGraph$new()
##     expect_error(gg$gg2js())

##     ## Check the empty edge case and no.y
##     nodes = c(GRanges("1", IRanges(1001,2000), "*"), GRanges("1", IRanges(2001,3000), "*"),
##               GRanges("1", IRanges(3001,4000), "*"), GRanges("1", IRanges(4001,5000), "*"),
##               GRanges("1", IRanges(5001,6000), "*"), GRanges("1", IRanges(6001,7000), "*"),
##               GRanges("1", IRanges(7001,8000), "*"), GRanges("1", IRanges(8001,9000), "*"),
##               GRanges("1", IRanges(9001,10000), "*"))

##     gg = gGraph$new(nodes = nodes, looseterm=F)
##     json = gg$gg2js(save=F, no.y=T)

##     expect_equal(nrow(json$connections), 0)
##     expect_equal(nrow(json$intervals), 9)
##     expect_false(json$settings$y_axis$visible)
    
##     ## Test saving and loading a file using jabba data - comparison file is visually pre checked
##     gg = gGraph$new(jabba = jab)
##     gg$gg2js("../../inst/extdata/data.with.cn.test.json")

##     json = fromJSON(system.file('extdata', 'data.with.cn.test.json', package="gGnome"))
##     json1 = fromJSON(system.file('extdata', 'data.with.cn.json', package="gGnome"))
    
##     expect_equal(json, json1)
## })



## ## ### some non-exported functions

## ## rev.comp = function(gr){
## ##     strmap = setNames(c("+", "-"), c("-", "+"))
## ##     if (!inherits(gr, "GRanges")){
## ##         stop("Error: Input must be GRanges.")
## ##     } else if (!all(strand(gr) %in% strmap)) {
## ##         stop("Error: Input must be all strand specific.")
## ##     }
## ##     return(rev(gr.flipstrand(gr)))
## ## }


## ## ##test_that('rev.comp works', {
## ## ##
## ## ##    expect_error(rev.comp())
## ## ##    expect_error(rev.comp(data.frame()))
## ## ##    gr3 = dt2gr(dt)
## ## ##    expect_error(rev.comp(gr3))   ## Error in rev.comp(gr3) : Input must be all strand specific.
## ## ##    expect_equal(width(rev.comp(gr)[1]), 4)
## ## ##    expect_equal(as.character(strand(rev.comp(gr)[1])), "+")
## ## ##    expect_equal(width(rev.comp(gr)[2]), 3)
## ## ##    expect_equal(as.character(strand(rev.comp(gr)[2])), "+")
## ## ##
## ## ##})


## ## capitalize = function(string, un = FALSE){
## ##     if (!un){
## ##         capped <- grep("^[^A-Z].*$", string, perl = TRUE)
## ##         substr(string[capped], 1, 1) <- toupper(substr(string[capped],1, 1))
## ##     } else{
## ##         capped <- grep("^[A-Z].*$", string, perl = TRUE)
## ##         substr(string[capped], 1, 1) <- tolower(substr(string[capped],1, 1))
## ##     }
## ##     return(string)
## ## }





## ## test_that('capitalize works', {

## ##     str1 = "Foo FOO"
## ##     str2 = "2Foo . $%@"
## ##     str3 = "foobar foo"
## ##     expect_match(capitalize(str1), "Foo FOO")
## ##     expect_match(capitalize(str1, un=TRUE), "foo FOO")
## ##     expect_equal(capitalize(str2), "2Foo . $%@")
## ##     expect_equal(capitalize(str2, un=TRUE), "2Foo . $%@")    ## probably should report this bug to r-lib
## ##     expect_match(capitalize(str3), "Foobar foo")
## ##     expect_match(capitalize(str3, un=TRUE), "foobar foo")

## ## })


  

## ## ul = function(x, n=6){
## ##     n = pmin(pmin(dim(x)), n)
## ##     return(x[1:n, 1:n])
## ## }



## ## test_that('ul works', {

## ##     A = matrix(  c(2, 4, 3, 1, 5, 7),  nrow=2, ncol=3, byrow = TRUE)    
## ##     expect_equal(as.integer(ul(A, n=0)), 2)   ### Is this expected behavior? 
## ##     expect_equal(as.integer(ul(A, n=1)), 2)
## ##     expect_equal(dim(ul(A, n=2))[1], 2)
## ##     expect_equal(dim(ul(A, n=2))[2], 2)
## ##     expect_equal(dim(ul(A, n=999))[1], 2)
## ##     expect_equal(dim(ul(A, n=9999))[2], 2)   ### Is this expected behavior? 

## ## })




## ## tile.name = function(x){
## ##     if (!inherits(x, "GRanges")){
## ##         stop("Only takes GRanges as input for now.")
## ##     }
## ##     hb = hydrogenBonds(segs = x)
## ##     if (hb[, any(is.na(from) | is.na(to))]){
## ##         stop("Not fully strand paired.")
## ##     }
## ##     hb.map = hb[, c(setNames(from, to), setNames(to, from))]
## ##     seg.name = ifelse(strand(x)=="+",
## ##                       as.character(seq_along(x)),
## ##                       paste0("-", hb.map[as.character(seq_along(x))]))
## ##     return(seg.name)
## ## }


## ## test_that('tile.name works', {

## ##     ## if (!inherits(x, "GRanges")){
## ##     expect_error(tile.name(GRangesList()))
## ##     expect_equal(length(tile.name(test_segs)), 10)

## ## })



## ## test_that('e2j works', {

## ##     expect_equal(length(e2j(test_segs, test_es)), 2)

## ## })




## ## test_that('etype works', {

## ##     expect_equal(dim(etype(test_segs, test_es))[1], 12)
## ##     expect_equal(dim(etype(test_segs, test_es))[2], 16) 

## ## })






## ## write.tab = function (x, ..., sep = "\t", quote = F, row.names = F)
## ## {
## ##     if (!is.data.frame(x)){
## ##         x = as.data.frame(x)
## ##     }
## ##     write.table(x, ..., sep = sep, quote = quote, row.names = row.names)
## ## }



## ## test_that('write.tab() works', {

## ##     expect_equal(length(write.tab(dt)), 0)  ### throws dt to STDOUT

## ## })







## ## test_that('get.ploidy works', {

## ##     expect_error(get.ploidy(GRangesList()))
## ##     expect_equal(get.ploidy(test_segs), 3)
## ##     ## if (length(cnix <- grep("CN", colnames(mcols(segs)), ignore.case = T)) == 

## ## })





## ## test_that('dedup() works', {

## ##     expect_equal(dedup(c(rep(2, 10.5), rep(3, 20)))[30], "3.20")

## ## })










## ## ## read_vcf()
## ## ## read_vcf = function(fn, gr = NULL, hg = 'hg19', geno = NULL, swap.header = NULL, verbose = FALSE, add.path = FALSE, tmp.dir = '~/temp/.tmpvcf', ...)
## ## ##test_that('read_vcf', {
## ## #    ## error
## ## #    expect_error(read_vcf('foobar'))
## ## #    ## default 
## ## #    expect_equal(length(read_vcf(somatic_vcf)), 60)
## ## #    expect_equal(length(seqnames(seqinfo(read_vcf(somatic_vcf)))), 84)
## ## #    ## gr  gr= GRanges('1:10075-10100')
## ## #    ## hg
## ## ##    expect_match(unique(as.data.frame(seqinfo(read_vcf(somatic_vcf, hg='hg12345')))$genome), 'hg12345')
## ## #    ## geno
## ## #    ## swap.header
## ## #    expect_equal(length(seqnames(seqinfo(read_vcf(somatic_vcf, swap.header='/Users/ebiederstedt/bamUtils/tests/testthat/new_header.vcf')))), 2)
## ## #    ## verbose
## ## #    expect_equal(length(read_vcf(somatic_vcf, verbose=TRUE)), 60)
## ## #    ## check 'if (!file.exists(swap.header))'
## ## #    expect_error(read_vcf(somatic_vcf, swap.header='foobar'))
## ## #
## ## #})



## ## test_that('chr2num works', {

## ##     expect_equal(as.logical(chr2num("ChrX")), NA)
## ##     expect_equal(chr2num("chrX"), 23)
## ##     expect_equal(chr2num("chrY"), 24)

## ## })



## ## test_that('affine.map works', {

## ##     expect_equal(affine.map(49), 0.5)

## ## })


## ## test_that('gr.flatmap works', {

## ##     expect_equal((gr.flatmap(example_genes, windows=GRanges('1:10000-20000'))$window.segs)$start, 1)
## ##     expect_equal((gr.flatmap(example_genes, windows=GRanges('1:10000-20000'))$window.segs)$end, 10001)
## ##     expect_equal(length((gr.flatmap(example_genes, windows=GRanges('1:10000-20000'))$window.segs)$grl.segs), 0)

## ## })




## ## ### XT's tests

## ## ##-------------------------------------------------------##
## ## test_that('constructors and essential functions', {
## ##     ## small example, nested tDUP
## ##     ## default
## ##     expect_equal(dim(etype(test_segs, test_es))[1], 12)
## ##     expect_equal(dim(etype(test_segs, test_es))[2], 16)
## ##     expect_equal(unique(as.integer(etype(test_segs, test_es)$toChr)), 5)
## ##     expect_equal(any(etype(test_segs, test_es)$fromLoose), FALSE)
## ##     expect_equal(any(etype(test_segs, test_es)$toLoose), FALSE)
## ##     expect_error(etype(GRangesList(), GRangesList()))  ## Error in etype(GRangesList(), GRangesList()) : Error:segs must be GRanges
## ##     expect_error(etype(GRanges(), GRangesList()))      ## Error in etype(GRanges(), GRangesList()) : Error:es must be data.frame
## ##     expect_error(etype(GRanges(), data.table()))       ## Error: 'from' & 'to' must be in es!
## ##     expect_error(gGraph$new(), NA)  ## test it works
## ##     foo = gGraph$new(segs=test_segs, es=test_es)
## ##     expect_equal(dim(foo$edges)[1], 12)
## ##     ## expect_equal(dim(foo$edges)[2], 16)
## ##     expect_equal(max((foo$edges)$cn), 3)
## ##     expect_equal(max((foo$edges)$fromStart), 18593415)
## ##     expect_equal(max((foo$edges)$fromEnd), 18793414)
## ##     foo = gGraph$new(segs=test_segs, es=test_es)
## ##     expect_equal(dim(foo$edges)[1], 12)
## ##     ## expect_equal(dim(foo$edges)[2], 16)
## ##     expect_equal(max((foo$edges)$cn), 3)
## ##     expect_equal(max((foo$edges)$fromStart), 18593415)
## ##     expect_equal(max((foo$edges)$fromEnd), 18793414)
## ##     expect_equal(dim((foo$nullGGraph())$edges)[1], 0)
## ##     expect_equal(dim((foo$nullGGraph())$edges)[2], 3)

## ## })

## ## test_that('karyograph', {

## ##     kag.tile = gGraph$new(tile = test_segs)
## ##     expect_true(inherits(kag.tile, "gGraph"))

## ## })




## ## ##-------------------------------------------------------##
## ## test_that('gread', {

## ##     jab = system.file('extdata', 'jabba.simple.rds', package="gGnome")
## ##     message("JaBbA result: ", jab)
## ##     prego = system.file('extdata', 'intervalFile.results', package='gGnome')
## ##     message("PREGO results: ", prego)
## ##     weaver = system.file('extdata', 'weaver', package='gGnome')
## ##     message("Weaver results: ", weaver)
## ##     expect_error(gread('no_file_here'))
## ##     jab_bgraph = gread(jab)
## ##     expect_true(is(jab_bgraph, "bGraph"))
## ##     ## preg_bgraph = gread(prego)
## ##     ## expect_true(is(preg_bgraph, "bGraph")) ### 'gGraph'
## ##     ## wv_bgraph = gread(weaver)   
## ##     ## expect_true(is(wv_bgraph, "bGraph"))  ### 'gGraph'
## ##     ## if (is.list(file)){
## ##     list_foo = gread(readRDS(system.file("extdata", "jabba.simple.rds", package="gGnome")))
## ##     expect_true(is(list_foo, 'bGraph'))
    
## ## })





## ## ## ##-------------------------------------------------------##
## ## ## test_that('gtf2json', {
## ## ##     expect_error(gread('no_file_here'))
## ## ##     expect_equal(gtf2json(system.file('extdata', 'test.gtf', package='gGnome')), "./gtf.json")
## ## ##     system(paste('rm', "./gtf.json"))
## ## ## })




## ## ## * could not find function "setxor"
## ## ##-------------------------------------------------------##

## ## test_that('setxor', {

## ##     A = c(1, 2, 3)
## ##     B = c(1, 4, 5)
## ##     expect_equal(setxor(A, B), c(2, 3, 4, 5))

## ## })





## ## ##-------------------------------------------------------##
## ## test_that('special ranges functions for skew-symmetric graph', {


## ##     segments = readRDS(jab)$segstats
## ##     junctions = readRDS(jab)$junctions
## ##     expect_equal(length(seg.fill(GRanges())), 0)
## ##     expect_equal(length(seg.fill(segments)), 2346)
## ##     ## check 'verbose'
## ##     expect_equal(length(seg.fill(segments, verbose=TRUE)), 2346)
## ##     expect_equal(length(seg.fill(segments %Q% (strand=="+"), verbose=TRUE)), 2346)
## ##     expect_equal(dim(hydrogenBonds(segments))[1], length(segments))
## ##     expect_equal(dim(hydrogenBonds(segments))[2], 3)
## ##     expect_equal(unique(hydrogenBonds(segments)$type), 'hydrogen')

## ## })


## ## ## Error: Test failed: 'gWalks'
## ## ## * length(gw <<- as(grl, "gWalks")) not equal to sum(values(grl)$cn > 0).
## ## ## 1/1 mismatches
## ## ## [1] 32 - 630 == -598
## ## ## * `bg <<- as(gw, "bGraph")` threw an error.
## ## ## Message: Error: Given edge data is not skew-symmetric!!!
## ## ## Class:   simpleError/error/condition
## ## ## * object 'bg' not found
## ## ## 1: expect_equal(length(bg$junctions), sum(values(junctions)$cn > 0)) at :12
## ## ## 2: quasi_label(enquo(object), label)
## ## ## 3: eval_bare(get_expr(quo), get_env(quo))
 
 
## ## test_that('gWalks reduce', {
## ##     ## Set up gWalks
## ##     grl = readRDS(jab.gw.grl)
## ##     gw = as(grl, "gWalks")
## ##     gw.dt = gr2dt(gw$nodes)

## ##     ## Testing reduce by copy
## ##     reduced = gw$reduce(mod=FALSE)$nodes

## ##     expect_equal(length(which(duplicated(reduced))), 0)
    
## ##     for(i in length(reduced)) {
## ##         gr = reduced[i]
## ##         expect_equal(sum(gw.dt[start == start(gr) & end == end(gr) & strand == as.character(strand(gr))][,cn]), gr$cn)
## ##     }


## ##     ## Testing reduce by reference
## ##     gw.copy = gw$clone()
## ##     gw.copy$reduce()

## ##     expect_equal(length(which(duplicated(reduced))), 0)
    
## ##     for(i in 1:length(gw.copy$nodes)) {
## ##         gr = gw.copy$nodes[i]
## ##         expect_equal(sum(gw.dt[start == start(gr) & end == end(gr) & strand == as.character(strand(gr))][,cn]), gr$cn)
## ##     }

## ##     ## Checking paths
## ##     for(i in 1:length(reduced$path)) {
## ##         expect_true(max(reduced$path[[i]]) <= length(reduced$nodes))
## ##     }

## ##     for(i in 1:length(gw.copy$path)) {
## ##         expect_true(max(gw.copy$path[[i]]) <= length(gw.copy$nodes))
## ##     }
## ## })



## ## ## I think downloading data without warnings is evil

## ## ## trying URL 'ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_27/GRCh37_mapping/gencode.v27lift37.basic.annotation.gff3.gz'
## ## ## Content type 'unknown' length 39550148 bytes (37.7 MB)
## ## ## ==================================================
## ## ## trying URL 'ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_27/GRCh37_mapping/gencode.v27lift37.basic.annotation.gff3.gz'
## ## ## Content type 'unknown' length 39550148 bytes (37.7 MB)
## ## ## ==================================================



## ## ##-------------------------------------------------------##
## ## ##test_that('fusions', {
## ## ##    juncs = system.file('extdata', 'testing_junctions.rds', package="gGnome")
## ## ##    message("Junctions for testing: ", juncs)
## ## ##    juncs = readRDS(juncs)

## ##     ## make sure the gene annotation can be loaded
## ## ##    expect_error(cds <<- read_gencode(type = "cds"), NA)
## ## ##    expect_error(fusions())
## ## ##    expect_error(fusions(junc = juncs, cds = cds), NA) ## no problem
## ## ##})

## ## ## ##-------------------------------------------------------##
## ## ## test_that('graph distance and proximity', {
## ## ##     query = readRDS()
## ## ##     expect_error()
## ## ## })

## ## ##-------------------------------------------------------##
## ## test_that('bGraph walk and walk2', {
## ##     expect_true(inherits(jab.gg <<- gread(jab), "bGraph"))
## ##     expect_true(inherits(subg <<- jab.gg$hood(jab.gg$junctions[[1]], 5e5), "bGraph"))
## ##     expect_true(inherits(subg.gw2 <<- subg$walk2(T, F), "gWalks")) ## default to CPLEX
## ##     expect_true(inherits(jab.gw2.gurobi <<- subg$walk2(T, F, gurobi=TRUE), "gWalks"))
## ##     ## expect_true(inherits(subg.gw <<- subg$walk(gurobi=TRUE), "gWalks"))
    
## ## })

## ## ##-------------------------------------------------------##
## ## ## test_that('able to make JSON output', {
## ## ##     expect_true(inherits(jab.gw <<- as(readRDS(jab.gw.grl), "gWalks"), "gWalks"))
## ## ##     expect_equal(jab.gw$json("testing_gw.json"), "testing_gw.json")
## ## ## })

