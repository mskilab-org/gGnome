
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






test_that('junctions works', {

    expect_error(junctions(data.frame()))
    expect_error(junctions(GRanges()))
    expect_equal(length(junctions(GRangesList())), 0)
    expect_equal(length(junctions(grl1)), 250)

})




test_that('ra.duplicated works', {

    expect_false(all(ra.duplicated(junctions(grl1))))
    ##expect_equal(ra.duplicated(GRangesList()), logical(0))

})



## gGraph = R6::R6Class("gGraph"

## initialize = function(tile=NULL, junctions=NULL, cn = FALSE, jabba=NULL,
##     weaver=NULL, prego=NULL, segs=NULL, es=NULL, ploidy=NULL, purity=NULL, regular=TRUE)


## segstats, edges, junctions, G, adj, A, parts, seqinfo, purity, ploidy, td, win, ig


test_that('gGraph works', {

    ggnew = gGraph$new()
    expect_true(is(ggnew, 'gGraph'))
    foobar = gGraph$new(segs = test_segs, es=test_es)
    expect_true(is(foobar, 'gGraph'))
    foojab = gGraph$new(jabba = jab, segs = test_segs, es=test_es)
    expect_true(is(foojab, 'gGraph'))
    fooweaver = gGraph$new(weaver=weaver, segs = test_segs, es=test_es)
    expect_true(is(fooweaver, 'gGraph'))
    fooprego = gGraph$new(prego=prego, segs = test_segs, es=test_es)
    expect_true(is(fooprego, 'gGraph'))
    foocn = gGraph$new(segs = test_segs, es=test_es, cn=TRUE)
    expect_true(is(foocn, 'gGraph'))
    fooregular = gGraph$new(segs = test_segs, es=test_es)
    expect_true(is(fooregular, 'gGraph'))
    ##
    ##
    added_junctions = foojab$addJuncs(readRDS(jab)$junc)
    expect_true(is(added_junctions, 'gGraph')) 

})


## FUCNTIONS: initialize, set.seqinfo, nullGGraph, simpleGraph, dipGraphd, addJuncs, addSegs, karyograph, simplify,
##     decouple, add, jabb2gg, wv2gg, pr2gg, print, plot, window, layout, summary, gg2td, son, html, gg2js, components, 
##     subgraph, filling, 

## ACTIVE BINDINGS: segstats, edges, junctions, G, adj, A, parts, seqinfo, purity, ploidy, td, win, ig

test_that('gGraph, empty constructor/set.seqinfo/length', {

    ## Test with unspecified genome - seqinfo(default values)
    gg = gGraph$new()
    expect_true(is(gg, 'gGraph'))
    
    ## Check the active bindings
    expect_equal(gg$length(), 0)
    expect_equal(length(gg), 0)
    expect_equal(length(gg$nodes), 0)
    expect_equal(dim(gg$edges)[1], 0)
    expect_equal(length(gg$junctions), 0)
    expect_error(gg$graph, NA)
    expect_equal(length(gg$adj), 0)
    expect_equal(gg$parts, NULL)
    expect_equal(length(gg$seqinfo), 0)
    expect_equal(gg$td, NULL)
    expect_equal(length(gg$win), 0)
    
    ## Test with specified genome - invalid - just check seqinfo (set.seqinfo w/invalid genome)
    gg = gGraph$new(genome = 'fake genome')
    expect_equal(gg$length(), 0)
    expect_equal(length(gg$seqinfo), 0)
    
    ## Test with specified genome - valid - just check seqinfo (set.seqinfo w/valid genome)
    gg = gGraph$new(genome = hg_seqlengths())
    expect_equal(gg$length(), 0)
    expect_equal(length(gg$seqinfo), 25)

    ## set.seqinfo, genome != NULL, gname != NULL
    gg$set.seqinfo(genome = hg_seqlengths(), gname = 'foobar')
    expect_equal(gg$length(), 0)
    expect_equal(length(gg$seqinfo), 25)
    expect_equal(unique(gg$seqinfo@genome), 'foobar')

    ## Setting from another gGraph object
    gg = gGraph$new(genome = gg)
    expect_equal(gg$length(), 0)
    expect_equal(length(gg$seqinfo), 25)
})


## FIXME: Definitely could add more tests to catch errors and things in constructor but it seems to be working
test_that('gGraph, Nodes and Edges Constructor/active bindings/looseNodes', {
    ## Make sure it add the right nodes and edges
    ## 1) edges = NULL
    ##     a) looseterm = FALSE
    ##     b) looseterm = TRUE
    ## 2) edges != NULL
    ##     a) looseterm = FALSE
    ##     b) looseterm = TRUE
    ## - check edges, nodes, seqinfo, looseNodes

    ## CASE 1a
    nodes = c(GRanges("1",IRanges(1,100),"*"), GRanges("1",IRanges(101,200),"*"),
              GRanges("1",IRanges(201,300),"*"), GRanges("1",IRanges(301,400),"*"),
              GRanges("1",IRanges(401,500),"*"))

    ## Set the seqinfo to make sure it carries over
    nodes = gUtils::gr.fix(nodes, hg_seqlengths())
    seq = hg_seqlengths()
    names(seq) = NULL

    ## Check the null case
    gg = gGraph$new(nodes = GRanges(), looseterm = FALSE)

    expect_equal(gg$nodes, GRanges())
    expect_equal(dim(gg$edges)[1], 0)
    expect_error(gg$graph, NA)
    expect_equal(length(seqinfo(gg)), 0)
    
    ## Check with looseterm = F
    gg = gGraph$new(nodes = nodes, looseterm = FALSE)

    expect_equal(sort(granges(gg$nodes)), sort(granges(nodes)))
    expect_equal(dim(gg$edges)[1], 0)
    expect_equal(length(gg$looseNodes()), 0)
    
    ## CASE 1b
    gg = gGraph$new(nodes = nodes, looseterm = TRUE)

    loosenodes = c(GRanges("1",IRanges(1,1),"*"), GRanges("1",IRanges(100,100),"*"),
                   GRanges("1",IRanges(101,101),"*"), GRanges("1",IRanges(200,200),"*"),
                   GRanges("1",IRanges(201,201),"*"), GRanges("1",IRanges(300,300),"*"),
                   GRanges("1",IRanges(301,301),"*"), GRanges("1",IRanges(400,400),"*"),
                   GRanges("1",IRanges(401,401),"*"), GRanges("1",IRanges(500,500),"*"))
    loosenodes = gUtils::gr.fix(loosenodes, hg_seqlengths())
    
    expect_equal(sort(granges(nodes)), sort(granges(gg$nodes)))
    expect_equal(dim(gg$edges)[1], 0)
    expect_equal(length(gg$looseNodes()), 10)
    expect_equal(sort(granges(gg$looseNodes())), sort(granges(loosenodes)))
    
    ## CASE 2a
    edges = data.table(n1 = c(3,2,4,1,5), n2 = c(3,4,2,5,1), n1.side = c(1,1,0,0,1), n2.side = c(0,0,0,1,0))

    gg = gGraph$new(nodes = nodes, edges = edges, looseterm = FALSE)

    expect_true(is(gg, 'gGraph'))
    
    ## Check that the correct things are stored in our gGraph
    expect_equal(sort(granges(nodes)), sort(granges(gg$nodes)))
    expect_equal(edges, gg$edges)
    expect_equal(length(gg), 5)
    expect_equal(length(gg$looseNodes()), 0)
    expect_equal(seqinfo(gg)@seqlengths, seq)
    
    ## CASE 2b
    gg = gGraph$new(nodes = nodes, edges = edges, looseterm = TRUE)

    loosenodes = c(GRanges("1", IRanges(100,100), "*"),
                   GRanges("1", IRanges(400,400), "*"), GRanges("1", IRanges(401,401), "*"))
    loosenodes = gUtils::gr.fix(loosenodes, hg_seqlengths())
    
    ## Check that the correct things are stored in our gGraph
    expect_equal(sort(granges(nodes)), sort(granges(gg$nodes)))
    expect_equal(edges, gg$edges)
    expect_equal(length(gg), 5)
    expect_equal(sort(granges(gg$looseNodes())), sort(granges(loosenodes)))
})


## FIXME: test the circular thing idk how to do that
test_that('gGraph, simpleGraph', {
    ggnew = gGraph$new()
    
    ## simpleGraph = function(genome = NULL, chr=FALSE, include.junk=FALSE)
    ## Testing simpleGraph(default values)
    gg = ggnew$simpleGraph()
    expect_true(is(gg, 'gGraph'))
    expect_equal(length(gg$nodes), 25)
    expect_equal(dim(gg$edges)[1], 0)
    expect_equal(length(gg$junctions), 0)
    expect_error(gg$graph, NA) ## check it works
    ##expect_equal(length(gg$adj), 2500) -- FIXME: Don't know what this does but it doesn't work
    ## expect_equal(gg$parts, NULL) -- FIXME: also broken, all don't know what it does
    expect_equal(length(gg$seqinfo), 25)
    expect_equal(length(gg$looseNodes()), 50)

    ## Don't know what the fuck is going on here
    expect_true(is(gg$td, 'gTrack'))
    expect_equal((gg$td)$ygap, 2)
    expect_match((gg$td)$name, 'CN')
    expect_equal(length(gg$win), 25)

    ## Testing simpleGraph with genome = Something
    gg = ggnew$simpleGraph(genome = hg_seqlengths())
    expect_equal(length(gg$nodes), 25)
    expect_equal(length(gg$looseNodes()), 50)
})


## TESTING TRIM
test_that('gGraph, trim', {
    ## 1) trim within single node
    ##   a) has subgraph
    ##   b) doesn't have subgraph
    ## 2) trim across multiple nodes
    ##   a) all have subgraphs
    ##   b) none have subgraphs
    ## 3) Trim is not continuous
    ## Want to reduplicate all of these after with mod = T and make sure the graphs are exactly the same
    ## FIXME: there are some random things like what if its not balanced and there is a copy number or it is a bg?
    ##           -- This just isn't good practice should be in bg class
    
    ## Build a gGraph
    gr = GRanges("1", IRanges(1001,2000), "+")
    gr = c(gr, GRanges("1", IRanges(2001,3000), "+"),
           GRanges("1", IRanges(3001,4000), "+"))
    gr = c(gr, gr.flipstrand(gr))
    es = data.table(to = c(2,3,4,5), from = c(1,2,5,6))
    
    graph = gGraph$new(segs = gr, es = es)

    ## Make sure trim was successful - trim on edges
    expect_equal(streduce(graph$nodes), streduce(gr))
     
    ## CASE 1b
    gr1 = GRanges("1", IRanges(1200,1500), "+")
    tmp = graph$trim(gr1)

    ## Make sure trim worked
    expect_equal(streduce(tmp$nodes), streduce(gr1))
    
    ## CASE 2b
    gr2 = GRanges("1", IRanges(1800,3500), "+")
    tmp2 = graph$trim(gr2)

    expect_equal(streduce(tmp2$nodes), streduce(gr2))
   
    ## Case 3b
    tmp3 = graph$trim(c(gr1,gr2))

    expect_equal(streduce(tmp3$nodes), streduce(c(gr1,gr2)))

    ## Subgraphs
    ##addGraphs = list(gGraph$new(tile = gr[1])$trim(gr[1]),
    ##                 gGraph$new(tile = gr[2])$trim(gr[2]),
    ##                 gGraph$new(tile = gr[3])$trim(gr[3]))
    ##graph$resetSubgraphs(subs = addGraphs)

    ## CASE 1a - We know trim works from above so check subgraphs
    ##tmp = graph$trim(gr1)

    ##expect_equal(tmp$subgraphs[[1]]$segstats, addGraphs[[1]]$trim(gr1)$segstats)

    ## CASE 2a
    ##tmp2 = graph$trim(gr2)

    ##for(i in 1:length(tmp2$subgraphs)) {
    ##    expect_equal(tmp2$subgraphs[[i]]$segstats, addGraphs[[i]]$trim(gr2)$segstats)
    ##}
    
    ## Case 3a
    ##gr3 = GRanges("1", IRanges(1800,2500), "+")
    ##tmp3 = graph$trim(c(gr1, gr3))

    ##expect_equal(tmp3$subgraphs[[1]]$segstats, addGraphs[[1]]$trim(gr1)$segstats)
    ##expect_equal(tmp3$subgraphs[[2]]$segstats, addGraphs[[1]]$trim(gr3)$segstats)
    ##expect_equal(tmp3$subgraphs[[3]]$segstats, addGraphs[[2]]$trim(gr3)$segstats)

    ##FIXME: Make some tests here about the mod thing
    
    ##FIXME: not checking edges at all but constructor should handle it
})


## Test the addSegs/makeSegs functionality
test_that('gGraph, addSegs', {
    ## TESTING NO INPUTTED SUBGRAPHS
    gr = GRanges("1", IRanges(1001,2000), "+")
    gr = c(gr, GRanges("1", IRanges(2001,3000), "+"),
           GRanges("1", IRanges(3001,4000), "+"))

    graph = gGraph$new(tile = gr)$trim(gr)

    ## Expect number of nodes to be just changes
    expect_equal(graph$length(), 3)

    ## Expect subgraphs to be half length of all nodes
    ##expect_equal(length(graph$subgraphs), graph$length())

    ## Expect there to only two of each subgraph
    ##subs_count = table(graph$nodes$subIndex)
    ##expect_equal(max(subs_count), 2)
    ##expect_equal(min(subs_count), 2)

    ## Expect all subgraphs to be length 0
    ##for (i in 1:length(graph$subgraphs)) {
    ##    expect_equal(graph$subgraphs[[1]]$length(), 0)
    ##}

    ## TESTING INPUTTED SUBGRAPHS
    ##subgraphs = list(gGraph$new(tile = gr[1])$trim(gr[1]),
    ##                 gGraph$new(tile = gr[2])$trim(gr[2]))
    
    ##graph = gGraph$new(tile = gr, subs = subgraphs)$trim(gr[1:3])

    ## Expect number of nodes to be changes
    ##expect_equal(graph$length(), 3)

    ## Expect subgraphs to be length of nodes (half total)
    ##expect_equal(length(graph$subgraphs), graph$length())

    ## Expect there to only two of each subgraph
    ##subs_count = table(graph$segstats$subIndex)
    ##expect_equal(max(subs_count), 2)
    ##expect_equal(min(subs_count), 2)

    ## Expcet that the two graphs we added are in the right spots
    ##expect_equal(graph$subgraphs[[ graph$segstats[2]$subIndex ]], subgraphs[[1]])
    ##expect_equal(graph$subgraphs[[ graph$segstats[3]$subIndex ]], subgraphs[[2]])
    
    ## TESTING WITH A SUBGRAPH BEING SPLIT FOR A NEW NDOE
    gr1 = GRanges("1", IRanges(1500,1700), "+")
    graph$addSegs(gr1)

    ## Expect the length to be 2 longer
    expect_equal(graph$length(), 5)

    ## Expect subgraphs to be half length of all nodes
    ##expect_equal(length(graph$subgraphs), graph$length())

    ## Expect there to only two of each subgraph
    ##subs_count = table(graph$nodes$subIndex)
    ##expect_equal(max(subs_count), 2)
    ##expect_equal(min(subs_count), 2)

    gr1 = gr.breaks(gr1, gr[1])
    
    ## Expect the nodes we split to have subgraphs on that range
    ##expect_equal(streduce(graph$subgraphs[[ graph$segstats[2]$subIndex ]]$segstats),
    ##             streduce(gr1[1]))
    ##expect_equal(streduce(graph$subgraphs[[ graph$segstats[3]$subIndex ]]$segstats),
    ##             streduce(gr1[2]))
    ##expect_equal(streduce(graph$subgraphs[[ graph$segstats[4]$subIndex ]]$segstats),
    ##             streduce(gr1[3]))
})


test_that('gGraph, mergeGraphs', {
    ## CASES
    ## 1) Self is null, adding null
    ## 2) Self is null, adding non-null ---- Could be caught with a clause a the beginning to save work
    ## 3) Self is non-null, adding null
    ## 4) Self is non-null, adding non-null
    ##    a) Neither have subgraphs
    ##    c) Both have subgraphs
    ##    b) self or added has subgraphs
    ## Not handling decoupling at all - this should just work though, I shouldn't need to test it as long as I fix it first

    
    gr = GRanges("1", IRanges(1001,2000), "+")
    gr = c(gr, GRanges("1", IRanges(2001,3000), "+"),GRanges("1", IRanges(3001,4000), "+"))

    graph = gGraph$new(tile = gr)$trim(gr)

    ## Case 1
    nullgg = gGraph$new()
    tmp = nullgg$mergeGraphs(nullgg, decouple=F)

    expect_equal(nullgg, tmp)

    ## Case 2
    tmp = nullgg$mergeGraphs(graph, decouple=F)

    expect_equal(graph, tmp)
    
    ## Case 3
    tmp = graph$mergeGraphs(nullgg, decouple=F)

    expect_equal(graph, tmp)

    ## Case 4a
    gr1 = GRanges("1", IRanges(4001,5000), "+")
    gr1 = c(gr1, GRanges("1", IRanges(5001,6000), "+"))

    graph1 = gGraph$new(tile = gr1)$trim(gr1)
    tmp = graph$mergeGraphs(graph1, decouple=F)

    ## FIXME: Check segs and es

    ## Check that all the subgraphs are null
    ##for(i in 1:length(tmp$subgraphs)) {
    ##    expect_equal(tmp$subgraphs[[i]]$length(), 0)
    ##}
    
    ## Case 4b - We know the main level works, just check subgraphs
    ##subgraphs = list(gGraph$new(tile = gr1[1])$trim(gr1[1]),
    ##                 gGraph$new(tile = gr1[2])$trim(gr1[2]))

    ##graph1$resetSubgraphs(subgraphs)
    tmp = graph$mergeGraphs(graph1, decouple=F)
    
    ## FIXME: insert subgraph test
    
    ## Case 4c
    ##subgraphs = list(gGraph$new(tile = gr[1])$trim(gr[1]),
    ##                 gGraph$new(tile = gr[2])$trim(gr[2]))

    ##graph$resetSubgraphs(subgraphs)
    tmp = graph$mergeGraphs(graph1, decouple=F)

    ## FIXME: insert subgraph test
    
})


test_that('gGraph, decouple', {
    
})


### some non-exported functions

rev.comp = function(gr){
    strmap = setNames(c("+", "-"), c("-", "+"))
    if (!inherits(gr, "GRanges")){
        stop("Error: Input must be GRanges.")
    } else if (!all(strand(gr) %in% strmap)) {
        stop("Error: Input must be all strand specific.")
    }
    return(rev(gr.flipstrand(gr)))
}


##test_that('rev.comp works', {
##
##    expect_error(rev.comp())
##    expect_error(rev.comp(data.frame()))
##    gr3 = dt2gr(dt)
##    expect_error(rev.comp(gr3))   ## Error in rev.comp(gr3) : Input must be all strand specific.
##    expect_equal(width(rev.comp(gr)[1]), 4)
##    expect_equal(as.character(strand(rev.comp(gr)[1])), "+")
##    expect_equal(width(rev.comp(gr)[2]), 3)
##    expect_equal(as.character(strand(rev.comp(gr)[2])), "+")
##
##})


capitalize = function(string, un = FALSE){
    if (!un){
        capped <- grep("^[^A-Z].*$", string, perl = TRUE)
        substr(string[capped], 1, 1) <- toupper(substr(string[capped],1, 1))
    } else{
        capped <- grep("^[A-Z].*$", string, perl = TRUE)
        substr(string[capped], 1, 1) <- tolower(substr(string[capped],1, 1))
    }
    return(string)
}





test_that('capitalize works', {

    str1 = "Foo FOO"
    str2 = "2Foo . $%@"
    str3 = "foobar foo"
    expect_match(capitalize(str1), "Foo FOO")
    expect_match(capitalize(str1, un=TRUE), "foo FOO")
    expect_equal(capitalize(str2), "2Foo . $%@")
    expect_equal(capitalize(str2, un=TRUE), "2Foo . $%@")    ## probably should report this bug to r-lib
    expect_match(capitalize(str3), "Foobar foo")
    expect_match(capitalize(str3, un=TRUE), "foobar foo")

})


  

ul = function(x, n=6){
    n = pmin(pmin(dim(x)), n)
    return(x[1:n, 1:n])
}



test_that('ul works', {

    A = matrix(  c(2, 4, 3, 1, 5, 7),  nrow=2, ncol=3, byrow = TRUE)    
    expect_equal(as.integer(ul(A, n=0)), 2)   ### Is this expected behavior? 
    expect_equal(as.integer(ul(A, n=1)), 2)
    expect_equal(dim(ul(A, n=2))[1], 2)
    expect_equal(dim(ul(A, n=2))[2], 2)
    expect_equal(dim(ul(A, n=999))[1], 2)
    expect_equal(dim(ul(A, n=9999))[2], 2)   ### Is this expected behavior? 

})




tile.name = function(x){
    if (!inherits(x, "GRanges")){
        stop("Only takes GRanges as input for now.")
    }
    hb = hydrogenBonds(segs = x)
    if (hb[, any(is.na(from) | is.na(to))]){
        stop("Not fully strand paired.")
    }
    hb.map = hb[, c(setNames(from, to), setNames(to, from))]
    seg.name = ifelse(strand(x)=="+",
                      as.character(seq_along(x)),
                      paste0("-", hb.map[as.character(seq_along(x))]))
    return(seg.name)
}


test_that('tile.name works', {

    ## if (!inherits(x, "GRanges")){
    expect_error(tile.name(GRangesList()))
    expect_equal(length(tile.name(test_segs)), 10)

})



test_that('e2j works', {

    expect_equal(length(e2j(test_segs, test_es)), 2)

})




test_that('etype works', {

    expect_equal(dim(etype(test_segs, test_es))[1], 12)
    expect_equal(dim(etype(test_segs, test_es))[2], 16) 

})






write.tab = function (x, ..., sep = "\t", quote = F, row.names = F)
{
    if (!is.data.frame(x)){
        x = as.data.frame(x)
    }
    write.table(x, ..., sep = sep, quote = quote, row.names = row.names)
}



test_that('write.tab() works', {

    expect_equal(length(write.tab(dt)), 0)  ### throws dt to STDOUT

})







test_that('get.ploidy works', {

    expect_error(get.ploidy(GRangesList()))
    expect_equal(get.ploidy(test_segs), 3)
    ## if (length(cnix <- grep("CN", colnames(mcols(segs)), ignore.case = T)) == 

})





test_that('dedup() works', {

    expect_equal(dedup(c(rep(2, 10.5), rep(3, 20)))[30], "3.20")

})










## read_vcf()
## read_vcf = function(fn, gr = NULL, hg = 'hg19', geno = NULL, swap.header = NULL, verbose = FALSE, add.path = FALSE, tmp.dir = '~/temp/.tmpvcf', ...)
##test_that('read_vcf', {
#    ## error
#    expect_error(read_vcf('foobar'))
#    ## default 
#    expect_equal(length(read_vcf(somatic_vcf)), 60)
#    expect_equal(length(seqnames(seqinfo(read_vcf(somatic_vcf)))), 84)
#    ## gr  gr= GRanges('1:10075-10100')
#    ## hg
##    expect_match(unique(as.data.frame(seqinfo(read_vcf(somatic_vcf, hg='hg12345')))$genome), 'hg12345')
#    ## geno
#    ## swap.header
#    expect_equal(length(seqnames(seqinfo(read_vcf(somatic_vcf, swap.header='/Users/ebiederstedt/bamUtils/tests/testthat/new_header.vcf')))), 2)
#    ## verbose
#    expect_equal(length(read_vcf(somatic_vcf, verbose=TRUE)), 60)
#    ## check 'if (!file.exists(swap.header))'
#    expect_error(read_vcf(somatic_vcf, swap.header='foobar'))
#
#})



test_that('chr2num works', {

    expect_equal(as.logical(chr2num("ChrX")), NA)
    expect_equal(chr2num("chrX"), 23)
    expect_equal(chr2num("chrY"), 24)

})



test_that('affine.map works', {

    expect_equal(affine.map(49), 0.5)

})


test_that('gr.flatmap works', {

    expect_equal((gr.flatmap(example_genes, windows=GRanges('1:10000-20000'))$window.segs)$start, 1)
    expect_equal((gr.flatmap(example_genes, windows=GRanges('1:10000-20000'))$window.segs)$end, 10001)
    expect_equal(length((gr.flatmap(example_genes, windows=GRanges('1:10000-20000'))$window.segs)$grl.segs), 0)

})




### XT's tests

##-------------------------------------------------------##
test_that('constructors and essential functions', {
    ## small example, nested tDUP
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
    ## expect_equal(dim(foo$edges)[2], 16)
    expect_equal(max((foo$edges)$cn), 3)
    expect_equal(max((foo$edges)$fromStart), 18593415)
    expect_equal(max((foo$edges)$fromEnd), 18793414)
    foo = gGraph$new(segs=test_segs, es=test_es)
    expect_equal(dim(foo$edges)[1], 12)
    ## expect_equal(dim(foo$edges)[2], 16)
    expect_equal(max((foo$edges)$cn), 3)
    expect_equal(max((foo$edges)$fromStart), 18593415)
    expect_equal(max((foo$edges)$fromEnd), 18793414)
    expect_equal(dim((foo$nullGGraph())$edges)[1], 0)
    expect_equal(dim((foo$nullGGraph())$edges)[2], 3)

})

test_that('karyograph', {

    kag.tile = gGraph$new(tile = test_segs)
    expect_true(inherits(kag.tile, "gGraph"))

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
    jab_bgraph = gread(jab)
    expect_true(is(jab_bgraph, "bGraph"))
    ## preg_bgraph = gread(prego)
    ## expect_true(is(preg_bgraph, "bGraph")) ### 'gGraph'
    ## wv_bgraph = gread(weaver)   
    ## expect_true(is(wv_bgraph, "bGraph"))  ### 'gGraph'
    ## if (is.list(file)){
    list_foo = gread(readRDS(system.file("extdata", "jabba.simple.rds", package="gGnome")))
    expect_true(is(list_foo, 'bGraph'))
    
})





## ##-------------------------------------------------------##
## test_that('gtf2json', {
##     expect_error(gread('no_file_here'))
##     expect_equal(gtf2json(system.file('extdata', 'test.gtf', package='gGnome')), "./gtf.json")
##     system(paste('rm', "./gtf.json"))
## })




## * could not find function "setxor"
##-------------------------------------------------------##

test_that('setxor', {

    A = c(1, 2, 3)
    B = c(1, 4, 5)
    expect_equal(setxor(A, B), c(2, 3, 4, 5))

})





##-------------------------------------------------------##
test_that('special ranges functions for skew-symmetric graph', {


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


## Error: Test failed: 'gWalks'
## * length(gw <<- as(grl, "gWalks")) not equal to sum(values(grl)$cn > 0).
## 1/1 mismatches
## [1] 32 - 630 == -598
## * `bg <<- as(gw, "bGraph")` threw an error.
## Message: Error: Given edge data is not skew-symmetric!!!
## Class:   simpleError/error/condition
## * object 'bg' not found
## 1: expect_equal(length(bg$junctions), sum(values(junctions)$cn > 0)) at :12
## 2: quasi_label(enquo(object), label)
## 3: eval_bare(get_expr(quo), get_env(quo))
 
 
test_that('gWalks reduce', {
    ## Set up gWalks
    grl = readRDS(jab.gw.grl)
    gw = as(grl, "gWalks")
    gw.dt = gr2dt(gw$nodes)

    ## Testing reduce by copy
    reduced = gw$reduce(mod=FALSE)$nodes

    expect_equal(length(which(duplicated(reduced))), 0)
    
    for(i in length(reduced)) {
        gr = reduced[i]
        expect_equal(sum(gw.dt[start == start(gr) & end == end(gr) & strand == as.character(strand(gr))][,cn]), gr$cn)
    }


    ## Testing reduce by reference
    gw.copy = gw$clone()
    gw.copy$reduce()

    expect_equal(length(which(duplicated(reduced))), 0)
    
    for(i in 1:length(gw.copy$nodes)) {
        gr = gw.copy$nodes[i]
        expect_equal(sum(gw.dt[start == start(gr) & end == end(gr) & strand == as.character(strand(gr))][,cn]), gr$cn)
    }

    ## Checking paths
    for(i in 1:length(reduced$path)) {
        expect_true(max(reduced$path[[i]]) <= length(reduced$nodes))
    }

    for(i in 1:length(gw.copy$path)) {
        expect_true(max(gw.copy$path[[i]]) <= length(gw.copy$nodes))
    }
})



## I think downloading data without warnings is evil

## trying URL 'ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_27/GRCh37_mapping/gencode.v27lift37.basic.annotation.gff3.gz'
## Content type 'unknown' length 39550148 bytes (37.7 MB)
## ==================================================
## trying URL 'ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_27/GRCh37_mapping/gencode.v27lift37.basic.annotation.gff3.gz'
## Content type 'unknown' length 39550148 bytes (37.7 MB)
## ==================================================



##-------------------------------------------------------##
##test_that('fusions', {
##    juncs = system.file('extdata', 'testing_junctions.rds', package="gGnome")
##    message("Junctions for testing: ", juncs)
##    juncs = readRDS(juncs)

    ## make sure the gene annotation can be loaded
##    expect_error(cds <<- read_gencode(type = "cds"), NA)
##    expect_error(fusions())
##    expect_error(fusions(junc = juncs, cds = cds), NA) ## no problem
##})

## ##-------------------------------------------------------##
## test_that('graph distance and proximity', {
##     query = readRDS()
##     expect_error()
## })

##-------------------------------------------------------##
test_that('bGraph walk and walk2', {
    expect_true(inherits(jab.gg <<- gread(jab), "bGraph"))
    expect_true(inherits(subg <<- jab.gg$hood(jab.gg$junctions[[1]], 5e5), "bGraph"))
    expect_true(inherits(subg.gw2 <<- subg$walk2(T, F), "gWalks")) ## default to CPLEX
    expect_true(inherits(jab.gw2.gurobi <<- subg$walk2(T, F, gurobi=TRUE), "gWalks"))
    ## expect_true(inherits(subg.gw <<- subg$walk(gurobi=TRUE), "gWalks"))
    
})

##-------------------------------------------------------##
## test_that('able to make JSON output', {
##     expect_true(inherits(jab.gw <<- as(readRDS(jab.gw.grl), "gWalks"), "gWalks"))
##     expect_equal(jab.gw$json("testing_gw.json"), "testing_gw.json")
## })
