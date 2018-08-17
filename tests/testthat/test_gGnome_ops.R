library(gGnome)
library(testthat)
library(gUtils)

svaba = system.file('extdata', "HCC1143.svaba.somatic.sv.vcf", package = "gGnome")
delly = system.file('extdata', "delly.final.vcf.gz", package = "gGnome")
novobreak = system.file('extdata', "novoBreak.pass.flt.vcf", package = "gGnome")

HGSL = c("1"=249250621, "2"=243199373, "3"=198022430, "4"=191154276, "5"=180915260, "6"=171115067, "7"=159138663, "X"=155270560, "8"=146364022, "9"=141213431, "10"=135534747, "11"=135006516, "12"=133851895, "13"=115169878, "14"=107349540, "15"=102531392, "16"=90354753, "17"=81195210, "18"=78077248, "20"=63025520, "Y"=59373566, "19"=59128983, "22"=51304566, "21"=48129895, "M"=16571)


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

remixt = system.file('extdata', 'remixt', package='gGnome')
message("remixt results: ", remixt)

genome = seqinfo(test_segs)

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
   expect_is(gGraph$new(remixt = remixt), "gGraph")
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
     gn$mark(col="purple")
     
     expect_is(gn, "gNode")
     expect_equal(length(gn), length(gg))
     expect_equal(gn$gr, gg$gr %Q% (strand == "+"))
     expect_equal(gn$id, (gg$gr %Q% (strand == "+"))$snode.id)
     expect_equal(gn$id, gn$sid)

     ##gNode ego
     gn=gNode$new(3, gg)
     expect_equal(gg[c(3:4)]$nodes$dt[, start], gn$ego(1)$dt[, start])
     expect_equal(gg[c(3, 4, 2)]$nodes$dt[, start], gn$ego(2)$dt[, start])

     ##loose.left, loose right
     gn=gNode$new(1, gg)
     expect_equal(gn$loose.left, FALSE)
     expect_equal(gn$loose.right, TRUE)
     expect_equal(gn$degree, 4)

     ##terminal, degrees
     expect_equal(gr2dt(gn$terminal)[, start], 100)
     expect_equal(gn$ldegree, 1)
     expect_equal(gn$rdegree, 0)     

     ##unioning gNodes
 ##    gn2=gNode$new(2, gg)    
 
    ## expect_equal(union(gn, gn2)$dt, gg$nodes[c(1:2)]$dt)
    ## gg=gGraph$new(nodes=nodes1, edges=edges)
    ## expect_error(union(gg$nodes, gn2))
     
})

test_that('gNode subsetting', {
    nodes1 = c(GRanges("1",IRanges(1,100),"*", cn=1), GRanges("1",IRanges(101,200),"*",cn=1),
               GRanges("1",IRanges(201,300),"*", cn=1), GRanges("1",IRanges(301,400),"*", cn=1),
               GRanges("1",IRanges(401,500),"*",cn=1))
    edges = data.table(n1 = c(3,2,4,1,3), n2 = c(3,4,2,5,4), n1.side = c(1,1,0,0,1), n2.side = c(0,0,0,1,0))        
    gg = gGraph$new(nodes = nodes1, edges = edges)      
    gn = gg$nodes    
    gn2= gNode$new(2, gg)
    gn3= gNode$new(1, gg)   
    
    ##Right and Left
    expect_equal(gn2$right$dt[, start], 301)
    expect_equal(gn2$left$dt[, start], 301)

    ##Quick subgraph test
    sub=gn$subgraph
    expect_equal(length(sub), length(gg))   
    ##Various overriden functions
    node1=gn[1]
    nodes2=gn[c(1:4)]
    expect_equal(length(c(gn2, gn2)), 2)
    expect_equal(length(setdiff(gn, gn)), 0)
    expect_identical(intersect(gn3, gg$nodes)$dt, gg$nodes[1]$dt)
    
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
##    expect_identical(gn[loose.left == FALSE]$dt, gn[c(1:4)]$dt) ## causing problems on travis only
#     expect_equal(gn[snode.id > 3], gn[4:5]) #### this syntax is not playing well with Travis
##    expect_identical(gn[start > 150 & end < 450 & loose.left == FALSE]$dt, gn[3:4]$dt) #### this syntax is not playing well with Travis
})

test_that('gEdge works',{
    nodes1 = c(GRanges("1",IRanges(1,100),"*"), GRanges("1",IRanges(101,200),"*"),
               GRanges("1",IRanges(201,300),"*"), GRanges("1",IRanges(301,400),"*"),
               GRanges("1",IRanges(401,500),"*"))
    edges = data.table(n1 = c(3,2,4,1,3), n2 = c(3,4,2,5,4), n1.side = c(1,1,0,0,1), n2.side = c(0,0,0,1,0))     
    gg = gGraph$new(nodes = nodes1, edges = edges)    
    ge=gEdge$new(1, gg)
    ge2=gEdge$new(2, gg)
    ge3=c(ge, ge2)
    
    ##Mark
    ge$mark(changed=TRUE)
    expect_equal(gg$edges[1]$dt[, changed], TRUE)
    
    ##Basic active fields
    expect_equal(ge$length, 1)     
    expect_equal(length(ge), 1)
    expect_equal(ge$left, gg$nodes[3])    
    expect_equal(ge$right, gg$nodes[3])
    
    ##Overriden functions   
    expect_equal(c(ge, ge2)$dt[, sedge.id], gg$edges[c(1, 2)]$dt[, sedge.id])    
    expect_equal(length(setdiff(ge, ge)), 0)
    expect_equal(intersect(c(ge, ge3), ge)$dt[, sedge.id], 1)          
    
    ##Miscellaneous functions    
    starts=gr2dt(ge$junctions$grl)[, start]
    expect_equal(gr2dt(ge$junctions$grl)[, start][1], 300)
    expect_equal(gr2dt(ge$junctions$grl)[, end][2], 201)
    expect_equal(class(ge$subgraph)[1], "gGraph")
    
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
    expect_equal(length(setdiff(jj, jj)), 0)
    expect_equal(length(intersect(jj, jj)),length(jj))
    
    expect_equal(length(jj), 500)
    expect_equal(unlist(jj$grl), unlist(juncs))
   
   ## expect_equal(jj$dt, as.data.table(values(juncs)))
    
    ## c() / +
    jj2 = c(jj, jj)
    expect_equal(length(jj2), length(jj)*2)   
    expect_equal(as.data.table(jj2$grl), as.data.table(c(jj$grl, jj$grl)))
  juncs.delly = Junction$new(delly)
  juncs.novobreak = Junction$new(novobreak)
  juncs.svaba = Junction$new(svaba)

  expect_true(length(juncs.delly)==210)
  expect_true(length(juncs.novobreak)==421)
  expect_true(length(juncs.svaba)==500)
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
  expect_equal(length(setdiff(jj, jj)), 0)
  expect_equal(length(intersect(jj, jj)),length(jj))
  
    expect_equal(length(jj), 500)
    expect_equal(unlist(jj$grl), unlist(juncs))   
    expect_equal(jj$dt, as.data.table(juncs))
  
  ## c() / +
  jj2 = c(jj, jj)
  expect_equal(length(jj2), length(jj)*2)
  expect_equal(as.matrix(as.data.table(jj2$grl)), as.matrix(as.data.table(c(jj$grl, jj$grl))))

  jj2 = jj + 100
  expect_true(all(all(width(jj2$grl)==101)))
  
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
    
    ##subset
    expect_equal(gr2dt(jj[1]$grl)[,group][1], 1)        

})

test_that('gGraph, empty constructor/length', {
  ## Test with unspecified genome - seqinfo(default values)
  gg = gGraph$new()
  expect_true(is(gg, 'gGraph'))
  
  ## Check the active bindings
  expect_equal(length(seqinfo(gg)), 0)
  
  ## Test with specified genome - valid - just check seqinfo (set.seqinfo w/valid genome)
  ##gg = gGraph$new(genome = HGSL)
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
  gr = gr.fix(gr, HGSL)
  
  es = data.table(n1 = c(1,2,3,4,4), n2 = c(2,3,4,1,3), n1.side = c(1,1,1,0,1), n2.side = c(0,0,0,1,0))
  
  graph = gGraph$new(nodes = gr, edges = es)

  ## CASE 1
  gr1 = GRanges("1", IRanges(1200,1500), "+")
  gr1 = gr.fix(gr1, HGSL)   
  tmp = graph$copy$trim(gr1)   
  expect_identical(as.data.table(streduce(tmp$nodes$gr)), as.data.table(streduce(gr1)))
  ##keep as is  expect_equal(granges(tmp$nodes$gr), granges(gr1))
  ##keep as is    expect_equal(length(gg$edges), 0)
  expect_equal(length(tmp), 1)
  
  ## CASE 2
  gr2 = GRanges("1", IRanges(2200,4500), "+")
  gr2 = gr.fix(gr2, HGSL)
  edges = data.table(n1 = c(1,2), n2 = c(2,3), n1.side = c(1,1), n2.side = c(0,0))
  
  tmp2 = graph$copy$trim(gr2)
  expect_identical(as.data.table(streduce(tmp2$nodes$gr)), as.data.table(streduce(gr2)))  
  es = tmp2$edges$dt[order(n1,n2,n1.side,n2.side),][, c("type") := NULL]    
  ## expect_equal(es, edges[order(n1,n2,n1.side,n2.side),])

   expect_equal(length(tmp2), 3)
   
    ## Case 3
    ranges(gr2) = IRanges(1800,2200)
    gr2 = c(gr2, GRanges("1", IRanges(2800,4500), "+"))
   edges = data.table(n1 = c(2,4,5,6), n2 = c(3,5,6,2), n1.side = c(1,1,1,0), n2.side = c(0,0,0,1))      
    graph$trim(c(gr1,gr2), mod=T)
    
    expect_equal(streduce(graph$nodes$gr), streduce(c(gr1,gr2)))
    
    es = graph$edges$dt[order(n1,n2,n1.side,n2.side),][, c("type") := NULL]
    ##  expect_equal(es, edges[order(n1,n2,n1.side,n2.side),])   
    expect_equal(length(graph), 6)

    
})

test_that('some public gGraph fields',{
     nodes1 = c(GRanges("1",IRanges(1,100),"*"), GRanges("1",IRanges(101,200),"*"),
                GRanges("1",IRanges(201,300),"*"), GRanges("1",IRanges(301,400),"*"),
                GRanges("1",IRanges(401,500),"*"))
     edges = data.table(n1 = c(3,2,4,1,3), n2 = c(3,4,2,5,4), n1.side = c(1,1,0,0,1), n2.side = c(0,0,0,1,0))
     
     ##gTrack
     gg = gGraph$new(nodes = nodes1, edges = edges)        
     expect_is(gg$gtrack(), "gTrack")
     expect_is(gg$gt, "gTrack")
     expect_equal(gg$gtrack()$ygap, 2)
     expect_equal(gg$gtrack()$name, "gGraph")

     #gwalks
     expect_is(gg$walks(), "gWalk")
     expect_equal(gg$walks()$dt[, wid], c(200, 300, 100))
 
     ##subgraph
     gr1=gg$nodes$gr[1]
     sub=gg$subgraph(gr1, 300)
     expect_equal(sub$dt[, start], gg[1:3]$dt[, start])    
     expect_equal(sub$dt[2, loose.right], TRUE)
     expect_equal(sub$dt[2, loose.left], TRUE)
     
     ##addJuncs
     graph=copy(gg)
     graph$addJuncs(graph$junctions)
     starts=data.table(r=duplicated(graph$junctions$dt[, start]))
     expect_equal(nrow(starts[r==TRUE,]), 13)
     ##mergeOverlaps    
    ## expect_identical(gg$mergeOverlaps(), gg)
   ##  nodes1 = c(GRanges("1",IRanges(1,100),"*"), GRanges("1",IRanges(50,200),"*"),
    #            GRanges("1",IRanges(201,300),"*"), GRanges("1",IRanges(250,400),"*"),
    ##            GRanges("1",IRanges(401,500),"*"))

  ##   gg=gGraph$new(nodes=nodes1, edges=edges)    
  ##     expect_equal(length(gg$mergeOverlaps()), 7)
  
  
  
})


test_that('gWalk works', {   

    ##create gWalk with null grl      
    nodes1 = c(GRanges("1",IRanges(1,100),"*"), GRanges("1",IRanges(101,200),"*"),
               GRanges("1",IRanges(201,300),"*"), GRanges("1",IRanges(301,400),"*"),
               GRanges("1",IRanges(401,500),"*"))
    edges = data.table(n1 = c(3,2,4,1,3), n2 = c(3,4,2,5,4), n1.side = c(1,1,0,0,1), n2.side = c(0,0,0,1,0))
    gg = gGraph$new(nodes = nodes1, edges = edges)
    
    col=data.table(x=1)
    gw=gWalk$new(snode.id=2,sedge.id=NULL, grl=NULL, graph=gg)
    expect_identical(gw$graph, gg)
    expect_equal(gw$length, 1)
    expect_identical(gw$nodes$dt, gg$nodes[2]$dt)
    ##empty gWalk
    empt=gWalk$new()
    expect_equal(length(empt), 0)

    ##set function
    expect_error(gw$set("walk.id"=1))
    gw$set(x=3)
    expect_equal(gw$dt[, x], 3)

    
    
    ##create gWalk with null sedge.id
   ## gw1=gWalk$new(snode.id=1, sedge.id=NULL, grl=NULL, graph=gg, meta=col)
   ## expect_equal(unlist(gw1$dt[, snode.id]), 1)

    ##create gWalk with null snode.id
   ## col=data.table(x=1:3)   
   ## gw4=gWalk$new(snode.id=NULL, sedge.id=c(1:3), grl=NULL, graph=gg, meta=col)
   ## expect_identical(gw4$dt[,walk.id], c(1:3))    
    
    ##subsetting
   ## gw2=gWalk$new(snode.id=c(1:3),sedge.id=NULL, grl=NULL, graph=gg, meta=col)
   ## expect_identical(gw[1]$dt, gw$dt)
    
    

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
##     nodes1 = gUtils::gr.fix(nodes1, HGSL)
##     seq = HGSL
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
##     loosenodes = gUtils::gr.fix(loosenodes, HGSL)

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
##     loosenodes = gUtils::gr.fix(loosenodes, HGSL)

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
##     gg = ggnew$simpleGraph(genome = HGSL)
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


test_that('gGraph, simplify', {
     nodes = c(GRanges("1", IRanges(1001,2000), "*"), GRanges("1", IRanges(2001,3000), "*"),
               GRanges("1", IRanges(3001,4000), "*"), GRanges("1", IRanges(4001,5000), "*"),
               GRanges("1", IRanges(5001,6000), "*"), GRanges("1", IRanges(6001,7000), "*"),
               GRanges("1", IRanges(7001,8000), "*"), GRanges("1", IRanges(8001,9000), "*"),
               GRanges("1", IRanges(9001,10000), "*"))
    
     edges = data.table(n1 = c(1,2,3,5,6,7,8,1,4,6,3),
                        n2 = c(2,3,4,6,7,8,9,4,8,6,2),
                        n1.side = c(1,1,1,1,1,1,1,1,1,0,1),
                        n2.side = c(0,0,0,0,0,0,0,0,0,1,0))         
     graph = gGraph$new(nodes = nodes, edges = edges)
     g=copy(graph)    
     graphSimple = graph$simplify()    
     ## Check to make sure the simplifed version correctly merged everything
     nodes = c(GRanges("1", IRanges(1001,2000), "*"), GRanges("1", IRanges(2001,4000), "*"),
               GRanges("1", IRanges(4001,5000), "*"),
               GRanges("1", IRanges(5001,6000), "*"), GRanges("1", IRanges(6001,7000), "*"),
               GRanges("1", IRanges(7001,8000), "*"), GRanges("1", IRanges(8001,10000), "*"))
     nodes = gr.fix(nodes, hg_seqlengths())    
     edges = data.table(n1 = c(1,2,2,3,1,4,5,6,5),
                        n2 = c(2,3,2,7,3,5,6,7,5),
                        n1.side = c(1,1,1,1,1,1,1,1,0),
                        n2.side = c(0,0,0,0,0,0,0,0,1))
     expect_equal(length(graphSimple), 7)
     expect_equal(length(graphSimple$edges), 9)
     dt1=gr2dt(nodes)
     dt2=gr2dt(graphSimple$nodes$gr)
     expect_equal(dt1[, start], dt2[, start])
     expect_equal(dt1[, end], dt2[, end])
     g$disjoin()     
    ## expect_equal(graphSimple$edges$dt[order(n1,n2,n1.side,n2.side),n1],
      ##            edges[order(n1,n2,n1.side,n2.side),n1])
   ##  expect_equal(graphSimple$edges$dt[order(n1,n2,n1.side,n2.side),n2],
     ##             edges[order(n1,n2,n1.side,n2.side),n2])    
})


test_that('%&% works for gNodes and gEdges', {
    ##gNode querying
    nodes = c(GRanges("1", IRanges(1001,2000), "*"), GRanges("1", IRanges(2001,3000), "*"),
               GRanges("1", IRanges(3001,4000), "*"), GRanges("1", IRanges(4001,5000), "*"),
               GRanges("1", IRanges(5001,6000), "*"), GRanges("1", IRanges(6001,7000), "*"),
               GRanges("1", IRanges(7001,8000), "*"), GRanges("1", IRanges(8001,9000), "*"),
               GRanges("1", IRanges(9001,10000), "*"))    
     edges = data.table(n1 = c(1,2,3,5,6,7,8,1,4,6,3),
                        n2 = c(2,3,4,6,7,8,9,4,8,6,2),
                        n1.side = c(1,1,1,1,1,1,1,1,1,0,1),
                        n2.side = c(0,0,0,0,0,0,0,0,0,1,0))     
     graph = gGraph$new(nodes = nodes, edges = edges)
     gr=GRanges("1", IRanges(3500, 8500))    
     gn=graph$nodes
     gr1=gn %&% gr
     gr1=gr1$gr
     expect_equal(graph$nodes$gr[3:8], gr1)

     ##gEdge querying
     ##1. gEdge with gEdge
     ge1=graph$edges[4:7]
     ge2=graph$edges[5:8]
     ge3=ge1 %&% ge2
     expect_equal(ge3, graph$edges[5:7])
     ##2. gEdge with Junction    
     ge4=ge1 %&% ge2$junctions
     expect_equal(ge4, graph$edges[5:7])
     
})

test_that('dist measures distances correctly, c.gGraph works', {
    ##dist
    nodes = c(GRanges("1", IRanges(1001,2000), "*"), GRanges("1", IRanges(2001,3000), "*"),
               GRanges("1", IRanges(3001,4000), "*"), GRanges("1", IRanges(4001,5000), "*"),
               GRanges("1", IRanges(5001,6000), "*"), GRanges("1", IRanges(6001,7000), "*"),
               GRanges("1", IRanges(7001,8000), "*"), GRanges("1", IRanges(8001,9000), "*"),
               GRanges("1", IRanges(9001,10000), "*"))    
     edges = data.table(n1 = c(1,2,3,5,6,7,8,1,4,6,3),
                        n2 = c(2,3,4,6,7,8,9,4,8,6,2),
                        n1.side = c(1,1,1,1,1,1,1,1,1,0,1),
                        n2.side = c(0,0,0,0,0,0,0,0,0,1,0))     
     graph = gGraph$new(nodes = nodes, edges = edges)    
     gr1=GRanges("1", IRanges(3500, 8500))
     graph$dist(gr1)    
     ##c.gGraph
     expect_equal(length(c(graph, graph)), 18)
     bg=c(graph, graph)
     expect_equal(bg$nodes$dt[1:9, start], graph$nodes$dt[, start])
     ##subsetting
     expect_equal(length(graph[loose.left==FALSE]), 7)
     expect_equal(length(graph[, n1.side=="right"]),9)
    
     ##     gr2
##     expect_equal(graph$dist(gr1), 0)
  ##   expect_equal(graph$
})



## test_that('gGraph, window', {
##     ## Cases
##     ## 1) There is not a cn
##     ## 2) different pad values

##     ## CASE 1: pad = 0
##     nodes = c(GRanges("1", IRanges(400,600), "*"), GRanges("1", IRanges(500,800), "*"), GRanges("1", IRanges(1000,1400), "*"))
##     graph = gGraph$new(nodes = nodes)

##     result = c(GRanges("1", IRanges(400,800), "*"), GRanges("1", IRanges(1000,1400), "*"))
##     result = gr.fix(result, HGSL)

##     expect_equal(graph$window(), graph$win)
##     expect_equal(graph$win, result)

##     ## Case 2: pad != 0

##     result = GRanges("1", IRanges(200,1600), "*")
##     result = gr.fix(result, HGSL)

##     expect_equal(graph$window(200), result)
## })


 test_that('gGraph$json', {

     ## Make sure it throws an error when the graph is empty
     gg = gGraph$new()
     expect_error(gg$gg2js())
     
     ## Check the empty edge case and no.y
     nodes = c(GRanges("1", IRanges(1001,2000), "*"), GRanges("1", IRanges(2001,3000), "*"),
               GRanges("1", IRanges(3001,4000), "*"), GRanges("1", IRanges(4001,5000), "*"),
               GRanges("1", IRanges(5001,6000), "*"), GRanges("1", IRanges(6001,7000), "*"),
               GRanges("1", IRanges(7001,8000), "*"), GRanges("1", IRanges(8001,9000), "*"),
               GRanges("1", IRanges(9001,10000), "*"))
     gg = gGraph$new(nodes = nodes)
##     json = gg$json(save=F, no.y=T)
     
  ##   expect_equal(nrow(json$connections), 0)
 ##    expect_equal(nrow(json$intervals), 9)
  ##   expect_false(json$settings$y_axis$visible)
     
    ## Test saving and loading a file using jabba data - comparison file is visually pre checked
     gg = gGraph$new(jabba = jab)
  ##   gg$json("../../inst/extdata/data.with.cn.test.json")
      
   ##  json = fromJSON(system.file('extdata', 'data.with.cn.test.json', package="gGnome"))
    ## json1 = fromJSON(system.file('extdata', 'data.with.cn.json', package="gGnome"))
     
    ## expect_equal(json, json1)
       })
      
      ##test_that('connect nodes makes appropriate edge', {    
      ##  nodes = c(GRanges("1", IRanges(1001,2000), "*"), GRanges("1", IRanges(2001,3000), "*"),
    ##           GRanges("1", IRanges(3001,4000), "*"), GRanges("1", IRanges(4001,5000), "*"),
      ##         GRanges("1", IRanges(5001,6000), "*"), GRanges("1", IRanges(6001,7000), "*"),
        ##       GRanges("1", IRanges(7001,8000), "*"), GRanges("1", IRanges(8001,9000), "*"),
          ##     GRanges("1", IRanges(9001,10000), "*"))    
##     edges = data.table(n1 = c(1,2,3,5,6,7,8,1,4,6,3),
  ##                      n2 = c(2,3,4,6,7,8,9,4,8,6,2),
    ##                    n1.side = c(1,1,1,1,1,1,1,1,1,0,1),
      ##                  n2.side = c(0,0,0,0,0,0,0,0,0,1,0))     
   ##
   ## graph = gGraph$new(nodes = nodes, edges = edges)    


    ##some errors
   ## expect_error(gGraph$connectNodes(13, 15))
   ## expect_error(gGraph$connectNodes(c(1, 8), 3))

    
##})



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

## ## ##-------------------------------------------------------##
## ## ## test_that('able to make JSON output', {
## ## ##     expect_true(inherits(jab.gw <<- as(readRDS(jab.gw.grl), "gWalks"), "gWalks"))
## ## ##     expect_equal(jab.gw$json("testing_gw.json"), "testing_gw.json")
## ## ## })

