library(testthat)
library(gUtils)
library(gTrack)

svaba = system.file('extdata', "HCC1143.svaba.somatic.sv.vcf", package = "gGnome")
delly = system.file('extdata', "delly.final.vcf.gz", package = "gGnome")
novobreak = system.file('extdata', "novoBreak.pass.flt.vcf", package = "gGnome")

HGSL = c("1"=249250621, "2"=243199373, "3"=198022430, "4"=191154276, "5"=180915260, "6"=171115067, "7"=159138663, "X"=155270560, "8"=146364022, "9"=141213431, "10"=135534747, "11"=135006516, "12"=133851895, "13"=115169878, "14"=107349540, "15"=102531392, "16"=90354753, "17"=81195210, "18"=78077248, "20"=63025520, "Y"=59373566, "19"=59128983, "22"=51304566, "21"=48129895, "M"=16571)


context('testing gGnome')

## message("Toy segments: ", system.file('extdata', 'testing.segs.rds', package="gGnome"))
## test_segs = readRDS(system.file('extdata', 'testing.segs.rds', package="gGnome"))

## message("Toy edges: ", system.file('extdata', 'testing.es.rds', package="gGnome"))
## test_es = readRDS(system.file('extdata', 'testing.es.rds', package="gGnome"))

## jab = system.file('extdata', 'jabba.simple.rds', package="gGnome")
## message("JaBbA result: ", jab)

## jab.gw.grl = system.file('extdata', 'gw.grl.rds', package = "gGnome")
## message("Example JaBbA genome walks: ", jab.gw.grl)

## prego = system.file('extdata', 'intervalFile.results', package='gGnome')
## message("PREGO results: ", prego)

## weaver = system.file('extdata', 'weaver', package='gGnome')
## message("Weaver results: ", weaver)

## remixt = system.file('extdata', 'remixt', package='gGnome')
## message("remixt results: ", remixt)

## genome = seqinfo(test_segs)

## ## 

## gr2 = GRanges(1, IRanges(c(1,9), c(6,14)), strand=c('+','-'), seqinfo=Seqinfo("1", 25), field=c(1,2))
## dt = data.table(seqnames=1, start=c(2,5,10), end=c(3,8,15))

## pr=gGraph$new(prego=prego)

## ## FIXME: REQUIREMENTS FOR THE BELOW TESTS
## ## TEST FOR GGRAPH DEFAULT CONSTRUCTOR
## ## TEST FOR QUERYLOOKUP


## test_that('Constructors', {       
##    expect_is(gGraph$new(jabba = jab), "gGraph")
##    expect_is(gGraph$new(weaver = weaver), "gGraph")
##    expect_is(pr, "gGraph")
##    expect_is(gGraph$new(remixt = remixt), "gGraph")
## })


## test_that('gNode Class Constructor/length, gGraph length/active $nodes', {
##      nodes1 = c(GRanges("1",IRanges(1,100),"*"), GRanges("1",IRanges(101,200),"*"),
##                 GRanges("1",IRanges(201,300),"*"), GRanges("1",IRanges(301,400),"*"),
##                 GRanges("1",IRanges(401,500),"*"))
##      edges = data.table(n1 = c(3,2,4,1,3), n2 = c(3,4,2,5,4), n1.side = c(1,1,0,0,1), n2.side = c(0,0,0,1,0))

##      gg = gGraph$new(nodes = nodes1, edges = edges)    
##      ## Testing some errors here    
##      expect_error(gNode$new(0, gg))
##      expect_error(gNode$new(100, gg))
##      expect_error(gNode$new(-30, gg))
##      expect_error(gNode$new(4, gGraph$new()))
##      expect_error(gNode$new(c(1,2,-10), gg))
     
##      ## Testing with no snode.id
##      gn = gNode$new(graph = gg)

##      expect_equal(length(gn), 0)
##      expect_identical(gn$graph, gg)
##      expect_equal(length(gn$sid), 0)
     
##      ## Testing building using active fields of gg
##      gn = gg$nodes
##      gn$mark(col="purple")
     
##      expect_is(gn, "gNode")
##      expect_equal(length(gn), length(gg))
##      expect_equal(gn$gr, gg$gr %Q% (strand == "+"))
##      expect_equal(gn$id, (gg$gr %Q% (strand == "+"))$snode.id)
##      expect_equal(gn$id, gn$sid)

##      ##gNode ego
##      gn=gNode$new(3, gg)
##      expect_equal(gg[c(3:4)]$nodes$dt[, start], gn$ego(1)$dt[, start])
##      expect_equal(gg[c(3, 4, 2)]$nodes$dt[, start], gn$ego(2)$dt[, start])

##      ##loose.left, loose right
##      gn=gNode$new(1, gg)
##      expect_equal(gn$loose.left, FALSE)
##      expect_equal(gn$loose.right, TRUE)
##      expect_equal(gn$degree, 1)

##      ##terminal, degrees
##      expect_equal(gr2dt(gn$terminal)[, start], 100)
##      expect_equal(gn$ldegree, 1)
##      expect_equal(gn$rdegree, 0)     

##      ##unioning gNodes
##  ##    gn2=gNode$new(2, gg)    
 
##     ## expect_equal(union(gn, gn2)$dt, gg$nodes[c(1:2)]$dt)
##     ## gg=gGraph$new(nodes=nodes1, edges=edges)
##     ## expect_error(union(gg$nodes, gn2))
     
## })

## test_that('gNode subsetting', {
##     nodes1 = c(GRanges("1",IRanges(1,100),"*", cn=1), GRanges("1",IRanges(101,200),"*",cn=1),
##                GRanges("1",IRanges(201,300),"*", cn=1), GRanges("1",IRanges(301,400),"*", cn=1),
##                GRanges("1",IRanges(401,500),"*",cn=1))
##     edges = data.table(n1 = c(3,2,4,1,3), n2 = c(3,4,2,5,4), n1.side = c(1,1,0,0,1), n2.side = c(0,0,0,1,0))        
##     gg = gGraph$new(nodes = nodes1, edges = edges)         
##     gn = gg$nodes    
##     gn2= gNode$new(2, gg)
##     gn3= gNode$new(1, gg)   
   
##     ##Right and Left
##     expect_equal(gn2$right$dt[, start], 301)
##     expect_equal(gn2$left$dt[, start], 301)

##     ##Quick subgraph test
##     sub=gn$subgraph
##     expect_equal(length(sub), length(gg))   
##     ##Various overriden functions
##     node1=gn[1]
##     nodes2=gn[c(1:4)]
##     expect_equal(length(c(gn2, gn2)), 2)
##     expect_equal(length(setdiff(gn, gn)), 0)
##     expect_identical(intersect(gn3, gg$nodes)$dt, gg$nodes[1]$dt)
    
##     ## Testing positive indicies
##     expect_equal(gn[1:5]$gr, gn$gr)
##     expect_equal(length(gn[1:5]), 5)
##     expect_identical(gn[1:5]$graph, gg)
    
##     expect_equal(gn[c(2,4)]$gr, gn$gr[c(2,4)])
##     expect_equal(length(gn[c(2,4)]), 2)
##     expect_equal(gn[c(2,4)]$id, c(2,4))
    
##     expect_equal(gn[1]$gr, gn$gr[1])
##     expect_equal(length(gn[1]), 1)
##     expect_equal(gn[1]$graph, gg)

##     expect_error(gn[6])
##     expect_error(gn[-10])
##     expect_error(gn[0])
    
##     ## Indexing via snode.id name
##     expect_equal(gn["4"]$gr, gg$gr[4])
##     expect_equal(gn["-1"]$gr, gg$gr[6])
##     expect_equal(gn[c("1","-3")]$gr, gg$gr[c(1,8)])

##     expect_error(gn["10"])
##     expect_error(gn["-6"])
    
##     ## Some complex indexing
##     expect_equal(gn["1"][-1], gn[-1])
##     expect_equal(gn[2:4]["-2"][-1], gn[2])
##     expect_equal(gn["4"][-1][-1], gn["4"])
##     expect_equal(gn[c(-1,-4)]$gr, gg$gr[c(6,9)])
   
##     ## Indexing via queries
## ##    expect_identical(gn[loose.left == FALSE]$dt, gn[c(1:4)]$dt) ## causing problems on travis only
## #     expect_equal(gn[snode.id > 3], gn[4:5]) #### this syntax is not playing well with Travis
## ##    expect_identical(gn[start > 150 & end < 450 & loose.left == FALSE]$dt, gn[3:4]$dt) #### this syntax is not playing well with Travis
## })

## test_that('gEdge works',{
##     nodes1 = c(GRanges("1",IRanges(1,100),"*"), GRanges("1",IRanges(101,200),"*"),
##                GRanges("1",IRanges(201,300),"*"), GRanges("1",IRanges(301,400),"*"),
##                GRanges("1",IRanges(401,500),"*"))
##     edges = data.table(n1 = c(3,2,4,1,3), n2 = c(3,4,2,5,4), n1.side = c(1,1,0,0,1), n2.side = c(0,0,0,1,0))     
##     gg = gGraph$new(nodes = nodes1, edges = edges)    
##     ge=gEdge$new(1, gg)
##     ge2=gEdge$new(2, gg)
##     ge3=c(ge, ge2)
    
##     ##Mark
##     ge$mark(changed=TRUE)
##     expect_equal(gg$edges[1]$dt[, changed], TRUE)
    
##     ##Basic active fields
##     expect_equal(ge$length, 1)     
##     expect_equal(length(ge), 1)
##     expect_equal(ge$left, gg$nodes[3])    
##     expect_equal(ge$right, gg$nodes[3])
    
##     ##Overriden functions   
##     expect_equal(c(ge, ge2)$dt[, sedge.id], gg$edges[c(1, 2)]$dt[, sedge.id])    
##     expect_equal(length(setdiff(ge, ge)), 0)
##     expect_equal(intersect(c(ge, ge3), ge)$dt[, sedge.id], 1)          
    
##     ##Miscellaneous functions    
##     starts=gr2dt(ge$junctions$grl)[, start]
##     expect_equal(gr2dt(ge$junctions$grl)[, start][1], 300)
##     expect_equal(gr2dt(ge$junctions$grl)[, end][2], 201)
##     expect_equal(class(ge$subgraph)[1], "gGraph")
    
## })



## test_that('Junction', {   
##     juncs = readRDS(jab)$junctions
##     badjuncs = GRangesList(GRanges("1", IRanges(1,100), "+"))   
##     ## Some errors
##     expect_error(Junction$new(GRanges("1", IRanges(1,100), "+")))
##     expect_error(Junction$new(badjuncs))   
##     ## Empty Junctions
##     jj = Junction$new(GRangesList())
##     expect_equal(length(jj), 0)
    
##     ## Build   
##     jj = Junction$new(juncs)
##     expect_equal(length(setdiff(jj, jj)), 0)
##     expect_equal(length(intersect(jj, jj)),length(jj))
    
##     expect_equal(length(jj), 500)
##     expect_equal(unlist(jj$grl), unlist(juncs))
   
##    ## expect_equal(jj$dt, as.data.table(values(juncs)))
    
##     ## c() / +
##     jj2 = c(jj, jj)
##     expect_equal(length(jj2), length(jj)*2)   
##     expect_equal(as.matrix(as.data.table(jj2$grl)), as.matrix(as.data.table(c(jj$grl, jj$grl))))
##   juncs.delly = Junction$new(delly)
##   juncs.novobreak = Junction$new(novobreak)
##   juncs.svaba = Junction$new(svaba)

##   expect_true(length(juncs.delly)==210)
##   expect_true(length(juncs.novobreak)==421)
##   expect_true(length(juncs.svaba)==500)
##   juncs = readRDS(jab)$junctions
##   badjuncs = GRangesList(GRanges("1", IRanges(1,100), "+"))   
##   ## Some errors
##   expect_error(Junction$new(GRanges("1", IRanges(1,100), "+")))
##   expect_error(Junction$new(badjuncs))   
##   ## Empty Junctions
##   jj = Junction$new(GRangesList())
##   expect_equal(length(jj), 0)
  
##   ## Build   
##   jj = Junction$new(juncs)
##   expect_equal(length(setdiff(jj, jj)), 0)
##   expect_equal(length(intersect(jj, jj)),length(jj))
  
##     expect_equal(length(jj), 500)
##     expect_equal(unlist(jj$grl), unlist(juncs))   
##     expect_equal(jj$dt, as.data.table(values(juncs)))
  
##   ## c() / +
##   jj2 = c(jj, jj)
##   expect_equal(length(jj2), length(jj)*2)
##   expect_equal(as.matrix(as.data.table(jj2$grl)), as.matrix(as.data.table(c(jj$grl, jj$grl))))

##   jj2 = jj + 100
##   expect_true(all(all(width(jj2$grl)==101)))
    
##   ## $breakpoints
##     bps = jj$breakpoints
##     bps = c(GenomicRanges::shift(bps %Q% (strand == "+"), 1), bps %Q% (strand == "-"))
##     expect_equal(length(findOverlaps(unlist(juncs), bps)), length(juncs)*2)
    
##     bps = jj2$breakpoints
##     bps = c(GenomicRanges::shift(bps %Q% (strand == "+"), 1), bps %Q% (strand == "-"))
##   expect_equal(length(findOverlaps(unlist(juncs), bps)), length(juncs)*2)   
##     ## $gGraph
##     gg = jj$graph
##     gg1 = gGraph$new(juncs = juncs)
    
##     expect_is(gg, "gGraph")
##     expect_equal(length(gg), length(gg1))
    
##     ##subset
##     expect_equal(gr2dt(jj[1]$grl)[,group][1], 1)        

## })

## test_that('gGraph, empty constructor/length', {
##   ## Test with unspecified genome - seqinfo(default values)
##   gg = gGraph$new()
##   expect_true(is(gg, 'gGraph'))
  
##   ## Check the active bindings
##   expect_equal(length(seqinfo(gg)), 0)
  
##   ## Test with specified genome - valid - just check seqinfo (set.seqinfo w/valid genome)
##   ##gg = gGraph$new(genome = HGSL)
##   ##expect_equal(length(gg$seqinfo), 25)
## })

## ## TESTING TRIM
## test_that('gGraph, trim', {
##   ## 1) trim within single node
##   ## 2) trim across multiple nodes
##   ## 3) Trim is not continuous
  
##   ## Build a gGraph
##   gr = c(GRanges("1", IRanges(1001,2000), "*"), GRanges("1", IRanges(2001,3000), "*"),
##          GRanges("1", IRanges(3001,4000), "*"), GRanges("1", IRanges(4001,5000), "*"))
##   gr = gr.fix(gr, HGSL)
  
##   es = data.table(n1 = c(1,2,3,4,4), n2 = c(2,3,4,1,3), n1.side = c(1,1,1,0,1), n2.side = c(0,0,0,1,0))
  
##   graph = gGraph$new(nodes = gr, edges = es)

##   ## CASE 1
##   gr1 = GRanges("1", IRanges(1200,1500), "+")
##   gr1 = gr.fix(gr1, HGSL)   
##   tmp = graph$copy$trim(gr1)   
##   expect_identical(as.data.table(streduce(tmp$nodes$gr)), as.data.table(streduce(gr1)))
##   ##keep as is  expect_equal(granges(tmp$nodes$gr), granges(gr1))
##   ##keep as is    expect_equal(length(gg$edges), 0)
##   expect_equal(length(tmp), 1)
  
##   ## CASE 2
##   gr2 = GRanges("1", IRanges(2200,4500), "+")
##   gr2 = gr.fix(gr2, HGSL)
##   edges = data.table(n1 = c(1,2), n2 = c(2,3), n1.side = c(1,1), n2.side = c(0,0))
  
##   tmp2 = graph$copy$trim(gr2)
##   expect_identical(as.data.table(streduce(tmp2$nodes$gr)), as.data.table(streduce(gr2)))  
##   es = tmp2$edges$dt[order(n1,n2,n1.side,n2.side),][, c("type") := NULL]    
##   ## expect_equal(es, edges[order(n1,n2,n1.side,n2.side),])

##    expect_equal(length(tmp2), 3)
   
##     ## Case 3
##     ranges(gr2) = IRanges(1800,2200)
##     gr2 = c(gr2, GRanges("1", IRanges(2800,4500), "+"))
##    edges = data.table(n1 = c(2,4,5,6), n2 = c(3,5,6,2), n1.side = c(1,1,1,0), n2.side = c(0,0,0,1))      
##     graph$trim(c(gr1,gr2))
    
##     expect_equal(streduce(graph$nodes$gr), streduce(c(gr1,gr2)))
    
##     es = graph$edges$dt[order(n1,n2,n1.side,n2.side),][, c("type") := NULL]
##     ##  expect_equal(es, edges[order(n1,n2,n1.side,n2.side),])   
##     expect_equal(length(graph), 6)

    
## })

## test_that('some public gGraph fields',{
##      nodes1 = c(GRanges("1",IRanges(1,100),"*"), GRanges("1",IRanges(101,200),"*"),
##                 GRanges("1",IRanges(201,300),"*"), GRanges("1",IRanges(301,400),"*"),
##                 GRanges("1",IRanges(401,500),"*"))
##      edges = data.table(n1 = c(3,2,4,1,3), n2 = c(3,4,2,5,4), n1.side = c(1,1,0,0,1), n2.side = c(0,0,0,1,0))
     
##      ##gTrack
##      gg = gGraph$new(nodes = nodes1, edges = edges)        
##      expect_is(gg$gtrack(), "gTrack")
##      expect_is(gg$gt, "gTrack")
##      expect_equal(gg$gtrack()$ygap, 2)
##      expect_equal(gg$gtrack()$name, "gGraph")

##      #gwalks
##      expect_is(gg$walks(), "gWalk")
##      expect_equal(gg$walks()$dt[, wid], c(200, 300, 100))
 
##      ##subgraph
##      gr1=gg$nodes$gr[1]
##      sub=gg$subgraph(gr1, 100)
##      expect_equal(sub$dt[c(1:2), start], c(1, 402))    
##      expect_equal(sub$dt[2, loose.right], FALSE)
##      expect_equal(sub$dt[2, loose.left], TRUE)
     
##      ## ##addJuncs
##      ## graph=copy(gg)
##      ## graph$addJuncs(graph$junctions)
##      ## starts=data.table(r=duplicated(as.data.table(unlist(graph$junctions$grl))[, start]))
##      ## expect_equal(nrow(starts[r==TRUE,]), 2)

##      ##clusters
##      gg=gGraph$new(nodes=nodes1, edges=edges)
##      gg$clusters()
##      expect_equal(gg$dt[c(1:5), cluster], c(1, 2, 2, 2, 1))

##      ##eclusters
##      pr$eclusters()
     
##      ##dim
##      expect_equal(dim(gg), c(5, 5))  
     
##      ##mergeOverlaps    
##     ## expect_identical(gg$mergeOverlaps(), gg)
##    ##  nodes1 = c(GRanges("1",IRanges(1,100),"*"), GRanges("1",IRanges(50,200),"*"),
##     #            GRanges("1",IRanges(201,300),"*"), GRanges("1",IRanges(250,400),"*"),
##     ##            GRanges("1",IRanges(401,500),"*"))

##   ##   gg=gGraph$new(nodes=nodes1, edges=edges)    
##   ##     expect_equal(length(gg$mergeOverlaps()), 7)
  
  
  
## })


## test_that('gWalk works', {   

##     ##create gWalk with null grl      
##     nodes1 = c(GRanges("1",IRanges(1,100),"*"), GRanges("1",IRanges(101,200),"*"),
##                GRanges("1",IRanges(201,300),"*"), GRanges("1",IRanges(301,400),"*"),
##                GRanges("1",IRanges(401,500),"*"))
##     edges = data.table(n1 = c(3,2,4,1,3), n2 = c(3,4,2,5,4), n1.side = c(1,1,0,0,1), n2.side = c(0,0,0,1,0))
##     gg = gGraph$new(nodes = nodes1, edges = edges)
##     grl=GRangesList(GRanges("1",IRanges(1,100),"*"), GRanges("1",IRanges(101,200),"*"),
##                GRanges("1",IRanges(201,300),"*"), GRanges("1",IRanges(301,400),"*"),
##                GRanges("1",IRanges(401,500),"*"))
    
##     col=data.table(x=1)
##     gw=gWalk$new(snode.id=2,sedge.id=NULL, grl=NULL, graph=gg)
##     expect_identical(gw$graph, gg)
##     expect_equal(gw$length, 1)
##     expect_identical(gw$nodes$dt, gg$nodes[2]$dt)
##     ##empty gWalk
##     empt=gWalk$new()
##     expect_equal(length(empt), 0)

##     ##set function
##     expect_error(gw$set("walk.id"=1))
##     gw$set(x=3)
##     expect_equal(gw$dt[, x], 3)

##     ##gWalk from grl
##     gw2=gWalk$new(grl=grl)
##     expect_is(gw2, "gWalk")

##     ##subsetting   
##     expect_equal(unlist(gw2[1:3]$dt[, snode.id]), c(-1, -2, -3))
##     expect_equal(gw2[walk.id==3]$dt[, name], "3")
    
##     ##dts
##     expect_equal(gw2$dts(makelists=FALSE), gw2$dts()[, snode.id:=NULL])

##     ##disjoin
##     gw3=gWalk$new(snode.id=c(1:4), graph=gg)
##     gr=GRanges("1", IRanges(50, 250), "+")
##     gw3$disjoin(gr=gr)
##     expect_equal(unlist(gw3$dt[1, snode.id]), c(1, 2))
##     expect_equal(unlist(gw3$dt[3, snode.id]), c(4, 5))

##     ##simplify
##     gw2$simplify()

##     ##gTrack
##     expect_is(gw2$gtrack(), "gTrack")
    
    
    
##     ##create gWalk with null sedge.id
##    ## gw1=gWalk$new(snode.id=1, sedge.id=NULL, grl=NULL, graph=gg, meta=col)
##    ## expect_equal(unlist(gw1$dt[, snode.id]), 1)

##     ##create gWalk with null snode.id
##    ## col=data.table(x=1:3)   
##    ## gw4=gWalk$new(snode.id=NULL, sedge.id=c(1:3), grl=NULL, graph=gg, meta=col)
##    ## expect_identical(gw4$dt[,walk.id], c(1:3))    
    
##     ##subsetting
##    ## gw2=gWalk$new(snode.id=c(1:3),sedge.id=NULL, grl=NULL, graph=gg, meta=col)
##    ## expect_identical(gw[1]$dt, gw$dt)
    
    

## })

## test_that('querying functions work', {   
    
##     ##gNode querying %&%
##     nodes = c(GRanges("1", IRanges(1001,2000), "*"), GRanges("1", IRanges(2001,3000), "*"),
##                GRanges("1", IRanges(3001,4000), "*"), GRanges("1", IRanges(4001,5000), "*"),
##                GRanges("1", IRanges(5001,6000), "*"), GRanges("1", IRanges(6001,7000), "*"),
##                GRanges("1", IRanges(7001,8000), "*"), GRanges("1", IRanges(8001,9000), "*"),
##                GRanges("1", IRanges(9001,10000), "*"))    
##      edges = data.table(n1 = c(1,2,3,5,6,7,8,1,4,6,3),
##                         n2 = c(2,3,4,6,7,8,9,4,8,6,2),
##                         n1.side = c(1,1,1,1,1,1,1,1,1,0,1),
##                         n2.side = c(0,0,0,0,0,0,0,0,0,1,0))     
##      graph = gGraph$new(nodes = nodes, edges = edges)     
##      gr=GRanges("1", IRanges(3500, 8500))
##      grl=GRangesList(gr)
##      gn=graph$nodes
##      gr1=gn %&% gr
##      gr1=gr1$gr
##      expect_equal(graph$nodes$gr[3:8], gr1)

##      ##gEdge querying %&%
##      ##1. gEdge with gEdge
##      ge1=graph$edges[4:7]
##      ge2=graph$edges[5:8]
##      ge3=ge1 %&% ge2
##      expect_equal(ge3, graph$edges[5:7])
##      ##2. gEdge with Junction    
##      ge4=ge1 %&% ge2$junctions
##      expect_equal(ge4, graph$edges[5:7])

##      ##gWalk querying %&%

##      ##junction querying %&%
     
##      ##gedge %^%
##      ##1. gEdge with gEdge
##      expect_equal(ge1 %^% ge2, c(FALSE, TRUE, TRUE, TRUE))
##      expect_equal(ge1  %^% grl, c(FALSE, FALSE, FALSE, FALSE))
##      ##2. gEdge with junction
##      expect_equal(ge1 %^% graph$junctions, c(TRUE, TRUE, TRUE, TRUE))
##      expect_equal(ge1 %^% ge2, ge1 %^% ge2$junctions)          
##      ##3. gEdge with gNode
##      expect_equal(ge1 %^% gn[7:8], c(FALSE, TRUE, TRUE, TRUE))
##      expect_equal(ge1[1] %^% gn[5], TRUE) 

##      ##gNode %^%
##      ##1. gNode with gNode
##      expect_equal(gn[1:4] %^% gn [4:9], c(FALSE, FALSE, FALSE, TRUE))     
##      ##2. gNode with gEdge
##      expect_equal(gn[4:6] %^% ge1, c(FALSE, TRUE, TRUE))
##      ##3. gNode with GRangesList
##      expect_equal(gn %^% gr, c(FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE))
     
##      ##junction %^%
##      expect_equal(ge1$junctions %^% ge2$junctions, ge1 %^% ge2) 
     
##      ##gWalk %^%
##      gw=gWalk$new(grl=grl)
##      ##gWalk with junction
##      expect_equal(gw %^% graph$junctions[7], FALSE)
##      expect_equal(gw %^% graph$junctions[4], TRUE)
##      ##gWalk with gEdge
##      expect_equal(gw %^% ge1, TRUE)
##      expect_equal(gw %^% ge2[3], FALSE)
##      ##gWalk with gNode
##      expect_equal(gw %^% gn[3], TRUE)
##      expect_equal(gw %^% gn[1], FALSE)
## })
## ## ## FIXME: Definitely could add more tests to catch errors and things in constructor but it seems to be working
## ## test_that('gGraph, gNodes and Edges Constructor/active bindings/loosegNodes', {
## ##     ## Make sure it add the right nodes and edges
## ##     ## 1) edges = NULL
## ##     ##     a) looseterm = FALSE
## ##     ##     b) looseterm = TRUE
## ##     ## 2) edges != NULL
## ##     ##     a) looseterm = FALSE
## ##     ##     b) looseterm = TRUE
## ##     ## - check edges, nodes, seqinfo, loosegNodes

## ##     ## CASE 1a
## ##     nodes1 = c(GRanges("1",IRanges(1,100),"*"), GRanges("1",IRanges(101,200),"*"),
## ##                GRanges("1",IRanges(201,300),"*"), GRanges("1",IRanges(301,400),"*"),
## ##                GRanges("1",IRanges(401,500),"*"))

## ##     ## Set the seqinfo to make sure it carries over
## ##     nodes1 = gUtils::gr.fix(nodes1, HGSL)
## ##     seq = HGSL
## ##     names(seq) = NULL

## ##     ## Check the null case
## ##     gg = gGraph$new(nodes = GRanges(), looseterm = FALSE)

## ##     expect_equal(gg$nodes, GRanges())
## ##     expect_equal(dim(gg$edges)[1], 0)
## ##     expect_error(gg$graph, NA)
## ##     expect_equal(length(seqinfo(gg)), 0)

## ##     ## Check with looseterm = F
## ##     gg = gGraph$new(nodes = nodes1, looseterm = FALSE)

## ##     expect_equal(sort(granges(gg$nodes)), sort(granges(nodes1)))
## ##     expect_equal(dim(gg$edges)[1], 0)
## ##     expect_equal(length(gg$loosegNodes()), 0)

## ##     ## CASE 1b
## ##     gg = gGraph$new(nodes = nodes1, looseterm = TRUE)

## ##     loosenodes = c(GRanges("1",IRanges(1,1),"*"), GRanges("1",IRanges(100,100),"*"),
## ##                    GRanges("1",IRanges(101,101),"*"), GRanges("1",IRanges(200,200),"*"),
## ##                    GRanges("1",IRanges(201,201),"*"), GRanges("1",IRanges(300,300),"*"),
## ##                    GRanges("1",IRanges(301,301),"*"), GRanges("1",IRanges(400,400),"*"),
## ##                    GRanges("1",IRanges(401,401),"*"), GRanges("1",IRanges(500,500),"*"))
## ##     loosenodes = gUtils::gr.fix(loosenodes, HGSL)

## ##     expect_equal(sort(granges(nodes1)), sort(granges(gg$nodes)))
## ##     expect_equal(dim(gg$edges)[1], 0)
## ##     expect_equal(length(gg$loosegNodes()), 10)
## ##     expect_equal(sort(granges(gg$loosegNodes())), sort(granges(loosenodes)))

## ##     ## CASE 2a
## ##     edges = data.table(n1 = c(3,2,4,1,3), n2 = c(3,4,2,5,4), n1.side = c(1,1,0,0,1), n2.side = c(0,0,0,1,0))

## ##     gg = gGraph$new(nodes = nodes1, edges = edges, looseterm = FALSE)

## ##     expect_true(is(gg, 'gGraph'))

## ##     ## Check that the correct things are stored in our gGraph
## ##     expect_equal(sort(granges(nodes1)), sort(granges(gg$nodes)))

## ##     es = gg$edges[order(n1,n2,n1.side,n2.side),][, c("type") := NULL]
## ##     expect_equal(es, edges[order(n1,n2,n1.side,n2.side),])
## ##     expect_equal(length(gg), 5)
## ##     expect_equal(length(gg$loosegNodes()), 0)
## ##     expect_equal(seqinfo(gg)@seqlengths, seq)

## ##     ## CASE 2b
## ##     gg = gGraph$new(nodes = nodes1, edges = edges, looseterm = TRUE)

## ##     loosenodes = c(GRanges("1", IRanges(100,100), "*"),
## ##                    GRanges("1", IRanges(400,400), "*"), GRanges("1", IRanges(401,401), "*"))
## ##     loosenodes = gUtils::gr.fix(loosenodes, HGSL)

## ##     ## Check that the correct things are stored in our gGraph
## ##     expect_equal(sort(granges(nodes1)), sort(granges(gg$nodes)))

## ##     es = gg$edges[order(n1,n2,n1.side,n2.side),][, c("type") := NULL]

## ##     expect_equal(es, edges[order(n1,n2,n1.side,n2.side),])
## ##     expect_equal(length(gg), 5)
## ##     expect_equal(sort(granges(gg$loosegNodes())), sort(granges(loosenodes)))

## ##     gg1 = gGraph$new(nodes = gg$nodes, edges = gg$edges, looseterm = TRUE)

## ##     expect_equal(sort(granges(gg1$nodes)), sort(granges(gg$nodes)))
## ##     expect_equal(gg1$edges[order(n1,n2,n1.side,n2.side),], gg$edges[order(n1,n2,n1.side,n2.side),])
## ##     expect_equal(length(gg), length(gg1))
## ##     expect_equal(sort(granges(gg$loosegNodes())), sort(granges(loosenodes)))

## ## })


## ## ## Just checks that the jabba input is correct processed not that the result is correct
## ## test_that('gGraph, jabba input', {
## ##     expect_error(gGraph$new(jabba = "badinput"))
## ##     badinput = list(1,2,3)
## ##     names(badinput) = c("segstats", "ab.edges", "purity")
## ##     expect_error(gGraph$new(jabba = badinput))

## ##     gg = gGraph$new(jabba = jab)
## ##     expect_true(is(gg, "gGraph"))

## ##     gg = gGraph$new(jabba = jab, regular=T)
## ##     gg1 = gGraph$new()$simpleGraph()

## ##     expect_equal(length(gr.findoverlaps(gg$nodes, gg1$nodes)), length(gg$nodes))
## ## })


## ## ## FIXME: test the circular thing idk how to do that
## ## test_that('gGraph, simpleGraph', {
## ##     ggnew = gGraph$new()

## ##     ## simpleGraph = function(genome = NULL, chr=FALSE, include.junk=FALSE)
## ##     ## Testing simpleGraph(default values)
## ##     gg = ggnew$simpleGraph()
## ##     expect_true(is(gg, 'gGraph'))
## ##     expect_equal(length(gg$nodes), 25)
## ##     expect_equal(dim(gg$edges)[1], 0)
## ##     expect_equal(length(gg$junctions), 0)
## ##     expect_error(gg$graph, NA) ## check it works
## ##     ##expect_equal(length(gg$adj), 2500) -- FIXME: Don't know what this does but it doesn't work
## ##     ## expect_equal(gg$parts, NULL) -- FIXME: also broken, all don't know what it does
## ##     expect_equal(length(gg$seqinfo), 25)
## ##     expect_equal(length(gg$loosegNodes()), 50)

## ##     ## Don't know what the fuck is going on here
## ##     expect_true(is(gg$td, 'gTrack'))
## ##     expect_equal((gg$td)$ygap, 2)
## ##     expect_match((gg$td)$name, 'CN')
## ##     expect_equal(length(gg$win), 25)

## ##     ## Testing simpleGraph with genome = Something
## ##     gg = ggnew$simpleGraph(genome = HGSL)
## ##     expect_equal(length(gg$nodes), 25)
## ##     expect_equal(length(gg$loosegNodes()), 50)
## ## })





## ## ## Test the addSegs/makeSegs functionality
## ## test_that('gGraph, addSegs/karyograph no juncs', {

## ##     expect_equal(gGraph$new()$nodes, gGraph$new(tile=NULL)$nodes)
## ##     expect_equal(gGraph$new()$edges, gGraph$new(tile=NULL)$edges)

## ##     ## Check again at addSegs level
## ##     gg = gGraph$new()
## ##     expect_error(gg$addSegs(tile=NULL))
## ##     expect_error(gg$addSegs(tile=5))

## ##     gr = GRanges("fakeChr", IRanges(1,10000), "*")
## ##     expect_error(gGraph$new(tile=gr))

## ##     gr = c(GRanges("1", IRanges(1001,2000), "*"), GRanges("1", IRanges(2001,3000), "*"),
## ##            GRanges("2", IRanges(3001,4000), "*"), GRanges("3", IRanges(5000,6000), "*"),
## ##            GRanges("3", IRanges(5500,6500), "*"))

## ##     expect_error(gg$addSegs(tile=gr))

## ##     ## Check building a gGraph from a tile, no cn
## ##     gg = gGraph$new(tile = gr)
## ##     gg1 = gGraph$new()$simpleGraph()
## ##     nodes = gr.breaks(gr,gg1$nodes)    

## ##     expect_equal(length(gg), 34)
## ##     expect_equal(length(gg$loosegNodes()),50) 
## ##     expect_equal(sort(granges(gg$nodes)), sort(granges(nodes)))

## ##     pairs = table(gg$nodes$pair.id)
## ##     expect_true(all(pairs == 1))

## ##     edges = data.table(n1 = c(1,2,3,5,6,8,9,10,11),
## ##                        n2 = c(2,3,4,6,7,9,10,11,12),
## ##                        n1.side = c(1,1,1,1,1,1,1,1,1),
## ##                        n2.side = c(0,0,0,0,0,0,0,0,0))

## ##     expect_equal(nrow(gg$edges), 9)
## ##     es = gg$edges[order(n1,n2,n1.side,n2.side),][, c("type") := NULL]
## ##     expect_equal(es, edges[order(n1,n2,n1.side,n2.side),][, c("type") := NULL])

## ##     ## Try again with cn

## ## })


## ## test_that('gGraph, mergeGraphs', {
## ##     ## CASES
## ##     ## 1) Self is null, adding null
## ##     ## 2) Self is null, adding non-null ---- Could be caught with a clause a the beginning to save work
## ##     ## 3) Self is non-null, adding null
## ##     ## 4) Self is non-null, adding non-null
## ##     ##     a) Neither have edges
## ##     ##     b) Either has edges but not both
## ##     ##     b) Elaborate both have edges

## ##     gr = GRanges("1",IRanges(1,1000),"*")

## ##     graph = gGraph$new(nodes = gr)

## ##     ## CASE 1
## ##     nullgg = gGraph$new()
## ##     tmp = nullgg$mergeGraphs(nullgg, decouple=F)

## ##     expect_equal(nullgg, tmp)

## ##     ## CASE 2
## ##     tmp = nullgg$mergeGraphs(graph, decouple=F)

## ##     expect_equal(graph, tmp)

## ##     ## CASE 3
## ##     tmp = graph$mergeGraphs(nullgg, decouple=F)

## ##     expect_equal(graph, tmp)

## ##     ## CASE 4a
## ##     gr1 = GRanges("1", IRanges(4001,5000), "*")

## ##     graph1 = gGraph$new(nodes = gr1)
## ##     tmp = graph$mergeGraphs(graph1, decouple=F)

## ##     ## Check nodes and edges
## ##     expect_equal(sort(granges(c(graph$nodes,graph1$nodes))), sort(granges(tmp$nodes)))
## ##     expect_equal(dim(tmp$edges)[1], 0)
## ##     expect_equal(length(tmp), 2)

## ##     ## CASE 4b
## ##     nodes = c(GRanges("1",IRanges(1,100),"*"), GRanges("1",IRanges(101,200),"*"),
## ##               GRanges("1",IRanges(201,300),"*"), GRanges("1",IRanges(301,400),"*"),
## ##               GRanges("1",IRanges(401,500),"*"))

## ##     nodes1 = c(GRanges("1",IRanges(150,250),"*"), GRanges("1",IRanges(401,500),"*"),
## ##                GRanges("1",IRanges(601,1000),"*"))

## ##     edges = data.table(n1 = c(3,2,4,1,5), n2 = c(3,4,2,5,1), n1.side = c(1,1,0,0,1), n2.side = c(0,0,0,1,0))
## ##     edges1 = data.table(n1 = c(2,1,3,2), n2 = c(3,2,1,2), n1.side = c(1,1,0,0), n2.side = c(0,0,1,1))

## ##     ## graph has edges and graph1 does not
## ##     graph = gGraph$new(nodes = nodes, edges = edges)
## ##     graph1 = gGraph$new(nodes = nodes1)
## ##     tmp = graph$mergeGraphs(graph1, decouple=F)

## ##     expect_equal(sort(granges(tmp$nodes)), sort(granges(c(graph$nodes, graph1$nodes))))

## ##     es = graph$edges[order(n1,n2,n1.side,n2.side),][, c("type") := NULL]
## ##     expect_equal(es, tmp$edges[order(n1,n2,n1.side,n2.side),][, c("type") := NULL])
## ##     expect_equal(length(tmp), 8)    

## ##     ## graph does not have edges and graph1 does
## ##     graph = gGraph$new(nodes = nodes)
## ##     graph1 = gGraph$new(nodes = nodes1, edges = edges1)
## ##     tmp = graph$mergeGraphs(graph1, decouple=F)

## ##     expect_equal(sort(granges(tmp$nodes)), sort(granges(c(graph$nodes, graph1$nodes))))

## ##     es = graph1$edges[order(n1,n2,n1.side,n2.side),][, c("type") := NULL]
## ##     expect_equal(es, tmp$edges[order(n1,n2,n1.side,n2.side),][, c("type") := NULL])
## ##     expect_equal(length(tmp), 8)

## ##     ## CASE 4c
## ##     graph = gGraph$new(nodes = nodes, edges = edges)
## ##     graph1 = gGraph$new(nodes = nodes1, edges = edges1)
## ##     tmp = graph$clone()
## ##     tmp$mergeGraphs(graph1, mod=T)

## ##     expect_equal(sort(granges(tmp$nodes)), sort(granges(c(graph$nodes, graph1$nodes))))

## ##     es = rbind(graph$edges, graph1$edges[, ":="(n1 = n1 + length(graph), n2 = n2 + length(graph))])
## ##     es = es[order(n1,n2,n1.side,n2.side),][, c("type") := NULL]
## ##     expect_equal(es, tmp$edges[order(n1,n2,n1.side,n2.side),][, c("type") := NULL])
## ##     expect_equal(length(tmp), 8)


## ##     ## TESTS FOR SUBGRAPHS
## ##     ## Check that all the subgraphs are null
## ##     ##for(i in 1:length(tmp$subgraphs)) {
## ##     ##    expect_equal(tmp$subgraphs[[i]]$length(), 0)
## ##     ##}

## ##     ## Case 4b - We know the main level works, just check subgraphs
## ##     ##subgraphs = list(gGraph$new(tile = gr1[1])$trim(gr1[1]),
## ##     ##                 gGraph$new(tile = gr1[2])$trim(gr1[2]))

## ##     ##graph1$resetSubgraphs(subgraphs)
## ##     ##tmp = graph$mergeGraphs(graph1, decouple=F)

## ##     ## FIXME: insert subgraph test

## ##     ## Case 4c
## ##     ##subgraphs = list(gGraph$new(tile = gr[1])$trim(gr[1]),
## ##     ##                 gGraph$new(tile = gr[2])$trim(gr[2]))

## ##     ##graph$resetSubgraphs(subgraphs)
## ##     ##tmp = graph$mergeGraphs(graph1, decouple=F)

## ##     ## FIXME: insert subgraph test

## ## })


## ## test_that('gGraph, mergeOverlaps', {

## ## })


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
})


test_that('gGnome tutorial', {
  pdf('test.pdf')
    ## SvAbA
  svaba = jJ(system.file('extdata', "HCC1143.svaba.somatic.sv.vcf", package = "gGnome"))
  expect_equal(length(svaba), 500)
  ## DELLY
  delly = jJ(system.file('extdata', "delly.final.vcf.gz", package = "gGnome"))
  expect_equal(length(delly), 210)
  ## novobreak
  novobreak = jJ(system.file('extdata', "novoBreak.pass.flt.vcf", package = "gGnome"))
  expect_equal(length(novobreak), 421)
  ## BEDPE
  bedpe = jJ(system.file('extdata', "junctions.bedpe", package = "gGnome"))
  expect_equal(length(bedpe), 83)

    ## can use both row and column subsetting on Junction metadata
  head(novobreak[1:2, 1:10])

  ## can use data.table style expressions on metadata to subset Junctions
  ## here, filter novobreak translocations with quality greater than 50
  novobreak[ALT == "<TRA>" & QUAL>50, 1:10][1:2, 1:5]

  ## subsetting SvAbA junctions with >5 bases of homologous sequence
  svaba[nchar(INSERTION)>10, ][1:2, 1:5]

  ## subsetting SVabA junctions with same sign and nearby breakpoints (i.e. small $span)
  svaba[svaba$sign>0 & svaba$span<1e5][1:2, 1:5]

  ## subsetting junctions with infinite span (ie different chromosome) and homology length >5
  delly[is.infinite(delly$span) & HOMLEN>5, ][1:2,1:5]

  ## subset svaba by those intersect with DELLY using gUtils subset %&% operator
  length(svaba %&% delly)
  
  ## increase the overlap substanntially by padding delly calls with 100bp (using + operator)
  length(svaba %&% (delly+100))
  
  ## basic set operations also work
  length(S4Vectors::setdiff(svaba, delly+100))
  
  length(S4Vectors::union(svaba, delly+100))
  
  ## gGraph from svaba input
  gg = gG(juncs = svaba)

  ## we use gTrack to plot the gTrack associated with this gGraph
  ## the second argument to gTrack plot is a string or GRanges representing the
  ## window to plot, the links argument enables drawing of junctions from GRangesList

  ## generate breaks using gUtils function to tile genome at evenly spaced 1MB intervals
  breaks = gr.tile(seqinfo(svaba), 1e6)

  ## gGraph from svaba input
  gg2 = gG(breaks = breaks, juncs = svaba)

  ## set gGraph metadata
  gg2$set(name = 'with breaks')

  ## compare graphs towards the beginning of chromosome 1
  ## (gTracks can be concatenated to plot multiple tracks)

  nodes = gr.tile(seqlengths(svaba)["1"], 1e7)

  ## generate 20 random edges (n1, n2, n1.side, n2.side)
  edges = data.table(
    n1 = sample(length(nodes), 20, replace = TRUE),
    n2 = sample(length(nodes), 20, replace = TRUE))
  edges[, n1.side := ifelse(runif(.N)>0.5, 'right', 'left')]
  edges[, n2.side := ifelse(runif(.N)>0.5, 'right', 'left')]

  gg3 = gG(nodes = nodes, edges = edges)

  pad = runif(length(nodes))*width(nodes)
  
  gg3 = gG(nodes = nodes + pad , edges = edges)
  
  ## PREGO is the original cancer SV graph caller from Oesper et al 2012
  gg.prego = gG(prego = system.file('extdata/hcc1954', 'prego', package='gGnome'))

  ## Weaver is from Li et al 2016
  gg.weaver = gG(weaver = system.file('extdata/hcc1954', 'weaver', package='gGnome'))

  ## RemiXt is from McPherson et al 2018
  gg.remixt = gG(remixt = system.file('extdata/hcc1954', 'remixt', package='gGnome'))

  ## JaBbA is from Imielinski Lab, Yao et al (in preparation)
  gg.jabba = gG(jabba = system.file('extdata/hcc1954', 'jabba.rds', package="gGnome"))

  plot(c(gg.prego$gt, gg.weaver$gt, gg.remixt$gt, gg.jabba$gt), '4')

  ## accessing the meta data features of a gGraph
  gg.remixt$meta

  ## not that the y.field points to a column of the node metadata, accessed
  ## by $nodes$dt
  gg.remixt$nodes$dt[1:2, ]

  ## setting the y.field to NULL (previously "cn")
  gg.remixt$set(y.field = NULL)

  ## now the gTrack when plotted on chromosome 4 will no longer plot 
  ## "cn", instead the nodes / segments will stack to stay out of each other's 
  ## way 

  ## returns gNode
  gg.jabba$nodes[1:2]

  ## node indices can be negative, in which case the node orientation is flipped
  ## (note difference from standard R subsetting syntax for negative indices)
  gg.jabba$nodes[-c(1:2)]

  ## returns GRanges
  gg.jabba$nodes[1:2]$gr

  ## returns data.table
  gg.jabba$nodes$dt[1:2]

  ## returns gEdge
  gg.jabba$edges[1:2]

  ## returns data.table
  gg.jabba$edges$dt[1:2]

  ## returns Junction
  gg.jabba$edges$junctions[1:2]

  ## select high copy gNode's in JaBbA object
  highcopy = gg.jabba$nodes[cn>100, ]

  ## use $dt gNode accessor to get the value of metadata associated with these nodes
  mean(highcopy$dt$cn)

  ## select edges associated with a long INSERTION character string
  biginsert = gg.jabba$edges[nchar(INSERTION)>20, ]

  ## subset ALT edges
  gg$edges[type == 'ALT']

  ## enumerate ALT edges classes
  table(gg$edges[type == 'ALT']$dt$class)

  ## subset INV-like edges
  gg$edges[class == 'INV-like']

  ## use the from= and to= arguments with signed node ids
  ## to query for edges connecting specifying node sets

  ## this is an "INV-like" ALT edge connecting the right side of 6 to the right side of 8
  gg.jabba$edges[from = 6, to = -8]

  ## this is another "INV-like" ALT edge connecting the left side of 6 to the left side of 8
  gg.jabba$edges[from = -6, to = 8]

  ## find the distributed of FILTER metadata among these junctions harboring an insertion
  table(biginsert$dt$FILTER)

  highcopy$left[cn<20]

  ## all of the nodes connected to a junction with a templated insertion
  biginsert$nodes

  ## the reference edge connected to the right of first node connected to the
  ## "biginsert" junction (i.e. the one with a templated insertion)
  biginsert$nodes[1]$eright[type == 'REF']


  ## note that these two expressions will give the same output
  gg.jabba$nodes[1]$right[1:2]

  gg.jabba$nodes[-1]$left[-c(2:1)]

  gencode = track.gencode(stack.gap = 1e5, cex.label = 0.8, height = 20)

  ## note that we use the $gtrack() method instead of the $gt active binding because
  ## we erased the $y.field metadata from gg.remixt above.
  ## not that the second argument tells gTrack to plot in the vicinity of the high copy nodes
  plot(c(gencode, gg.remixt$gtrack(y.field = 'cn'), gg.jabba$gt), highcopy$gr+1e5)

  ## select ReMiXT edges overlapping JaBba edges with long insertions
  eli.remixt = gg.remixt$edges %&% biginsert

  ## select ReMiXT nodes overlapping JaBba nodes connected to JaBbA edges with long insertions
  nli.remixt = gg.remixt$nodes %&% biginsert$nodes

  eli.remixt$mark(col = 'blue')

  ## set the metadata column "col" of remixt nodes overlapping long insert jabba edges to the value "green"
  nli.remixt$mark(col = 'green')

  ## we can also mark the analogous regions in the JaBbA model
  biginsert$mark(col = 'blue') ## marking gEdge
  biginsert$nodes$mark(col = 'green') ## marking gNode associated with these gEdges

  ## these metadata values will be interpreted by gTrack as segment and connection colors
  ## we plot both jabba and remixt objects (which have been changed by the above commands)
  ## near the vicinity of 5 of the biginsert junctions
  #plot(c(gg.remixt$gtrack(y.field = 'cn'), gg.jabba$gt), unlist(biginsert$grl[1:2])[, c()]+1e6)

  ## gGraph uses the bracket syntax for subsetting, where the indices before the comma
  ## corresponds to nodes and after the comma correspond to edges
  ## as with gNode and gEdge, both integers and metadata expressions will work 
  ggs1 = gg.jabba[cn>100, ]
  ggs1$set(name = 'gGraph\nsubset')

  ## this syntax is equivalent to the above, but uses the `gNode` subgraph command
  ggs2 = gg.jabba$nodes[cn>100]$subgraph
  ggs2$set(name = 'gNode\nsubgraph')

  ## here instead of subsetting on nodes, we trim the graph around a set of `GRanges`
  ## note that trim works in place, so if we want another copy we use $copy to clone the
  ## gGraph object
  ggs3 = gg.jabba$copy$trim(highcopy$gr+1e5)
  ggs3$set(name = 'trimmed\nJaBbA')

  ## since this function uses GRanges as input, we can apply it to the ReMiXT graph
  ## as well. 
  ggs4 = gg.remixt$copy$trim(highcopy$gr+1e5)
  ggs4$set(name = 'trimmed\nRemiXT', y.field = 'cn')

#  plot(c(gencode, ggs1$gt, ggs2$gt, ggs3$gt, ggs4$gt), highcopy$gr+2e5)

  ## define a simple window on chromosome 1
  win = GRanges('1:1-1e7')

  ## create a simpel graph with 3MB bins
  tiles1 = gr.tile(win, 3e6);
  gg1 = gG(breaks = tiles1, meta = data.table(name = 'gg1'))

  ## create a second graph tiling the window with 2 MB bins
  tiles2 = gr.tile(win, 2e6);
  gg2 = gG(breaks = tiles2, meta = data.table(name = 'gg2'))

  ## this gGraph metadata will tell gTrack to plot the node.id with each node
  gg1$set(gr.labelfield = 'node.id')
  gg2$set(gr.labelfield = 'node.id')

  ## plot these two simple graphs
#  plot(c(gg1$gt, gg2$gt), win)

  ## concatenate gg1 and gg2 
  gg3 = c(gg1, gg2)
  gg3$set(name = 'c(gg1, gg2)', height = 20)

  ## disjoin gg3 collapses the graphs into each other
  ## by taking the disjoin bins of any overlapping nodes
  gg3d = gg3$copy$disjoin()
  gg3d$set(name = 'disjoined')

  ## simplify collapses reference adjacent nodes that lack
  ## an intervening ALT junction or loose end.
  gg3ds = gg3d$copy$simplify()
  gg3ds$set(name = 'simplified')

  ## reduce is equivalent to a disjoin followed by a simplify
  gg3r = gg3$copy$reduce()
  gg3r$set(name = 'reduced')

  ## plot
  plot(c(gg3$gt, gg3d$gt, gg3ds$gt, gg3r$gt), win)

  ## randomly sample 4 width 1 GRanges representing SNV
  snv = gr.sample(win, 4, wid = 1)

  ## disjoin with gr= argument breaks the graph at these SNV
  gg1d = gg1$copy$disjoin(gr = snv)

  ## plot results
  plot(gg1d$gt, win)

  ## new edges are specified as data.table with n1, n2, n1.side and n2.side
  gg1d$add(edges = data.table(n1 = 3, n1.side = 'left', n2 = 7, n2.side = 'right'))

  ## plot 
#  plot(gg1d$gt, win)


  ## connect syntax specifies edges as pairs of "signed" node ids
  ## this means that the edge leaves the <right> side of the first signed node
  ## and enters the <left> side of the second signed node

  ## thus here we create an edge leaving the right side of 5 and entering the right side of 8
  ## (i.e. the left side of -8)
  gg1d$connect(5, -8)

  ## this connects the right side of 3 and the left side of 9
  gg1d$connect(3, 9)

  ## plot 
#  plot(gg1d$gt, win)


  ## adding a junction to a copy of gg3
  gg3j = gg3$copy$add(junctions = svaba[7])
  gg3j$set(name = 'add junction')

  ## note that we have instantiated 4 separate edges, connecting all
  ## eligible breakpoint pairs on our input graph
  gg3j$edges[type == 'ALT', ]

  ## alternatively let's create a disjoint graph containing the breakpoints of this junction
  bp = unlist(svaba[7]$grl)

  ## this uses an alternative syntax of disjoin with collapse = FALSE flag (ie where we only
  ## do a partial disjoin by chopping up reference nodes without collapsing overlapping nodes)
  gg3d = gg3$copy$disjoin(gr = bp, collapse = FALSE)
  gg3d$set(name = 'disjoin w bp')

  ## now we can add an ALT edge to just one of the 4 pairs of breakpoints
  gg3de = gg3d$copy$connect(18,-14)
  gg3de$set(name = 'connect')

  ## plot results
  plot(c(gg3j$gt, gg3d$gt, gg3de$gt), win)

  ## copy gg2
  gg2$set(name = 'original')
  gg2c = gg2$copy

  ## replaces current copy of the first SNV with three separate copies
  ## i.e. representing different variants
  gg2c$rep(2, 3)

  ## rep adds a metadata field "parent.node.id" to the graph
  ## which allows us to track the original node.id prior to replication
  ## we set gr.labelfield here to plot the parent.node.id instead of the node.id
  ## to see this correspondence
  gg2c$set(name = 'replicate', gr.labelfield = 'parent.node.id')

  ##
##  plot(c(gg2$gt, gg2c$gt), win)


  ## n1 is the third copy of the node previously known as 2
  n1 = gg2c$nodes[parent.node.id == 2][1]

  ## N2 is the node previously known as 4
  n2 = gg2c$nodes[parent.node.id == 4]

  ## retrieve the path from these two nodes and replicate it in the graph
  p = gg2c$paths(n1,n2)
  gg2c$rep(p, 2)

  ## replot with current node.id
  gg2c$set(gr.labelfield = 'node.id')
##  plot(c(gg2c$gt), win)

  ## we can now use $connect with gNode arguments to connect these two alleles
  gg2c$connect(5, 10)

  ## now let's mark the (new) shortest path between nodes 1 and 2
  p = gg2c$paths(1, 2)
  p$mark(col = 'pink')

##  plot(c(gg2c$gt), win)

  ## label weakly connected components in the jabba graph
  ## the $cluster node metadata field stores the output of running this method
  gg.jabba$clusters('weak')

  ## inspecting these shows that most nodes are part of a large weakly-connected component
  ## and the remaining are part of 1-node clusters
  sort(table(gg.jabba$nodes$dt$cluster))

  ## analysis of strongly connected components reveals some more structure
  ## we peek at one of these clusters, marking its nodes blue
  gg.jabba$clusters('strong')

  gg.jabba$nodes[cluster == 240]$mark(col = 'blue')

  ## then plotting shows an interesting amplicon
##  plot(gg.jabba$gt, gg.jabba$nodes[cluster == 240]$gr+1e5)


  ## we first select nodes that are 1 Mbp in width, then compute clusters
  gg.jabba$nodes[width<1e6]$clusters('weak')

  ## note that this syntax still sets the $clusters metadata field of the original 
  ## graph, giving any >1 Mbp nodes a cluster ID of NA
  table(is.na(gg.jabba$nodes$dt$cluster), gg.jabba$nodes$dt$width>1e6)

  ## we peek at one of these interesting clusters, marking it with a blue color
  gg.jabba$nodes[cluster == 111]$mark(col = 'green')

  ## interestingly, we have re-discovered the ERBB2 BFB-driven amplification highlighted above
  gg.jabba$set(height = 30)
##  plot(c(gencode, gg.jabba$gt), gg.jabba$nodes[cluster == 111]$gr[, c()]+1e5)


  ## this will populate the ALT edges of the gGraph with metadata fields $ecluster, $ecycle, and $epath
  ## where $ecluster is the concatenation of $ecycle and $epath labels
  gg.jabba$eclusters()

  ## paths are labeled by a "p" prefix, and cycles labeled by a "c" prefix
  ## here we see a multi-junction cluster p52 with 6 edges
  sort(table(gg.jabba$edges$dt$ecluster))

  ## we can mark these edges (and their associated nodes) with a special color
  gg.jabba$edges[ecluster == "p52"]$mark(col = 'purple')
  gg.jabba$edges[ecluster == "p52"]$nodes$mark(col = 'purple')

  ## here, the edges and nodes of the cluster that we have discovered
  ## are highlighted in purple 
##  plot(c(gg.jabba$gt), unlist(gg.jabba$edges[ecluster == "p52"]$grl)[, c()]+1e4)

  ## path between nodes 1 and 10
  p1 = gg.jabba$paths(1, 1000)
  p1

  ## $dt accessor gives us walk metadata, which we can set with $set method
  ## by default it contains the $length which is the number of nodes in the walk, the width which is the
  ## total base pairs contain in the walk
  p1$set(name = 'my first walk')
  p1$dt

  ## like the gGraph, a gWalk contains $nodes and $edges accessors which correspond to all the nodes that
  ## and edges that contribute to that walk
  p1$nodes

  p1$edges

  ## these edges and nodes reference the gGraph from which they were derived
  ## that gGraph can be accessed via the $graph accessor
  identical(p1$graph, gg.jabba)

  ## we can view the nodes of the gWalk as GRangesList
  p1$grl


  ## sign only matters if we use ignore.strand = FALSE
  p2 = gg.jabba$paths(1, -1000)

  ## so p2 will be virtually identical to p1
  identical(p1$grl, p2$grl)

  ## if we compute path in a strand specific manner, then we may get a different result
  gg.jabba$paths(1, -1000, ignore.strand = FALSE)

  ## this long and winding path involves a completely separate chromosome, and may represent
  ## an interesting long range allele in this cancer genome.
  ## like with gNode and gEdge we can "mark" the nodes and edges of the gGraph that comprise
  ## this gWalk
  p1$mark(col = 'pink')
  gg.jabba$set(gr.labelfield = 'node.id')

  ## we can generate a gTrack from this via the $gt method and the GRanges comprising the  genomic "footprint"
  ## of this gWalk using the $footprint accessor
##  plot(c(gg.jabba$gt, p1$gt), p1$footprint+1e6)

  ## the paths method is vectorized, therefore we can provide a vector of sources and sinks
  ## as well as gNode arguments for either

  ## all low copy nodes on chromosome 1
  n1 = gg.jabba$nodes[seqnames == 2 & cn<=2]

  ## all low copy nodes on chromosome 21
  n2 = gg.jabba$nodes[seqnames == 21 & cn<=2]

  ## paths between low copy nodes on chromosome 11 and 17
  p4 = gg.jabba$paths(n1, n2)

  ## paths is vectorized so we can subset using integer indices
  p4[1:2]

  ## reverse complement gWalks using negative indices
  ## note the signs of the intervals in the $gr metadata string
  p4[-1]

  ## can subset gWalks using metadata
  ## e.g. we canlook for longer walks, eg those shorter than 10MB
  p5 = p4[wid<10e6]

  ## we can mark the nodes and edges of gWalk just as we would for a gNode or gEDge
  ## this marks the original graph
  p5$mark(col = 'purple')

  ## plot
 # plot(c(gg.jabba$gt, p5$gt), p5$footprint+1e6)

  ## this expression subset our graph to the nodes and edges
  ## that contribute to edge cluster "p52" (see previous section on clusters
  ## and communities)
  gg.sub = gg.jabba[, ecluster == 'p52']


  ## walk decompositiion of this small subgraph
  ## generates 28 possible linear alleles
  walks = gg.sub$walks()

  ## we can order these walks based on their width and choose the longest
  walks = walks[rev(order(wid))][1]

  ## and plot with the walk track up top (with nodes already marked purple from before)
#  plot(c(gg.jabba$gt, walks$gt), walks$footprint+1e4)


  ## retrieve signed node ids associated with a gWalk
  nid = p5$snode.id

  ## retrieve signed edge ids associated with a gWalk
  eid = p5$sedge.id

  ## you can instantiate a gWalk from signed node ids
  gW(snode.id = nid, graph = p5$graph)

  ## or from signed edge ids
  gW(sedge.id = eid, graph = p5$graph)

  ## not every node id or edge id sequence will work
  ## for example the reverse node sequence won't (necessarily) be in the graph
  revnid = list(rev(nid[[1]]))

  ## this will error out
  expect_error(gW(snode.id = revnid, graph = p5$graph))

  ## however the reverse complement (reverse and multiply by -1) of a legal
  ## sequence will always work
  rcnid = list(-rev(nid[[1]]))
  gW(snode.id = rcnid, graph = p5$graph)

  ## if we use drop = TRUE on a list of node ids we won't error out 
  ## if some of the walks are illegal, just return a gWalk whose length is shorter
  ## than the input
  nid2 = c(nid, revnid, rcnid)

  ## the result here is length 2, though the input is length 3
  ## this is because revnid is "ignored"
  p6 = gW(snode.id = nid2, graph = p5$graph, drop = TRUE)



  ## we can instantiate from GRangesList
  ## to demo we extract the grl from the walk above
  grl = p6$grl

  ## this command will thread these provided grl onto the existing graph
  p7 = gW(grl = grl, graph = p5$graph)

  ## reset the colors in our graph
  p7$mark(col = 'gray')

  ## let's say we chop up i.e. hypersegment the p5 graph
  ## the above instantiation will still work .. the output
  ## gWalk will however be "chopped up" to be compatible with the
  ## chopped up graph

  ## create 500kb tiles on p5's genome
  tiles = gr.tile(seqinfo(p5), 5e5);

  ## disjoin a copy p5$graph by these tiles
  ggd = p5$graph$copy$disjoin(gr = tiles)

  ## the disjoint graph has many more nodes and edges because we have added reference edges
  ## at every tile breakpoint
  dim(p5$graph)
  dim(ggd)

  ## however instantiation from the grl will still work
  p7d = gW(grl = grl, graph = ggd)
  p7d$graph$set(name = 'chopt')

  ## plotting will help visualize the differences
  ## you can see that the top version of the walk and the top version of
  ## the graph is more "chopped up"
  plot(c(p5$graph$gt, ggd$gt, p7[1]$gt, p7d[1]$gt), p7d$footprint + 1e5)



  ## let's create a new grl concatenating the original and "chopped" up grl
  grl2 = unname(grl.bind(p7$grl, p7d$grl))

  ## by default, disjoin = FALSE, and thus will not collapse the graphs corresponding to the inputted walks
  ## each walk will create a separate (linear) subgraph
  p8 = gW(grl = grl2)
  p8$nodes$mark(col = 'gray') ## reset node color
  p8$graph$set(name = 'non-dis')

  ## disjoin = TRUE will create a single disjoint gGraph that results from "collapsing" the walks represented
  ## by the input GLR
  p9 = gW(grl = grl2, disjoin = TRUE)
  p9$graph$set(name = 'disjoint')

  ## plotting these with the original graph to visualize
  ## you can see the "induced subgraph" for p8 and p9 only spans
  ## the footprint of these walks

  plot(c(ggd$gt, p8$graph$gt, p8$gt, p9$graph$gt, p9$gt), p9$footprint + 1e5)

  ## we can use simplify to "unchop" p8
  ## note that this will not collapse the disjoint paths, only remove reference edges,
  ## the resulting graph will continue to have four separate components representing each
  ## "haplotype"
  p8s = p8$copy$simplify()
  p8s$graph$set(name = 'simp', border = 'black')

  ## similarly we can use disjoin on the non-disjoint walks instantiated above
  ## this will collapse the graph to a non-overlapping set of nodes 
  p8d = p8$copy$disjoin()
  p8d$graph$set(name = 'disj', border = 'black')

  ## visualizing the results of the graphs
  gt = c(p8$graph$gt, p8$gt, p8s$graph$gt, p8s$gt, p8d$graph$gt, p8d$gt)
  gt$name = paste(gt$name, c('gG', 'gW'))
  plot(gt, p8$footprint + 1e5)


  ## revisiting walks traversing ecluster p52
  gg.sub = gg.jabba[, ecluster == 'p52']

  ## walk decompositiion of this small subgraph
  ## generates 28 possible linear alleles
  walks = gg.sub$walks()

  ## now we can use eval to annotate walks with how many ALT junctions they contain
  ## ALT junctions
  ## the expression evaluates the edge metadata field type and returns a scalar result,
  ## one for each walk
  numalt = walks$eval(sum(type == 'ALT'))

  ## we can set a new column in the walks metadata to this result
  walks$set(nalt = numalt)

  ## we use eval to identify the number of short intervals contained in this walk
  ## width is a node metadata 
  walks$set(nshort = walks$eval(sum(width<1e4)))

  ## by default eval tries to evalute the expression on nodes and then on edges
  ## if nodes and edges share some metadata field then we may want to specify
  ## exactly which data type we want eval to run on

  ## first let's use $mark to add a metadata field "type" to the nodes of this walk (edges
  ## by default already has a metadata field "type")
  walks$nodes$mark(type = 'bla')

  ## now if we rerun the above expression for numalt, it will give us a new result
  ## this is because the expression is successfully evaluated on the nodes metadata field
  ## "type"
  identical(walks$eval(sum(type == 'ALT')), numalt)

  ## if we specify edge= argument to $eval then we will get the old result
  ## i.e. forcing evaluation on the edge metadata
  identical(walks$eval(edge = sum(type == 'ALT')), numalt)

  ## and if we use force nodes evaluation with node=, we will again get a non-identical result
  identical(walks$eval(node = sum(type == 'ALT')), numalt)


  ## we need a GENCODE (style) object either as a GRanges (cached as an RDS on mskilab.com)
  ## or directly from GENCODE (https://www.gencodegenes.org/)
  gff = readRDS(gzcon(url('http://mskilab.com/gGnome/hg19/gencode.v19.annotation.gtf.gr.rds')))

  ## we are looking for any fusions connecting the genes CNOT6, ASAP1, and EXT1 using the
  ## genes= argument
  ## (we know there are complex fusions here because we've run a previous genome wide analysis,
  ## i.e. without setting the "genes =" argument, which discovered complex in.frame fusions in these genes)
  fus = fusions(gg.jabba, gff, genes = c('CNOT6', 'ASAP1', 'EXT1'))
  length(fus)

  ## fusions will output many "near duplicates" which just represent various combinations
  ## of near equivalent transcripts, we can filter these down using gWalk operations
  ufus = fus[in.frame == TRUE][!duplicated(genes)]

  ## there are 5 unique gene in-frame gene combinations, we plot the first
  ## connecting ASAP1 to CNOT6 with a chunk of intergenic genome in between

  ## ufus[1] connects the first 20 amino acids of ASAP1 to the downstream 400+
  ## amino acids of EXT1
  ufus[1]$dt$gene.pc

  ## this walk has 6 aberrant junctions, as shown by $numab metadata
  ufus[1]$dt$numab

  ## indeed that is verified by this expression
  length(ufus[1]$edges[type == 'ALT'])

  ## here we plot the walk on top of the JaBbA-derived  gGraph, which you will notice
  ## has been "chopped up" to include features of relevant genes. 
  plot(c(gencode, ufus$graph$gt, ufus[1]$gt), ufus[1]$footprint+1e4)

  ufus = fus[frame.rescue == TRUE]

  ## In this fusion model, a frame-shifted chunk of NSD1 spans 35 amino acids 
  ## and has been essentially inserted into the middle of an unrearranged
  ## ASAP1 transcript. 
  ufus[1]$dt$gene.pc

  ## there are 4 unique gene in-frame gene combinations, we plot the first
  ## connecting ASAP1 to CNOT6 with a chunk of intergenic genome in between
  ## here we plot the walk on top of the JaBbA-derived  gGraph, which you will notice
  ## has been "chopped up" to include features of relevant genes. 
  plot(c(gencode, ufus$graph$gt, ufus[1]$gt), ufus[1]$footprint+1e4)


  ## load Hnisz et al 2013 superenhancers mapped to hg19
  se = readRDS(gzcon(url('http://mskilab.com/gGnome/hg19/Hnisz2013.se.rds')))

  ## many of these are redundant / overlapping so we will reduce them with some padding
  ## to reduce computation
  ser = reduce(se)

  ## read gff (if did not do it above)
  ## gff = readRDS(gzcon(url('http://mskilab.com/gGnome/hg19/gencode.v19.annotation.gtf.gr.rds')))

  genes = gff %Q% (type == 'gene' & gene_name %in% c('TERT', 'BRD9'))

  ## useful (optional) params to tune for performance include "chunksize" and "mc.cores"
  ## which will chunk up and parallelize the path search, in this case 1000
  ## intervals at a time across 5 cores, can also monitor progress with verbose = TRUE
  px = proximity(gg.jabba, ser, genes[, 'gene_name'], chunksize = 2000, mc.cores = 5)

  ## peek at the first proximity, we can see the reldist, altdist, refdist
  ## and additional metadata features inherited from the genes object
  px[1]

  ## make a gTrack for the super-enhancers, coloring by tissue
  gt.se = gTrack(se, gr.colorfield = 'tissue', name = 'SupEnh')

  ## plot the first super-enhancer connecting to BRD9
  px[1]$mark(col = 'purple')

  plot(c(gencode, gt.se, px$graph$gt, px[1]$gt), px[1]$footprint+1e5)

  ## use $eval to count ALT junctions for each walk
  px$set(numalt = px$eval(sum(type == 'ALT')))

  ## let's look for a superenhancer connecting to TERT
  this.px = px[numalt>2 & refdist == Inf & gene_name == 'TERT']

  expect_equal(length(this.px), 56)

  ## mark it up
  this.px[1]$mark(col = 'purple')
  plot(c(gencode, gt.se, this.px$graph$gt, this.px[1]$gt), this.px[1]$footprint+1e5)
  dev.off()
})
