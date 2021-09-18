library(testthat)
library(gUtils)
library(gTrack)

HGSL = c("1"=249250621, "2"=243199373, "3"=198022430, "4"=191154276, "5"=180915260, "6"=171115067, "7"=159138663, "X"=155270560, "8"=146364022, "9"=141213431, "10"=135534747, "11"=135006516, "12"=133851895, "13"=115169878, "14"=107349540, "15"=102531392, "16"=90354753, "17"=81195210, "18"=78077248, "20"=63025520, "Y"=59373566, "19"=59128983, "22"=51304566, "21"=48129895, "M"=16571)

context('testing gGnome tutorial')

test_that('gGnome tutorial', {
  setDTthreads(1)
  pdf('test.pdf')
    ## ----cache=FALSE,echo=FALSE,results="hide", message = FALSE, warning = FALSE----------
    knitr::opts_chunk$set(collapse = TRUE, fig.width = 8, fig.height = 8, message = FALSE, warning = FALSE)
    library(gTrack)
    library(rtracklayer)
    library(kableExtra)    
    library(magrittr)
    library(tidyr)


    ## ----cache=FALSE,message=FALSE, warning=FALSE, include = FALSE------------------------
    ## Load the package
    library(gGnome)


    ## ---- junctions, cache=FALSE,warning=FALSE--------------------------------------------
    ## SvAbA
    svaba = jJ(system.file('extdata', "HCC1143.svaba.somatic.sv.vcf", package = "gGnome"))

    ## DELLY
    delly = jJ(system.file('extdata', "delly.final.vcf.gz", package = "gGnome"))

    ## novobreak
    novobreak = jJ(system.file('extdata', "novoBreak.pass.flt.vcf", package = "gGnome"))

    ## BEDPE
    bedpe = jJ(system.file('extdata', "junctions.bedpe", package = "gGnome"))


    ## ----cache=FALSE,warning=FALSE, collapse = TRUE---------------------------------------
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


    ## ----cache=FALSE,warning=FALSE--------------------------------------------------------
    svaba$grl[1:2]


    ## ----cache=FALSE,warning=FALSE--------------------------------------------------------

    ## subset svaba by those intersect with DELLY using gUtils subset %&% operator
    length(svaba %&% delly)

    ## increase the overlap substanntially by padding delly calls with 100bp (using + operator)
    length(svaba %&% (delly+100))

    ## basic set operations also work
    length(setdiff(svaba, delly+100))

    ## length(union(svaba, delly+100))



    ## ----cache=FALSE,warning=FALSE--------------------------------------------------------
    ## any names work in the arguments to merge, these will be reflected in metadata as a $seen.by column
    ## (using padding of 1kb and c() to remove existing metadata)
    res = merge(svaba = svaba[, c()], delly = delly[, c()], 
                novo = novobreak[, c()], anynameworks = bedpe[,c()], pad = 1e3)
    head(res)

    ## here we can use the $dt (data.table) accessor quickly make an UpSetR plot
    library(UpSetR)

    ## munge res$dt into upset friendly format
    df = as.data.frame(sign(as.matrix(res$dt[,.(seen.by.anynameworks, seen.by.delly, seen.by.novo, seen.by.svaba)])))
    upset(df)


    ## ----cache=FALSE,warning=FALSE--------------------------------------------------------
    ### gGraph from svaba input
    gg = gG(juncs = svaba)



    ## ----cache=FALSE,warning=FALSE--------------------------------------------------------
    ## we use gTrack to plot the gTrack associated with this gGraph
    ## the second argument to gTrack plot is a string or GRanges representing the
    ## window to plot, the links argument enables drawing of junctions from GRangesList
    plot(gg$gt, '1', links = svaba$grl)


    ## ----cache=FALSE,warning=FALSE--------------------------------------------------------
    ### generate breaks using gUtils function to tile genome at evenly spaced 1MB intervals
    breaks = gr.tile(seqinfo(svaba), 1e6)

    ### gGraph from svaba input
    gg2 = gG(breaks = breaks, juncs = svaba)

    ### set gGraph metadata
    gg2$set(name = 'with breaks')

    ### compare graphs towards the beginning of chromosome 1
    ### (gTracks can be concatenated to plot multiple tracks)
    plot(c(gg$gt, gg2$gt), '1:1-5e7', links = svaba$grl)



    ## ----cache=FALSE,warning=FALSE--------------------------------------------------------
    ## tiling of chromosome 1 
    nodes = gr.tile(seqlengths(svaba)["1"], 1e7)

    ## generate 20 random edges (n1, n2, n1.side, n2.side)
    edges = data.table(
               n1 = sample(length(nodes), 20, replace = TRUE),
               n2 = sample(length(nodes), 20, replace = TRUE))
    edges[, n1.side := ifelse(runif(.N)>0.5, 'right', 'left')]
    edges[, n2.side := ifelse(runif(.N)>0.5, 'right', 'left')]

    gg3 = gG(nodes = nodes, edges = edges)

    plot(gg3$gt, '1')



    ## ----cache=FALSE,warning=FALSE--------------------------------------------------------
    pad = runif(length(nodes))*width(nodes)

    gg3 = gG(nodes = nodes + pad , edges = edges)

    plot(gg3$gt, '1')


    ## ----cache=FALSE,warning=FALSE, fig.height=15-----------------------------------------
    ## PREGO is the original cancer SV graph caller from Oesper et al 2012
    gg.prego = gG(prego = system.file('extdata/hcc1954', 'prego', package='gGnome'))

    ## Weaver is from Li et al 2016
    gg.weaver = gG(weaver = system.file('extdata/hcc1954', 'weaver', package='gGnome'))

    ## RemiXt is from McPherson et al 2018
    gg.remixt = gG(remixt = system.file('extdata/hcc1954', 'remixt', package='gGnome'))

    ## JaBbA is from Imielinski Lab, Yao et al (in preparation)
    gg.jabba = gG(jabba = system.file('extdata/hcc1954', 'jabba.rds', package="gGnome"))

    plot(c(gg.prego$gt, gg.weaver$gt, gg.remixt$gt, gg.jabba$gt), '4')


    ## ----cache=FALSE,warning=FALSE, fig.height=10-----------------------------------------

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
    plot(gg.remixt$gt, '4')


    ## ----cache=FALSE,warning=FALSE, fig.height=6------------------------------------------
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



    ## ----cache=FALSE,warning=FALSE--------------------------------------------------------
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



    ## ----cache=FALSE,warning=FALSE--------------------------------------------------------
    ## all of the low copy gNodes connected to the left of highcopy intervals
    highcopy$left[cn<20]

    ## all of the nodes connected to a junction with a templated insertion
    biginsert$nodes

    ## the reference edge connected to the right of first node connected to the
    ## "biginsert" junction (i.e. the one with a templated insertion)
    biginsert$nodes[1]$eright[type == 'REF']


    ## note that these two expressions will give the same output
    gg.jabba$nodes[1]$right[1:2]

    gg.jabba$nodes[-1]$left[-c(2:1)]



    ## ----cache=FALSE,warning=FALSE, fig.height = 12---------------------------------------
    ## track.gencode, pulls hg19 GENCODE by default, modulates the stacking and font sizes
    gencode = track.gencode(stack.gap = 1e5, cex.label = 0.8, height = 20)

    ## note that we use the $gtrack() method instead of the $gt active binding because
    ## we erased the $y.field metadata from gg.remixt above.
    ## not that the second argument tells gTrack to plot in the vicinity of the high copy nodes
    plot(c(gencode, gg.remixt$gtrack(y.field = 'cn'), gg.jabba$gt), highcopy$gr+1e5)


    ## ----cache=FALSE,warning=FALSE--------------------------------------------------------
    ## select ReMiXT edges overlapping JaBba edges with long insertions
    eli.remixt = gg.remixt$edges %&% biginsert

    ## select ReMiXT nodes overlapping JaBba nodes connected to JaBbA edges with long insertions
    nli.remixt = gg.remixt$nodes %&% biginsert$nodes



    ## ----cache=FALSE,warning=FALSE, fig.height=8------------------------------------------
    ## set the metadata column "col" of remixt edges corresponding to long insertion jabba nodes to the value "blue"
    eli.remixt$mark(col = 'blue')

    ## set the metadata column "col" of remixt nodes overlapping long insert jabba edges to the value "green"
    nli.remixt$mark(col = 'green')

    ## we can also mark the analogous regions in the JaBbA model
    biginsert$mark(col = 'blue') ## marking gEdge
    biginsert$nodes$mark(col = 'green') ## marking gNode associated with these gEdges

    ## these metadata values will be interpreted by gTrack as segment and connection colors
    ## we plot both jabba and remixt objects (which have been changed by the above commands)
    ## near the vicinity of 5 of the biginsert junctions
    plot(c(gg.remixt$gtrack(y.field = 'cn'), gg.jabba$gt), unlist(biginsert$grl[1:2])[, c()]+1e6)


    ## ----cache=FALSE,warning=FALSE, fig.height=18-----------------------------------------
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

    plot(c(gencode, ggs1$gt, ggs2$gt, ggs3$gt, ggs4$gt), highcopy$gr+2e5)



    ## ----cache=FALSE,warning=FALSE, fig.height=10-----------------------------------------
    ## define subgraph 100kbp around our "high copy region"
    ## we copy so we keep our gg.jabba intact
    ggs = gg.jabba$copy$subgraph(highcopy$gr, d = 1e5)

    ## adjust GENCODE track to keep things pretty
    gencode = track.gencode(stack.gap = 2e5, cex.label = 0.5, height = 10, name = 'GENCODE')

    ## we could plot ggs or just plot gg.jabba and use the ggs$footprint to guide us
    ## we set the upper y axis limit y1 to 20 to visualize low level copy number changes 
    plot(c(gencode, gg.jabba$gt), ggs$footprint+1e5, y1 = 20)



    ## ----cache=FALSE,warning=FALSE, fig.height=6------------------------------------------
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
    plot(c(gg1$gt, gg2$gt), win)


    ## ----cache=FALSE,warning=FALSE, fig.height=10-----------------------------------------
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



    ## ----cache=FALSE,warning=FALSE, fig.height=6------------------------------------------
    ## randomly sample 4 width 1 GRanges representing SNV
    snv = gr.sample(win, 4, wid = 1)

    ## disjoin with gr= argument breaks the graph at these SNV
    gg1d = gg1$copy$disjoin(gr = snv)

    ## plot results
    plot(gg1d$gt, win)



    ## ----cache=FALSE,warning=FALSE, fig.height=6------------------------------------------
    ## new edges are specified as data.table with n1, n2, n1.side and n2.side
    gg1d$add(edges = data.table(n1 = 3, n1.side = 'left', n2 = 7, n2.side = 'right'))

    ## plot 
    plot(gg1d$gt, win)



    ## ----cache=FALSE,warning=FALSE, fig.height=6------------------------------------------

    ## connect syntax specifies edges as pairs of "signed" node ids
    ## this means that the edge leaves the <right> side of the first signed node
    ## and enters the <left> side of the second signed node

    ## thus here we create an edge leaving the right side of 5 and entering the right side of 8
    ## (i.e. the left side of -8)
    gg1d$connect(5, -8)

    ## this connects the right side of 3 and the left side of 9
    gg1d$connect(3, 9)

    ## plot 
    plot(gg1d$gt, win)


    ## ----cache=FALSE,warning=FALSE, fig.height=16-----------------------------------------

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


    ## ----cache=FALSE,warning=FALSE, fig.height=6------------------------------------------

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
    plot(c(gg2$gt, gg2c$gt), win)


    ## ----cache=FALSE,warning=FALSE, fig.height=6------------------------------------------
    ## n1 is the third copy of the node previously known as 2
    n1 = gg2c$nodes[parent.node.id == 2][1]

    ## N2 is the node previously known as 4
    n2 = gg2c$nodes[parent.node.id == 4]

    ## retrieve the path from these two nodes and replicate it in the graph
    p = gg2c$paths(n1,n2)
    gg2c$rep(p, 2)

    ## replot with current node.id
    gg2c$set(gr.labelfield = 'node.id')
    plot(c(gg2c$gt), win)


    ## ----cache=FALSE,warning=FALSE, fig.height=6------------------------------------------
    ## we can now use $connect with gNode arguments to connect these two alleles
    gg2c$connect(5, 10)

    ## now let's mark the (new) shortest path between nodes 1 and 2
    p = gg2c$paths(1, 2)
    p$mark(col = 'pink')

    plot(c(gg2c$gt), win)


    ## ----cache=FALSE,warning=FALSE, fig.height=8------------------------------------------
    ## label weakly connected components in the jabba graph
    ## the $cluster node metadata field stores the output of running this method
    gg.jabba$clusters(mode = 'weak')

    ## inspecting these shows that most nodes are part of a large weakly-connected component
    ## and the remaining are part of 1-node clusters
    sort(table(gg.jabba$nodes$dt$cluster))

    ## analysis of strongly connected components reveals some more structure
    ## we peek at one of these clusters, marking its nodes blue
    gg.jabba$clusters(mode = 'strong')

    gg.jabba$nodes[cluster == 240]$mark(col = 'blue')

    ## then plotting shows an interesting amplicon
    plot(gg.jabba$gt, gg.jabba$nodes[cluster == 240]$gr+1e5)


    ## ----cache=FALSE,warning=FALSE, fig.height=8------------------------------------------
    ## we first select nodes that are 1 Mbp in width, then compute clusters
    gg.jabba$nodes[width<1e6]$clusters('weak')

    ## note that this syntax still sets the $clusters metadata field of the original 
    ## graph, giving any >1 Mbp nodes a cluster ID of NA
    table(is.na(gg.jabba$nodes$dt$cluster), gg.jabba$nodes$dt$width>1e6)

    ## we peek at one of these interesting clusters, marking it with a blue color
    gg.jabba$nodes[cluster == 111]$mark(col = 'green')

    ## interestingly, we have re-discovered the ERBB2 BFB-driven amplification highlighted above
    gg.jabba$set(height = 30)
    plot(c(gencode, gg.jabba$gt), gg.jabba$nodes[cluster == 111]$gr[, c()]+1e5)


    ## ----cache=FALSE,warning=FALSE, fig.height=10-----------------------------------------
    ## this will populate the ALT edges of the gGraph with metadata fields $ecluster, $ecycle, and $epath
    ## where $ecluster is the concatenation of $ecycle and $epath labels
    gg.jabba$eclusters()

    ## paths are labeled by a "p" prefix, and cycles labeled by a "c" prefix
    ## here we see a multi-junction cluster p52 with 6 edges
    sort(table(gg.jabba$edges$dt$ecluster))

    ## we can mark these edges (and their associated nodes) with a special color
    gg.jabba$edges[ecluster == 41]$mark(col = 'purple')
    gg.jabba$edges[ecluster == 41]$nodes$mark(col = 'purple')

    ## here, the edges and nodes of the cluster that we have discovered
    ## are highlighted in purple 
    plot(c(gg.jabba$gt), gg.jabba$edges[ecluster == 41]$nodes$gr %>% streduce(1e5))


    ## ----cache=FALSE,warning=FALSE, fig.height=10-----------------------------------------

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
    plot(c(gg.jabba$gt, p1$gt), p1$footprint+1e6)


    ## ----cache=FALSE,warning=FALSE, fig.height=10-----------------------------------------
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
    plot(c(gg.jabba$gt, p5$gt), p5$footprint+1e6)



    ## ----cache=FALSE,warning=FALSE, fig.height=10-----------------------------------------
    ## this expression subset our graph to the nodes and edges
    ## that contribute to edge cluster 41 (see previous section on clusters
    ## and communities)
    gg.sub = gg.jabba[, ecluster == 41]

    ## walk decompositiion of this small subgraph
    ## generates 28 possible linear alleles
    walks = gg.sub$walks()

    ## we can choose the longest walk (most nodes traversed)
    walks = walks[which.max(length)]

    ## and plot with the walk track up top (with nodes already marked purple from before)
    plot(c(gg.jabba$gt, walks$gtrack(name = "Longest walk")), walks$footprint+1e4)


    ## ----cache=FALSE,warning=FALSE, fig.height=10, error = TRUE---------------------------
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


    ## ----cache=FALSE,warning=FALSE, fig.height=10-----------------------------------------

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



    ## ----cache=FALSE,warning=FALSE, fig.height=15-----------------------------------------

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



    ## ----cache=FALSE,warning=FALSE, fig.height=15-----------------------------------------
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


    ## ----cache=FALSE,warning=FALSE, fig.height=12-----------------------------------------
    ## revisiting walks traversing ecluster 41
    gg.sub = gg.jabba[, ecluster == 41]

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


    ## ----cache=FALSE,warning=FALSE, fig.height=10-----------------------------------------

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


    ## ----cache=FALSE,warning=FALSE, fig.height=10-----------------------------------------
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



    ## ----cache=FALSE,warning=FALSE, fig.height=10-----------------------------------------
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
    px = proximity(gg.jabba, ser, genes[, 'gene_name'], chunksize = 2000, mc.cores = 1)

    ## peek at the first proximity, we can see the reldist, altdist, refdist
    ## and additional metadata features inherited from the genes object
    px[1]

    ## make a gTrack for the super-enhancers, coloring by tissue
    gt.se = gTrack(se, gr.colorfield = 'tissue', name = 'SupEnh')

    ## plot the first super-enhancer connecting to BRD9
    px[1]$mark(col = 'purple')

    plot(c(gencode, gt.se, px$graph$gt, px[1]$gt), px[1]$footprint+1e5)


    ## ----cache=FALSE,warning=FALSE, fig.height=15-----------------------------------------
    ## use $eval to count ALT junctions for each walk
    px$set(numalt = px$eval(edge = sum(type == 'ALT')))

    ## let's look for a superenhancer connecting to TERT
    this.px = px[numalt>2 & refdist == Inf & gene_name == 'TERT']

    ## check out the first proximity
    this.px[1]

    ## mark it up
    this.px[1]$mark(col = 'purple')

    plot(c(gencode, gt.se, this.px$graph$gt, this.px[1]$gt), this.px[1]$footprint+1e5)



    ## ---- events, warning = FALSE, fig.width = 8, collapse = TRUE, warning = FALSE, results = "markup", message = FALSE----
    ## load the graph for HCC1954
    hcc1954 = gG(jabba = system.file("extdata", "hcc1954", "jabba.rds", package = "gGnome"))

    ## Identify all supported SV event types
    hcc1954 = events(hcc1954, verbose = FALSE)

    ## Summary of identified events
    hcc1954$meta$event[, table(type)]

    ## plot the locus of a BFB event
    plot(hcc1954$gt, hcc1954[bfb>0]$footprint + 1e6); title("BFB in HCC1954")

    ## for the following examples load the CCLE models
    ccle = dir(system.file("extdata", package = "gGnome"), ".+jabba.simple.rds", full = TRUE)
    names(ccle) = gsub(".*gGnome/.*extdata/(.*)\\.jabba\\.simple\\.rds$", "\\1", ccle)


    ## ---- deletions and rigma, warning = FALSE, fig.width = 8, collapse = TRUE, cache = FALSE, warning = FALSE, results = "hide", message = FALSE----
    hcc1954 = del(hcc1954)
    plot(hcc1954$gt, hcc1954$meta$rigma$footprint %>% GRanges %>% streduce(5e5)); title("Rigma in HCC1954")


    ## ---- duplications and pyrgo, warning = FALSE, cache = FALSE, fig.width = 8, collapse = TRUE, cache = FALSE, warning = FALSE, results = "hide", message = FALSE----
    mfe280 = gG(jabba = ccle["MFE_280"])
    mfe280 = dup(mfe280)
    plot(mfe280$gt, mfe280$meta$pyrgo$footprint %>% head(3) %>% GRanges %>% streduce(5e5)); title("Pyrgos in MFE-280")


    ## ---- chromothripsis, warning = FALSE, cache = FALSE, fig.width = 8, collapse = TRUE, cache = FALSE, warning = FALSE, results = "hide", message = FALSE----
    h2081 = gG(jabba = ccle["NCI_H2081"])
    h2081 = chromothripsis(h2081)
    plot(h2081$gt, streduce(h2081$gr %Q% which(chromothripsis>0), 1e6)); title("Chromothripsis in NCI-H2081")


    ## ---- amplicon, warning=FALSE, fig.width = 8, collapse = TRUE, cache = FALSE, warning = FALSE, results = "hide", message = FALSE----
    h526 = gG(jabba = ccle["NCI_H526"])
    h526 = amp(h526)
    plot(h526$gt, streduce(h526$gr %Q% which(tyfonas>0), 1e6))

    hara = gG(jabba = ccle["HARA"])
    hara = amp(hara)
    plot(hara$gt, streduce(hara$gr %Q% which(bfb>0), 1e6))

    hcc827 = gG(jabba = ccle["HCC827"])
    hcc827 = amp(hcc827)
    plot(hcc827$gt, streduce(hcc827$gr %Q% which(dm>0), 1e6))


    ## ---- cp and tic, warning = FALSE, cache = FALSE, fig.width = 8, collapse = TRUE, warning = FALSE, echo = TRUE, results = "hide", message = FALSE----
    h2228 = gG(jabba = ccle["NCI_H2228"])
    h2228 = chromoplexy(h2228)
    plot(h2228$gt, h2228$edges[which(chromoplexy>0)]$shadow %>% streduce(5e6)); title("Chromoplexy in NCI-H2228")

    jhos2 = gG(jabba = ccle["JHOS_2"])
    jhos2 = tic(jhos2)
    plot(jhos2$gt, jhos2$gr %Q% which(tic %in% head(sort(unique(tic)), 3)) %>% streduce(5e4)); title("TICs in JHOS-2")


    ## ---- others, warning = FALSE, cache = FALSE, fig.width = 8, collapse = TRUE, warning = FALSE, results = "hide", message = FALSE----
    hcc1954 = simple(hcc1954)
    plot(hcc1954$gt, hcc1954$edges[grepl("^INV[0-9]+$", simple)]$shadow %>% streduce(1e5)); title("Inversion in HCC1954")
    plot(hcc1954$gt, hcc1954$edges[grepl("^INVDUP[0-9]+$", simple)]$shadow %>% streduce(1e5)); title("Inverted duplication in HCC1954")

    h526 = simple(h526)
    plot(h526$gt, h526$edges[grepl("TRA", simple)]$shadow %>% streduce(1e5))


    ## ----cache=FALSE,warning=FALSE--------------------------------------------------------
    saveRDS(hcc1954, 'hcc1954.rds')

    saveRDS(h526, 'h526.rds')


    ## ----cache=FALSE,warning=FALSE--------------------------------------------------------
    hcc1954.cov.file = system.file('extdata/hcc1954', 'hcc1954.chr8.cov.rds', package="gGnome")


    ## ----cache=FALSE,warning=FALSE--------------------------------------------------------
    jsdt = data.table(sample = c('hcc1954', 'h526'), graph = c('hcc1954.rds', 'h526.rds'), coverage = c(hcc1954.cov.file, ''))


    ## ----cache=FALSE,warning=FALSE--------------------------------------------------------
    gGnome.js(data = jsdt,
                          outdir = './demo_gGnome.js',
                          cov.field = 'reads')


    ## ----echo=FALSE, out.width='100%'-----------------------------------------------------
    knitr::include_graphics(system.file('extdata/tutorial', 'gGnome.js.screenshot.png', package="gGnome"))


    ## ---- json, echo = TRUE, eval = FALSE-------------------------------------------------
    ## gg$json(filename = "gGnome.js/json/[sample_name].json")


    ## ---- cov2csv, echo = TRUE, eval = FALSE----------------------------------------------
    ## cov2csv(hcc1954.cov.file, field = 'reads', output_file = 'demo_gGnome.js/scatterPlot/hcc1954.csv')

})
