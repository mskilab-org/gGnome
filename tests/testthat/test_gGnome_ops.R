
library(gGnome)
library(testthat)
library(gUtils)

options(gGnome.verbose=TRUE)


context('testing gGnome')


message("Toy segments: ", system.file('extdata', 'testing.segs.rds', package="gGnome"))
test_segs = readRDS(system.file('extdata', 'testing.segs.rds', package="gGnome"))

message("Toy edges: ", system.file('extdata', 'testing.es.rds', package="gGnome"))
test_es = readRDS(system.file('extdata', 'testing.es.rds', package="gGnome"))

jab = system.file('extdata', 'jabba.simple.rds', package="gGnome")
message("JaBbA result: ", jab)

prego = system.file('extdata', 'intervalFile.results', package='gGnome')
message("PREGO results: ", prego)

weaver = system.file('extdata', 'weaver', package='gGnome')
message("Weaver results: ", weaver)

svabavcf = system.file('extdata', 'HCC1143.svaba.somatic.sv.vcf', package='gGnome')
message("SvABA results: ", svabavcf)

novobreakvcf = system.file('extdata', 'novoBreak.pass.flt.vcf', package='gGnome')
message("Novobreaks results: ", novobreakvcf)

dellyvcf = system.file('extdata', 'delly.final.vcf.gz', package='gGnome')
message("Delly results: ", dellyvcf)

lumpyvcf = system.file('extdata', 'filter.PE2.SR2.sv.vcf.gz', package='gGnome')
message("Lumpy results: ", lumpyvcf)



jabba = readRDS(jab)

junctions = jabba$junctions
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
    expect_error(seqinfo(ggnew), NA) ## check works
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
    fooregular = gGraph$new(segs = test_segs, es=test_es, regular=FALSE)
    expect_true(is(fooregular, 'gGraph'))
    foopurityploidy = gGraph$new(segs = test_segs, es=test_es, ploidy=3334, purity=233432)
    expect_true(is(foopurityploidy, 'gGraph'))  
    ##
    ##
    added_junctions = foojab$addJuncs(readRDS(jab)$junc)
    expect_true(is(added_junctions, 'gGraph')) 

})



## FUCNTIONS: initialize, set.seqinfo, nullGGraph, simpleGraph, dipGraph, addJuncs, addSegs, karyograph, simplify,
##     decouple, add, jabb2gg, wv2gg, pr2gg, print, plot, window, layout, summary, gg2td, son, html, gg2js, components, 
##     subgraph, filling, 

## ACTIVE BINDINGS: segstats, edges, junctions, G, adj, A, parts, seqinfo, purity, ploidy, td, win, ig

test_that('gGraph works, default', {
    
    options(gGnome.verbose=TRUE)
    ggnew = gGraph$new()
    expect_true(is(ggnew, 'gGraph'))
    ## ACCESS ACTIVE BINDINGS
    expect_equal(length(ggnew$segstats), 0)
    expect_equal(dim(ggnew$edges)[1], 0)
    expect_equal(length(ggnew$junctions), 0)
    expect_error(ggnew$G, NA)  ## check it works; IGRAPH 84fc0c4 D--- 0 0 -- + edges from 84fc0c4:
    expect_equal(length(ggnew$adj), 0)
    expect_equal(length(ggnew$A), 0)
    expect_equal(ggnew$parts, NULL)
    expect_equal(length(ggnew$seqinfo), 0)
    expect_equal(ggnew$purity, NULL)
    expect_equal(ggnew$ploidy, NULL)
    expect_equal(ggnew$td, NULL)
    expect_equal(length(ggnew$win), 0)
    expect_equal(ggnew$ig, NULL)
    ## FUNCTIONS
    ## set.seqinfo = function(genome=NULL, gname=NULL, drop=FALSE)
    ggnew_setseq = ggnew$set.seqinfo()
    expect_true(is(ggnew_setseq, 'gGraph'))
    expect_equal(length(ggnew_setseq$segstats), 0)
    expect_equal(dim(ggnew_setseq$edges)[1], 0)
    expect_equal(length(ggnew_setseq$junctions), 0)
    expect_error(ggnew_setseq$G, NA)  ## check it works
    expect_equal(length(ggnew_setseq$adj), 0)
    expect_equal(length(ggnew_setseq$A), 0)
    expect_equal(ggnew_setseq$parts, NULL)
    expect_equal(length(ggnew_setseq$seqinfo), 0)
    expect_equal(ggnew_setseq$purity, NULL)
    expect_equal(ggnew_setseq$ploidy, NULL)
    expect_equal(ggnew_setseq$td, NULL)
    expect_equal(length(ggnew_setseq$win), 0)
    expect_equal(ggnew_setseq$ig, NULL)
    ## set.seqinfo, drop = TRUE
    ggnew_setseq_drop = ggnew$set.seqinfo(gname = 'foobar', drop=TRUE)
    expect_true(is(ggnew_setseq_drop, 'gGraph'))
    expect_equal(length(ggnew_setseq_drop$segstats), 0)
    expect_equal(dim(ggnew_setseq_drop$edges)[1], 0)
    expect_equal(length(ggnew_setseq_drop$junctions), 0)
    expect_error(ggnew_setseq_drop$G, NA)   ## check it works
    expect_equal(length(ggnew_setseq_drop$adj), 0)
    expect_equal(length(ggnew_setseq_drop$A), 0)
    expect_equal(ggnew_setseq_drop$parts, NULL)
    expect_equal(length(ggnew_setseq_drop$seqinfo), 0)
    expect_equal(ggnew_setseq_drop$purity, NULL)
    expect_equal(ggnew_setseq_drop$ploidy, NULL)
    expect_equal(ggnew_setseq_drop$td, NULL)
    expect_equal(length(ggnew_setseq_drop$win), 0)
    expect_equal(ggnew_setseq_drop$ig, NULL)
    ## set.seqinfo, genome != NULL, gname != NULL
    ggnew_setseq_hg = ggnew$set.seqinfo(genome = hg_seqlengths(), gname = 'foobar', drop = TRUE)
    expect_true(is(ggnew_setseq_hg, 'gGraph'))
    expect_equal(length(ggnew_setseq_hg$segstats), 0)
    expect_equal(dim(ggnew_setseq_hg$edges)[1], 0)
    expect_equal(length(ggnew_setseq_hg$junctions), 0)
    expect_error(ggnew_setseq_hg$G, NA) ## check it works
    expect_equal(length(ggnew_setseq_hg$adj), 0)
    expect_equal(length(ggnew_setseq_hg$A), 0)
    expect_equal(ggnew_setseq_hg$parts, NULL)
    expect_equal(length(ggnew_setseq_hg$seqinfo), 25)   ### checks!
    expect_equal(ggnew_setseq_hg$purity, NULL)
    expect_equal(ggnew_setseq_hg$ploidy, NULL)
    expect_equal(ggnew_setseq_hg$td, NULL)
    expect_equal(length(ggnew_setseq_hg$win), 0)
    expect_equal(ggnew_setseq_hg$ig, NULL)
    ##
    ## nullGraph = function(regular=TRUE, genome=NULL)
    ggnew_setseq_nullGraph = ggnew$nullGGraph()
    expect_true(is(ggnew_setseq_nullGraph, 'gGraph'))
    expect_equal(length(ggnew_setseq_nullGraph$segstats), 0)
    expect_equal(dim(ggnew_setseq_nullGraph$edges)[1], 0)
    expect_equal(length(ggnew_setseq_nullGraph$junctions), 0)
    expect_error(ggnew_setseq_nullGraph$G, NA) ## check it works
    expect_equal(length(ggnew_setseq_nullGraph$adj), 0)
    expect_equal(length(ggnew_setseq_nullGraph$A), 0)
    expect_equal(ggnew_setseq_nullGraph$parts, NULL)
    expect_equal(length(ggnew_setseq_nullGraph$seqinfo), 25)   ### checks! "null" means there is no node, you can still have a "space" of possible values when the set is empty
    expect_equal(ggnew_setseq_nullGraph$purity, NULL)
    expect_equal(ggnew_setseq_nullGraph$ploidy, NULL)
    expect_equal(ggnew_setseq_nullGraph$td, NULL)
    expect_equal(length(ggnew_setseq_nullGraph$win), 0)
    expect_equal(ggnew_setseq_nullGraph$ig, NULL)
    ##
    ## simpleGraph = function(genome = NULL, chr=FALSE, include.junk=FALSE, ploidy = NULL)
    ggnew_setseq_simpleGraph = ggnew$simpleGraph()
    expect_true(is(ggnew_setseq_simpleGraph, 'gGraph'))
    expect_equal(length(ggnew_setseq_simpleGraph$segstats), 50)
    expect_equal(dim(ggnew_setseq_simpleGraph$edges)[1], 0)
    expect_equal(length(ggnew_setseq_simpleGraph$junctions), 0)
    expect_error(ggnew_setseq_simpleGraph$G, NA) ## check it works
    expect_equal(length(ggnew_setseq_simpleGraph$adj), 2500)
    expect_equal(length(ggnew_setseq_simpleGraph$A), 2500)
    ## expect_equal(ggnew_setseq_simpleGraph$parts, NULL)
    expect_equal(length(ggnew_setseq_simpleGraph$seqinfo), 25)   ### checks! "null" means there is no node, you can still have a "space" of possible values when the set is empty
    expect_equal(ggnew_setseq_simpleGraph$purity, NULL)
    expect_equal(ggnew_setseq_simpleGraph$ploidy, NULL)
    expect_true(is(ggnew_setseq_simpleGraph$td, 'gTrack'))
    expect_equal((ggnew_setseq_simpleGraph$td)$ygap, 2)
    expect_match((ggnew_setseq_simpleGraph$td)$name, 'CN')
    expect_equal(length(ggnew_setseq_simpleGraph$win), 25)
    ## STILL ERROR
    ## >  ggnew_setseq_simpleGraph$ig
    ##Error in log(private$segs$cn, 1.4) : 
    ##  non-numeric argument to mathematical function
    ##In addition: Warning messages:
    ##1: replacing previous import ‘VariantAnnotation::select’ by ‘plotly::select’ when loading ‘skitools’ 
    ##2: replacing previous import ‘ggplot2::last_plot’ by ‘plotly::last_plot’ when loading ‘skitools’ 
    ## dipGraph = function(genome = NULL, chr=FALSE, include.junk=FALSE)
    ggnew_dd = ggnew$dipGraph()
    expect_true(is(ggnew_dd, 'gGraph'))
    expect_equal(length(ggnew_dd$segstats), 50)
    expect_equal(dim(ggnew_dd$edges)[1], 0)
    expect_equal(length(ggnew_dd$junctions), 0)
    expect_error(ggnew_dd$G, NA) ## check it works
    expect_equal(length(ggnew_dd$adj), 2500)
    expect_equal(length(ggnew_dd$A), 2500)
    ## > ggnew_dd$ig
    ## > ggnew_dd$parts
    expect_equal(length(ggnew_dd$seqinfo), 25)   
    expect_equal(ggnew_dd$purity, NULL)
    expect_equal(ggnew_dd$ploidy, 2)  ## checks!
    expect_true(is(ggnew_dd$td, 'gTrack'))
    expect_equal((ggnew_dd$td)$ygap, 2)
    expect_match((ggnew_dd$td)$name, 'CN')
    expect_equal(length(ggnew_dd$win), 25)
    ##
    ## dipGraph
    ggnew_dd_junk = ggnew$dipGraph(genome = hg_seqlengths(), chr=TRUE, include.junk=TRUE)
    expect_true(is(ggnew_dd_junk, 'gGraph'))
    expect_equal(length(ggnew_dd_junk$segstats), 50)
    expect_equal(dim(ggnew_dd_junk$edges)[1], 0)
    expect_equal(length(ggnew_dd_junk$junctions), 0)
    expect_error(ggnew_dd_junk$G, NA) ## check it works
    expect_equal(length(ggnew_dd_junk$adj), 2500)
    expect_equal(length(ggnew_dd_junk$A), 2500)
    ## expect_equal(ggnew_dd$parts, NULL)
    expect_equal(length(ggnew_dd_junk$seqinfo), 25)   
    expect_equal(ggnew_dd_junk$purity, NULL)
    expect_equal(ggnew_dd_junk$ploidy, 2)  ## checks!
    expect_true(is(ggnew_dd_junk$td, 'gTrack'))
    expect_equal((ggnew_dd_junk$td)$ygap, 2)
    expect_match((ggnew_dd_junk$td)$name, 'CN')
    expect_equal(length(ggnew_dd_junk$win), 25)
    ##
    ## addJuncs = function(junc, cn=TRUE)
    expect_error(ggnew$addJuncs())
    added_juncs = ggnew$addJuncs(junc = junctions)
    expect_true(is(added_juncs, 'gGraph'))
    expect_equal(length(added_juncs$segstats), 1112)
    expect_equal(dim(added_juncs$edges)[1], 1596)
    expect_equal(dim(added_juncs$edges)[2], 15)
    expect_equal(length(added_juncs$junctions), 267)
    expect_error(added_juncs$G, NA) ## check it works
    expect_equal(length(added_juncs$adj), 1236544)
    expect_equal(length(added_juncs$A),  1236544)
    ## expect_equal(added_juncs$parts, NULL)
    expect_equal(length(added_juncs$seqinfo), 25)   
    expect_equal(added_juncs$purity, NULL)
    expect_equal(added_juncs$ploidy, NULL)  ## checks!
    expect_true(is(added_juncs$td, 'gTrack'))
    expect_equal((added_juncs$td)$ygap, 2)
    expect_match((added_juncs$td)$name, 'CN')
    expect_equal(length(added_juncs$win), 25)
    ## addJuncs, cn = FALSE
    ## added_juncs_cnfalse = ggnew$addJuncs(junc = junctions, cn = FALSE) ERROR
    ##
    ## addSegs = function(tile)
    added_segs = ggnew$addSegs(tile = test_segs)
    expect_true(is(added_segs, 'gGraph'))
    expect_equal(length(added_segs$segstats), 1124)
    expect_equal(dim(added_segs$edges)[1], 1074)
    expect_equal(dim(added_segs$edges)[2], 14)
    expect_equal(length(added_segs$junctions), 0)
    expect_error(added_segs$G, NA) ## check it works
    expect_equal(length(added_segs$adj), 1263376)
    expect_equal(length(added_segs$A),  1263376)
    ## expect_equal(added_segs$parts, NULL)
    expect_equal(length(added_segs$seqinfo), 25)   
    expect_equal(added_segs$purity, NULL)
    expect_equal(added_segs$ploidy, NULL)  ## checks!
    expect_true(is(added_segs$td, 'gTrack'))
    expect_equal((added_segs$td)$ygap, 2)
    expect_match((added_segs$td)$name, 'CN')
    expect_equal(length(added_segs$win), 25)
    ## karyograph = function(tile=NULL, juncs=NULL, cn=FALSE, regular=FALSE)
    ## default
    ggnew_karyograph = ggnew$karyograph()
    expect_true(is(ggnew_karyograph, 'gGraph'))
    expect_equal(length(ggnew_karyograph$segstats), 1124)
    expect_equal(dim(ggnew_karyograph$edges)[1], 1074)
    expect_equal(dim(ggnew_karyograph$edges)[2], 23)
    expect_equal(length(ggnew_karyograph$junctions), 0)
    expect_error(ggnew_karyograph$G, NA) ## check it works
    expect_equal(length(ggnew_karyograph$adj), 1263376)
    expect_equal(length(ggnew_karyograph$A),  1263376)
    ## expect_equal(length((ggnew_karyograph$parts)$membership), 50)
    ##expect_equal(length((ggnew_karyograph$parts)$csize), 25)
    ##expect_equal(length(ggnew_karyograph$seqinfo), 25)   
    expect_equal(ggnew_karyograph$purity, NULL)
    expect_equal(ggnew_karyograph$ploidy, NULL)  ## checks!
    expect_true(is(ggnew_karyograph$td, 'gTrack'))
    expect_equal((ggnew_karyograph$td)$ygap, 2)
    expect_match((ggnew_karyograph$td)$name, 'CN')
    expect_equal(length(ggnew_karyograph$win), 25)
    ##
    ## simplify = function(mod=TRUE)
    ggnew_simplify = ggnew$simplify()
    expect_true(is(ggnew_simplify, 'gGraph'))
    expect_equal(length(ggnew_simplify$segstats), 50)
    expect_equal(dim(ggnew_simplify$edges)[1], 0)
    expect_equal(dim(ggnew_simplify$edges)[2], 11)
    expect_equal(length(ggnew_simplify$junctions), 0)
    expect_error(ggnew_simplify$G, NA) ## check it works
    expect_equal(length(ggnew_simplify$adj), 2500)
    expect_equal(length(ggnew_simplify$A),  2500)
    ##expect_equal(length((ggnew_simplify$parts)$membership), 50)
    ##expect_equal(length((ggnew_simplify$parts)$csize), 25)
    expect_equal(length(ggnew_simplify$seqinfo), 25)   
    expect_equal(ggnew_simplify$purity, NULL)
    expect_equal(ggnew_simplify$ploidy, NULL)  ## checks!
    expect_true(is(ggnew_simplify$td, 'gTrack'))
    expect_equal((ggnew_simplify$td)$ygap, 2)
    expect_match((ggnew_simplify$td)$name, 'CN')
    expect_equal(length(ggnew_simplify$win), 25)
    ## simplify = function(mod=FALSE)
    ggnew_simplify = ggnew$simplify(mod=FALSE)
    expect_true(is(ggnew_simplify, 'gGraph'))
    expect_equal(length(ggnew_simplify$segstats), 50)
    expect_equal(dim(ggnew_simplify$edges)[1], 0)
    expect_equal(dim(ggnew_simplify$edges)[2], 11)
    expect_equal(length(ggnew_simplify$junctions), 0)
    expect_error(ggnew_simplify$G, NA) ## check it works
    expect_equal(length(ggnew_simplify$adj), 2500)
    expect_equal(length(ggnew_simplify$A),  2500)
    ##expect_equal(length((ggnew_simplify$parts)$membership), 50)
    ##expect_equal(length((ggnew_simplify$parts)$csize), 25)
    expect_equal(length(ggnew_simplify$seqinfo), 25)   
    expect_equal(ggnew_simplify$purity, NULL)
    expect_equal(ggnew_simplify$ploidy, NULL)  ## checks!
    expect_true(is(ggnew_simplify$td, 'gTrack'))
    expect_equal((ggnew_simplify$td)$ygap, 2)
    expect_match((ggnew_simplify$td)$name, 'CN')
    expect_equal(length(ggnew_simplify$win), 25)
    ## decouple = function(mod=TRUE)
    ggnew_decouple = ggnew$decouple()
    expect_true(is(ggnew_decouple, 'gGraph'))
    expect_equal(length(ggnew_decouple$segstats), 50)
    expect_equal(dim(ggnew_decouple$edges)[1], 0)
    expect_equal(dim(ggnew_decouple$edges)[2], 19)
    expect_equal(length(ggnew_decouple$junctions), 0)
    expect_error(ggnew_decouple$G, NA) ## check it works
    expect_equal(length(ggnew_decouple$adj), 2500)
    expect_equal(length(ggnew_decouple$A),  2500)
    ##expect_equal(length((ggnew_decouple$parts)$membership), 50)
    ##expect_equal(length((ggnew_decouple$parts)$csize), 25)
    expect_equal(length(ggnew_decouple$seqinfo), 25)   
    expect_equal(ggnew_decouple$purity, NULL)
    expect_equal(ggnew_decouple$ploidy, NULL)  ## checks!
    expect_true(is(ggnew_decouple$td, 'gTrack'))
    expect_equal((ggnew_decouple$td)$ygap, 2)
    expect_match((ggnew_decouple$td)$name, 'CN')
    expect_equal(length(ggnew_decouple$win), 25)
    ##
    ggnew_decouple = ggnew$decouple(mod=FALSE)
    expect_true(is(ggnew_decouple, 'gGraph'))
    expect_equal(length(ggnew_decouple$segstats), 50)
    expect_equal(dim(ggnew_decouple$edges)[1], 0)
    expect_equal(dim(ggnew_decouple$edges)[2], 19)
    expect_equal(length(ggnew_decouple$junctions), 0)
    expect_error(ggnew_decouple$G, NA) ## check it works
    expect_equal(length(ggnew_decouple$adj), 2500)
    expect_equal(length(ggnew_decouple$A),  2500)
    ##expect_equal(length((ggnew_decouple$parts)$membership), 50)
    ##expect_equal(length((ggnew_decouple$parts)$csize), 25)
    expect_equal(length(ggnew_decouple$seqinfo), 25)   
    expect_equal(ggnew_decouple$purity, NULL)
    expect_equal(ggnew_decouple$ploidy, NULL)  ## checks!
    expect_true(is(ggnew_decouple$td, 'gTrack'))
    expect_equal((ggnew_decouple$td)$ygap, 2)
    expect_match((ggnew_decouple$td)$name, 'CN')
    expect_equal(length(ggnew_decouple$win), 25)
    ## 
    ## add = function(gg, mod=FALSE)
    added = ggnew$add(gg=ggnew_karyograph, mod=TRUE)
    expect_true(is(added, 'gGraph'))
    expect_equal(length(added$segstats), 100)
    expect_equal(dim(added$edges)[1], 0)
    expect_equal(dim(added$edges)[2], 10)
    expect_equal(length(added$junctions), 0)
    expect_error(added$G, NA) ## check it works
    expect_equal(length(added$adj), 10000)
    expect_equal(length(added$A),  10000)
    ##expect_true(is(added$parts, 'list'))
    ##expect_equal(length(added$parts), 3)
    ##expect_equal(length(added$parts$membership), 100)
    ##expect_equal(length((added$parts)$csize), 25)
    expect_equal(length(added$seqinfo), 25)   
    expect_equal(added$purity, NULL)
    expect_equal(added$ploidy, NULL)  ## checks!
    expect_true(is(added$td, 'gTrack'))
    expect_equal((added$td)$ygap, 2)
    expect_match((added$td)$name, 'CN')
    expect_equal(length(added$win), 25)
    ## add() several times
    five_adds = ggnew$add(gg=ggnew_karyograph)$add(gg=ggnew_karyograph)$add(gg=ggnew_karyograph)$add(gg=ggnew_karyograph)$add(gg=ggnew_karyograph)
    expect_true(is(five_adds, 'gGraph'))
    expect_equal(length(five_adds$segstats), 600)
    expect_equal(dim(five_adds$edges)[1], 0)
    expect_equal(dim(five_adds$edges)[2], 10)
    expect_equal(length(five_adds$junctions), 0)
    expect_error(five_adds$G, NA) ## check it works
    expect_equal(length(five_adds$adj), 360000)
    expect_equal(length(five_adds$A), 360000)
    ##expect_true(is(five_adds$parts, 'list'))
    ##expect_equal(length(five_adds$parts), 3)
    ##expect_equal(length(five_adds$parts$membership), 300)
    ##expect_equal(length((five_adds$parts)$csize), 25)
    expect_equal(length(five_adds$seqinfo), 25)   
    expect_equal(five_adds$purity, NULL)
    expect_equal(five_adds$ploidy, NULL)  ## checks!
    expect_true(is(five_adds$td, 'gTrack'))
    expect_equal((five_adds$td)$ygap, 2)
    expect_match((five_adds$td)$name, 'CN')
    expect_equal(length(five_adds$win), 25) 
    ##
    ## jabb2gg()
    jabbd = ggnew$jab2gg(jabba)
    expect_true(is(jabbd, 'gGraph'))
    expect_equal(length(jabbd$segstats), 2346)
    expect_equal(dim(jabbd$edges)[1], 2714)
    expect_equal(dim(jabbd$edges)[2], 12)
    expect_equal(length(jabbd$junctions), 0)
    expect_error(jabbd$G, NA) ## check it works
    ## ERROR expect_equal(length(jabbd$adj),  90000)
    ## ERROR expect_equal(length(jabbd$A),   90000)
    ##expect_true(is(jabbd$parts, 'list'))
    ##expect_equal(length(jabbd$parts), 3)
    ##expect_equal(length(jabbd$parts$membership), 300)
    ##expect_equal(length((jabbd$parts)$csize), 25)
    expect_equal(length(jabbd$seqinfo), 85)   
    expect_equal(jabbd$purity, 0.96)
    expect_equal(round(jabbd$ploidy, 2), 3.85)  
    expect_true(is(jabbd$td, 'gTrack'))
    expect_equal((jabbd$td)$ygap, 2)
    expect_match((jabbd$td)$name, 'CN')
    expect_equal(length(jabbd$win), 85) 
    ##
    ## seqlengths(jabba$segstats) = NA 
    ## should trip 'if (all(is.na(s.sl <- seqlengths(jabba$segstats)))){'
    jabba_seqna = jabba
    seqlengths(jabba_seqna$segstats) = NA 
    jabbd_seqna = ggnew$jab2gg(jabba_seqna)
    expect_true(is(jabbd_seqna, 'gGraph'))
    expect_equal(length(jabbd_seqna$segstats), 2346)
    expect_equal(dim(jabbd_seqna$edges)[1], 2714)
    expect_equal(dim(jabbd_seqna$edges)[2], 12)
    expect_equal(length(jabbd_seqna$junctions), 0)
    expect_error(jabbd_seqna$G, NA) ## check it works
    ## ERROR expect_equal(length(jabbd$adj),  90000)
    ## ERROR expect_equal(length(jabbd$A),   90000)
    ##expect_true(is(jabbd$parts, 'list'))
    ##expect_equal(length(jabbd$parts), 3)
    ##expect_equal(length(jabbd$parts$membership), 300)
    ##expect_equal(length((jabbd$parts)$csize), 25)
    expect_equal(length(jabbd_seqna$seqinfo), 85)   
    expect_equal(jabbd_seqna$purity, 0.96)
    expect_equal(round(jabbd_seqna$ploidy, 2), 3.85)  
    expect_true(is(jabbd_seqna$td, 'gTrack'))
    expect_equal((jabbd_seqna$td)$ygap, 2)
    expect_match((jabbd_seqna$td)$name, 'CN')
    expect_equal(length(jabbd_seqna$win), 85) 
    ##
    ## wv2gg()
    weavd = ggnew$wv2gg(weaver)
    expect_true(is(weavd, 'gGraph'))
    expect_equal(length(weavd$segstats), 16184)
    expect_equal(dim(weavd$edges)[1], 16974)
    expect_equal(dim(weavd$edges)[2], 15)
    expect_equal(length(weavd$junctions), 420)
    expect_error(weavd$G, NA) ## check it works
    ## ERROR expect_equal(length(weavd$adj),  90000)
    ## ERROR expect_equal(length(weavd$A),   90000)
    ##expect_true(is(weavd$parts, 'list'))
    ##expect_equal(length(weavd$parts), 3)
    ##expect_equal(length(weavd$parts$membership), 300)
    ##expect_equal(length((weavd$parts)$csize), 25)
    expect_equal(length(weavd$seqinfo), 25)   
    expect_equal(weavd$purity, NULL)
    expect_equal(round(weavd$ploidy, 2), 8.54)  
    expect_true(is(weavd$td, 'gTrack'))
    expect_equal((weavd$td)$ygap, 2)
    expect_match((weavd$td)$name, 'CN')
    expect_equal(length(weavd$win), 25) 
    ## pr2gg()
    pregod = ggnew$pr2gg(prego)
    expect_true(is(pregod, 'gGraph'))
    expect_equal(length(pregod$segstats), 1208)
    expect_equal(dim(pregod$edges)[1], 1380)
    expect_equal(dim(pregod$edges)[2], 21)
    expect_equal(length(pregod$junctions), 420)
    expect_error(pregod$G, NA) ## check it works
    ## ERROR expect_equal(length(weavd$adj),  90000)
    ## ERROR expect_equal(length(weavd$A),   90000)
    ##expect_true(is(weavd$parts, 'list'))
    ##expect_equal(length(weavd$parts), 3)
    ##expect_equal(length(weavd$parts$membership), 300)
    ##expect_equal(length((weavd$parts)$csize), 25)
    expect_equal(length(pregod$seqinfo), 84)   
    expect_equal(pregod$purity, 1)
    expect_equal(round(pregod$ploidy, 2), 2.1)  
    expect_true(is(pregod$td, 'gTrack'))
    expect_equal((pregod$td)$ygap, 2)
    expect_match((pregod$td)$name, 'CN')
    expect_equal(length(pregod$win), 24) 
    ##
    ## print() to STDOUT
    expect_error(ggnew$print(), NA)
    ## plot() nothing is returned here, so let's do this for now:
    expect_error(ggnew$plot(), NA)
    expect_error(ggnew$plot(colorful=TRUE), NA)
    ## window()
    expect_equal(length(ggnew$window()), 24)
    ## 
    ## > ggnew$layout()
    ## Error in log(private$segs$cn, 1.4) : 
    ## summary()
    expect_true(is.character(ggnew$summary()))
    ## length()
    expect_equal(ggnew$length(), NULL)
    ##
    expect_true(is(ggnew$gg2td(), 'gTrack'))
    ## JSON
    ## > ggnew$json()
    expect_error(ggnew$json(), NA)
    ## HTML
    ## > ggnew$html()
    expect_error(ggnew$html()) ## Error in ggnew$html() : Get from https://github.com/mskilab/gGnome.js
    ## No gGnome.js repository found on your system.
    ## Error in ggnew$html() : Get from https://github.com/mskilab/gGnome.js
    ## gg2j()
    ## > ggnew$gg2js()
    expect_error(ggnew$gg2js(), NA)
    ##
    ## component()
    component = ggnew$components(mc.cores=2)
    expect_true(is(component, 'list'))
    expect_equal(length(component), 10)
    expect_equal(length(component$segstats), 0)
    expect_equal(dim(component$edges)[1], NULL)
    expect_equal(dim(component$edges)[2], NULL)
    expect_equal(length(component$junctions), 0)
    expect_equal(component$G, NULL) 
    expect_equal(length(component$adj), 0)
    expect_equal(length(component$A),  0)
    expect_equal(length(component$parts), 0)
    expect_equal(length(component$seqinfo), 0)   
    expect_equal(component$purity, NULL)
    expect_equal(component$ploidy, NULL)  ## checks!
    expect_equal(component$td, NULL)
    expect_equal(length(component$win), 0) 
    ##
    ## subgraph = function(v=numeric(0), na.rm=T, mod=T)
    ## default 
    subgraphed = ggnew$subgraph()
    expect_true(is(subgraphed, 'gGraph'))
    expect_equal(length(subgraphed), 10)
    expect_equal(length(subgraphed$segstats), 1208)
    expect_equal(dim(subgraphed$edges)[1], 1380)
    expect_equal(dim(subgraphed$edges)[2], 29)
    expect_equal(length(subgraphed$junctions), 420)
    expect_error(subgraphed$G, NA) 
    expect_equal(length(subgraphed$adj), 1459264)
    expect_equal(length(subgraphed$A),  1459264)
    ##expect_equal(length(subgraphed$parts), 0)
    expect_equal(length(subgraphed$seqinfo), 84)   
    expect_equal(subgraphed$purity, 1)
    expect_equal(round(subgraphed$ploidy, 2), 2.1)  ## checks!
    expect_true(is(subgraphed$td, 'gTrack'))
    expect_equal(length(subgraphed$win), 24) 
    ##
    ## vertices5K = ggnew$subgraph(v=5000)
    ##
    filled = ggnew$fillin()
    expect_true(is(filled, 'gGraph'))
    expect_equal(length(filled$segstats), 1208)
    expect_equal(dim(filled$edges)[1], 1380)
    expect_equal(dim(filled$edges)[2],  29)
    expect_equal(length(filled$junctions), 420)
    expect_error(filled$G, NA) 
    expect_equal(length(filled$adj), 1459264)
    expect_equal(length(filled$A),  1459264)
    ##expect_equal(length(filled$parts), 3)
    expect_equal(length(filled$seqinfo), 84)   
    expect_equal(filled$purity, 1)
    expect_equal(round(filled$ploidy, 2), 2.1)  ## checks!
    expect_true(is(filled$td, 'gTrack'))
    expect_equal(length(filled$win), 24)
    ## 
    ## trim = function(gr=NULL, mod=FALSE)
    ## default
    trimmed = ggnew$trim()
    expect_true(is(trimmed, 'gGraph'))
    expect_equal(length(trimmed$segstats), 1208)
    expect_equal(dim(trimmed$edges)[1], 1380)
    expect_equal(dim(trimmed$edges)[2], 29)
    expect_equal(length(trimmed$junctions), 420)
    expect_error(trimmed$G, NA) 
    expect_equal(length(trimmed$adj), 1459264)
    expect_equal(length(trimmed$A),  1459264)
    ##expect_equal(length(trimmed$parts), 3)
    expect_equal(length(trimmed$seqinfo), 84)    
    expect_equal(trimmed$purity,  1)
    expect_equal(round(trimmed$ploidy, 2), 2.1)  ## checks!
    expect_true(is(trimmed$td, 'gTrack'))
    expect_equal(length(trimmed$win), 24)
    ## 
    ## trimmed_mod = ggnew$trim(gr=gr2, mod=TRUE) 
    ## 
    ## gotg
    ## default
    gotg = ggnew$get.g() 
    expect_true(is(gotg, 'gGraph'))
    expect_equal(length(gotg$segstats), 1208)
    expect_equal(dim(gotg$edges)[1], 1380)
    expect_equal(dim(gotg$edges)[2], 29)
    expect_equal(length(gotg$junctions), 420)
    expect_error(gotg$G, NA) 
    expect_equal(length(gotg$adj), 1459264)
    expect_equal(length(gotg$A), 1459264)
    ##expect_equal(length(gotg$parts), 3)
    expect_equal(length(gotg$seqinfo), 84)   
    expect_equal(gotg$purity, 1)
    expect_equal(round(gotg$ploidy, 2), 2.1)   ## checks!
    expect_true(is(gotg$td, 'gTrack'))
    expect_equal(length(gotg$win), 24)
    ## get.g = function(force=FALSE)
    gotg_forced = ggnew$get.g(force=TRUE) 
    expect_true(is(gotg_forced, 'gGraph'))
    expect_equal(length(gotg_forced$segstats), 1208)
    expect_equal(dim(gotg_forced$edges)[1], 1380)
    expect_equal(dim(gotg_forced$edges)[2], 29)
    expect_equal(length(gotg_forced$junctions), 420)
    expect_error(gotg_forced$G, NA) 
    expect_equal(length(gotg_forced$adj), 1459264)
    expect_equal(length(gotg_forced$A), 1459264)
    ##expect_equal(length(gotg_forced$parts), 3)
    expect_equal(length(gotg_forced$seqinfo), 84)   
    expect_equal(gotg_forced$purity, 1)
    expect_equal(round(gotg_forced$ploidy, 2), 2.1)   ## checks!
    expect_true(is(gotg_forced$td, 'gTrack'))
    expect_equal(length(gotg_forced$win), 24)
    ## with added
    gotgadd = added$get.g(force=TRUE) 
    expect_true(is(gotgadd, 'gGraph'))
    expect_equal(length(gotgadd$segstats), 1208)
    expect_equal(dim(gotgadd$edges)[1], 1380)
    expect_equal(dim(gotgadd$edges)[2], 29)
    expect_equal(length(gotgadd$junctions), 420)
    expect_error(gotgadd$G, NA) 
    expect_equal(length(gotgadd$adj), 1459264)
    expect_equal(length(gotgadd$A),  1459264)
    ##expect_equal(length(gotgadd$parts), 3)  
    ##expect_equal(length((gotgadd$parts)$membership), 1208)
    ##expect_equal(length((gotgadd$parts)$csize), 25)
    ##expect_equal(length((gotgadd$parts)$no), 1)  ## 25
    expect_equal(length(gotgadd$seqinfo), 84)   
    expect_equal(gotgadd$purity, 1)
    expect_equal(round(gotgadd$ploidy, 2), 2.1) 
    expect_true(is(gotgadd$td, 'gTrack'))
    expect_equal((gotgadd$td)$ygap, 2)
    expect_match((gotgadd$td)$name, 'CN')
    expect_equal(length(gotgadd$win), 24) 
    ##  hood = function(win, d=NULL, k=NULL, pad=0, bagel=FALSE, ignore.strand=T, verbose=FALSE)  
    ##gr2_win = ggnew$hood(win=gr2)
    expect_error(ggnew$hood(win=grl1)) ## Error in .local(x, y, ...) : setdiff() between a GRanges and a GRangesList object is not supported
    ##grl2_win = ggnew$hood(win=grl.unlist(grl2))
    ##hgseq_win = ggnew$hood(win = si2gr(hg_seqlengths()))
    ##expect_true(is(hgseq_win, 'gGraph'))
    ##expect_equal(length(hgseq_win$segstats), 0)
    ##expect_equal(dim(hgseq_win$edges)[1], 0)
    ##expect_equal(dim(hgseq_win$edges)[2], 19)
    ##expect_equal(length(hgseq_win$junctions), 0)
    ##expect_error(hgseq_win$G, NA) 
    ##expect_equal(length(hgseq_win$adj), 0)
    ##expect_equal(length(hgseq_win$A),  0)
    ##expect_equal(length(hgseq_win$parts), 0)   
    ##expect_equal(hgseq_win$purity, NULL)
    ##expect_equal(hgseq_win$ploidy, NULL)  ## checks!
    ##expect_equal(length(hgseq_win$td), 0)
    ##expect_equal(length(hgseq_win$win), 0) 
    ##
    ## dist
    distanced1 = ggnew$dist(GRanges('1:5500-6000'), GRanges('1:5000-5500'))
    expect_equal(as.numeric(distanced1), 0)
    ## ERROR distanced2 = ggnew$dist(GRanges('1:5500-6000'), GRanges('1:15000-15500'))
    distanced_diffchroms = ggnew$dist(GRanges('2:5500-6000'), GRanges('3:5000-5500'))
    expect_equal(as.numeric(distanced_diffchroms), Inf)
    ##
    ## fillup
    filledup = ggnew$fillup()
    expect_true(is(filledup, 'gGraph'))
    expect_equal(length(filledup$segstats), 1208)
    expect_equal(dim(filledup$edges)[1], 1380)
    expect_equal(dim(filledup$edges)[2], 29)
    expect_equal(length(filledup$junctions), 420)
    expect_error(filledup$G, NA) 
    expect_equal(length(filledup$adj), 1459264)
    expect_equal(length(filledup$A),  1459264)
    ## expect_equal(length(filledup$parts), 3)  
    ##expect_equal(length((filledup$parts)$membership), 100)
    ##expect_equal(length((filledup$parts)$csize), 25)
    ##expect_equal(length((filledup$parts)$no), 1)  ## 25
    expect_equal(length(filledup$seqinfo), 84)   
    expect_equal(filledup$purity, 1)
    expect_equal(round(filledup$ploidy, 2), 2.1)  ## checks!
    expect_true(is(filledup$td, 'gTrack'))
    expect_equal((filledup$td)$ygap, 2)
    expect_match((filledup$td)$name, 'CN')
    expect_equal(length(filledup$win), 24) 
    ## isBalance
    expect_true(ggnew$isBalance())
    ## get.loose
    expect_equal(length(ggnew$get.loose()), 0)

})


### try gGraph with inputs above
### gGraph$new(tile = test_segs)


## FUCNTIONS: initialize, set.seqinfo, nullGGraph, simpleGraph, dipGraph, addJuncs, addSegs, karyograph, simplify,
##     decouple, add, jabb2gg, wv2gg, pr2gg, print, plot, window, layout, summary, gg2td, son, html, gg2js, components, 
##     subgraph, filling, 


## segstats, edges, grl, td, path, values


test_that('check gGraph w/ inputs works', {

    options(gGnome.verbose=TRUE)
    gg = gGraph$new(segs = test_segs, es=test_es, cn=TRUE)
    expect_true(is(gg, 'gGraph'))
    expect_equal(length(gg$segstats), 10)
    expect_equal(dim(gg$edges)[1], 12)
    expect_equal(length(gg$junctions), 2)
    expect_error(gg$G, NA)  ## IGRAPH 8d66213 D--- 10 12 --
    expect_equal(length(gg$adj), 100)
    expect_equal(length(gg$A), 100)   
    expect_equal(length((gg$parts)$membership), 10)
    expect_equal(length(gg$seqinfo), 93)
    expect_equal(gg$purity, NULL)
    expect_equal(gg$ploidy, 3)
    expect_true(is(gg$td, 'gTrack'))
    expect_equal(width(gg$win), 1000000)
    ## expect_equal(gg$ig, NULL)
    ## 
    ## set.seqinfo
    setseq = gg$set.seqinfo()
    expect_true(is(setseq, 'gGraph'))
    expect_equal(length(setseq$segstats), 10)
    expect_equal(dim(setseq$edges)[1], 12)
    expect_equal(dim(setseq$edges)[2], 29)
    expect_equal(length(setseq$junctions), 2)
    expect_error(setseq$G, NA)  ## IGRAPH 8d66213 D--- 10 12 --
    expect_equal(length(setseq$adj), 100)
    expect_equal(length(setseq$A), 100)   
    expect_equal(length((setseq$parts)$membership), 10)
    expect_equal(length(setseq$seqinfo), 93)
    expect_equal(setseq$purity, NULL)
    expect_equal(setseq$ploidy, 3)
    expect_true(is(setseq$td, 'gTrack'))
    expect_equal(width(setseq$win), 1000000)
    ## expect_equal(setseq$ig, NULL)   
    ## gg_setseq_hg = gg$set.seqinfo(genome = hg_seqlengths(), gname = 'foobar', drop = TRUE)
    gg_setseq_hg = gg$set.seqinfo(genome = hg_seqlengths(), gname = 'foobar', drop = TRUE)
    expect_true(is(gg_setseq_hg, 'gGraph'))
    expect_equal(length(gg_setseq_hg$segstats), 10)
    expect_equal(dim(gg_setseq_hg$edges)[1], 12)
    expect_equal(dim(gg_setseq_hg$edges)[2], 29)
    expect_equal(length(gg_setseq_hg$junctions), 2)
    expect_error(gg_setseq_hg$G, NA)  ## IGRAPH 8d66213 D--- 10 12 --
    expect_equal(length(gg_setseq_hg$adj), 100)
    expect_equal(length(gg_setseq_hg$A), 100)   
    expect_equal(length((gg_setseq_hg$parts)$membership), 10)
    expect_equal(length(gg_setseq_hg$seqinfo), 25)   ### CHANGE happened here
    expect_equal(gg_setseq_hg$purity, NULL)
    expect_equal(gg_setseq_hg$ploidy, 3)
    expect_true(is(gg_setseq_hg$td, 'gTrack'))
    expect_equal(width(gg_setseq_hg$win), 1000000)
    ###expect_equal(gg_setseq_hg$ig, NULL)   

})




test_that('gGraph works, JaBbA input', {
    
    options(gGnome.verbose=TRUE)
    ggjab = gGraph$new(jabba=jabba, cn=TRUE)
    expect_true(is(ggjab, 'gGraph'))
    ## ACCESS ACTIVE BINDINGS
    expect_equal(length(ggjab$segstats), 2346)
    expect_equal(dim(ggjab$edges)[1], 2714)
    expect_equal(dim(ggjab$edges)[2], 12)
    expect_equal(length(ggjab$junctions), 269)
    expect_error(ggjab$G, NA)  ## check it works
    expect_equal(length(ggjab$adj), 5503716)
    expect_equal(length(ggjab$A), 5503716)
    expect_equal(length(ggjab$parts), 3)
    expect_equal(length((ggjab$parts)$csize), 63)
    expect_equal(length(ggjab$seqinfo), 85)
    expect_equal(ggjab$purity, 0.96)
    expect_equal(round(ggjab$ploidy, 2), 3.85)
    expect_true(is(ggjab$td, 'gTrack'))
    expect_equal(length(ggjab$win), 85)
    ### SLOW ### expect_equal(ggjab$ig, NULL)
    ## FUNCTIONS
    ## set.seqinfo = function(genome=NULL, gname=NULL, drop=FALSE)
    ggjab_setseq = ggjab$set.seqinfo()
    expect_true(is(ggjab_setseq, 'gGraph'))
    expect_equal(length(ggjab_setseq$segstats), 2346)
    expect_equal(dim(ggjab_setseq$edges)[1], 2714)
    expect_equal(dim(ggjab_setseq$edges)[2], 12)
    expect_equal(length(ggjab_setseq$junctions), 269)
    expect_error(ggjab_setseq$G, NA)  ## check it works
    expect_equal(length(ggjab_setseq$adj), 5503716)
    expect_equal(length(ggjab_setseq$A), 5503716)
    expect_equal(length((ggjab_setseq$parts)$membership), 2346)
    expect_equal(length(ggjab_setseq$seqinfo), 85)
    expect_equal(ggjab_setseq$purity, 0.96)
    expect_equal(round(ggjab_setseq$ploidy, 2), 3.85)
    expect_true(is(ggjab_setseq$td, 'gTrack'))
    expect_equal(length(ggjab_setseq$win), 85)
    ### SLOW ### expect_equal(ggnew_setseq$ig, NULL)
    ## set.seqinfo, drop = TRUE
    ggjab_setseq_drop = ggjab$set.seqinfo(gname = 'foobar', drop=TRUE)
    expect_true(is(ggjab_setseq_drop, 'gGraph'))
    expect_equal(length(ggjab_setseq_drop$segstats), 2346)
    expect_equal(dim(ggjab_setseq_drop$edges)[1], 2714)
    expect_equal(dim(ggjab_setseq_drop$edges)[2], 12)
    expect_equal(length(ggjab_setseq_drop$junctions), 269)
    expect_error(ggjab_setseq_drop$G, NA)   ## check it works
    expect_equal(length(ggjab_setseq_drop$adj), 5503716)
    expect_equal(length(ggjab_setseq_drop$A), 5503716)
    expect_equal((ggjab_setseq_drop$parts)$no, 63)
    expect_equal(length(ggjab_setseq_drop$seqinfo), 85)
    expect_equal(ggjab_setseq_drop$purity,  0.96)
    expect_equal(round(ggjab_setseq_drop$ploidy, 2), 3.85)
    expect_true(is(ggjab_setseq_drop$td, 'gTrack'))
    expect_equal(length(ggjab_setseq_drop$win), 85)
    ## SLOW expect_equal(ggjab_setseq_drop$ig, NULL)
    ## set.seqinfo, genome != NULL, gname != NULL
    ##ggjab = gGraph$new(jabba=jabba, cn=TRUE)
    ##ggjab_setseq_hg = ggjab$set.seqinfo(genome = hg_seqlengths(), gname = 'foobar', drop = TRUE)
    ##expect_true(is(ggjab_setseq_hg, 'gGraph'))
    ##expect_equal(length(ggjab_setseq_hg$segstats), 0)
    ##expect_equal(dim(ggjab_setseq_hg$edges)[1], 0)
    ##expect_equal(length(ggjab_setseq_hg$junctions), 0)
    ##expect_error(ggjab_setseq_hg$G, NA) ## check it works
    ##expect_equal(length(ggjab_setseq_hg$adj), 0)
    ##expect_equal(length(ggjab_setseq_hg$A), 0)
    ##expect_equal(ggjab_setseq_hg$parts, NULL)
    ##expect_equal(length(ggjab_setseq_hg$seqinfo), 25)   ### checks!
    ##expect_equal(ggjab_setseq_hg$purity, NULL)
    ##expect_equal(ggjab_setseq_hg$ploidy, NULL)
    ##expect_equal(ggjab_setseq_hg$td, NULL)
    ##expect_equal(length(ggjab_setseq_hg$win), 0)
    ### SLOW ### expect_equal(ggjab_setseq_hg$ig, NULL)
    ##
    ## nullGraph = function(regular=TRUE, genome=NULL)
    ggjab_setseq_nullGraph = ggjab$nullGGraph()
    expect_true(is(ggjab_setseq_nullGraph, 'gGraph'))
    expect_equal(length(ggjab_setseq_nullGraph$segstats), 0)
    expect_equal(dim(ggjab_setseq_nullGraph$edges)[1], 0)
    expect_equal(length(ggjab_setseq_nullGraph$junctions), 0)
    expect_error(ggjab_setseq_nullGraph$G, NA) ## check it works
    expect_equal(length(ggjab_setseq_nullGraph$adj), 0)
    expect_equal(length(ggjab_setseq_nullGraph$A), 0)
    expect_equal(ggjab_setseq_nullGraph$parts, NULL)
    expect_equal(length(ggjab_setseq_nullGraph$seqinfo), 85)
    expect_equal(ggjab_setseq_nullGraph$purity, 0.96)
    expect_equal(round(ggjab_setseq_nullGraph$ploidy, 2), 3.85)
    expect_equal(ggjab_setseq_nullGraph$td, NULL)
    expect_equal(length(ggjab_setseq_nullGraph$win), 0)
    ### SLOW ### expect_equal(ggjab_setseq_nullGraph$ig, NULL)
    ##
    ## simpleGraph = function(genome = NULL, chr=FALSE, include.junk=FALSE, ploidy = NULL)
    ggjab_setseq_simpleGraph = ggjab$simpleGraph()
    expect_true(is(ggjab_setseq_simpleGraph, 'gGraph'))
    expect_equal(length(ggjab_setseq_simpleGraph$segstats), 50)
    expect_equal(dim(ggjab_setseq_simpleGraph$edges)[1], 0)
    expect_equal(length(ggjab_setseq_simpleGraph$junctions), 0)
    expect_error(ggjab_setseq_simpleGraph$G, NA) ## check it works
    expect_equal(length(ggjab_setseq_simpleGraph$adj), 2500)
    expect_equal(length(ggjab_setseq_simpleGraph$A), 2500)
    ## expect_equal(ggnew_setseq_simpleGraph$parts, NULL)
    expect_equal(length(ggjab_setseq_simpleGraph$seqinfo), 25)   ### checks! "null" means there is no node, you can still have a "space" of possible values when the set is empty
    expect_equal(ggjab_setseq_simpleGraph$purity, NULL)
    expect_equal(ggjab_setseq_simpleGraph$ploidy, NULL)
    expect_true(is(ggjab_setseq_simpleGraph$td, 'gTrack'))
    expect_equal((ggjab_setseq_simpleGraph$td)$ygap, 2)
    expect_match((ggjab_setseq_simpleGraph$td)$name, 'CN')
    expect_equal(length(ggjab_setseq_simpleGraph$win), 25)
    ## $ig
    ## dipGraph = function(genome = NULL, chr=FALSE, include.junk=FALSE)
    ggjab_dd = ggjab$dipGraph()
    expect_true(is(ggjab_dd, 'gGraph'))
    expect_equal(length(ggjab_dd$segstats), 50)
    expect_equal(dim(ggjab_dd$edges)[1], 0)
    expect_equal(length(ggjab_dd$junctions), 0)
    expect_error(ggjab_dd$G, NA) ## check it works
    expect_equal(length(ggjab_dd$adj), 2500)
    expect_equal(length(ggjab_dd$A), 2500)
    ## > ggjab_dd$ig
    ## > ggjab_dd$parts
    expect_equal(length(ggjab_dd$seqinfo), 25)
    expect_equal(ggjab_dd$purity, NULL)
    expect_equal(ggjab_dd$ploidy, 2)  ## checks!
    expect_true(is(ggjab_dd$td, 'gTrack'))
    expect_equal((ggjab_dd$td)$ygap, 2)
    expect_match((ggjab_dd$td)$name, 'CN')
    expect_equal(length(ggjab_dd$win), 25)
    ##
    ## dipGraph
    ggjab_dd_junk = ggjab$dipGraph(genome = hg_seqlengths(), chr=TRUE, include.junk=TRUE)
    expect_true(is(ggjab_dd_junk, 'gGraph'))
    expect_equal(length(ggjab_dd_junk$segstats), 50)
    expect_equal(dim(ggjab_dd_junk$edges)[1], 0)
    expect_equal(length(ggjab_dd_junk$junctions), 0)
    expect_error(ggjab_dd_junk$G, NA) ## check it works
    expect_equal(length(ggjab_dd_junk$adj), 2500)
    expect_equal(length(ggjab_dd_junk$A), 2500)
    ## expect_equal(ggjab_dd$parts, NULL)
    expect_equal(length(ggjab_dd_junk$seqinfo), 25)
    expect_equal(ggjab_dd_junk$purity, NULL)
    expect_equal(ggjab_dd_junk$ploidy, 2)  ## checks!
    expect_true(is(ggjab_dd_junk$td, 'gTrack'))
    expect_equal((ggjab_dd_junk$td)$ygap, 2)
    expect_match((ggjab_dd_junk$td)$name, 'CN')
    expect_equal(length(ggjab_dd_junk$win), 25)
    ##
    ## addJuncs = function(junc, cn=TRUE)
    expect_error(ggjab$addJuncs())
    added_juncs = ggjab$addJuncs(junc = junctions)
    expect_true(is(added_juncs, 'gGraph'))
    expect_equal(length(added_juncs$segstats), 1112)
    expect_equal(dim(added_juncs$edges)[1], 1596)
    expect_equal(dim(added_juncs$edges)[2], 15)
    expect_equal(length(added_juncs$junctions), 267)
    expect_error(added_juncs$G, NA) ## check it works
    expect_equal(length(added_juncs$adj), 1236544)
    expect_equal(length(added_juncs$A),  1236544)
    ## expect_equal(added_juncs$parts, NULL)
    expect_equal(length(added_juncs$seqinfo), 25)
    expect_equal(added_juncs$purity, NULL)
    expect_equal(added_juncs$ploidy, NULL)  ## checks!
    expect_true(is(added_juncs$td, 'gTrack'))
    expect_equal((added_juncs$td)$ygap, 2)
    expect_match((added_juncs$td)$name, 'CN')
    expect_equal(length(added_juncs$win), 25)
    ## addJuncs, cn = FALSE
    ## added_juncs_cnfalse = ggnew$addJuncs(junc = junctions, cn = FALSE) ERROR
    ##
    ## addSegs = function(tile)
    added_segs = ggjab$addSegs(tile = test_segs)
    expect_true(is(added_segs, 'gGraph'))
    expect_equal(length(added_segs$segstats), 1124)
    expect_equal(dim(added_segs$edges)[1], 1074)
    expect_equal(dim(added_segs$edges)[2], 14)
    expect_equal(length(added_segs$junctions), 0)
    expect_error(added_segs$G, NA) ## check it works
    expect_equal(length(added_segs$adj), 1263376)
    expect_equal(length(added_segs$A),  1263376)
    ## expect_equal(added_segs$parts, NULL)
    expect_equal(length(added_segs$seqinfo), 25)
    expect_equal(added_segs$purity, NULL)
    expect_equal(added_segs$ploidy, NULL)  ## checks!
    expect_true(is(added_segs$td, 'gTrack'))
    expect_equal((added_segs$td)$ygap, 2)
    expect_match((added_segs$td)$name, 'CN')
    expect_equal(length(added_segs$win), 25)
    ## karyograph = function(tile=NULL, juncs=NULL, cn=FALSE, regular=FALSE)
    ## default
    ggjab_karyograph = ggjab$karyograph()
    expect_true(is(ggjab_karyograph, 'gGraph'))
    expect_equal(length(ggjab_karyograph$segstats), 1124)
    expect_equal(dim(ggjab_karyograph$edges)[1], 1074)
    expect_equal(dim(ggjab_karyograph$edges)[2], 23)
    expect_equal(length(ggjab_karyograph$junctions), 0)
    expect_error(ggjab_karyograph$G, NA) ## check it works
    expect_equal(length(ggjab_karyograph$adj), 1263376)
    expect_equal(length(ggjab_karyograph$A),  1263376)
    ## expect_equal(length((ggjab_karyograph$parts)$membership), 50)
    ## expect_equal(length((ggjab_karyograph$parts)$csize), 25)
    ## expect_equal(length(ggjab_karyograph$seqinfo), 25)
    expect_equal(ggjab_karyograph$purity, NULL)
    expect_equal(ggjab_karyograph$ploidy, NULL)  ## checks!
    expect_true(is(ggjab_karyograph$td, 'gTrack'))
    expect_equal((ggjab_karyograph$td)$ygap, 2)
    expect_match((ggjab_karyograph$td)$name, 'CN')
    expect_equal(length(ggjab_karyograph$win), 25)
    ##
    ## simplify = function(mod=TRUE)
    ggjab_simplify = ggjab$simplify()
    expect_true(is(ggjab_simplify, 'gGraph'))
    expect_equal(length(ggjab_simplify$segstats), 50)
    expect_equal(dim(ggjab_simplify$edges)[1], 0)
    expect_equal(dim(ggjab_simplify$edges)[2], 11)
    expect_equal(length(ggjab_simplify$junctions), 0)
    expect_error(ggjab_simplify$G, NA) ## check it works
    expect_equal(length(ggjab_simplify$adj), 2500)
    expect_equal(length(ggjab_simplify$A),  2500)
    ##expect_equal(length((ggjab_simplify$parts)$membership), 50)
    ##expect_equal(length((ggjab_simplify$parts)$csize), 25)
    expect_equal(length(ggjab_simplify$seqinfo), 25)
    expect_equal(ggjab_simplify$purity, NULL)
    expect_equal(ggjab_simplify$ploidy, NULL)  ## checks!
    expect_true(is(ggjab_simplify$td, 'gTrack'))
    expect_equal((ggjab_simplify$td)$ygap, 2)
    expect_match((ggjab_simplify$td)$name, 'CN')
    expect_equal(length(ggjab_simplify$win), 25)
    ## simplify = function(mod=FALSE)
    ggjab_simplify = ggjab$simplify(mod=FALSE)
    expect_true(is(ggjab_simplify, 'gGraph'))
    expect_equal(length(ggjab_simplify$segstats), 50)
    expect_equal(dim(ggjab_simplify$edges)[1], 0)
    expect_equal(dim(ggjab_simplify$edges)[2], 11)
    expect_equal(length(ggjab_simplify$junctions), 0)
    expect_error(ggjab_simplify$G, NA) ## check it works
    expect_equal(length(ggjab_simplify$adj), 2500)
    expect_equal(length(ggjab_simplify$A),  2500)
    ##expect_equal(length((ggjab_simplify$parts)$membership), 50)
    ##expect_equal(length((ggjab_simplify$parts)$csize), 25)
    expect_equal(length(ggjab_simplify$seqinfo), 25)
    expect_equal(ggjab_simplify$purity, NULL)
    expect_equal(ggjab_simplify$ploidy, NULL)  ## checks!
    expect_true(is(ggjab_simplify$td, 'gTrack'))
    expect_equal((ggjab_simplify$td)$ygap, 2)
    expect_match((ggjab_simplify$td)$name, 'CN')
    expect_equal(length(ggjab_simplify$win), 25)
    ##
    ## decouple = function(mod=TRUE)
    ggjab_decouple = ggjab$decouple()
    expect_true(is(ggjab_decouple, 'gGraph'))
    expect_equal(length(ggjab_decouple$segstats), 50)
    expect_equal(dim(ggjab_decouple$edges)[1], 0)
    expect_equal(dim(ggjab_decouple$edges)[2], 19)
    expect_equal(length(ggjab_decouple$junctions), 0)
    expect_error(ggjab_decouple$G, NA) ## check it works
    expect_equal(length(ggjab_decouple$adj), 2500)
    expect_equal(length(ggjab_decouple$A),  2500)
    ##expect_equal(length((ggjab_decouple$parts)$membership), 50)
    ##expect_equal(length((ggjab_decouple$parts)$csize), 25)
    expect_equal(length(ggjab_decouple$seqinfo), 25)
    expect_equal(ggjab_decouple$purity, NULL)
    expect_equal(ggjab_decouple$ploidy, NULL)  ## checks!
    expect_true(is(ggjab_decouple$td, 'gTrack'))
    expect_equal((ggjab_decouple$td)$ygap, 2)
    expect_match((ggjab_decouple$td)$name, 'CN')
    expect_equal(length(ggjab_decouple$win), 25)
    ##
    ggjab_decouple = ggjab$decouple(mod=FALSE)
    expect_true(is(ggjab_decouple, 'gGraph'))
    expect_equal(length(ggjab_decouple$segstats), 50)
    expect_equal(dim(ggjab_decouple$edges)[1], 0)
    expect_equal(dim(ggjab_decouple$edges)[2], 19)
    expect_equal(length(ggjab_decouple$junctions), 0)
    expect_error(ggjab_decouple$G, NA) ## check it works
    expect_equal(length(ggjab_decouple$adj), 2500)
    expect_equal(length(ggjab_decouple$A),  2500)
    ##expect_equal(length((ggjab_decouple$parts)$membership), 50)
    ##expect_equal(length((ggjab_decouple$parts)$csize), 25)
    expect_equal(length(ggjab_decouple$seqinfo), 25)
    expect_equal(ggjab_decouple$purity, NULL)
    expect_equal(ggjab_decouple$ploidy, NULL)  ## checks!
    expect_true(is(ggjab_decouple$td, 'gTrack'))
    expect_equal((ggjab_decouple$td)$ygap, 2)
    expect_match((ggjab_decouple$td)$name, 'CN')
    expect_equal(length(ggjab_decouple$win), 25)
    ##
    ## add = function(gg, mod=FALSE)
    added = ggjab$add(gg=ggjab_karyograph, mod=TRUE)
    expect_true(is(added, 'gGraph'))
    expect_equal(length(added$segstats), 100)
    expect_equal(dim(added$edges)[1], 0)
    expect_equal(dim(added$edges)[2], 10)
    expect_equal(length(added$junctions), 0)
    expect_error(added$G, NA) ## check it works
    expect_equal(length(added$adj), 10000)
    expect_equal(length(added$A),  10000)
    ##expect_true(is(added$parts, 'list'))
    ##expect_equal(length(added$parts), 3)
    ##expect_equal(length(added$parts$membership), 100)
    ##expect_equal(length((added$parts)$csize), 25)
    expect_equal(length(added$seqinfo), 25)
    expect_equal(added$purity, NULL)
    expect_equal(added$ploidy, NULL)  ## checks!
    expect_true(is(added$td, 'gTrack'))
    expect_equal((added$td)$ygap, 2)
    expect_match((added$td)$name, 'CN')
    expect_equal(length(added$win), 25)
    ## add() several times
    five_adds = ggjab$add(gg=ggjab_karyograph)$add(gg=ggjab_karyograph)$add(gg=ggjab_karyograph)$add(gg=ggjab_karyograph)$add(gg=ggjab_karyograph)
    expect_true(is(five_adds, 'gGraph'))
    expect_equal(length(five_adds$segstats), 600)
    expect_equal(dim(five_adds$edges)[1], 0)
    expect_equal(dim(five_adds$edges)[2], 10)
    expect_equal(length(five_adds$junctions), 0)
    expect_error(five_adds$G, NA) ## check it works
    expect_equal(length(five_adds$adj), 360000)
    expect_equal(length(five_adds$A), 360000)
    ##expect_true(is(five_adds$parts, 'list'))
    ##expect_equal(length(five_adds$parts), 3)
    ##expect_equal(length(five_adds$parts$membership), 300)
    ##expect_equal(length((five_adds$parts)$csize), 25)
    expect_equal(length(five_adds$seqinfo), 25)
    expect_equal(five_adds$purity, NULL)
    expect_equal(five_adds$ploidy, NULL)  ## checks!
    expect_true(is(five_adds$td, 'gTrack'))
    expect_equal((five_adds$td)$ygap, 2)
    expect_match((five_adds$td)$name, 'CN')
    expect_equal(length(five_adds$win), 25)
    ##
    ##
    ## wv2gg()
    weavd = ggjab$wv2gg(weaver)
    expect_true(is(weavd, 'gGraph'))
    expect_equal(length(weavd$segstats), 16184)
    expect_equal(dim(weavd$edges)[1], 16974)
    expect_equal(dim(weavd$edges)[2], 15)
    expect_equal(length(weavd$junctions), 420)
    expect_error(weavd$G, NA) ## check it works
    ## ERROR expect_equal(length(weavd$adj),  90000)
    ## ERROR expect_equal(length(weavd$A),   90000)
    ##expect_true(is(weavd$parts, 'list'))
    ##expect_equal(length(weavd$parts), 3)
    ##expect_equal(length(weavd$parts$membership), 300)
    ##expect_equal(length((weavd$parts)$csize), 25)
    expect_equal(length(weavd$seqinfo), 25)
    expect_equal(weavd$purity, NULL)
    expect_equal(round(weavd$ploidy, 2), 8.54)
    expect_true(is(weavd$td, 'gTrack'))
    expect_equal((weavd$td)$ygap, 2)
    expect_match((weavd$td)$name, 'CN')
    expect_equal(length(weavd$win), 25)
    ## pr2gg()
    pregod = ggjab$pr2gg(prego)
    expect_true(is(pregod, 'gGraph'))
    expect_equal(length(pregod$segstats), 1208)
    expect_equal(dim(pregod$edges)[1], 1380)
    expect_equal(dim(pregod$edges)[2], 21)
    expect_equal(length(pregod$junctions), 420)
    expect_error(pregod$G, NA) ## check it works
    ## ERROR expect_equal(length(weavd$adj),  90000)
    ## ERROR expect_equal(length(weavd$A),   90000)
    ##expect_true(is(weavd$parts, 'list'))
    ##expect_equal(length(weavd$parts), 3)
    ##expect_equal(length(weavd$parts$membership), 300)
    ##expect_equal(length((weavd$parts)$csize), 25)
    expect_equal(length(pregod$seqinfo), 84)
    expect_equal(pregod$purity, 1)
    expect_equal(round(pregod$ploidy, 2), 2.1)
    expect_true(is(pregod$td, 'gTrack'))
    expect_equal((pregod$td)$ygap, 2)
    expect_match((pregod$td)$name, 'CN')
    expect_equal(length(pregod$win), 24)
    ##
    ## print() to STDOUT
    ggjab = gGraph$new(jabba=jabba, cn=TRUE)
    expect_error(ggjab$print(), NA) ## check it runs
    ## plot() nothing is returned here, so let's do this for now:
    expect_error(ggjab$plot(), NA) ## check it runs
    expect_error(ggjab$plot(colorful=TRUE), NA)
    ## window()
    expect_equal(length(ggjab$window()), 85)
    ##
    ## > ggnew$layout()
    ## Error in log(private$segs$cn, 1.4) :
    ## summary()
    expect_true(is.character(ggjab$summary()))
    ## length()
    ### ERROR  expect_equal(ggjab$length(), NULL)
    ##
    expect_true(is(ggjab$gg2td(), 'gTrack'))
    ## JSON
    ## > ggjab$json()
    expect_error(ggjab$json(), NA)
    ## HTML
    ## > ggjab$html()
    expect_error(ggjab$html()) ## Error in ggjab$html() : Get from https://github.com/mskilab/gGnome.js
    ## No gGnome.js repository found on your system.
    ## Error in ggnew$html() : Get from https://github.com/mskilab/gGnome.js
    ## gg2j()
    ## > ggjab$gg2js()
    expect_error(ggjab$gg2js(), NA)
    ##
    ## component()
    ### ggjab = gGraph$new(jabba=jabba, cn=TRUE)      ############################################################################
    ## component = ggjab$components(mc.cores=2)
    ## expect_true(is(component, 'list'))
    ## expect_equal(length(component), 63)
    ## expect_equal(length(component$segstats), 0)
    ## expect_equal(dim(component$edges)[1], NULL)
    ## expect_equal(dim(component$edges)[2], NULL)
    ## expect_equal(length(component$junctions), 0)
    ## expect_equal(component$G, NULL)
    ## expect_equal(length(component$adj), 0)
    ## expect_equal(length(component$A),  0)
    ## expect_equal(length(component$parts), 0)
    ## expect_equal(length(component$seqinfo), 0)
    ## expect_equal(component$purity, NULL)
    ## expect_equal(component$ploidy, NULL)  ## checks!
    ## expect_equal(component$td, NULL)
    ## expect_equal(length(component$win), 0)
    ##
    ## subgraph = function(v=numeric(0), na.rm=T, mod=T)
    ## default
    subgraphed = ggjab$subgraph()
    expect_true(is(subgraphed, 'gGraph'))
    ### ERROR Error in add_edges(stickyG, t(as.matrix(hB[, .(from, to)]))) :  expect_equal(length(subgraphed), 63)
    expect_equal(length(subgraphed$segstats), 2346)
    expect_equal(dim(subgraphed$edges)[1], 2714)
    expect_equal(dim(subgraphed$edges)[2], 12)
    expect_equal(length(subgraphed$junctions), 269)
    expect_error(subgraphed$G, NA)
    expect_equal(length(subgraphed$adj), 5503716)
    expect_equal(length(subgraphed$A),  5503716)
    expect_equal(length(subgraphed$parts), 3)
    expect_equal(length(subgraphed$seqinfo), 85)
    expect_equal(subgraphed$purity, 0.96)
    expect_equal(round(subgraphed$ploidy, 2), 3.85)  ## checks!
    expect_true(is(subgraphed$td, 'gTrack'))
    expect_equal(length(subgraphed$win), 85)
    ##
    ## vertices5K = ggjab$subgraph(v=5000)
    ##
    filled = ggjab$fillin()
    expect_true(is(filled, 'gGraph'))
    expect_equal(length(filled$segstats), 2346)
    expect_equal(dim(filled$edges)[1], 2714)
    expect_equal(dim(filled$edges)[2],  12)
    expect_equal(length(filled$junctions), 269)
    expect_error(filled$G, NA)
    expect_equal(length(filled$adj), 5503716)
    expect_equal(length(filled$A),  5503716)
    ##expect_equal(length(filled$parts), 3)
    expect_equal(length(filled$seqinfo), 85)
    expect_equal(filled$purity, 0.96)
    expect_equal(round(filled$ploidy, 2), 3.85)  ## checks!
    expect_true(is(filled$td, 'gTrack'))
    expect_equal(length(filled$win), 85)
    ##
    ## trim = function(gr=NULL, mod=FALSE)
    ## default
    trimmed = ggjab$trim()
    expect_true(is(trimmed, 'gGraph'))
    expect_equal(length(trimmed$segstats), 2346)
    expect_equal(dim(trimmed$edges)[1], 2714)
    expect_equal(dim(trimmed$edges)[2], 12)
    expect_equal(length(trimmed$junctions), 269)
    expect_error(trimmed$G, NA)
    expect_equal(length(trimmed$adj), 5503716)
    expect_equal(length(trimmed$A),  5503716)
    ##expect_equal(length(trimmed$parts), 0)
    expect_equal(length(trimmed$seqinfo), 85)
    expect_equal(trimmed$purity,  0.96)
    expect_equal(round(trimmed$ploidy, 2), 3.85)  ## checks!
    expect_true(is(trimmed$td, 'gTrack'))
    expect_equal(length(trimmed$win), 85)
    ##
    ## trimmed_mod = ggjab$trim(gr=gr2, mod=TRUE)
    ##  ############################################################################
    ## gotg
    ## default
    gotg = ggjab$get.g()
    expect_true(is(gotg, 'gGraph'))
    expect_equal(length(gotg$segstats), 2346)
    expect_equal(dim(gotg$edges)[1], 2714)
    expect_equal(dim(gotg$edges)[2], 12)
    expect_equal(length(gotg$junctions), 269)
    expect_error(gotg$G, NA)
    expect_equal(length(gotg$adj), 5503716)
    expect_equal(length(gotg$A), 5503716)
    ##expect_equal(length(gotg$parts), 0)
    expect_equal(length(gotg$seqinfo), 85)
    expect_equal(gotg$purity, 0.96)
    expect_equal(round(gotg$ploidy, 2), 3.85)   ## checks!
    expect_true(is(gotg$td, 'gTrack'))
    expect_equal(length(gotg$win), 85)
    ## get.g = function(force=FALSE)
    gotg_forced = ggjab$get.g(force=TRUE)
    expect_true(is(gotg_forced, 'gGraph'))
    expect_equal(length(gotg_forced$segstats), 2346)
    expect_equal(dim(gotg_forced$edges)[1], 2714)
    expect_equal(dim(gotg_forced$edges)[2], 12)
    expect_equal(length(gotg_forced$junctions), 269)
    expect_error(gotg_forced$G, NA)
    expect_equal(length(gotg_forced$adj), 5503716)
    expect_equal(length(gotg_forced$A), 5503716)
    ##expect_equal(length(gotg_forced$parts), 3)
    expect_equal(length(gotg_forced$seqinfo), 85)
    expect_equal(gotg_forced$purity, 0.96)
    expect_equal(round(gotg_forced$ploidy, 2), 3.85)   ## checks!
    expect_true(is(gotg_forced$td, 'gTrack'))
    expect_equal(length(gotg_forced$win), 85)
    ## with added
    gotgadd = added$get.g(force=TRUE)
    expect_true(is(gotgadd, 'gGraph'))
    expect_equal(length(gotgadd$segstats), 1208)
    expect_equal(dim(gotgadd$edges)[1], 1380)
    expect_equal(dim(gotgadd$edges)[2], 29)
    expect_equal(length(gotgadd$junctions), 420)
    expect_error(gotgadd$G, NA) ## check works
    expect_equal(length(gotgadd$adj), 1459264)
    expect_equal(length(gotgadd$A),  1459264)
    ##expect_equal(length(gotgadd$parts), 3)
    ##expect_equal(length((gotgadd$parts)$membership), 1208)
    ##expect_equal(length((gotgadd$parts)$csize), 25)
    ##expect_equal(length((gotgadd$parts)$no), 1)  ## 25
    expect_equal(length(gotgadd$seqinfo), 84)
    expect_equal(gotgadd$purity, 1)
    expect_equal(round(gotgadd$ploidy, 2), 2.1)
    expect_true(is(gotgadd$td, 'gTrack'))
    expect_equal((gotgadd$td)$ygap, 2)
    expect_match((gotgadd$td)$name, 'CN')
    expect_equal(length(gotgadd$win), 24)
    ##  hood = function(win, d=NULL, k=NULL, pad=0, bagel=FALSE, ignore.strand=T, verbose=FALSE)
    ##gr2_win = ggnew$hood(win=gr2)
    expect_error(ggjab$hood(win=grl1)) ## Error in .local(x, y, ...) : setdiff() between a GRanges and a GRangesList object is not supported
    ##grl2_win = ggnew$hood(win=grl.unlist(grl2))
    ##hgseq_win = ggnew$hood(win = si2gr(hg_seqlengths()))
    ##expect_true(is(hgseq_win, 'gGraph'))
    ##expect_equal(length(hgseq_win$segstats), 0)
    ##expect_equal(dim(hgseq_win$edges)[1], 0)
    ##expect_equal(dim(hgseq_win$edges)[2], 19)
    ##expect_equal(length(hgseq_win$junctions), 0)
    ##expect_error(hgseq_win$G, NA)
    ##expect_equal(length(hgseq_win$adj), 0)
    ##expect_equal(length(hgseq_win$A),  0)
    ##expect_equal(length(hgseq_win$parts), 0)
    ##expect_equal(hgseq_win$purity, NULL)
    ##expect_equal(hgseq_win$ploidy, NULL)  ## checks!
    ##expect_equal(length(hgseq_win$td), 0)
    ##expect_equal(length(hgseq_win$win), 0)
    ##
    ## dist
    ## ERROR    distanced1 = ggjab$dist(GRanges('1:5500-6000'), GRanges('1:5000-5500'))
    ## ERROR    expect_equal(as.numeric(distanced1), 0)
    ## ERROR distanced2 = ggjab$dist(GRanges('1:5500-6000'), GRanges('1:15000-15500'))
    ## ERROR    distanced_diffchroms = ggjab$dist(GRanges('2:5500-6000'), GRanges('3:5000-5500'))
    ## ERROR    expect_equal(as.numeric(distanced_diffchroms), Inf)
    ##
    ## fillup
    ## ERROR    filledupjab = ggjab$fillup()
    ## ERROR    expect_true(is(filledupjab, 'gGraph'))
    ## ERROR    expect_equal(length(filledupjab$segstats), 1208)
    ## ERROR    expect_equal(dim(filledupjab$edges)[1], 1380)
    ## ERROR    expect_equal(dim(filledupjab$edges)[2], 29)
    ## ERROR    expect_equal(length(filledupjab$junctions), 420)
    ## ERROR    expect_error(filledupjab$G, NA)
    ## ERROR    expect_equal(length(filledupjab$adj), 1459264)
    ## ERROR    expect_equal(length(filledupjab$A),  1459264)
    ## expect_equal(length(filledupjab$parts), 3)
    ##expect_equal(length((filledupjab$parts)$membership), 100)
    ##expect_equal(length((filledupjab$parts)$csize), 25)
    ##expect_equal(length((filledupjab$parts)$no), 1)  ## 25
    ## ERROR    expect_equal(length(filledupjab$seqinfo), 84)
    ## ERROR    expect_equal(filledupjab$purity, 1)
    ## ERROR    expect_equal(round(filledupjab$ploidy, 2), 2.1)  ## checks!
    ## ERROR    expect_true(is(filledupjab$td, 'gTrack'))
    ## ERROR    expect_equal((filledupjab$td)$ygap, 2)
    ## ERROR    expect_match((filledupjab$td)$name, 'CN')
    ## ERROR    expect_equal(length(filledupjab$win), 24)
    ## isBalance
    ## ERROR    expect_true(ggjab$isBalance())
    ## get.loose
    ## ERROR    expect_equal(length(ggjab$get.loose()), 0)
    
})





### bGraph, default
### bGraph inherits all active bindings from gGraph
## METHODS
## print
## subgraph
## walk
## walk2
##




test_that('check bGraph works, default', {
    
    options(gGnome.verbose=TRUE)
    bgnew = bGraph$new()
    expect_true(is(bgnew, 'bGraph'))
    expect_equal(length(bgnew$segstats), 0)
    expect_equal(dim(bgnew$edges)[1], 0)
    expect_equal(dim(bgnew$edges)[2], 3)
    expect_equal(length(bgnew$junctions), 0)
    expect_error(bgnew$G, NA)
    expect_equal(length(bgnew$adj), 0)
    expect_equal(length(bgnew$A),  0)
    ##expect_equal(length(bgnew$parts), 0)
    expect_equal(length(bgnew$seqinfo), 0)
    expect_equal(bgnew$purity, NULL)
    expect_equal(bgnew$ploidy, NULL)  ## checks!
    expect_equal(bgnew$td, NULL)
    expect_equal(length(bgnew$win), 0)
    ##
    ##
    ##
    ## print()
    ## variable outputs 'numeric(0)'
    expect_error(bgnew$print(), NA) ## check it works
    ##
    ## subgraph()
    subgraphed = bgnew$subgraph()
    expect_true(is(subgraphed, 'bGraph'))
    expect_equal(length(subgraphed$segstats), 0)
    expect_equal(dim(subgraphed$edges)[1], 0)
    expect_equal(dim(subgraphed$edges)[2], 3)
    expect_equal(length(subgraphed$junctions), 0)
    expect_error(subgraphed$G, NA)
    expect_equal(length(subgraphed$adj), 0)
    expect_equal(length(subgraphed$A),  0)
    ##expect_equal(length(bgnew$parts), 0)
    expect_equal(length(subgraphed$seqinfo), 0)
    expect_equal(subgraphed$purity, NULL)
    expect_equal(subgraphed$ploidy, NULL)  ## checks!
    expect_equal(subgraphed$td, NULL)
    expect_equal(length(subgraphed$win), 0)
    ##
    ##
    ## walk()
    expect_equal(bgnew$walk(), NULL)
    ##
    ## walk2()
    expect_equal(bgnew$walk2(), NULL)
    
})




###
### gWalks
### METHODS:
### nullGWalks
### pairup
### set.seg.cn
### gw2gg
### gw2grl
### grl2gw
### gw2td
### json
### gw2js
### v2e
### p2e
### epath
### simplify
### reduce
### subset
### print
### len
### metaCols
### window
### plot
### isStrandPaired
### rpaths
### label
###
### ACTIVE BINDINGS:
### segstats
### edges
### grl
### td
### path
### values




test_that('check gWalks works, default', {
    newgwalk = gWalks$new()
    expect_equal(newgwalk$segstats, NULL)
    ### edges
    ### grl
    ### td
    expect_equal(newgwalk$path, NULL)
    expect_equal(newgwalk$value, NULL)
    ##
    ## nullGWalks
    nullgw = newgwalk$nullGWalks()
    expect_equal(nullgw$segstats, NULL)
    ### edges
    ### grl
    ### td
    expect_equal(nullgw$path, NULL)
    expect_equal(nullgw$value, NULL)
    ##
    ## pairup
    ## ERROR
    ### pairedup = newgwalk$pairup()
    setsegcnd = nullgw$set.seg.cn()
    expect_equal(length(setsegcnd$segstats), 1)
    ### edges
    ### grl
    ### td
    expect_equal(setsegcnd$path, NULL)
    expect_equal(setsegcnd$value, NULL)
    ##
    ## gw2gg
    ## gw2ggfoo = newgwalk$gw2gg()
    ## Error in do.call(c, new.paths) : second argument must be a list
    ##
    ## gw2grl
    ## gw2grlfoo = newgwalk$gw2grl()
    ##
    ## grl2gw
    grl2gwfoo = newgwalk$grl2gw(grl=grl2)
    expect_equal(length(grl2gwfoo$segstats), 502)
    expect_equal(dim(grl2gwfoo$edges)[1], 251)
    expect_equal(dim(grl2gwfoo$edges)[2], 16)
    expect_true(is(grl2gwfoo$grl, 'GRangesList'))
    expect_equal(length(grl2gwfoo$grl), 251)
    expect_true(is(grl2gwfoo$td, 'gTrack'))
    expect_equal(length(grl2gwfoo$path), 251)
    expect_equal(grl2gwfoo$value, NULL)
    ##
    ## gw2td
    expect_true(is(newgwalk$gw2td(), 'gTrack'))
    ##
    ## json
    jsonfoo = newgwalk$json()
    expect_match(jsonfoo, 'data.json')
    ##
    ## gw2js
    ## gw2jsfoo = newgwalk$gw2js()
    ## expect_equal(length(newgwalk$gw2js()), 0)
    ##
    ## ERRORs
    ### v2e
    ### p2e
    ##
    ## epath
    epathfoo = newgwalk$epath()
    expect_equal(length(epathfoo), 502)
    ##
    ##
    ## > simplified = newgwalk$simplify()
    ## Error in FUN(X[[i]], ...) : object 'cn' not found
    ##
    ## reduced = newgwalk$reduce()
    ## Error in do.call(c, new.paths) : second argument must be a list
    ##
    ## subset()
    ##subsetted = newgwalk$subset(GRanges())
    ##expect_equal(subsetted$segstats, NULL)
    ### edges
    ### grl
    ### td
    ##expect_equal(subsetted$path, NULL)
    ##expect_equal(subsetted$values, NULL)
    ### path
    ### values
    ##
    ## print
    expect_equal(newgwalk$print(), NULL)
    ##
    ## len
    expect_equal(newgwalk$len(), 502)
    ##
    ## metaCols
    expect_equal(dim(newgwalk$metaCols())[1], 502)
    expect_equal(dim(newgwalk$metaCols())[2], 3)
    ##
    ## > newgwalk$window()
    ## Error in (function (classes, fdef, mtable)  :
    ##   unable to find an inherited method for function ‘seqinfo’ for signature ‘"NULL"’
    ## > newgwalk$plot()
    ## Error in (function (classes, fdef, mtable)  :
    ##   unable to find an inherited method for function ‘seqinfo’ for signature ‘"NULL"’
    ## ERROR
    ##  newgwalk$isStrandPaired()
    ##
    
})








test_that('rev.comp works', {
    
    expect_error(rev.comp())
    expect_error(rev.comp(data.frame()))
    gr3 = dt2gr(dt)  ## not stranded
    expect_error(rev.comp(gr3))   ## Error in rev.comp(gr3) : Input must be all strand specific.
    ##expect_equal(as.integer(width(rev.comp(gr2[1]))), 6)
    ##expect_equal(as.character(strand(rev.comp(gr)[1])), "+")
    ##expect_equal(width(rev.comp(gr)[2]), 3)
    ##expect_equal(as.character(strand(rev.comp(gr)[2])), "+")
    
})





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





test_that('ul works', {
    
    A = matrix(c(2, 4, 3, 1, 5, 7),  nrow=2, ncol=3, byrow = TRUE)
    expect_equal(ul(A, n=0), NULL)   ### Is this expected behavior?
    expect_equal(as.integer(ul(A, n=1)), 2)
    expect_equal(dim(ul(A, n=2))[1], 2)
    expect_equal(dim(ul(A, n=2))[2], 2)
    expect_equal(dim(ul(A, n=999))[1], 2)
    expect_equal(dim(ul(A, n=9999))[2], 2)   ### Is this expected behavior?
    
})






test_that('e2j works', {
    
    expect_equal(length(e2j(test_segs, test_es)), 2)
    
})




test_that('etype works', {
    
    expect_equal(dim(etype(test_segs, test_es))[1], 12)
    expect_equal(dim(etype(test_segs, test_es))[2], 16)
    
})







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
    ##expect_equal(dim(foo$edges)[1], 12)
    ##expect_equal(dim(foo$edges)[2], 16)
    expect_equal(max((foo$edges)$cn), 3)
    expect_equal(max((foo$edges)$fromStart), 18593415)
    expect_equal(max((foo$edges)$fromEnd), 18793414)
    foo = gGraph$new(segs=test_segs, es=test_es)
    ##expect_equal(dim(foo$edges)[1], 12)
    ## expect_equal(dim(foo$edges)[2], 16)
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




test_that('test grl.match() ', {
    
    foo = grl.match(grl1, grl2)
    expect_equal(dim(foo)[1], 502)
    expect_equal(dim(foo)[2], 8)
    
})



##-------------------------------------------------------##
test_that('testing ra_breaks', {
    
    breaks = ra_breaks(svabavcf)
    expect_true(is(breaks, 'GRangesList'))
    expect_equal(length(breaks), 494)
    novobreaks = ra_breaks(novobreakvcf)
    expect_true(is(novobreaks, 'GRangesList'))
    expect_equal(length(novobreaks), 383)
    dellybreaks = ra_breaks(dellyvcf)
    expect_true(is(dellybreaks, 'GRangesList'))
    expect_equal(length(dellybreaks), 210)
    lumpybreaks = ra_breaks(lumpyvcf)
    expect_true(is(lumpybreaks, 'GRangesList'))
    expect_equal(length(lumpybreaks), 2115)
    
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




##-------------------------------------------------------##
##test_that('gWalks', {
##
##    jab = system.file('extdata', 'jabba.simple.rds', package="gGnome")
##    message("JaBbA result: ", jab)
##    segments = readRDS(jab)$segstats
##    junctions = readRDS(jab)$junctions
##    grl = system.file("extdata", "gw.grl.rds", package="gGnome")
##    message("Walks for testing:", grl)
##    grl = readRDS(grl)
##    expect_equal(length(gw <<- as(grl, "gWalks")), sum(values(grl)$cn>0))
##    expect_error(bg <<- as(gw, "bGraph"), NA)
##    expect_equal(length(bg$junctions), sum(values(junctions)$cn>0))
##    expect_true(inherits(gw.simp <<- gw$simplify(mod=FALSE), "gWalks"))
##    expect_error(bg.simp <<- as(gw.simp, "bGraph"), NA)
##    expect_error(bg.dc <<- bg.simp$decouple(mod=FALSE), NA)
##    ## expect_equal(length(bg.dc$junctions), length(bg$junctions)) 291>269
##    ## why does simplifying gwalks then decouple create more junctions????
##
##})


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





