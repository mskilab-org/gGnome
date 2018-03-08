
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
    ## ggnew$ig
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

})




## segstats, edges, grl, td, path, values




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

##-------------------------------------------------------##
test_that('bGraph walk and walk2', {
    expect_true(inherits(jab.gg <<- gread(jab), "bGraph"))
    expect_true(inherits(subg <<- jab.gg$hood(jab.gg$junctions[[1]], 5e5), "bGraph"))
    expect_true(inherits(subg.gw2 <<- subg$walk2(T, F), "gWalks")) ## default to CPLEX
    expect_true(inherits(jab.gw2.gurobi <<- subg$walk2(T, F, gurobi=TRUE), "gWalks"))
    ## expect_true(inherits(subg.gw <<- subg$walk(gurobi=TRUE), "gWalks"))
})
