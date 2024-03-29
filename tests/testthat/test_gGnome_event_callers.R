library(testthat)
library(gUtils)
library(gTrack)

gg.jabba = gG(jabba = system.file('extdata/hcc1954', 'jabba.rds', package="gGnome"))

test_that('event caller', {
  ## run events caller with qrp
############# not running this currently since it is erroring out (See: https://app.travis-ci.com/github/mskilab/gGnome/builds/237955222 )
  gg.events = events(gg.jabba, QRP = TRUE)
  expect_is(gg.events$meta$events, 'data.table')
  expect_true(any(grepl('qrp', gg.events$meta$events$type)))

  ## run events caller without qrp
  gg.events.no.qrp = events(gg.jabba, QRP = FALSE)
  expect_is(gg.events.no.qrp$meta$events, 'data.table')
  expect_false(any(grepl('qrp', gg.events.no.qrp$meta$events$type)))

  # test events.to.gr
  egr = events.to.gr(gg.events)
  expect_error(events.to.gr(gg.jabba))
  expect_error(events.to.gr(list()))
  empt = gG()
  empt$set(events = data.table())
  expect_warning(events.to.gr(empt))
})

## test gGnome amp function
test_that('amp caller', {
    setDTthreads(1)

    ## test amp without NOS
    gg.amp.vanilla = amp(gg = gg.jabba, mark.nos = FALSE)
    expect_false("nos" %in% colnames(gg.amp.vanilla$nodes$dt)) ## no amps should be foundb
    expect_false("nos" %in% gg.amp.vanilla$meta$amp[, type])

    ## test amp with NOS
    gg.amp.nos = amp(gg = gg.jabba, mark.nos = TRUE, min.nodes = 0, min.jun = 0) ## don't filter...?
    expect_true("nos" %in% colnames(gg.amp.nos$nodes$dt))
    expect_true("nos" %in% gg.amp.nos$meta$amp[, type])
})
    
test_that('tic', {
    setDTthreads(1)
    ccle = dir(system.file("extdata", package = "gGnome"), ".+jabba.simple.rds", full = TRUE)
    names(ccle) = gsub(".*gGnome/.*extdata/(.*)\\.jabba\\.simple\\.rds$", "\\1", ccle)
    jhos2 = gG(jabba = ccle["JHOS_2"])
    jhos2 = tic(jhos2)

    expect_true('tic' %in% names(jhos2$meta))
})


test_that('chromothripsis', {
    setDTthreads(1)
    ccle = dir(system.file("extdata", package = "gGnome"), ".+jabba.simple.rds", full = TRUE)
    names(ccle) = gsub(".*gGnome/.*extdata/(.*)\\.jabba\\.simple\\.rds$", "\\1", ccle)
    h2081 = gG(jabba = ccle["NCI_H2081"])
    h2081 = chromothripsis(h2081)

    expect_true('chromothripsis' %in% names(h2081$meta))
})

test_that('dup/pyrgo', {
    setDTthreads(1)
    ccle = dir(system.file("extdata", package = "gGnome"), ".+jabba.simple.rds", full = TRUE)
    names(ccle) = gsub(".*gGnome/.*extdata/(.*)\\.jabba\\.simple\\.rds$", "\\1", ccle)
    mfe280 = gG(jabba = ccle["MFE_280"])
    mfe280_dup = dup(mfe280)
    expect_true(mfe280_dup$meta$pyrgo[,.N] == 8)
})

test_that('microhomology', {
    
  # make simple graph
  nodes1 = c(GRanges("1",IRanges(1,100),"*"), GRanges("1",IRanges(101,200),"*"),
             GRanges("1",IRanges(201,300),"*"), GRanges("1",IRanges(301,400),"*"),
             GRanges("1",IRanges(401,500),"*"))
  edges = data.table(n1 = c(3,2,4,1,3), n2 = c(3,4,2,5,4), n1.side = c(1,1,0,0,1), n2.side = c(0,0,0,1,0))
  gg = gGraph$new(nodes = nodes1, edges = edges)    

  # get short fasta
  fa = system.file("extdata", 'microhomology.test.fasta', package = "gGnome")
  m = microhomology(gg, fa)

  # make sure mh annotations were added
  expect_true(length(setdiff(c('mh5', 'mh10', 'mh50', 'mh100'), names(m$edges$dt))) == 0)
  expect_true(is.numeric(m$edges[type == 'ALT']$dt$mh100) & all(!is.na(m$edges[type == 'ALT']$dt$mh100)))

  # test with Junction input
  m = microhomology(gg$junctions, fa)
  # make sure mh annotations were added
  expect_true(length(setdiff(c('mh5', 'mh10', 'mh50', 'mh100'), names(m$dt))) == 0)
  expect_true(is.numeric(m[type == 'ALT']$dt$mh100) & all(!is.na(m[type == 'ALT']$dt$mh100)))

  expect_error(microhomology(gg$nodes, fa))
  expect_true(length(microhomology(jJ(), fa)) == 0)

  expect_true(length(microhomology(jJ(), fa)) == 0)
  gg_ref_only = gGraph$new(nodes = nodes1, edges = edges[5])
  m = microhomology(gg_ref_only, fa)
  expect_true(all(is.na(m$edges$mh100)))

  # bad seqnames - seqnames don't match between gGraph and FASTA
  nodes1 = c(GRanges("2",IRanges(1,100),"*"), GRanges("2",IRanges(101,200),"*"),
             GRanges("2",IRanges(201,300),"*"), GRanges("2",IRanges(301,400),"*"),
             GRanges("2",IRanges(401,500),"*"))
  edges = data.table(n1 = c(3,2,4,1,3), n2 = c(3,4,2,5,4), n1.side = c(1,1,0,0,1), n2.side = c(0,0,0,1,0))
  gg = gGraph$new(nodes = nodes1, edges = edges)    
  expect_error(microhomology(gg, fa))

  ## test microhomology with prefix only
  fa = rtracklayer::import(system.file("extdata", "microhomology", "ex1", "ref.fasta", package = "gGnome"), format = "fasta")
  jj = readRDS(system.file("extdata", "microhomology", "ex1", "jj.grl.rds", package = "gGnome"))
  gg = gG(junctions = jj)

  hgg = gGnome::microhomology(gg,
                              hg = fa,
                              prefix_only = TRUE,
                              pad = c(30, 40, 50),
                              ignore_missing = TRUE)
  
  expect_true(hgg$edges$dt[type == "ALT", mh30 == 3])
  expect_true(hgg$edges$dt[type == "ALT", mh40 == 3])
  expect_true(hgg$edges$dt[type == "ALT", mh50 == 3])

  ## ensure that this still works with junction input
  hjj = gGnome::microhomology(jJ(jj),
                              hg = fa,
                              prefix_only = TRUE,
                              pad = c(30, 40, 50),
                              ignore_missing = TRUE)
  expect_true(hjj$dt[, mh30 == 3])
  expect_true(hjj$dt[, mh40 == 3])
  expect_true(hjj$dt[, mh50 == 3])


  fa = rtracklayer::import(system.file("extdata", "microhomology", "ex2", "ref.fasta", package = "gGnome"), format = "fasta")
  jj = readRDS(system.file("extdata", "microhomology", "ex2", "jj.grl.rds", package = "gGnome"))
  
  gg = gG(junctions = jj)
  hgg = gGnome::microhomology(gg,
                              hg = fa,
                              prefix_only = TRUE,
                              pad = c(20, 40, 60),
                              ignore_missing = TRUE)
  expect_true(hgg$edges$dt[type == "ALT", mh20 == 3])
  expect_true(hgg$edges$dt[type == "ALT", mh40 == 3])
  expect_true(hgg$edges$dt[type == "ALT", mh60 == 3])
  
})
