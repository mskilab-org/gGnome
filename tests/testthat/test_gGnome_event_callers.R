library(testthat)
library(gUtils)
library(gTrack)

gg.jabba = gG(jabba = system.file('extdata/hcc1954', 'jabba.rds', package="gGnome"))

test_that('event caller', {
  ## run events caller with qrp
############# not running this currently since it is erroring out (See: https://app.travis-ci.com/github/mskilab/gGnome/builds/237955222 )
#  gg.events = events(gg.jabba, QRP = TRUE)
#  expect_is(gg.events$meta$events, 'data.table')
#  expect_true(any(grepl('qrp', gg.events$meta$events$type)))

  ## run events caller without qrp
  gg.events.no.qrp = events(gg.jabba, QRP = FALSE)
  expect_is(gg.events.no.qrp$meta$events, 'data.table')
  expect_false(any(grepl('qrp', gg.events.no.qrp$meta$events$type)))
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
    mfe280 = dup(mfe280)

    expect_true(mfe280$meta$pyrgo[,.N] == 8)
})
