library(testthat)
library(gUtils)
library(gTrack)

gg.jabba = gG(jabba = system.file('extdata/hcc1954', 'jabba.rds', package="gGnome"))

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
    

