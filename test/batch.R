setwd("~/gitLcl/gGnome"); load_all()

jab = readRDS("data/jabba/jabba.simple.rds")
win = readRDS("data/win_17.21.rds")

chr_17.21 = si2gr(GENOME)[c(17,21)]

g1=gGraph$new(jabba = jab)
g1
comp1 = g1$components()
g2=comp1[[2]]

g1win = g1$trim(win)
g1hood = g1$hood(win, d=1e5, pad=1e3)

g1$dist(win[1], win[2])
prox1 = g1$proximity(query, subject)

g0 = gGraph$new(tile=win, junctions = g1$junctions[1:5])
