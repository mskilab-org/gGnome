## cheating:
library(gTrack)

setwd("~/gitLcl/gGnome"); load_all()

jab = readRDS("inst/extdata/jabba.simple.rds")
win = readRDS("inst/extdata/win_17.21.rds")

chr_17.21 = si2gr(GENOME)[c(17,21)]

## create g1 from HCC1143 JaBbA result
g1=gGraph$new(jabba = jab)
g1
comp1 = components(g1)
g1$parts

## select the second component
g1.2=comp1[[2]]
plot(g1.2)

## walk this easy subgraph
g1.2$isJunctionBalanced()
b1.2 = bGraph$new(g1.2)
b1.2
w1.2 = b1.2$walk()

## pick a neighborhood
g1win = g1$trim(win)
g1hood = g1$hood(win, d=1e5, pad=1e3)

## see instantly
plot(g1win)
plot(g1hood)

## dist and prox still need more design
g1$dist(win[1], win[2])
prox1 = g1$proximity(query, subject)

## recreate classic Campbell example
r1 = gr.rand(w=5e6, GENOME)
s5 =
