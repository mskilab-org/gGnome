library(gUtils)
library(gTrack)

library(devtools)
load_all("~/gitLcl/gGnome/")

## ============= Example 1a: nested tandem duplication ============= ##
## get a 1Mb gr
set.seed(123)
gr = gr.rand(2e5, GENOME)
strand(gr) = "+"

## equal length 5 parts
segs = rep(gr, 5)
for (i in 2:5){
    segs[i] = segs[1] %+% ((i-1)*2e5)
}
win = streduce(segs+5e4)

## create the correct paths resulting from this evolutionary history
## ABCDE
## ABCDBCDE
## ABCDBDBCDE
paths = list(c(1,2,3,4,2,4,2,3,4,5),
             c(1,2,3,4,5))

load_all()
## the corresponding gWalks
grl = GRangesList(lapply(paths,
             function(x) segs[x]))
gw = as(grl, "gWalks")
gg = as(gw, "bGraph")
gw2 = gg$walk2(F, F)

gw.td = gw$td
gw.td$name = "i:gWalks"
gg.td = gg$td
gg.td$name = "gGraph"
gw2.td = gw2$td
gw2.td$name = "o:gWalks"

gTrack::plot(c(gw.td, gg.td, gw2.td), win, links=juncs)

## ============= Example 1b: classic BFB cycle ============= ##
## LATER

## ============= Example 2: analyzing one case with three tools  ============= ##
wv = system.file("extdata", "weaver", package="gGnome")
pg = system.file("extdata", "intervalFile.results", package="gGnome")
jab = system.file("extdata", "jabba.simple.rds", package="gGnome")

load_all()

j = gread(jab)
w = gread(wv)
p = gread(pg)

j.td = j$td; j.td$name = "JaBbA"
w.td = w$td; w.td$name = "Weaver"
p.td = p$td; p.td$name = "PREGO"

cool.win = GRanges("11:68793238-70619925")
cool.hood = j$hood(cool.win, 1e6)
cool.hood = as(cool.hood, "bGraph")

cool.walk = cool.hood$walk2(F,F)
cool.walk.td = cool.walk$td; cool.walk.td$name = "gwalks"
## NOTE: solve the walk problem here!!!
cool.hood.td = cool.hood$td; cool.hood.td$name = "CN"
cool.win.to.see = cool.walk$window(pad=2e5)
plot(c(cool.hood.td, cool.walk.td), cool.win.to.see)

## ============= Example 3: heterogeneous graph features among pan-cancer cohort  ============= ##
