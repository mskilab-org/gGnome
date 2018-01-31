library(devtools,quietly = T)
library(gUtils,quietly = T)
library(gTrack,quietly = T) ## what's wrong with gTrack???
source("~/links/gTrack.R")

devtools::load_all("~/git/gGnome/", quiet=TRUE)

## get a 1Mb gr
set.seed(123)
gr = gUtils::gr.rand(2e5, GENOME)
strand(gr) = "+"

## equal length 5 parts
segs = rep(gr, 5)
for (i in 2:5){
    segs[i] = segs[1] %+% ((i-1)*2e5)
}
cat("These are the five segments we created.")
segs

win = streduce(segs+5e4)

## create the correct paths resulting from this evolutionary history
## ABCDE
## ABCDBCDE
## ABCDBDBCDE
paths = list(c(1,2,3,4,2,4,2,3,4,5),
             c(1,2,3,4,5))

## step 1: create the contig as gWalks
grl = GRangesList(lapply(paths,
             function(x) segs[x]))
gw = as(grl, "gWalks")
cat('This is the gWalks object we created:')
print(gw)
plot(gw)

## step 2: reduce down to graph
gg = as(gw, "bGraph")
cat('This is the bGraph object we created:')
gg
plot(gg)

## step 3: walk the graph again
gw2 = gg$walk2(verbose=FALSE, grl=FALSE)

gw.td = gw$td
gw.td$name = "i:gWalks"
gg.td = gg$td
gg.td$name = "gGraph"
gw2.td = gw2$td
gw2.td$name = "o:gWalks"

juncs = gg$junctions

gTrack::plot(c(gw.td, gg.td, gw2.td), win, links=juncs)

## step 4 (optional)
gw.simp = gw$simplify()
gg.simp = as(gw.simp, "bGraph")
gw.simp2 = gg.simp$walk2(F, F)

gw.simp.td = gw.simp$td; gw.simp.td$name = "i:simple gw"
gg.simp.td = gg.simp$td; gg.simp.td$name = "phased gg"
gw.simp2.td = gw.simp2$td; gw.simp2.td$name = "o:simple gw"

plot(c(gw.simp.td, gg.simp.td, gw.simp2.td), win, links=juncs)

## the data comes with the package
wv = system.file("extdata", "weaver", package="gGnome")
pg = system.file("extdata", "intervalFile.results", package="gGnome")
jab = system.file("extdata", "jabba.simple.rds", package="gGnome")

message("'gread' will determine what kind of data is this.")
j = gread(jab)
w = gread(wv)
p = gread(pg)

j.td = j$td; j.td$name = "JaBbA"
w.td = w$td; w.td$name = "Weaver"
p.td = p$td; p.td$name = "PREGO"

## plot one chr
## validated by coverage data
plot(c(p.td, w.td, j.td), "11:5e7-8e7")

##donotrun({
##    j$json("gGraph.json")
##    j.gw = j$walk2(verbose=F, grl=F)
##    j.gw$gw2json()
##})

cool.win = GRanges("11:68793238-70619925")
cool.hood = j$hood(cool.win, 1e6)
plot(cool.hood)

cool.hood = as(cool.hood, "bGraph")
cool.walk = cool.hood$walk2(F,F)
cool.walk.td = cool.walk$td; cool.walk.td$name = "gwalks"
cool.hood.td = cool.hood$td; cool.hood.td$name = "CN"
cool.win.to.see = cool.walk$window(pad=2e5)
plot(c(cool.hood.td, cool.walk.td), cool.win.to.see)


bfb.segs = gw$segstats
## hb = hydrogenBonds(segs = bfb.segs)
## hb.map = hb[, c(setNames(from, to), setNames(to, from))]
## seg.name = ifelse(strand(bfb.segs)=="+",
##                   as.character(seq_along(bfb.segs)),
##                   paste0("-", hb.map[as.character(seq_along(bfb.segs))]))
load_all()
names(bfb.segs) = tile.name(bfb.segs)
bfb.paths = list(c('1', '2', '3', '4', '5','-5', '-4', '-3',
              '3', '4', '-4', '-3', '3', '-3', '3', '4',
              '-4', '-3', '3', '4', '5', '-5','-4', '-3', '-2'),
              1:5)
bfb.grl = GRangesList(lapply(bfb.paths, function(x) bfb.segs[x]))
bfb.gw = as(bfb.grl, "gWalks")
bfb.gg = as(bfb.gw, "bGraph")
bfb.gw2 = bfb.gg$walk()

fbf(bfb.gg)
