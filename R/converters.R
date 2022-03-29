# assigning this operator since sometimes R tries to use the ggplot2 operator instead of the gUtils one (depending on the order of libraries loaded in a session)
`%+%` = gUtils::`%+%`

#' @name breakgraph
#' @title breakgraph
#'
#' @description
#' Builds a gGraph by breaking the reference genome at the points specified in tile
#' Treats tile as nothing but breakpoints, no metadata, nothing special
#' Juncs can be either a GRangesList of junctions in proper format or a Junction Object
#' If tile has length 0, this function creates a simple graph (1 node per chromosome, each one full length)
#' If genome is specified, it will try to use that genome instead - if it is invalid it wil use the default#' 
#' unlists a list of vectors, matrices, data.tables into a data.table indexed by the list id
#'
#' @param tile GRanges of tiles
#' @param juncs Junction object or grl coercible to Junctions object
#' @param genome seqinfo or seqlengths
#' @return list with gr and edges which can be input into standard gGnome constructor
#' @author Marcin Imielinski, Joe DeRose, Xiaotong Yao
#' @keywords internal
#' @noRd 
breakgraph = function(breaks = NULL,
                      juncs = NULL,
                      genome = NULL)
{  
  ## Make sure user entered some input
  if (is.null(breaks) & is.null(juncs)) {
    stop("Cannot have both breaks and juncs be NULL, must be some input")
  }

  ## If juncs is not Junction Object, try to convert it
  if(!is.null(juncs) && !inherits(juncs, "Junction")) {
    juncs = tryCatch(Junction$new(juncs),
                     error = function(e) {
                       NULL
                     })
    if (is.null(juncs)) {
      stop("Input is not a valid junctions set.")
    }
  }

  ## Validate breaks as GRanges
  if (!is.null(breaks) && !inherits(breaks, "GRanges")){
    breaks = tryCatch(GRanges(breaks),
                    error=function(e){
                      NULL
                    })
    if (is.null(breaks)){
      stop("Input cannot be converted into a GRanges object.")
    }
  }

  ## If the user provided a genome, check if it is valid. If it isn't, set to NULL
  if (!is.null(genome)) {
    tmp = tryCatch(si2gr(genome), error=function(e) NULL)
    if (is.null(tmp)) {
      genome = NULL
      warning('Provided genome not coercible to seqinfo')
      }
    else
    {
      genome = seqinfo(tmp)
    }
  }

  if (is.null(genome))
  {
    if (!is.null(breaks)){
        genome = seqinfo(breaks)
        }
    else if (!is.null(juncs)){
        genome = seqinfo(juncs$grl)
        }
    else{
      stop('Provided genome not coercible to Seqinfo object and no breaks or junc provided.  Must provide either breaks, junc, or genome argument to this constructor')
    }
    }


  ## Build a GRanges from the default genome
  ## FIXME: if using misc chr, definitely need to specify genome                           
  nodes = gr.stripstrand(si2gr(genome))
  edges = data.table(n1 = numeric(0), n2 = numeric(0), n1.side = numeric(0), n2.side = numeric(0), type = character(0))
  
  ## Break the genome based on whether there is breaks or juncs
  if (!is.null(breaks) && length(breaks) > 0) {
    tmp.br = gr.stripstrand(breaks)[, c()]
    tmp.br$break.id = 1:length(tmp.br)
    nodes = gr.disjoin(grbind(tmp.br, nodes))
  }

  if (!is.null(juncs) && length(juncs) > 0) {
    ## collapse node with positive and  
    juncsGR = gr.start(gr.fix(grl.unlist(juncs$grl), nodes), ignore.strand = FALSE) ## grl.ix keeps track of junction id
    nodes = gr.fix(nodes, juncsGR)

    ## keep track of nmeta to paste back later
    nodes$nid = 1:length(nodes)
    nmeta = as.data.table(values(nodes))

    ## disjoin nodes and bps
    bps = juncsGR[, c('grl.ix')] ## trim bps to first base
    dnodes = sort(gr.disjoin(grbind(nodes, gr.stripstrand(bps))))
    bpov = dnodes %*% bps ## merge back with bps to figure out where to trim dnodes
    strand(bpov) = strand(bps)[bpov$subject.id]

    ## mark whether bps have a negative strand or positive strand junction attaching to them
    ## (each bp should have at least one has.neg or has.pos = TRUE)
    has = gr2dt(bpov)[, .(has.pos = any(strand=='+'), has.neg = any(strand=='-')), keyby = query.id][.(1:length(dnodes)), ]

    ## merge this info back to dnodes, not that non bp intervals will have has.neg and has.pos = NA
    dnodes = merge(gr2dt(dnodes)[, query.id := 1:.N], has, by = 'query.id', all = TRUE)
    setkey(dnodes, query.id)

    dnodes[, is.bp := !is.na(has.neg)]
    dnodes[, has.neg := ifelse(is.na(has.neg), FALSE, has.neg)]
    dnodes[, has.pos := ifelse(is.na(has.pos), FALSE, has.pos)]

    ## merge dnodes with next interval if that interval
    ## precedes a has.pos = FALSE
    ## merge dnodes with previous interval if that interval
    ## follows a has.neg = FALSE
    ## the following grouping will encode this principle
    ## i.e. we always increment the group counter from node to node except when
    ## current node has.pos = FALSE and previous node has.neg = FALSE
    ## and at least one of the nodes is a bp node

    dnodes[, increment := !((is.bp | c(FALSE, is.bp[-.N])) & c(TRUE, !has.neg[-.N]) & !has.pos), by = seqnames]
    dnodes[, group := cumsum(sign(increment)), by = seqnames]

    ## collapse dnodes with same group in same seqnames
    ## but don't merge across nids
    nodes = dt2gr(dnodes[, .(start = start[1], end = end[.N]), by = .(seqnames, group, nid)], seqlengths = seqlengths(nodes))

    ## paste back metadata
    values(nodes) = nmeta[nodes$nid, ]

    nodes.left = gr.start(nodes); nodes.left$side = 'left'; 
    nodes.right = gr.end(nodes); nodes.right$side = 'right'
    nodes.right$nid = nodes.left$nid = 1:length(nodes)
    node.ends = unique(grbind(nodes.left, nodes.right))
    ov = gr2dt(juncsGR[, c('grl.ix', 'grl.iix')] %*% node.ends)[order(grl.iix), ]

    ## Map the Junctions to an edge table
    ## Mapping is as follows in terms of junctions a -> b
    ##         a- => leaves/enters left side of base a
    ##         a+ => leaves/enters right side of base a
    ##  a- <-> a+ => leaves left side of base a, enters right side of base a (NOT THE BASE NEXT TO a)
    ov[, strand := as.character(strand(juncsGR))[query.id]]
    setkeyv(ov, c('grl.ix', 'grl.iix'))

    ### ov should have exactly one row per grl.ix grl.iix combo
    new.edges = 
      merge(
        ov[.(1:length(juncs), rep(1, length(juncs))),
           .(grl.ix = 1:length(juncs),
             n1 = node.ends$nid[subject.id],
             n1.side = ifelse(strand == '-', 1, 0))],
        ov[.(1:length(juncs), rep(2, length(juncs))),
           .(
             grl.ix = 1:length(juncs),
             n2 = node.ends$nid[subject.id],
             n2.side = ifelse(strand == '-', 1, 0))],
        by = 'grl.ix'
      )[, type := 'ALT']

    ## new.edges = ov[, .(
    ##   n1 = node.ends$nid[subject.id[1]],
    ##   n2 = node.ends$nid[subject.id[2]],
    ##   n1.side = ifelse(strand[1] == '-', 1, 0),
    ##   n2.side = ifelse(strand[2] == "-", 1, 0),
    ##   type = "ALT"
    ##   ), keyby = grl.ix]

    if (length(setdiff(1:length(juncs), new.edges$grl.ix))>0)
      stop('Error in breakgraph generation - some junctions failed to be incorporated')

    ## Remove junctions that aren't in the genome selected
    ## this will have n1 or n2 NA
    new.edges = new.edges[!(is.na(n1) | is.na(n2)), ]

    ## reconcile new.edges with metadata
    if (!is.null(juncsGR))
      if (ncol(values(juncs$grl))>0)
        new.edges = cbind(new.edges, as.data.table(values(juncs$grl)[new.edges$grl.ix, , drop = FALSE]))
      
      edges = rbind(edges, new.edges, fill = TRUE)
  }
    
  ## now make reference edges
  nodesdt = as.data.table(nodes[, c()])
  nodesdt[, node.id := 1:.N]
  nodesdt[, endnext := end +1]

  ref.edges = merge(nodesdt, nodesdt, by.x = c('seqnames', 'endnext'),
                    by.y = c('seqnames', 'start'),
                    allow.cartesian = TRUE)
  if (nrow(ref.edges)>0)
  {

    ref.edges = ref.edges[, .(n1 = node.id.x, n1.side = 1, n2 = node.id.y, n2.side = 0)]
    ref.edges$type = "REF"
    edges = rbind(edges, ref.edges, fill = TRUE)

  }
  nodes$qid = NULL
  names(nodes) = NULL

  ## let nodes inherit break metadata using results of disjoin above
  if (!is.null(breaks) && !is.null(nodes$break.id)){
    if (ncol(values(breaks))>0){
      values(nodes) = values(breaks)[nodes$break.id, ]
    }
  }
  NONO.FIELDS = c('from', 'to')
  if (any(nono.ix <- names(edges) %in% NONO.FIELDS))
  {
    warning(paste('removing reserved edge metadata fields from edges:', paste(NONO.FIELDS, collapse = ',')))
    edges = edges[, !nono.ix, with = FALSE]
  }

  NONO.FIELDS = c('node.id', 'snode.id', 'index')
  if (any(nono.ix <- names(values(nodes)) %in% NONO.FIELDS))
  {
    warning(paste('removing reserved edge metadata fields from nodes:', paste(NONO.FIELDS, collapse = ',')))
    nodes = nodes[, !nono.ix]
  }

  return(list(nodes = nodes, edges = edges))
}


#' @name rck2gg
#' @title rck2gg
#'
#' @description
#' Constructor for creating gGraph from RCK output file
#'
#' @param rck.dirname (character) directory name containing RCK outputs
#' @param simplify (logical) merge adjacent regions with same total CN? default FALSE
#' @param haploid (logical) create total CN (unphased) graph? default TRUE. FALSE NOT IMPLEMENTED YET.
#' @param prefix (character) prefix in the RCK output files (namely {prefix}rck.scnt.tsv {prefix}rck.acnt.tsv)
#' 
#' @return list of gr and edges that can be input into standard gGraph constructor
#' @author Marcin Imielinski, Zi-Ning Choo, Xiaotong Yao
#' @keywords internal
#' @noRd
rck2gg = function(rck.dirname, haploid = TRUE, simplify = TRUE, prefix = '')
{
    if (!dir.exists(rck.dirname)) {
        stop("Input RCK directory not found")
    }
    scnt.fname = file.path(rck.dirname, paste0(prefix, "rck.scnt.tsv"))
    acnt.fname = file.path(rck.dirname, paste0(prefix, "rck.acnt.tsv"))
    if (!file.exists(scnt.fname) | !file.exists(acnt.fname)) {
        stop("Required output files rck.scnt.tsv and rck.acnt.tsv cannot be located")
    }

    ## read segment copy numbers
    segs.dt = fread(scnt.fname)

    ## extract total CN from allelic CN if haploid == TRUE
    seg.ptn = "cn=\\{'c1': \\{'A': ([0-9]+), 'B': ([0-9]+)\\}\\}"
    if (all( grep(seg.ptn, segs.dt$extra[1], value = FALSE) == 0)) {
        stop("rck.scnt.tsv not properly formatted")
    }

    segs.dt[, A := gsub(seg.ptn, "\\1", extra) %>% as.numeric]
    segs.dt[, B := gsub(seg.ptn, "\\2", extra) %>% as.numeric]
    segs.dt[, total := A + B]

    ## read adjacency copy numbers
    adjs.dt = fread(acnt.fname)

    ## check valid extra field
    adjs.ptn = "aid=\\w+;cn=\\{'c1': \\{'AA': ([0-9]+), 'AB': ([0-9]+), 'BA': ([0-9]+), 'BB': ([0-9]+)}};at=\\w+"
    if (all( grep(adjs.ptn, adjs.dt$extra[1], value = FALSE) == 0)) {
        stop("rck.acnt.tsv not properly formatted")
    }
    
    adjs.dt[, ":="(AA = gsub(adjs.ptn, "\\1", extra) %>% as.numeric, ## CN of junction endpoints
                   AB = gsub(adjs.ptn, "\\2", extra) %>% as.numeric,
                   BA = gsub(adjs.ptn, "\\3", extra) %>% as.numeric,
                   BB = gsub(adjs.ptn, "\\4", extra) %>% as.numeric)]

    ## total CN
    adjs.dt[, total := AA + AB + BA + BB]

    ## get n1 and n2 sides
    adjs.dt[, n1.side := ifelse(strand1 == "+", "right", "left")]
    adjs.dt[, n2.side := ifelse(strand2 == "+", "right", "left")]

    ## only ALT junctions with nonzero total CN
    ## alt.dt = adjs.dt[grep("^[0-9]+$", aid, value = FALSE),][total > 0,] ## if not start with R
    adjs.dt = adjs.dt[grepl("^[0-9]+$", aid) | (total > 0)]
    adjs.dt[, type := ifelse(grepl("^[0-9]+$", aid), "ALT", "REF")]

    if (haploid) {
        nodes.gr = dt2gr(segs.dt[, .(seqnames = chr, start, end, cn = A + B)])

        ## prepare edge data table
        edges.dt = adjs.dt[, .(chr1, coord1, chr2, coord2, n1.side, n2.side, cn = AA + AB + BA + BB, type)]

        ## n1 coordinates
        edges.n1 = GRanges(seqnames = edges.dt$chr1, ranges = IRanges(start = edges.dt$coord1, width = 1))

        ## n2 coordinates
        edges.n2 = GRanges(seqnames = edges.dt$chr2, ranges = IRanges(start = edges.dt$coord2, width = 1))

        ## add corresponding nodes
        n1.mt = gr.match(edges.n1, nodes.gr)
        n2.mt = gr.match(edges.n2, nodes.gr)

        edges.dt[, n1 := n1.mt]
        edges.dt[, n2 := n2.mt]

        nodes.gr$ywid = 0.8
    } else {
        nodes.gr = dt2gr(
            rbind(segs.dt[, .(seqnames = chr, start, end, cn = A, total = A + B, haplotype = "A")],
                  segs.dt[, .(seqnames = chr, start, end, cn = B, total = A + B, haplotype = "B")])
        )

        ## prepare edge data table
        edges.dt = rbind(
            adjs.dt[, .(chr1, coord1, chr2, coord2, n1.side, n2.side, type,
                        n1.haplotype = "A", n2.haplotype = "A", cn = AA, total)],
            adjs.dt[, .(chr1, coord1, chr2, coord2, n1.side, n2.side, type,
                        n1.haplotype = "A", n2.haplotype = "B", cn = AB, total)],
            adjs.dt[, .(chr1, coord1, chr2, coord2, n1.side, n2.side, type,
                        n1.haplotype = "B", n2.haplotype = "A", cn = BA, total)],
            adjs.dt[, .(chr1, coord1, chr2, coord2, n1.side, n2.side, type,
                        n1.haplotype = "B", n2.haplotype = "B", cn = BB, total)]
        )

        ## n1 coordinates
        edges.n1 = GRanges(seqnames = edges.dt$chr1,
                           ranges = IRanges(start = edges.dt$coord1, width = 1),
                           haplotype = edges.dt$n1.haplotype)

        ## n2 coordinates
        edges.n2 = GRanges(seqnames = edges.dt$chr2,
                           ranges = IRanges(start = edges.dt$coord2, width = 1),
                           haplotype = edges.dt$n2.haplotype)

        ## add corresponding nodes
        n1.mt = gr.match(edges.n1, nodes.gr, by = "haplotype")
        n2.mt = gr.match(edges.n2, nodes.gr, by = "haplotype")

        edges.dt[, n1 := n1.mt]
        edges.dt[, n2 := n2.mt]

        ## formatting
        nodes.gr$col = ifelse(nodes.gr$haplotype == "A", alpha("red", 0.5), alpha("blue", 0.5))
        nodes.gr$ywid = 0.8

        edges.dt[cn == 0, col := alpha("gray", 0.01)]
    }
    if (simplify) {
        edges.dt = edges.dt[cn > 0]
        nodes.gr = inferLoose(nodes.gr, edges.dt)
    }
    return(list(nodes = nodes.gr, edges = edges.dt))
}

    
#' @name pr2gg
#' @title pr2gg
#'
#' @description
#' ## Generates a gGraph from a prego output file
#'
#' @return list of gr and edges that can be input into standard gGraph constructor
#' @author Marcin Imielinski, Joe DeRose, Xiaotong Yao
#' @keywords internal
#' @noRd 
pr2gg = function(fn, simplify = TRUE)
{

  if (file.info(fn)$isdir)
    fn = paste0(fn, '/.')
  
  res.tmp = readLines(paste(dirname(fn), 'intervalFile.results', sep = '/'))
  chrm.map.fn = paste(dirname(fn), "chrm.map.tsv", sep = '/')

  if (file.exists(chrm.map.fn)){
    chrm.map = fread(chrm.map.fn)[,setNames(V1, V2)]
  }

  res = structure(lapply(split(res.tmp, cumsum(grepl("edges", res.tmp))),
                           function(x) {
                           rd = read.delim(textConnection(x),
                                           strings = F,
                                           skip = 1,
                                           header = F,
                                           col.names = c("node1", "chr1",
                                                         "pos1", "node2",
                                                         "chr2", "pos2", "cn"))
                           if (exists("chrm.map")){
                             rd$chr1 = chrm.map[rd$chr1]
                             rd$chr2 = chrm.map[rd$chr2]
                           }
                           else {
                               rd = rd[which(rd$chr1 %in% as.character(1:24) &
                                             rd$chr2 %in% as.character(1:24)),]
                               rd$chr1 = gsub("24", "Y", gsub("23","X",rd$chr1))
                               rd$chr2 = gsub("24", "Y", gsub("23","X",rd$chr2))
                           }
                           
                           return(rd)
                           }),
                    names = gsub(":", "", grep("edges", res.tmp, value = T)))
    res[[1]]$tag = paste0(res[[1]]$node1, ":", res[[1]]$node2)


  ## turn into our nodes
  nodes = GRanges(res[[1]]$chr1,
                  IRanges(res[[1]]$pos1,
                          res[[1]]$pos2),
                  strand = "+",
                  cn = res[[1]]$cn,
                  left.tag = res[[1]]$node1,
                  right.tag = res[[1]]$node2)
  

  edges = rbind(as.data.table(res[[2]])[, type := 'REF'],
                as.data.table(res[[3]])[, type := 'ALT'])


  if (nrow(edges)>0)
    {
      edges[, n1.left := match(node1, nodes$left.tag)]
      edges[, n1.right := match(node1, nodes$right.tag)]
      edges[, n2.left := match(node2, nodes$left.tag)]
      edges[, n2.right := match(node2, nodes$right.tag)]
      
      edges[, n1 := ifelse(is.na(n1.left), n1.right, n1.left)]
      edges[, n1.side := ifelse(is.na(n1.left), 'right', 'left')]
      edges[, n2 := ifelse(is.na(n2.left), n2.right, n2.left)]
      edges[, n2.side := ifelse(is.na(n2.left), 'right', 'left')]
      
      if (simplify)
      {
        edges = edges[cn>0, ]
      }
    }

  nodes = inferLoose(nodes, edges)

  return(list(nodes = gr.fix(nodes), edges = edges))
}


#' @name jab2gg
#' @title jab2gg
#' @description
#' 
#' Constructor for creating a gGraph object from a jabba output ## Generates a gGraph from a prego output file
#'
#' @return list of gr and edges that can be input into standard gGraph constructor
#' @author Marcin Imielinski, Joe DeRose, Xiaotong Yao
#' @keywords internal
#' @noRd 
jab2gg = function(jabba)
{
  ## Validate our input

  if (is.list(jabba)) {
    if (!all(is.element(c("segstats", "adj",
                          "purity", "ploidy"),
                        names(jabba)))){
      stop("The input is not a JaBbA output.")
    }
  } else if (is.character(jabba) && grepl(".rds$", jabba)){
    if (file.exists(jabba)){
      jabba = readRDS(jabba)
    } else {
      stop("JaBbA file not found")
    }
  }
  else if (inherits(jabba, 'gGraph'))
  {
    return(list(nodes = jabba$nodes$gr, edges = jabba$edges$dt))
  }
  else {
    stop("Error loading jabba object from provided .rds path or object: please check input")
  }

  ## second round check .. just in case .rds file had gGraph object inside it
  if (inherits(jabba, 'gGraph'))
  {
      ## don't discard purity/ploidy metadata if included
      return(list(nodes = jabba$nodes$gr,
                  edges = jabba$edges$dt,
                  purity = jabba$meta$purity,
                  ploidy = jabba$meta$ploidy))
  }

  if (is.null(jabba$segstats$loose))
    jabba$segstats$loose = FALSE

  if (is.null(jabba$segstats$cn))
    jabba$segstats$cn = NA

  afields = c('cn', 'type', 'parent')
  if (!is.null(jabba$asegstats) && inherits(jabba$asegstats, 'GRanges') && length(jabba$asegstats) == 2 * length(jabba$segstats) && length(setdiff(afields, names(mcols(jabba$asegstats)))) == 0){
      snodes = jabba$segstats
      aseg.dt = gr2dt(jabba$asegstats[, afields])
      aseg.dt.dcast = dcast.data.table(aseg.dt, parent ~ type, value.var = 'cn')
      setkey(aseg.dt.dcast, 'parent')
      snodes$cn.low = aseg.dt.dcast$low
      snodes$cn.high = aseg.dt.dcast$high
      snodes = snodes %Q% (loose == FALSE)
  } else {
      snodes = jabba$segstats %Q% (loose == FALSE)
  }
      
  snodes$index = 1:length(snodes)
  snodes$snode.id = ifelse(as.logical(strand(snodes)=='+'), 1, -1) * gr.match(snodes, unique(gr.stripstrand(snodes)))

  if (length(snodes)==0)
    return(gG(genome = seqinfo(segs)))


  sedges = spmelt(jabba$adj[jabba$segstats$loose == FALSE, jabba$segstats$loose == FALSE])
  
  nodes = snodes %Q% (strand == '+')

  if (!is.null(nodes$eslack.in) & !is.null(nodes$eslack.out))
  {
    nodes$loose.left = nodes$eslack.in>0
    nodes$loose.right = nodes$eslack.out>0      
    nodes$loose.cn.left = nodes$eslack.in
    nodes$loose.cn.right = nodes$eslack.out
#    nodes$loose = nodes$loose.left | nodes$loose.right
  }

    ## don't do this o/w edgesdt won't have n1, n2, n1.side, n2.side
  ## if (nrow(sedges)==0)
  ##   gG(nodes = nodes)

  setnames(sedges, c('from', 'to', 'cn'))
  setkeyv(sedges, c('from', 'to'))
  sedges[, type := 'REF']

  if (nrow(jabba$ab.edges)>0)
  {
    ab.edges = as.data.table(rbind(jabba$ab.edges[, 1:2, '+'],
                                   jabba$ab.edges[, 1:2, '-']))
    ab.edges[, jid := rep(1:nrow(jabba$ab.edges),2)]
    ab.edges = ab.edges[!is.na(from) & !is.na(to), ]
    if (nrow(ab.edges)>0)
    {
      if (ncol(values(jabba$junctions))>0)
      {
        ab.edges = as.data.table(cbind(ab.edges, values(jabba$junctions)[ab.edges$jid, ]))
        if (!is.null(ab.edges$cn))
          ab.edges[, cn := NULL] ## confilct with the cn inferred from adj
      }

      sedges[.(ab.edges$from, ab.edges$to), type := 'ALT']
      ab.cols = c('from', 'to', setdiff(names(ab.edges), c('sedge.id', names(sedges))))
      sedges = merge(sedges, ab.edges[, ab.cols, with = FALSE], by = c('from', 'to'), all.x = TRUE)
    }
  }

  ## rescue any hom-del ref edge
  if (!is.null(jabba$edges$cn))
  {
    if (any(data.table(jabba$edges)[type=="reference", cn==0])){
      sedges = rbind(sedges,
                     data.table(jabba$edges)[cn==0 & type=="reference",
                                             .(from, to, type="REF", cn)],
                     fill=T)
    }
  }
  edges = convertEdges(snodes, sedges, metacols = TRUE)

  return(list(nodes = nodes,
              edges = edges,
              purity = jabba$purity,
              ploidy = jabba$ploidy
              ))
  ## return(list(nodes = nodes[, intersect(c('cn', 'loose.left', 'loose.right', 'loose.cn.left', 'loose.cn.right'), names(values(nodes)))],
  ##             edges = edges))
}


#' @name wv2gg
#' @title wv2gg
#' @description
#' 
#' Constructor for creating a gGraph object from Weaver output
#'
#' @param weaver directory containing SV_CN_PHASE and REGION_CN_PHASE files
#' @return list of gr and edges that can be input into standard gGraph constructor
#' @author Marcin Imielinski, Joe DeRose, Xiaotong Yao
#' @keywords internal
#' @noRd 
wv2gg = function(weaver, simplify = TRUE)
{
  if (!file.info(weaver)$isdir)
    weaver = dirname(weaver)

  if (!dir.exists(weaver)){
    stop("Error: Invalid input weaver directory!")
  }
  
  if (!all(is.element(c("SV_CN_PHASE", "REGION_CN_PHASE"), dir(weaver))) ){
    stop('Error: Need "SV_CN_PHASE" and "REGION_CN_PHASE".')
  }
  
  region = data.table(read.delim(
    paste(weaver, "REGION_CN_PHASE", sep="/"),
    header = FALSE, sep = "\t"))
  
  sv.fn = paste(weaver, "SV_CN_PHASE", sep="/")
  if (file.size(sv.fn)>0){
    sv = data.table(read.delim(sv.fn, header = FALSE, sep = "\t"))
    names(sv) = c("chr1", "pos1", "side1", "allele1",
                  "chr2", "pos2", "side2", "allele2",
                  "cn", "unknown1", "unknown2", "timing", "class")[1:ncol(sv)]
  }
  else {
    sv = NULL
  }
  
  ## define the columns
  names(region) = c("seqnames", "start", "end", "acn", "bcn")
  region[, cn := acn + bcn]
  ## names(snp) = c("seqnames", "pos", "ref", "alt", "acn", "bcn")
  ## Xiaotong remove this on Apr 23
  ## don't need to do anything to the segment locations
  ## region$start = region$start+1 ## start coordinates appear to be 0 centric
  ## region$end = region$end+2 ## end coordinates appear to be "left-centric"
  ss = dt2gr(region)
  ss = gr.fix(ss)
  
  ## get junctions
  ## ALERT: in the file, +/- means right/left end of a segment
  ## exactly reverse of what we define a junction
  strmap = setNames(c("+", "-"), c("-", "+"))
  ## sv.select = sv[!is.na(allele1) & !is.na(allele2)]
  junc = NULL
  if (!is.null(sv)){
    sv = sv[which(allele1 !=0 & allele2 !=0), ]
    if (simplify)
      {
        sv = sv[cn>0, ]
      }

    if (nrow(sv)>0)
    {
        bps = grbind(
              dt2gr(
                  sv[, .(seqnames = chr1,
                         start = pos1,
                         end = pos1,
                         jix=.I, ii = 1,
                         strand = strmap[side1])],
                  seqlengths = seqlengths(ss)),
              dt2gr(
                  sv[, .(seqnames = chr2,
                         start = pos2,
                         end = pos2,
                         jix=.I, ii = 2,
                         strand = strmap[side2])],
                  seqlengths = seqlengths(ss))
        )

        ## Xiaotong remove this on Apr 23
        ## don't need to nudge anymore since we keep the segment faithful
        ## bps = grbind(
        ##   dt2gr(
        ##     sv[, .(seqnames = chr1,
        ##            start = ifelse(side1=="-", pos1+1, pos1+2),
        ##            end = ifelse(side1=="-", pos1+1, pos1+2),
        ##            jix=.I, ii = 1,
        ##            strand = strmap[side1])], seqlengths = seqlengths(ss)),
        ##   dt2gr(
        ##     sv[, .(seqnames = chr2,
        ##            start = ifelse(side2=="-", pos2+1, pos2+2),
        ##            end = ifelse(side2=="-", pos2+1, pos2+2),
        ##            jix=.I, ii = 2,
        ##            strand = strmap[side2])], seqlengths = seqlengths(ss))
        ## )


        ## sanity check, all raw.bp at this point should
        ## locate at left/right boundary of segements
        ss.ends = c(gr.start(ss), gr.end(ss))
        if (any(!bps %^% ss.ends)){
          warning("Eligible SVs not matching segment ends!")
        }

        ## create junctions
        junc = grl.pivot(split(bps, bps$ii))
        toget = intersect(c("allele1", "allele2", "cn", "unknown1", "unknown2", "timing", "class"), colnames(sv))
        values(junc) = sv[, toget, with=F]
      }
  }

  return(breakgraph(breaks = ss, juncs = junc))
}

#' @name remixt2gg
#' @title remixt2gg
#' @description
#' 
#' Constructor for creating a gGraph object from remiXT output
#'
#' @return list of gr and edges that can be input into standard gGraph constructor
#' @author Marcin Imielinski, Joe DeRose, Xiaotong Yao
#' @keywords internal
#' @noRd 
remixt2gg= function(remixt, simplify = TRUE)
{   
  if (!dir.exists(remixt)){
    stop("Input ReMixT directory not found.")
  } else if (length(rmt.out <- dir(remixt, "cn.tsv$|brk.tsv$", full.names=TRUE)) != 2){
    stop("Required output files cn.tsv$ and brk.tsv$ cannot be located.")
  }
  
  rmt.seg = fread(grep("cn.tsv", rmt.out, value=TRUE))
  rmt.seg[, ":="(start = data.table::shift(end)+1)]
  rmt.seg[, start := ifelse(start > end, 1, start)]
  rmt.seg[is.na(start), start:=1]
  rmt.seg[, cn := major_1 + minor_1]
  setnames(rmt.seg, 'chromosome', 'seqnames')
  rmt.tile = dt2gr(rmt.seg)
  rmt.bks = fread(grep("brk.tsv", rmt.out, value=TRUE))
  if (nrow(rmt.bks)>0){
    strmap = setNames(c("+", "-"), c("-", "+"))
    rmt.bks[, cn := cn_1] ## only consider major clone right now
    ## add 1 to "+" positions since our breakgraph convention is to add + junctions to the left of the specified breakpoint
    ## while the remixt convention is to add them to the base immediately following the breakpoint
    if (simplify)
    {
      rmt.bks = rmt.bks[cn>0, ]
    }
    bp1 = dt2gr(rmt.bks[, .(seqnames=chromosome_1,
                            start=position_1 + sign(strmap[strand_1]=='+'),
                            end=position_1 + sign(strmap[strand_1]=='+'),
                            strand=strmap[strand_1])], seqlengths = seqlengths(rmt.tile))
    bp2 = dt2gr(rmt.bks[, .(seqnames=chromosome_2,
                            start=position_2 + sign(strmap[strand_2]=='+'),
                            end=position_2 + sign(strmap[strand_2]=='+'),
                            strand=strmap[strand_2])], seqlengths = seqlengths(rmt.tile))
    juncs = grl.pivot(GRangesList(list(bp1, bp2)))
    values(juncs) = rmt.bks[, .(prediction_id, cn, cn_0, cn_1, cn_2, n_1, side_1, n_2, side_2)]
  } else {
    juncs = NULL
  }

  return(breakgraph(breaks = rmt.tile, juncs = juncs))
}


#' @name read_vcf
#' @title read_vcf: utility function to read VCF into GRanges object
#'
#' @name read_vcf
#' @importFrom VariantAnnotation readVcf
#' @keywords internal
#' @noRd
read_vcf = function (fn, gr = NULL, hg = "hg19", geno = NULL, swap.header = NULL,
                     verbose = FALSE, add.path = FALSE, tmp.dir = "~/temp/.tmpvcf",
                     ...)
{
    in.fn = fn
    if (verbose){
        cat("Loading", fn, "\n")}
    if (!is.null(gr)) {
        tmp.slice.fn = paste(tmp.dir, "/vcf_tmp", gsub("0\\.",
                                                       "", as.character(runif(1))), ".vcf", sep = "")
        cmd = sprintf("bcftools view %s %s > %s", fn, paste(gr.string(gr.stripstrand(gr)),
                                                            collapse = " "), tmp.slice.fn)
        if (verbose){
            cat("Running", cmd, "\n")
        }
        system(cmd)
        fn = tmp.slice.fn
    }
    if (!is.null(swap.header)) {
        if (!file.exists(swap.header)){
            stop(sprintf("Swap header file %s does not exist\n",
                         swap.header))
        }
        system(paste("mkdir -p", tmp.dir))
        tmp.name = paste(tmp.dir, "/vcf_tmp", gsub("0\\.", "",
                                                   as.character(runif(1))), ".vcf", sep = "")
        if (grepl("gz$", fn)){
            system(sprintf("zcat %s | grep '^[^#]' > %s.body",
                           fn, tmp.name))
        } else system(sprintf("grep '^[^#]' %s > %s.body", fn,
                            tmp.name))
        if (grepl("gz$", swap.header)){
            system(sprintf("zcat %s | grep '^[#]' > %s.header",
                           swap.header, tmp.name))
        } else{
            system(sprintf("grep '^[#]' %s > %s.header", swap.header,
                            tmp.name))
        }
        system(sprintf("cat %s.header %s.body > %s", tmp.name,
                       tmp.name, tmp.name))
        vcf = readVcf(tmp.name, hg, ...)
        system(sprintf("rm %s %s.body %s.header", tmp.name, tmp.name,
                       tmp.name))
    } else{
        vcf = readVcf(fn, hg, ...)
    }
    out = granges(vcf)
    if (!is.null(values(out))){
        values(out) = cbind(values(out), info(vcf))
    }
    else values(out) = info(vcf)
    if (!is.null(geno)) {
        if (geno){
            for (g in names(geno(vcf))) {
                geno = names(geno(vcf))
                warning(sprintf("Loading all geno field:\n\t%s",
                                paste(geno, collapse = ",")))
            }
            }
        gt = NULL
        for (g in geno) {
            m = as.data.frame(geno(vcf)[[g]])
            names(m) = paste(g, names(m), sep = "_")
            if (is.null(gt)){
                gt = m
            } else{
                gt = cbind(gt, m)
            }
        }
        values(out) = cbind(values(out), as(gt, "DataFrame"))
    }
    if (!is.null(gr)){
        system(paste("rm", tmp.slice.fn))
    }
    if (add.path){
        values(out)$path = in.fn
    }
    return(out)
}

#' @name read.juncs
#' @title read.juncs: parse junction data from various common formats
#'
#' @description Parsing various formats of structural variation data into junctions.
#'
#' @usage read.juncs(rafile,
#' keep.features = T,
#' seqlengths = NULL,
#' chr.convert = T,
#' geno=NULL,
#' flipstrand = FALSE,
#' swap.header = NULL,
#' breakpointer = FALSE,
#' seqlevels = NULL,
#' force.bnd = FALSE,
#' skip = NA)
#'
#' @param rafile path to the junctions file. See details for the compatible formats.
#' @param keep.features \code{logical}, if TRUE preserve meta data from the input
#' @param seqlengths a named \code{numeric} vector containing reference contig lengths
#' @param chr.convert \code{logical}, if TRUE strip "chr" prefix from contig names
#' @param geno \code{logical}, whether to parse the 'geno' fields of VCF
#' @param flipstrand \code{logical}, if TRUE will flip breakpoint strand
#' @param swap.header path to the alternative VCF header file
#' @param breakpointer \code{logical}, if TRUE will parse as breakpointer output
#' @param seqlevels vector for renaming the chromosomes
#' @param force.bnd if TRUE overwrite all junction "type" to "BND"
#' @param skip \code{numeric} lines to skip
#'
#' @details
#' A junction is a unordered pair of strand-specific genomic locations (breakpoints). Within a given
#' reference genome coordinate system, we call the direction in which coordinates increase "+". A breakpoint
#' is a width 1 (\code{start==end})genomic range with \code{strand} specified, and "+" means the side with larger
#' coordinate is fused with the other breakpoint in a junction.
#'
#' \code{rafile} must be one of the following formats:
#' 1) Some VCF (variant call format). We currently support the VCF output
#' from a number of structural variation detection methods, namely
#' SvABA (https://github.com/walaj/svaba),
#' DELLY (https://github.com/dellytools/delly),
#' LUMPY (https://github.com/arq5x/lumpy-sv),
#' novoBreak (https://sourceforge.net/projects/novobreak/). In theory,
#' VCF defined with BND style should be compatible but be cautious
#' when using the output from other methods since
#' no universal data definition is adopted by the community yet.
#' 2) BEDPE (http://bedtools.readthedocs.io/en/latest/content/general-usage.html#bedpe-format)
#' 3) Textual output from Breakpointer
#' (http://archive.broadinstitute.org/cancer/cga/breakpointer)
#' 4) R serialized object storing junctions (.rds)
#'
#' @section Warning:
#' We assume the orientation definition in the input is consistent with ours. Check with
#' the documentation of your respective method to make sure. If the contrary, use
#' \code{flipstrand=TRUE} to reconcile.
#'
#' @return a \code{GRangesList} of the junctions
#'
#' @keywords internal
#' @noRd
#' @importFrom VariantAnnotation readVcf info
#' @import data.table
read.juncs = function(rafile,
                     keep.features = T,
                     seqlengths = NULL,
                     chr.convert = T,
                     geno=NULL,
                     flipstrand = FALSE,
                     swap.header = NULL,
                     breakpointer = FALSE,
                     seqlevels = NULL,
                     force.bnd = FALSE,
                     skip = NA,
                     verbose = FALSE, 
                     get.loose = FALSE){
    if (is.na(rafile)){
        return(NULL)
    }
    ## if TRUE will return a list with fields $junctions and $loose.ends

    if (is.character(rafile)){
        if (grepl('.rds$', rafile)){
            ra = readRDS(rafile)
            ## validity check written for "junctions" class
            return(Junction$new(ra))
        } else if (grepl('(.bedpe$)', rafile)){
            ra.path = rafile
            cols = c('chr1', 'start1', 'end1', 'chr2', 'start2', 'end2', 'name', 'score', 'str1', 'str2')

            f = file(ra.path, open = "rb")
            headers = character(0)
            thisline = readLines(f, 1)
            while (grepl("^((#)|(chrom)|(chr))", thisline)) {
                headers = c(headers, thisline)
                thisline = readLines(f, 1)
            }
            ln = sum(length(headers), length(thisline))
            while (length(thisline) > 0) {
                ## thisline = readBin(f, "raw", n = 50000)
                ## sum(thisline == as.raw(10L))
                thisline = readLines(f, n = 50000)
                ln = length(thisline) + ln
            }
            lastheader = tail(headers, 1)
            ## ln = readLines(ra.path)
            if (is.na(skip)){
                ## nh = min(c(Inf, which(!grepl('^((#)|(chrom)|(chr))', ln))))-1
                nh = length(headers)
                ## if (is.infinite(nh)){
                ##     nh = 1
                ## }
            } else{
                nh = skip
            }

            if ( (ln-nh) <=0) {
                ## if (get.loose){
                ##     return(list(junctions = GRangesList(GRanges(seqlengths = seqlengths))[c()], loose.ends = GRanges(seqlengths = seqlengths)))
                ## }
                ## else{
                return(GRangesList(GRanges(seqlengths = seqlengths))[c()])
                ## }
            }

            if (nh ==0) {
                rafile = fread(rafile, header = FALSE)
            } else {

                if (nh == 1) {
                    header_arg = TRUE
                    skip_arg = 0
                    bedhead = NULL
                } else if (nh > 1) {
                    header_arg = F
                    skip_arg = nh
                    bedhead = gsub("^#", "", unlist(strsplit(lastheader, "\t|,")))
                }

                rafile = tryCatch(fread(ra.path, header = header_arg, skip = skip_arg), error = function(e) NULL)
                if (is.null(rafile)){
                    rafile = tryCatch(fread(ra.path, header = header_arg, skip = skip_arg, sep = '\t'), error = function(e) NULL)
                }

                if (is.null(rafile)){
                    rafile = tryCatch(fread(ra.path, header = header_arg, skip = skip_arg, sep = ','), error = function(e) NULL)
                }

                if (is.null(rafile)){
                    stop('Error reading bedpe')
                }

                if (!is.null(bedhead) && identical(length(bedhead), ncol(rafile))) {
                    colnames(rafile) = bedhead
                }
            }

            if (nrow(rafile)==0)
                return(GRangesList())
            ## this is not robust enough! there might be mismatching colnames
            setnames(rafile, 1:length(cols), cols)
            rafile[, str1 := ifelse(str1 %in% c('+', '-'), str1, '*')]
            rafile[, str2 := ifelse(str2 %in% c('+', '-'), str2, '*')]
        } else if (grepl('(vcf$)|(vcf.gz$)|(vcf.bgz$)', rafile)){
            vcf = VariantAnnotation::readVcf(rafile)

            ## vgr = rowData(vcf) ## parse BND format
            vgr = read_vcf(rafile, swap.header = swap.header, geno=geno)

            mc = data.table(as.data.frame(mcols(vgr)))

            if (!('SVTYPE' %in% colnames(mc))) {
                warning('Vcf not in proper format.  Is this a rearrangement vcf?')
                return(GRangesList());
            }

            if (any(w.0 <- (width(vgr)<1))){
                warning("Some breakpoint width==0.")
                ## right bound smaller coor
                ## and there's no negative width GR allowed
                bpid = names(vgr)
                names(vgr) = NULL ## for some reason the below lines doesn't like names sometimes
                vgr[which(w.0)] = GenomicRanges::shift(gr.start(vgr[which(w.0)]), -1)
                names(vgr) = bpid
            }

            ## BND format doesn't have duplicated rownames
            if (any(duplicated(names(vgr)))){
                names(vgr) = NULL
            } 

            ## no events
            if (length(vgr) == 0){
                return (GRangesList())
            }

            ## local function that turns old VCF to BND
            .vcf2bnd = function(vgr){
                if (!"END" %in% colnames(values(vgr))){
                    stop("Non BND SV should have the second breakpoint coor in END columns!")
                }

                if (!"CHR2" %in% colnames(values(vgr)) | any(is.na(vgr$CHR2))){
                    vgr$CHR2 = as.character(seqnames(vgr))
                }

                bp2 = data.table(as.data.frame(mcols(vgr)))
                bp2[, ":="(seqnames=CHR2, start=as.numeric(END), end=as.numeric(END))]
                slbp2 = bp2[, pmax(1, end), by = seqnames][, structure(V1, names = seqnames)]
                bp2.gr = dt2gr(bp2, seqlengths = slbp2)
                mcols(bp2.gr) = mcols(vgr)

                if (!is.null(names(vgr)) & !anyDuplicated(names(vgr))){
                    jid = names(vgr)
                } else {
                    jid = seq_along(vgr)
                }
                names(vgr) = paste(paste0("exp", jid), "1", sep=":")
                names(bp2.gr) = paste(paste0("exp", jid), "2", sep=":")

                nm = c(names(vgr), names(bp2.gr))
                vgr = resize(grbind(vgr, bp2.gr), 1)
                names(vgr) = nm

                if (all(grepl("[_:][12]$",names(vgr)))){
                    ## row naming same with Snowman
                    nm <- vgr$MATEID <- names(vgr)
                    ix <- grepl("1$",nm)
                    vgr$MATEID[ix] = gsub("(.*?)(1)$", "\\12", nm[ix])
                    vgr$MATEID[!ix] = gsub("(.*?)(2)$", "\\11", nm[!ix])
                    vgr$SVTYPE="BND"
                }
                return(vgr)
            }            

            ## GRIDSS FIX?
            if ("PARID" %in% colnames(mcols(vgr))) {
                vgr$MATEID = vgr$PARID
            }

            ## TODO: Delly and Novobreak
            ## fix mateids if not included
            ## if ("EVENT" %in% colnames(mc) && any(grepl("gridss", names(vgr)))){
            ##     if (verbose){
            ##         message("Recognized GRIDSS junctions")
            ##     }
            ##     if (!get.loose){
            ##         if (verbose){
            ##             message("Ignoring single breakends")
            ##             paired.ix = grep("[oh]$", names(vgr))
            ##             vgr = vgr[paired.ix]
            ##             mc  = mc[paired.ix]
            ##         }
            ##         ## GRIDSS has event id in "EVENT" column
            ##         ## GRIDSS naming of breakends ends with "o", "h", "b" (single, unpaired)                    
            ##         ematch = data.table(
            ##             ev = as.character(mc$EVENT),
            ##             nms = names(vgr)
            ##         )
            ##         ematch[, ":="(mateid = paste0(ev, ifelse(grepl("o$", nms), "h", "o")))]
            ##         ematch[, ":="(mateix = match(nms, mateid))]
            ##         if (len(mism.ix <- which(is.na(ematch$mateix))) > 0){
            ##             warning("Found ", len(mism.ix), " unpaired breakends, ignoring")
            ##             paired.ix = setdiff(seq_len(nrow(ematch)), mism.ix)
            ##             vgr = vgr[paired.ix]
            ##             mc  = mc[paired.ix]
            ##             ematch = ematch[paired.ix]
            ##         }
            ##         values(vgr)$MATEID = ematch$mateid
            ##     } else {
            ##         stop("Hasn't implemented single breakend parsing for GRIDSS!")
            ##     }
            ## } else
            if (!"MATEID" %in% colnames(mcols(vgr))) {
                ## TODO: don't assume every row is a different junction
                ## Novobreak, I'm looking at you.
                ## now delly...
                ## if SVTYPE is BND but no MATEID, don't pretend to be
                if (length(fake.bix <- which(values(vgr)$SVTYPE=="BND"))!=0){
                    values(vgr)$SVTYPE[fake.bix] = "TRA" ## values(vgr[fake.bix])$SVTYPE = "TRA"
                }

                ## add row names just like Snowman
                if (all(names(vgr)=="N" | ## Novobreak
                        is.null(names(vgr)) |
                        all(grepl("^DEL|DUP|INV|BND", names(vgr)))) ## Delly
                    ){
                    ## otherwise if all "N", as Novobreak
                    ## or starts with DEL|DUP|INV|BND, as Delly
                    ## expand and match MATEID
                    vgr=.vcf2bnd(vgr)
                }
            } else if (any(is.na(mid <- as.character(vgr$MATEID)))){
                ## like Lumpy, the BND rows are real BND but blended with non-BND rows
                ## treat them separately
                if (is.null(vgr$CHR2)){
                    vgr$CHR2 = as.character(NA)
                }

                names(vgr) = gsub("_", ":", names(vgr))
                vgr$MATEID = sapply(vgr$MATEID, function(x) gsub("_", ":", x))

                values(vgr) = data.table(as.data.frame(values(vgr)))

                ## break up the two junctions in one INV line!
                if ("STRANDS" %in% colnames(mc) & any(ns <- sapply(vgr$STRANDS, length)>1)){
                    ## first fix format errors, two strand given, but not comma separeted
                    ## so you'd have taken them as single
                    if (any(fuix <- sapply(vgr[which(!ns)]$STRANDS, stringr::str_count, ":")>1)){
                        which(!ns)[fuix] -> tofix
                        vgr$STRANDS[tofix] = lapply(vgr$STRANDS[tofix],
                                                    function(x){
                                                        strsplit(gsub("(\\d)([\\+\\-])", "\\1,\\2", x), ",")[[1]]
                                                    })
                        ns[tofix] = TRUE
                    }

                    ## for the one line two junction cases
                    ## split into two lines
                    vgr.double = vgr[which(ns)]
                    j1 = j2 = vgr.double
                    st1 = lapply(vgr.double$STRANDS, function(x)x[1])
                    st2 = lapply(vgr.double$STRANDS, function(x)x[2])
                    j1$STRANDS = st1
                    j2$STRANDS = st2
                    vgr.double = c(j1, j2)
                    names(vgr.double) = dedup(names(vgr.double))
                    vgr = c(vgr[which(!ns)], vgr.double)
                }
              
              mid <- as.logical(sapply(vgr$MATEID, length))
              vgr$loose.end = FALSE
              vgr.bnd = vgr[which(mid)]
              vgr.nonbnd = vgr[which(!mid)]

              if (length(vgr.nonbnd))
              {
                if (any(naix <- is.na(vgr.nonbnd$END)))
                  {
                    vgr.nonbnd$END[naix] = -1
                    vgr.nonbnd$loose.end[naix] = TRUE
                  }

                vgr.nonbnd = .vcf2bnd(vgr.nonbnd)
              }
              
                mc.bnd = data.table(as.data.frame(values(vgr.bnd)))
                mc.nonbnd = data.table(as.data.frame(values(vgr.nonbnd)))
                mc.bnd$MATEID = as.character(mc.bnd$MATEID)

                vgr = c(vgr.bnd[,c()], vgr.nonbnd[,c()])
                values(vgr) = rbind(mc.bnd, mc.nonbnd, fill = TRUE)
            }

            ## sanity check
            if (!any(c("MATEID", "SVTYPE") %in% colnames(mcols(vgr)))){
                stop("MATEID or SVTYPE not included. Required")
            }

            vgr$mateid = vgr$MATEID
            ## what's this???
            vgr$svtype = vgr$SVTYPE

            if (!is.null(info(vcf)$SCTG)){
                vgr$SCTG = info(vcf)$SCTG
            }

            if (force.bnd){
                vgr$svtype = "BND"
            }

            if (sum(vgr$svtype == 'BND')==0){
                warning('Vcf not in proper format.  Will treat rearrangements as if in BND format')
            }

            if (!all(vgr$svtype == 'BND')){
                warning(sprintf('%s rows of vcf do not have svtype BND, treat them as non-BND!',
                                sum(vgr$svtype != 'BND')))

            }

            bix = which(vgr$svtype == "BND")
            vgr = vgr[bix]
            alt <- sapply(vgr$ALT, function(x) x[1])

            ## Determine each junction's orientation
            if ("CT" %in% colnames(mcols(vgr))){
                if (verbose)
                {
                    message("CT INFO field found.")
                }
                if ("SVLEN" %in% colnames(values(vgr))){
                    ## proceed as Novobreak
                    ## ALERT: overwrite its orientation!!!!
                    del.ix = which(vgr$SVTYPE=="DEL")
                    dup.ix = which(vgr$SVTYPE=="DUP")
                    vgr$CT[del.ix] = "3to5"
                    vgr$CT[dup.ix] = "5to3"
                }

                ## also, Delly is like this
                ori = strsplit(vgr$CT, "to")
                iid = sapply(strsplit(names(vgr), ":"), function(x)as.numeric(x[2]))
                orimap = setNames(c("+", "-"), c("5", "3"))
                strd = orimap[sapply(seq_along(ori), function(i) ori[[i]][iid[i]])]
                strand(vgr) = strd
                vgr.pair1 = vgr[which(iid==1)]
                vgr.pair2 = vgr[which(iid==2)]
            } else if ("STRANDS" %in% colnames(mcols(vgr))){
                ## TODO!!!!!!!!!!!!!!!
                ## sort by name, record bp1 or bp2
                if (verbose)
                {
                    message("STRANDS INFO field found.")
                }
                iid = sapply(strsplit(names(vgr), ":"), function(x)as.numeric(x[2]))
                vgr$iid = iid
                vgr = vgr[order(names(vgr))]
                iid = vgr$iid

                ## get orientations
                ori = strsplit(substr(unlist(vgr$STRANDS), 1, 2), character(0))
                orimap = setNames(c("+", "-"), c("-", "+"))

                ## map strands
                strd = orimap[sapply(seq_along(ori), function(i) ori[[i]][iid[i]])]
                strand(vgr) = strd

                vgr.pair1 = vgr[which(iid==1)]
                vgr.pair2 = vgr[which(iid==2)]
            } else if (any(grepl("\\[|\\]", alt))){
                if (verbose)
                {
                    message("ALT field format like BND")
                }
                ## proceed as Snowman
                vgr$first = !grepl('^(\\]|\\[)', alt) ## ? is this row the "first breakend" in the ALT string (i.e. does the ALT string not begin with a bracket)
                vgr$right = grepl('\\[', alt) ## ? are the (sharp ends) of the brackets facing right or left
                vgr$coord = as.character(paste(seqnames(vgr), ':', start(vgr), sep = ''))
                vgr$mcoord = as.character(gsub('.*(\\[|\\])(.*\\:.*)(\\[|\\]).*', '\\2', alt))
                vgr$mcoord = gsub('chr', '', vgr$mcoord)

                ## add extra genotype fields to vgr
                if (all(is.na(vgr$mateid))){
                    if (!is.null(names(vgr)) & !any(duplicated(names(vgr)))){
                        warning('MATEID tag missing, guessing BND partner by parsing names of vgr')
                        vgr$mateid = paste(gsub('::\\d$', '', names(vgr)),
                        (sapply(strsplit(names(vgr), '\\:\\:'), function(x) as.numeric(x[length(x)])))%%2 + 1, sep = '::')
                    }
                    else if (!is.null(vgr$SCTG))
                    {
                        warning('MATEID tag missing, guessing BND partner from coordinates and SCTG')
                        ucoord = unique(c(vgr$coord, vgr$mcoord))
                        vgr$mateid = paste(vgr$SCTG, vgr$mcoord, sep = '_')
                        
                        if (any(duplicated(vgr$mateid)))
                        {
                            warning('DOUBLE WARNING! inferred mateids not unique, check VCF')
                            bix = bix[!duplicated(vgr$mateid)]
                            vgr = vgr[!duplicated(vgr$mateid)]
                        }
                    }
                    else{
                        stop('Error: MATEID tag missing')
                    }
                }
                
                vgr$mix = as.numeric(match(vgr$mateid, names(vgr)))

                pix = which(!is.na(vgr$mix))

                vgr.pair = vgr[pix]

                if (length(vgr.pair)==0){
                    stop('Error: No mates found despite nonzero number of BND rows in VCF')
                }

                vgr.pair$mix = match(vgr.pair$mix, pix)

                vix = which(1:length(vgr.pair)<vgr.pair$mix)
                vgr.pair1 = vgr.pair[vix]
                vgr.pair2 = vgr.pair[vgr.pair1$mix]

                ## now need to reorient pairs so that the breakend strands are pointing away from the breakpoint

                ## if "first" and "right" then we set this entry "-" and the second entry "+"
                tmpix = vgr.pair1$first & vgr.pair1$right
                if (any(tmpix)){
                    strand(vgr.pair1)[tmpix] = '-'
                    strand(vgr.pair2)[tmpix] = '+'
                }

                ## if "first" and "left" then "-", "-"
                tmpix = vgr.pair1$first & !vgr.pair1$right
                if (any(tmpix)){
                    strand(vgr.pair1)[tmpix] = '-'
                    strand(vgr.pair2)[tmpix] = '-'
                }

                ## if "second" and "left" then "+", "-"
                tmpix = !vgr.pair1$first & !vgr.pair1$right
                if (any(tmpix)){
                    strand(vgr.pair1)[tmpix] = '+'
                    strand(vgr.pair2)[tmpix] = '-'
                }

                ## if "second" and "right" then "+", "+"
                tmpix = !vgr.pair1$first & vgr.pair1$right
                if (any(tmpix)){
                    strand(vgr.pair1)[tmpix] = '+'
                    strand(vgr.pair2)[tmpix] = '+'
                }

                pos1 = as.logical(strand(vgr.pair1)=='+') ## positive strand junctions shift left by one (i.e. so that they refer to the base preceding the break for these junctions
                if (any(pos1)){
                    start(vgr.pair1)[pos1] = start(vgr.pair1)[pos1]-1
                    end(vgr.pair1)[pos1] = end(vgr.pair1)[pos1]-1
                }

                pos2 = as.logical(strand(vgr.pair2)=='+') ## positive strand junctions shift left by one (i.e. so that they refer to the base preceding the break for these junctions
                if (any(pos2)){
                    start(vgr.pair2)[pos2] = start(vgr.pair2)[pos2]-1
                    end(vgr.pair2)[pos2] = end(vgr.pair2)[pos2]-1
                }
            }

            ra = grl.pivot(GRangesList(vgr.pair1[, c()], vgr.pair2[, c()]))

            ## ALERT: vgr has already been subsetted to only include BND rows
            ## bix is the original indices, so NOT compatible!
            ## this.inf = values(vgr)[bix[pix[vix]], ]
            if (exists("pix") & exists("vix")){
                this.inf = values(vgr)[pix[vix], ]
            }
            if (exists("iid")){
                this.inf = values(vgr[which(iid==1)])
            }

            if (is.null(this.inf$POS)){
                this.inf = cbind(data.frame(POS = ''), this.inf)
            }
            if (is.null(this.inf$CHROM)){
                this.inf = cbind(data.frame(CHROM = ''), this.inf)
            }

            if (is.null(this.inf$MATL)){
                this.inf = cbind(data.frame(MALT = ''), this.inf)
            }

            this.inf$CHROM = seqnames(vgr.pair1)
            this.inf$POS = start(vgr.pair1)
            this.inf$MATECHROM = seqnames(vgr.pair2)
            this.inf$MATEPOS = start(vgr.pair2)
            this.inf$MALT = vgr.pair2$AL

            ## NOT SURE WHY BROKEN
            ## tmp = tryCatch(cbind(values(vgr)[bix[pix[vix]],], this.inf), error = function(e) NULL)
            ## if (!is.null(tmp))
            ##     values(ra) = tmp
            ## else
            ##     values(ra) = cbind(vcf@fixed[bix[pix[vix]],], this.inf)

            values(ra) = this.inf

            if (is.null(values(ra)$TIER)){
                ## baseline tiering of PASS vs non PASS variants
                ## ALERT: mind the naming convention by diff programs
                ## TODO: make sure it is compatible with Delly, Novobreak, Meerkat
                ## Snowman/SvABA uses "PASS"
                ## Lumpy/Speedseq uses "."
                values(ra)$tier = ifelse(values(ra)$FILTER %in% c(".", "PASS"), 2, 3)
            } else {
                values(ra)$tier = values(ra)$TIER
            }

            ## ra = ra.dedup(ra)
            if (!get.loose | is.null(vgr$mix)){
                return(ra)
            } else {
                npix = is.na(vgr$mix)
                vgr.loose = vgr[npix, c()] ## these are possible "loose ends" that we will add to the segmentation

                ## NOT SURE WHY BROKEN
                tmp =  tryCatch( values(vgr)[bix[npix], ],
                                error = function(e) NULL)
                if (!is.null(tmp)){
                    values(vgr.loose) = tmp
                } else{
                    values(vgr.loose) = cbind(vcf@fixed[bix[npix], ], info(vcf)[bix[npix], ])
                }

                return(list(junctions = ra, loose.ends = vgr.loose))
            }
        }
        else
      {
        stop('Unrecognized file extension: currently accepted are .rds, .bedpe, .vcf, .vcf.gz, vcf.bgz')
      }
        ## else {
        ##     rafile = read.delim(rafile)
        ## }
    }


    if (is.data.table(rafile)){
        rafile = as.data.frame(rafile)
    }

    if (nrow(rafile)==0){
        out = GRangesList()
        values(out) = rafile
        return(out)
    }
    
    ## flip breaks so that they are pointing away from junction
    if (flipstrand) {
        rafile$str1 = ifelse(rafile$strand1 == '+', '-', '+')
        rafile$str2 = ifelse(rafile$strand2 == '+', '-', '+')
    }

    if (!is.null(seqlevels)) ## convert seqlevels from notation in tab delim file to actual
    {
        rafile$chr1 = seqlevels[rafile$chr1]
        rafile$chr2 = seqlevels[rafile$chr2]
    }


    if (is.null(rafile$str1)){
        rafile$str1 = rafile$strand1
    }

    if (is.null(rafile$str2)){
        rafile$str2 = rafile$strand2
    }

    if (!is.null(rafile$pos1) & !is.null(rafile$pos2)){
        if (breakpointer){
            rafile$pos1 = rafile$T_BPpos1
            rafile$pos2 = rafile$T_BPpos2
        }

        if (!is.numeric(rafile$pos1)){
            rafile$pos1 = as.numeric(rafile$pos1)
        }

        if (!is.numeric(rafile$pos2)){
            rafile$pos2 = as.numeric(rafile$pos2)
        }

        ## clean the parenthesis from the string

        rafile$str1 <- gsub('[()]', '', rafile$str1)
        rafile$str2 <- gsub('[()]', '', rafile$str2)

        ## goal is to make the ends point <away> from the junction where - is left and + is right
        if (is.character(rafile$str1) | is.factor(rafile$str1)){
            rafile$str1 = gsub('0', '-', gsub('1', '+', gsub('\\-', '1', gsub('\\+', '0', rafile$str1))))
        }

        if (is.character(rafile$str2) | is.factor(rafile$str2)){
          rafile$str2 = gsub('0', '-', gsub('1', '+', gsub('\\-', '1', gsub('\\+', '0', rafile$str2))))
        }


        if (is.numeric(rafile$str1)){
            rafile$str1 = ifelse(rafile$str1>0, '+', '-')
        }

        if (is.numeric(rafile$str2)){
            rafile$str2 = ifelse(rafile$str2>0, '+', '-')
        }

        rafile$rowid = 1:nrow(rafile)

        bad.ix = is.na(rafile$chr1) | is.na(rafile$chr2) | is.na(rafile$pos1) | is.na(rafile$pos2) | is.na(rafile$str1) | is.na(rafile$str2) | rafile$str1 == '*'| rafile$str2 == '*' | rafile$pos1<0 | rafile$pos2<0

        rafile = rafile[which(!bad.ix), ]

        if (nrow(rafile)==0){
            return(GRanges())
        }

        seg = rbind(data.frame(chr = rafile$chr1, pos1 = rafile$pos1, pos2 = rafile$pos1, strand = rafile$str1, ra.index = rafile$rowid, ra.which = 1, stringsAsFactors = F),
                    data.frame(chr = rafile$chr2, pos1 = rafile$pos2, pos2 = rafile$pos2, strand = rafile$str2, ra.index = rafile$rowid, ra.which = 2, stringsAsFactors = F))

        if (chr.convert){
            seg$chr = gsub('chr', '', gsub('25', 'M', gsub('24', 'Y', gsub('23', 'X', seg$chr))))
        }

        out = seg2gr(seg, seqlengths = seqlengths)[, c('ra.index', 'ra.which')];
        out = split(out, out$ra.index)
    } else if (!is.null(rafile$start1) & !is.null(rafile$start2) & !is.null(rafile$end1) & !is.null(rafile$end2)){
        ra1 = gr.flipstrand(GRanges(rafile$chr1, IRanges(rafile$start1, rafile$end1), strand = rafile$str1))
        ra2 = gr.flipstrand(GRanges(rafile$chr2, IRanges(rafile$start2, rafile$end2), strand = rafile$str2))
        out = grl.pivot(GRangesList(ra1, ra2))
    }

    if (keep.features){
        values(out) = rafile[, ]
    }

    ## if (!is.null(pad)){
    ##     out = ra.dedup(out, pad = pad)
    ## }

    if (!get.loose){
        return(out)
    } else{
        return(list(junctions = out, loose.ends = GRanges()))
    }

    return(Junction$new(out))
}



#' @name karyotype
#' @title karyotype
#' @description
#'
#' returns gWalk (if karyo arg is not NULL) or gGraph of
#' cytoBands given chrom.sizes file, with built in colormap
#' for disjoining with other gGraphs and visualizing cytobands
#' in the context of a gGRaph
#'
#' @param karyo karyotype string to generate a gWalk of alleles representing karyotype
#' @param cytoband path or URL to UCSC style cytoband file
#' @param ... Additional arguments sent to the \code{gTrack} constructor
#' @export
#' @author Marcin Imielinski
karyotype = function(karyo = NULL, cytoband = NULL, ... )
{
  if (!is.null(karyo))
    stop('karyotype string TBD, please leave NULL for now')

  if (is.null(cytoband))
    cytoband = system.file("extdata", "hg19.cytoband.txt", package = 'gGnome')

  ucsc.bands = fread(cytoband)
  setnames(ucsc.bands, c('seqnames', 'start', 'end', 'name', 'stain'))

  ucsc.bands[, seqnames := gsub('chr', '', seqnames)]
  sl = ucsc.bands[, max(end), by = seqnames][order(suppressWarnings(as.numeric(seqnames)), seqnames), structure(V1, names = seqnames)]  
  ucsc.bands = dt2gr(ucsc.bands, seqlengths = sl)
  gg = gG(breaks = ucsc.bands)
  
  gg$set(colormaps = list(stain = c('gneg' = 'white', 'gpos25' = 'gray25', 'gpos50' = 'gray50', 'gpos75'= 'gray75', 'gpos100' = 'black', 'acen' = 'red', 'gvar' = 'pink', 'stalk' = 'blue')))
  gg$set(border = 'black', ...)

  return(gg)
}




#' @name pairNodesAndEdges
#' @title pairNodesAndEdges
#' @description
#' 
#' Adds the proper snode.id and index to nodes, adds proper sedge.id/edge.id to edges
#' Returns list(nodes, edges) with updated fields and miscellaneous metadata removed
#' USE THIS FUNCTION WITH CAUTION - gGraphFromNodes will require that sedge.id/edge.id is removed before running (FIXME: don't make it do that)
#'
#' @author Joe DeRose
#' @keywords internal
#' @noRd
pairNodesAndEdges = function(nodes, edges)
{                            
  if('loose' %in% names(values(nodes))) {
    not.loose = which(!nodes$loose)
  } else {
    not.loose = seq_along(nodes)
  }

  edges = as.data.table(edges)

  ix = not.loose[match(gr.stripstrand(nodes[not.loose]), gr.stripstrand(nodes[not.loose]))]
  nodes$snode.id = NA
  nodes$snode.id[not.loose] = ifelse(as.logical(strand(nodes[not.loose]) == '+'), ix, -ix)
  nodes$index = 1:length(nodes)
  
  edges[, tag := paste(nodes$snode.id[to], nodes$snode.id[from])]
  edges[, rtag := paste(-nodes$snode.id[from], -nodes$snode.id[to])]
  edges[, id := 1:.N]
  edges[, rid := match(tag, rtag)]                                 
  edges[, edge.id := igraph::clusters(igraph::graph.edgelist(cbind(id, rid)), 'weak')$membership]
  edges[, sedge.id := ifelse(duplicated(edge.id), -edge.id, edge.id)]
  tmp.edges = edges[, .(from = not.loose[from], to = not.loose[to],
                    cn = if("cn" %in% names(edges)) cn,
                    type = if("type" %in% names(edges)) type,
                    edge.id, sedge.id)]
  extracols = setdiff(colnames(edges), colnames(tmp.edges))
  if (length(extracols)>0)
    {
      edges = cbind(tmp.edges, edges[, extracols, with = FALSE])
    }

  nodes = nodes[not.loose]
  edges = edges[!is.na(to) & !is.na(from), ]

  return(list(nodes, edges))
}


#' @name inferLoose
#' @description
#' Given nodes and edges in unsigned / n1 / n2, n1.side / n2.side format
#' and $cn field on both nodes and edges infer loose ends by comparing
#' the node and edge cns.
#' 
#' @param nodes unstranded GRanges
#' @param edges edges dt such in the format of the input to gG
#' @param force whether to fill in the 'cn' field in edges when not found
#' @return GRanges with logical columns $loose.left and $loose.right computed
inferLoose = function(nodes, edges, force = TRUE)
{
  nodes = gr.stripstrand(nodes)
  nodes.out = nodes
  nodes$cn.left = nodes$cn.right = 0;

  if (is.null(nodes$cn))
    stop('cn field must be present in nodes')

  if (any(nodes$cn<0, na.rm = TRUE) | any((nodes$cn %% 1)!=0, na.rm = TRUE)){
    stop('cn must be non-negative integers')
  }    

  if (nrow(edges)>0)
  {
      if (is.null(edges$cn)){
          if (force){
              warning("'cn' field not found in edges, force into zeroes")
              edges[, cn := 0]
          } else {
              stop('cn field must be present in edges')
          }
      }

    if (any((edges$cn<0) | ((edges$cn %% 1)!=0), na.rm = TRUE))
      stop('cn must be non-negative integers')

    if (!is.character(edges$n1.side))
      {
        edges[, n1.side := ifelse(n1.side==1, 'right', 'left')]
        edges[, n2.side := ifelse(n2.side==1, 'right', 'left')]
      }

    .loose = function(nodes, edges)
    {
      left.esum = merge(edges[n1.side == 'left', sum(cn, na.rm = TRUE), keyby = n1],
                        edges[n2.side == 'left', sum(cn, na.rm = TRUE), keyby = n2], by.x = "n1", by.y = 'n2', all = TRUE)
      left.esum[, cn := rowSums(cbind(V1.x, V1.y), na.rm = TRUE)]
      setkey(left.esum, n1)
      
      right.esum = merge(edges[n1.side == 'right', sum(cn, na.rm = TRUE), keyby = n1],
                         edges[n2.side == 'right', sum(cn, na.rm = TRUE), keyby = n2], by.x = "n1", by.y = 'n2', all = TRUE)
      right.esum[, cn := rowSums(cbind(V1.x, V1.y), na.rm = TRUE)]
      setkey(right.esum, n1)  
      nodes$cn.left = pmax(0, left.esum[.(1:length(nodes)), cn], na.rm = TRUE)
      nodes$cn.right = pmax(0, right.esum[.(1:length(nodes)), cn], na.rm = TRUE)
      return(nodes)
    }

    nodes = .loose(nodes, edges)

    if (any(ix <- edges[, is.na(cn) & type == 'REF']))
    {
      ## infer reference cns using missing side cn
      missing.cn = rbind(data.table(id = 1:length(nodes), side = "left", missing = nodes$cn - nodes$cn.left),
                      data.table(id = 1:length(nodes), side = "right", missing = nodes$cn - nodes$cn.right))
      setkeyv(missing.cn, c('id', 'side'))
      edges[ix, cn := as.integer(pmin(missing.cn[.(n1, n1.side), missing], missing.cn[.(n2, n2.side), missing]))]
      nodes = .loose(nodes, edges)
    }

  }
    
  nodes.out$loose.left = nodes$cn != nodes$cn.left
  nodes.out$loose.right = nodes$cn != nodes$cn.right

  return(nodes.out)
}


#' @name haplograph
#' @title haplograph
#' @description
#'
#' Haplograph creates a graph from a set of walks each which is joined
#' to a reference graph backbone
#' i.e. each haplotype is a bubble on the original reference graph
#'
#' Haplotype ends also branch to each other, as long as the termini
#' of both haplotype's end nodes are contained in the other end node.
#' 
#' @param walks gWalk or GRangesList (with seqinfo fully populated)
#' @export
haplograph = function(walks, breaks = NULL)
{
  if (inherits(walks, 'gWalk'))
    walks = walks$grl

  if (is.null(breaks))
    breaks = grl.unlist(walks$grl)

  ## rebuild walks (otherwise inherit giant graph)
  walks = gW(grl = walks)

  sources = sapply(walks$snode.id, head, 1)
  sinks = sapply(walks$snode.id, tail, 1)

  ## mark sources and sinks so we can keep track of them in the final graph
  walks$graph$nodes$mark(is.source = FALSE, is.sink = FALSE)
  walks$graph$nodes[sources]$mark(is.source = TRUE, walk.id = 1:length(sources))
  walks$graph$nodes[sinks]$mark(is.sink= TRUE, walk.id = 1:length(sinks))

  walksd = walks$copy$disjoin(gr = grbind(breaks, grl.unlist(walks$grl)), collapse = FALSE)

  ## recompute sources and sinks to get walk "termini" ie the width ranges at the tips
  ## of each walk
  sources = sapply(walksd$snode.id, head, 1)
  sinks = sapply(walksd$snode.id, tail, 1)

  ## retrieve granges of walk sources and sinks to figure out breaks
  grs = walksd$graph$nodes[sources]$gr %>% gr.start(ignore.strand = FALSE)
  grs$is.source = TRUE; grs$is.sink = FALSE
  gre = walksd$graph$nodes[sinks]$gr %>% gr.end(ignore.strand = FALSE)
  gre$is.source = FALSE; gre$is.sink = TRUE

  ## need to figure out which end segments the sinks and sources match
  ## and then only include the walk pairs with mutual matches

  ## break negative stranded starts and
  ## positive stranded ends one base to the right
  breaks = grbind(breaks %>% disjoin,
    c(grs %>% GenomicRanges::shift(as.integer(sign(strand(grs)=='-'))),
      gre %>% GenomicRanges::shift(as.integer(sign(strand(gre)=='+')))))
  width(breaks) = 0
 
  ## make wild type graph using breaks associated with these starts
  gd = gG(breaks = breaks[, c()])    

  ## combine the graphs
  gn = c(wt = gd, variant = walksd$graph)

  ## now we need to suture the walks with each other and the wt graphs
  ## since we defined grs and gre on the walksd graph will need to figure
  ## out the ids in the combined graph, and then figure out the new
  ## edges to add via $connect

  #############
  ## WT graph suturing
  #############
  
  ## for each start and end node find its reference neighbor, which for 
  ## a positive start node is the right side of the -1 base node
  ## a negative start node is the right side of the +1 base node
  ## a positive end node is the left side of the +1 base
  ## a negative end node is the right side of the -1 base
  grs$rn.node.id.wt = gr.match(grs %>% GenomicRanges::shift(as.integer(sign((strand(grs)=='-') - 0.5))), gd$nodes$gr)
  grs$rn.side = ifelse(strand(grs)=='+', 'right', 'left')
  gre$rn.node.id.wt = gr.match(gre %>% GenomicRanges::shift(as.integer(sign((strand(gre)=='+') - 0.5))), gd$nodes$gr)
  gre$rn.side = ifelse(strand(gre)=='-', 'right', 'left')

  ## the edges will leave each of the
  ## positive start nodes on the left side
  ## negative start nodes on the right side
  ## positive end nodes on the right side
  ## negative end nodes on the left side
  grs$side = ifelse(strand(grs)=='+', 'left', 'right')
  gre$side = ifelse(strand(gre)=='-', 'left', 'right')

  ## map these nodes ids to the current graph gn
  grs$rn.node.id.new = gn$nodes[parent.graph == 'wt']$dt[match(grs$rn.node.id.wt, og.node.id), node.id]
  gre$rn.node.id.new = gn$nodes[parent.graph == 'wt']$dt[match(gre$rn.node.id.wt, og.node.id), node.id]

  ## match starts to appropriate variant nodes using both coordinate and og.node.id
  ## (can't just use og.node.id since the starts may have been split  
  grs$node.id.new = gn$nodes[parent.graph == 'variant']$dt[match(grs$node.id, og.node.id), node.id]
  gre$node.id.new = gn$nodes[parent.graph == 'variant']$dt[match(gre$node.id, og.node.id), node.id]

  ## create new edge data.table from grs and gre

  edges = (grbind(grs, gre) %>% as.data.table)[, .(n1 = node.id.new, n1.side = side,
                                                   n2 = rn.node.id.new, n2.side = rn.side)][!is.na(n1) & !is.na(n2), ]

  ## add these edges
  gn$connect(n1 = edges$n1, n2 = edges$n2, n1.side = edges$n1.side, n2.side = edges$n2.side, type = 'REF', meta = data.table(stype = rep('V-WT', nrow(edges))))

  #############
  ## variant-variant suturing
  #############

  ## now we also suture an end of terminal (ie source or sink) nodes in every variant walk i to the
  ## terminal end of any walk j whose end overlaps the same terminus of of walk j

  ## overlap grs and gre with termini
  termini = walks$nodes[is.source | is.sink]
  termini = gn$nodes[is.source | is.sink]

  grsov = (grs[, c('walk.id', 'node.id.new', 'is.source', 'is.sink', 'side', 'rn.side')] %>% GenomicRanges::shift(as.integer(sign((strand(grs)=='-') - 0.5)))) %*% termini$gr[, c('is.source', 'is.sink', 'node.id', 'walk.id')]
  greov = (gre[, c('walk.id', 'node.id.new', 'is.source', 'is.sink', 'side', 'rn.side')] %>% GenomicRanges::shift(as.integer(sign((strand(gre)=='+') - 0.5)))) %*% termini$gr[, c('is.source', 'is.sink', 'node.id', 'walk.id')]

  names(values(grsov)) = dedup(names(values(grsov)))
  names(values(greov)) = dedup(names(values(greov)))

  ## note: side and rn.side will remain as above

  ## now restrict to mutual matches ie distinct walk "end" pairs i j
  ## where for example the <end> of the x of i intersects the y of j
  ## AND the <end> of the y of j intersects the x of i
  ## where x, y \in {source, sink} 
  ovs = grbind(grsov, greov) %>% gr2dt

  if (nrow(ovs))
  {
    ovs = ovs[walk.id != walk.id.2, ]      
    

    ## tag just let's us match {i x } <-> {j y} "mates"
    ## the tag is done so that i is first if i<j so the
    ## both rows receive the same tag and can be grouped
    ovs[, tag := ifelse(walk.id< walk.id.2,
                        paste(is.sink, walk.id, is.sink.2, walk.id.2),
                        paste(is.sink.2, walk.id.2, is.sink, walk.id))]


    ## find i j mates ie those i x that match j y
    ## we only need to keep 1 of the 2 possible reference edges since they
    ## are equivalent (hence the !duplicated)

    ovs[, count := .N, by = tag] ## count should be only 1 or 2

    if (ovs[, !all(count %in% c(1,2))])
      stop('Something wrong with variant-variant suturing')

    ovs = ovs[count==2, ][!duplicated(tag), ]



    ## add these edges
    gn$connect(n1 = ovs$node.id.new, n2 = ovs$node.id, n1.side = ovs$side, n2.side = ovs$rn.side, type = 'REF', meta = data.table(stype = rep('V-V', nrow(ovs))))
  }
  
  ## remove loose ends at all starts and ends
  ## note: since we are using signed nodes then "loose.left" and "loose.right"
  ## is guaranteed  be oriented in the proper orientation (where 5' is left
  ## and 3' is right)
  tmp = gn$nodes[sign(grs$snode.id)*grs$node.id.new]
  tmp$loose.left = FALSE

  tmp = gn$nodes[sign(gre$snode.id)*gre$node.id.new]
  tmp$loose.right = FALSE

  return(gn)
}

#' @name cougar2gg
#' @title cougar2gg
#' @description
#'
#' Parse CouGaR results into a gGraph
#' 
#' @param cougar directory containing CouGaR results
#' @export
cougar2gg = function(cougar){
    ## "Convert the cougar output directory to gGraph."
    if (!dir.exists(cougar)){
        stop("Error: invalid input CouGaR directory!")
    }

    if (!dir.exists(paste(cougar, 'solve',sep = '/'))){
        stop("No CouGaR solutions found in the input directory!")
    }

    .parsesol = function(this.sol)
    {
        ## verbose = getOption("gGnome.verbose")
        tmp = unlist(.parseparens(this.sol[2]))
        tmp2 = as.data.table(matrix(tmp[nchar(stringr::str_trim(tmp))>0], ncol = 3, byrow = TRUE))
        segs = cbind(
            as.data.table(matrix(unlist(strsplit(tmp2$V1, ' ')), ncol = 2, byrow = TRUE))[, .(seqnames = V1, start = V2)],
            data.table(end = as.numeric(sapply(strsplit(tmp2$V2, ' '), '[', 2)), strand = '*'),
            as.data.table(matrix(unlist(strsplit(stringr::str_trim(tmp2$V3), ' ')),
                                 ncol = 4, byrow = TRUE))[, .(type = V1, cn = as.numeric(V2), ncov = V3, tcov  = V4)])
        segs = suppressWarnings(dt2gr(segs))

        ## any aberrant edges?
        tmp3 = unlist(.parseparens(this.sol[3]))
        if (length(tmp3)>0){
            tmp4 = as.data.table(matrix(tmp3[nchar(stringr::str_trim(tmp3))>0], ncol = 3, byrow = TRUE))
            abadj = cbind(
                as.data.table(matrix(unlist(strsplit(tmp4$V1, ' ')), ncol = 2, byrow = TRUE))[, .(seqnames1 = V1, pos1 = V2)],
                as.data.table(matrix(unlist(strsplit(tmp4$V2, ' ')), ncol = 2, byrow = TRUE))[, .(seqnames2 = V1, pos2 = V2)],
                as.data.table(matrix(unlist(strsplit(stringr::str_trim(tmp4$V3), ' ')),
                                     ncol = 4, byrow = TRUE))[
                  , .(type = V1, cn = as.numeric(V2), ncov = V3, tcov  = V4)]
            )
            ## decide if pos1 and pos2 are ordered
            abadj[seqnames1==seqnames2, ordered := pos1<pos2]
            chr2num = setNames(1:24, c(1:22, "X", "Y"))
            abadj[seqnames1!=seqnames2, ordered := chr2num[as.character(seqnames1)]<chr2num[as.character(seqnames2)]]

            ## set up orientation
            abadj[type==3 & ordered, ":="(strand1 = "+", strand2 = "+")]
            abadj[type==3 & !ordered, ":="(strand1 = "-", strand2 = "-")]
            abadj[type==2 & ordered, ":="(strand1 = "-", strand2 = "-")]
            abadj[type==2 & !ordered, ":="(strand1 = "+", strand2 = "+")]
            abadj[type==0 & ordered, ":="(strand1 = "-", strand2 = "+")]
            abadj[type==0 & !ordered, ":="(strand1 = "+", strand2 = "-")]
            abadj[type==1 & ordered, ":="(strand1 = "+", strand2 = "-")]
            abadj[type==1 & !ordered, ":="(strand1 = "-", strand2 = "+")]

            ## move 

            ## convert to junctions
            ## TODO: make it not depend on hg19!!!!!!!!
            jmd = abadj[, .(cougar_type = type, cn, ncov, tcov, ordered)]
            bp1 = gUtils::dt2gr(abadj[, .(
                seqnames = seqnames1, start = as.numeric(pos1), end = as.numeric(pos1), strand = strand1)],
                seqlengths = hg_seqlengths(chr = FALSE)[1:24])
            bp2 = gUtils::dt2gr(abadj[, .(
                seqnames = seqnames2, start = as.numeric(pos2), end = as.numeric(pos2), strand = strand2)],
                seqlengths = hg_seqlengths(chr = FALSE)[1:24])
            juncs = grl.pivot(GRangesList(bp1, bp2))
            values(juncs) = jmd
            
        } else {
            juncs = GRangesList()
        }
        ## segs$id = seq_along(segs)
        ## gg = gG(breaks = segs)
        ## gg = gG(breaks = segs, juncs = juncs)

        ## ====== NEW approach ====== ##
        ## using breaks and juncs to instantiate the graph, like `wv2gg`
        

        ## ====== OLD approach ====== ##
        ## nodes = c(segs, gr.flipstrand(segs))
        ## nodes$nid = ifelse(as.logical(strand(nodes) == '+'), 1, -1)*nodes$id
        ## nodes$ix = seq_along(nodes)
        ## nodes$rix = match(-nodes$nid, nodes$nid)
        ## adj = array(0, dim = rep(length(nodes),2))
        ## adj = sparseMatrix(length(nodes),length(nodes), x = 0)

        ## tmp = unlist(.parseparens(this.sol[3]))
        ## if (length(tmp)>0) ## are there any somatic edges?
        ## {
        ##     tmp2 = as.data.table(matrix(tmp[nchar(stringr::str_trim(tmp))>0], ncol = 3, byrow = TRUE))
        ##     abadj = cbind(
        ##         as.data.table(matrix(unlist(strsplit(tmp2$V1, ' ')), ncol = 2, byrow = TRUE))[, .(seqnames1 = V1, pos1 = V2)],
        ##         as.data.table(matrix(unlist(strsplit(tmp2$V2, ' ')), ncol = 2, byrow = TRUE))[, .(seqnames2 = V1, pos2 = V2)],
        ##         as.data.table(matrix(unlist(strsplit(stringr::str_trim(tmp2$V3), ' ')),
        ##                              ncol = 4, byrow = TRUE))[, .(type = V1, cn = as.numeric(V2), ncov = V3, tcov  = V4)]
        ##     )
        ##     abadj$strand1 = ifelse(abadj$type %in% c(0,2), '+', '-')
        ##     abadj$strand2 = ifelse(abadj$type %in% c(0,3), '+', '-')
            
        ##     abadj$start.match1 = match(abadj[, paste(seqnames1, pos1)], paste(seqnames(segs), start(segs)))
        ##     abadj$end.match1 = match(abadj[, paste(seqnames1, pos1)], paste(seqnames(segs), end(segs)))
        ##     abadj$start.match2 = match(abadj[, paste(seqnames2, pos2)], paste(seqnames(segs), start(segs)))
        ##     abadj$end.match2 = match(abadj[, paste(seqnames2, pos2)], paste(seqnames(segs), end(segs)))

        ##     ## if strand1 == '+' then end match
        ##     ## if strand1 == '-' then start match
        ##     ## if strand2 == '+' then start match
        ##     ## if strand2 == '-' then end match
            
        ##     abadj[, match1 := ifelse(strand1 == '+', end.match1, -start.match1)]
        ##     abadj[, match2 := ifelse(strand2 == '+', start.match2, -end.match2)]

            
        ##     abadj[, nmatch1 := match(match1, nodes$nid)]
        ##     abadj[, nmatch2 := match(match2, nodes$nid)]

        ##     abadj[, nmatch1r := match(-match1, nodes$nid)]
        ##     abadj[, nmatch2r := match(-match2, nodes$nid)]
            
        ##     adj[cbind(abadj$nmatch1, abadj$nmatch2)] = abadj$cn
        ##     adj[cbind(abadj$nmatch2r, abadj$nmatch1r)] = abadj$cn

        ##     abadj[, ":="(n1 = abs(match1), n1.side = ifelse(match1>0, "right", "left"),
        ##                  n2 = abs(match2), n2.side = ifelse(match2>0, "right", "left"))]
        ## }

        ## ## how many node copies are unaccounted for by aberrant edges on left and right
        ## node.diff.in = nodes$cn - colSums(adj)
        ## node.diff.out = nodes$cn - rowSums(adj)

        ## norm.adj = as.data.table(cbind(seq_along(segs), match(gr.end(segs), gr.start(segs))))[!is.na(V2), ]
        ## norm.adj = rbind(norm.adj, norm.adj[, .(V2 = -V1, V1 = -V2)])[, nid1 := match(V1, nodes$nid)][, nid2 := match(V2, nodes$nid)]

        ## ## now add non-aberrant edge copy numbers that are the minimum of the unaccounted
        ## ## for copy number going <out> of the source node and going <in> to the sink node
        ## adj[as.matrix(norm.adj[, .(nid1, nid2)])] =
        ##     pmin(node.diff.out[norm.adj[, nid1]], node.diff.in[norm.adj[, nid2]])

        ## nodes$eslack.in = nodes$cn - colSums(adj)
        ## nodes$eslack.out = nodes$cn - rowSums(adj)

        ## if (sum(adj!=0)>0)
        ## {
        ##     if (!identical(adj[which(adj>0)], adj[as.matrix(as.data.table(which(adj!=0, arr.ind = TRUE))[, .(row = nodes$rix[col], col = nodes$rix[row])])]))
        ##     {
        ##         stop('reciprocality violated')
        ##     }
        ## }
        ## end(nodes) = end(nodes)-1
        end(segs) = end(segs)-1 ## TODO figure out how to match segment ends with junction bps
        return(list(breaks = segs, juncs = juncs)) ## return(list(nodes, as(adj, 'Matrix')))
    }

    .parseparens = function(str)
    {
        cmd = gsub(',$', '', gsub(',\\)', ')', gsub('\\)', '),', gsub('\\(', 'list(',
                                                                      gsub('([^\\(^\\[^\\]^\\)]+)', '"\\1",', perl = TRUE, gsub('\\]', ')', gsub('\\[', '\\(', str)))))))
        eval(parse(text = cmd))
    }

    sols.fn = dir(dir(paste(cougar, 'solve',sep = '/'), full = TRUE)[1], '^g_', full = TRUE)
    sols.fn = sols.fn[which(!grepl("svg", sols.fn))]
    sols = lapply(
        sols.fn,
        ## dir(dir(paste(cougar, 'solve',sep = '/'), full = TRUE)[1], '^g_', full = TRUE),
        readLines)
    
    ## if (length(sols)==0){
    ##     return(self$nullGGraph())
    ## }

    ## parse cougar graphs
    graphs = lapply(seq_along(sols), function(i){x = sols[[i]]; .parsesol(x)})

    ## concatenate nodes and block diagonal bind adjacency matrices
    segs = do.call('grbind', lapply(graphs, '[[', "breaks"))## segs = do.call('c', lapply(graphs, '[[', 1))
    juncs = do.call('grl.bind', lapply(graphs, '[[', "juncs"))
    ## segs$id = paste(rep(seq_along(graphs), sapply(lapply(graphs, '[[', 1), length)), segs$id, sep = '.')
    ## segs$nid = paste(rep(seq_along(graphs), sapply(lapply(graphs, '[[', 1), length)), segs$nid, sep = '.')
    ## segs$ix = paste(rep(seq_along(graphs), sapply(lapply(graphs, '[[', 1), length)), segs$ix, sep = '.')
    ## segs$rix = paste(rep(seq_along(graphs), sapply(lapply(graphs, '[[', 1), length)), segs$rix, sep = '.')
    ## segs$rix = match(segs$rix, segs$ix)
    ## segs$ix = seq_along(segs)
    ## adj = do.call('bdiag', lapply(graphs, '[[', 2))

    ## final double check for identicality
    ## if (!(identical(adj[which(adj>0)], adj[as.matrix(as.data.table(which(adj!=0, arr.ind = TRUE))[, .(row = segs$rix[col], col = segs$rix[row])])])))
    ## {
    ##     stop('Reciprocality check failed!')
    ## }

    ## gg = gGraph$new(segs = segs, es = adj)$fillin()$decouple()
    return(list(breaks = segs, juncs = juncs))
}


#' @name alignments2gg
#' @title alignments2gg
#'
#' @description
#' Builds a gGraph representing a set of alignments provided as GRanges in "SAM" format
#' i.e. with $cigar, $flag, $qname e.g. output of read.bam
#'
#' The gGraph represents all the implicit nodes and edges in an "end to end" walk of all
#' sequences in the query (specified  $qname and $flag fields - i.e. specifiying R1 and R2 in
#' paired end alignments) through the alignment specified in the alignment record.
#'  
#' The outputted graph can be further walked to exhaustively or greedily identify linear alignments
#' to the reference.
#' 
#' @param tile GRanges of tiles
#' @param juncs Junction object or grl coercible to Junctions object
#' @param genome seqinfo or seqlengths
#' @return list with gr and edges which can be input into standard gGnome constructor
#' @author Marcin Imielinski, Joe DeRose, Xiaotong Yao
#' @keywords internal
#' @noRd 
alignments2gg = function(alignment, verbose = TRUE)
{

  if (inherits(alignment, 'GRangesList') | inherits(alignment, 'CompressedGRangesList')){
      alignment = grl.unlist(alignment)
  }
  if (!inherits(alignment, 'GRanges') || !all(c('qname', 'cigar', 'flag') %in%  names(values(alignment))))
    stop('alignment input must be GRanges with fields $qname $cigar and $flag')

  if (verbose)
    message('making cgChain')

  cg = gChain::cgChain(alignment)

  if (verbose)
    message('disjoining query ranges and lifting nodes to reference')

  lgr = gChain::links(cg)$x
  verboten = c("seqnames", "ranges",
    "strand", "seqlevels", "seqlengths", "isCircular", "start", "end",
    "width", "element")
  values(lgr) = cbind(values(lgr), values(cg)[, setdiff(names(values(cg)), verboten)])
  grc = gr.disjoin(grbind(lgr, si2gr(gChain::links(cg)$x)))
  grc$qname = seqnames(grc)
  gwc = gW(grl = split(grc, seqnames(grc)))

  nodes = gwc$graph$nodes
  grr = gChain::lift(cg, nodes$gr)
  grr$insertion = FALSE

  ## add a pad either to the right or left (basically, there should always be mapped sequence on one side on an insertion ..
  ## otherwise there is no alignment (ie pure insertions means no alignment)
  ix = setdiff(nodes$gr$node.id, grr$node.id)
  if (length(ix))
  {
    insertions = nodes[ix]$gr
    start(insertions) = ifelse(start(insertions)>1, start(insertions)-1, start(insertions))
    end(insertions) = ifelse(start(insertions)== 1 & end(insertions) < seqlengths(insertions)[as.character(seqnames(insertions))],
                             end(insertions)+1, end(insertions))
    insertions$insertion = TRUE
    grr = grbind(grr, gChain::lift(cg, insertions)) ## add the lifted insertions to the pile of intervals
  }

  ugrr = unique(gr.stripstrand(grr))

  ## there may be dups here if say the lift aligns the contig to both the negative and positive side
  ## of a contig
  grr$ugrr.id = match(gr.stripstrand(grr), ugrr)
  grr$grr.id = 1:length(grr)

  if (any(ugrr$insertion))
  {
    width(ugrr[ugrr$insertion]) = 0
  }

  ## find insertions ie nodes that did not survive the lift
  edges = gwc$graph$edges$dt

  if (verbose)
    message('lifting edges to reference')


  ## to lift to genome cordinates, merge old edges with new ids (will duplicate edges across multimaps)
  edges.new = edges %>% merge(gr2dt(grr), by.x = 'n1', by.y = 'node.id', allow.cartesian = TRUE) %>% merge( gr2dt(grr), by.x = 'n2', by.y = 'node.id', allow.cartesian = TRUE)
  edges.new[, n1 := ugrr.id.x] ## we map to the 
  edges.new[, n2 := ugrr.id.y]

  ## flip sides for nodes that are flipped (i.e. negative strand) during lift
  .flip = function(x) c(left = 'right', right = 'left')[x]

  edges.new$n1.side = ifelse(strand(grr)[edges.new$grr.id.x] == '-', .flip(edges.new$n1.side), edges.new$n1.side)
  edges.new$n2.side = ifelse(strand(grr)[edges.new$grr.id.y] == '-', .flip(edges.new$n2.side), edges.new$n2.side)

  ## now just need to replace any edges to and from an insertion
  ugrr$loose.left = ugrr$loose.right = NULL

  if (verbose)
    message('building graph')
  
  return(list(nodes = ugrr, edges = edges.new[, .(n1, n1.side, n2, n2.side)]))
}
