
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
    juncsGR = gr.fix(grl.unlist(juncs$grl), nodes) ## grl.ix keeps track of junction id
    nodes = gr.fix(nodes, juncsGR)

    ## keep track of nmeta to paste back later
    nodes$nid = 1:length(nodes)
    nmeta = as.data.table(values(nodes))

    ## disjoin nodes and bps
    bps = gr.start(juncsGR[, c('grl.ix')], ignore.strand = FALSE) ## trim bps to first base
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

    new.edges = ov[, .(
      n1 = node.ends$nid[subject.id[1]],
      n2 = node.ends$nid[subject.id[2]],
      n1.side = ifelse(strand[1] == '-', 1, 0),
      n2.side = ifelse(strand[2] == "-", 1, 0),
      type = "ALT"
      ), keyby = grl.ix]

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
    return(list(nodes = jabba$nodes$gr, edges = jabba$edges$dt))
  }

  snodes = jabba$segstats %Q% (loose == FALSE)
  snodes$index = 1:length(snodes)
  snodes$snode.id = ifelse(as.logical(strand(snodes)=='+'), 1, -1) * gr.match(snodes, unique(gr.stripstrand(snodes)))

  if (length(snodes)==0)
    return(gG(genome = seqinfo(segs)))

  sedges = spmelt(jabba$adj[jabba$segstats$loose == FALSE, jabba$segstats$loose == FALSE])
  
  nodes = snodes %Q% (strand == '+')
  nodes$loose.left = nodes$eslack.in>0
  nodes$loose.right = nodes$eslack.out>0

  nodes$loose.cn.left = nodes$eslack.in
  nodes$loose.cn.right = nodes$eslack.out
  
  if (nrow(sedges)==0)
    gG(nodes = nodes)

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
          ab.edges[, cn := NULL] ## confilct with the cn inferred from adj
      }

      sedges[.(ab.edges$from, ab.edges$to), type := 'ALT']
      sedges = merge(sedges, ab.edges, by = c('from', 'to'), all.x = TRUE)
    }
  }

    ## rescue any hom-del ref edge
    if (any(data.table(jabba$edges)[type=="reference", cn==0])){
        sedges = rbind(sedges,
                       data.table(jabba$edges)[cn==0 & type=="reference",
                                               .(from, to, type="REF", cn)],
                       fill=T)
    }

  edges = convertEdges(snodes, sedges, metacols = TRUE)

  return(list(nodes = nodes[, c('cn', 'loose.left', 'loose.right', 'loose.cn.left', 'loose.cn.right')],
              edges = edges))
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

            ln = readLines(ra.path)
            if (is.na(skip)){
                nh = min(c(Inf, which(!grepl('^((#)|(chrom))', ln))))-1
                if (is.infinite(nh)){
                    nh = 1
                }
            } else{
                nh = skip
            }

            if ((length(ln)-nh)==0){
                ## if (get.loose){
                ##     return(list(junctions = GRangesList(GRanges(seqlengths = seqlengths))[c()], loose.ends = GRanges(seqlengths = seqlengths)))
                ## }
                ## else{
                return(GRangesList(GRanges(seqlengths = seqlengths))[c()])
                ## }
            }

            if (nh ==0){
                rafile = fread(rafile, header = FALSE)
            } else {

                rafile = tryCatch(fread(ra.path, header = FALSE, skip = nh), error = function(e) NULL)
                if (is.null(rafile)){
                    rafile = tryCatch(fread(ra.path, header = FALSE, skip = nh, sep = '\t'), error = function(e) NULL)
                }

                if (is.null(rafile)){
                    rafile = tryCatch(fread(ra.path, header = FALSE, skip = nh, sep = ','), error = function(e) NULL)
                }

                if (is.null(rafile)){
                    stop('Error reading bedpe')
                }
            }
            ## this is not robust enough! there might be mismatching colnames
            setnames(rafile, 1:length(cols), cols)
            rafile[, str1 := ifelse(str1 %in% c('+', '-'), str1, '*')]
            rafile[, str2 := ifelse(str2 %in% c('+', '-'), str2, '*')]
        } else if (grepl('(vcf$)|(vcf.gz$)', rafile)){
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
                vgr[which(w.0)] = GenomicRanges::shift(gr.start(vgr[which(w.0)]), -1)
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
                bp2.gr = dt2gr(bp2)
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

            ## TODO: Delly and Novobreak
            ## fix mateids if not included
            if (!"MATEID" %in% colnames(mcols(vgr))) {
                ## TODO: don't assume every row is a different junction
                ## Novobreak, I'm looking at you.
                ## now delly...
                ## if SVTYPE is BND but no MATEID, don't pretend to be
                if (length(fake.bix <- which(values(vgr)$SVTYPE=="BND"))!=0){
                    values(vgr[fake.bix])$SVTYPE = "TRA"
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
                vgr.bnd = vgr[which(mid)]
                vgr.nonbnd = vgr[which(!mid)]

                vgr.nonbnd = .vcf2bnd(vgr.nonbnd)

                mc.bnd = data.table(as.data.frame(values(vgr.bnd)))
                mc.nonbnd = data.table(as.data.frame(values(vgr.nonbnd)))
                mc.bnd$MATEID = as.character(mc.bnd$MATEID)

                vgr = c(vgr.bnd[,c()], vgr.nonbnd[,c()])
                values(vgr) = rbind(mc.bnd, mc.nonbnd)
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
            }
            else if (any(grepl("\\[|\\]", alt))){
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
            }
            else{
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
        } else{
            rafile = read.delim(rafile)
        }
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
    chrom.sizes = system.file("extdata", "hg19.cytoband.txt", package = 'gGnome')

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
  nodes.out = nodes
  nodes$cn.left = nodes$cn.right = 0;

  if (is.null(nodes$cn))
    stop('cn field must be present in nodes')

  if (any(nodes$cn<0 || (nodes$cn %% 1)!=0, na.rm = TRUE))
    stop('cn must be non-negative integers')

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

    if (any(edges$cn<0 || (edges$cn %% 1)!=0, na.rm = TRUE))
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
