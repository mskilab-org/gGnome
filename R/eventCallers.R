#' @name proximity
#' @export
#' @rdname internal
#' proximity
#'
#' Takes a set of n "query" elements (GRangse object, e.g. genes) and determines their proximity to m "subject" elements
#' (GRanges object, e.g. regulatory elements) subject to set of rearrangement adjacencies (GRangesList with width 1 range pairs)
#'
#' This analysis makes the (pretty liberal) assumption that all pairs of adjacencies that can be linked on a gGraph path are in
#' cis (i.e. share a chromosome) in the tumor genome.
#'
#' Each output proximity is a gWalk that connects query-subject on the genome
#' described by gGraph gg.  Each gWalk is  annotated by the metadata of the
#' corresponding query-subject GRanges pair as well as fields "altdist" and "refdist"
#' specifying the "alternate and "reference" gGraph distance of the
#' query-subject pair.  The gWalk metadata field "reldist" specifies the relative
#' distance (i.e. ratio of altdist to refdist) for that walk. 
#' 
#' @param gg gGraph of the "alternate genome"
#' @param query GRanges of "intervals of interest" eg regulatory elements
#' @param subject GRanges of "intervals of interest" eg genes
#' @param ref gGraph of the "reference genome", by derigma is the reference genome but can be any gGraph
#' @param ignore.strand whether to ignore strand of input GRanges
#' @param verbose logical flag
#' @param mc.cores how many cores to use for the path exploration step or if chunksize is provided, across chunks (default 1)
#' @param chunksize chunks to split subject and query into to minimize memory usage, if mc.cores>1 then each chunk will be allotted a core
#' @param max.dist maximum genomic distance to store and compute (1MB by default) should the maximum distance at which biological interactions may occur
#' 
#' @return gWalk object each representing a proximity
proximity = function(gg,
                     query,
                     subject,
                     reduce = TRUE,
                     ignore.strand = TRUE,
                     verbose = F,
                     mc.cores = 1,
                     strict.ref = FALSE,
                     chunksize = NULL,
                     max.dist = 1e6)
{
    if (length(gg)==0)
        return(gW(graph = gg))
    
    if (!ignore.strand)
        stop('strand-aware proximity is TBD')
    
    if (length(query)==0 | length(subject)==0)
        return(gW(graph = gg))

    if (!is.null(chunksize)) ## recursively run proximity on subchunks of query and subject
    {
        numqchunks = ceiling(length(query)/chunksize)
        numschunks = ceiling(length(subject)/chunksize)
        qmap = data.table(qid = 1:length(query), chunkid = rep(1:numqchunks, each = chunksize)[1:length(query)])
        smap = data.table(sid = 1:length(subject), chunkid = rep(1:numschunks, each = chunksize)[1:length(subject)])

        setkey(qmap, chunkid)
        setkey(smap, chunkid)

        queue = expand.grid(qchunkid = unique(qmap$chunkid), schunkid = unique(smap$chunkid))

        ## we have to return grl so we can concatenate since each px object will be on a different graph
        mc.cores2 = ceiling(mc.cores/nrow(queue))

        grls = mclapply(1:nrow(queue), function(i)
        {
            if (verbose)
                message(sprintf('proximity queue %s of %s total', i, nrow(queue)))
            qchunkid = queue$qchunkid[i]
            schunkid = queue$schunkid[i]
            qix = qmap[.(qchunkid), qid]
            six = smap[.(schunkid), sid]
            px = proximity(gg, query[qix], subject[six], reduce = reduce, ignore.strand = ignore.strand,
                           verbose = verbose, mc.cores = mc.cores2, strict.ref = strict.ref, max.dist = max.dist)

            if (length(px)==0)
                return(GRangesList())
            
            ## mod qid and sid to "global" values      
            px$set(qid = qix[px$dt$qid], sid = six[px$dt$sid])
            return(px$grl)
        }, mc.cores = mc.cores)

        ls = sapply(grls, length)
        grls = grls[ls>0]

        if (length(grls)==0)
            return(gW(graph = gg))

        grl = do.call(grl.bind, grls)

        ## then we rethread on the original graph attaching metadata
        px = gW(grl = grl, graph = gg)

        px = px[order(dist), ]
        return(px)
    }

    if (is.null(names(query)))
        names(query) = 1:length(query)

    if (is.null(names(subject)))
        names(subject) = 1:length(subject)

    ra = gg$edges[type == 'ALT', ]$junctions$grl

    query.og = query
    subject.og = subject

    query$id = 1:length(query)
    subject$id = 1:length(subject)

    qix.filt = 1:length(query)
    six.filt = 1:length(subject)

    if (!is.infinite(max.dist))
    {
        qix.filt = gr.in(query, unlist(ra)+max.dist) ## to save time, filter only query ranges that are "close" to RA's
        six.filt = gr.in(subject, unlist(ra)+max.dist) ## to save time, filter only query ranges that are "close" to RA's
    }

    query = query[qix.filt] 
    subject = subject[six.filt]

    if (verbose)
        message(length(query), ' query and ', length(subject), ' subject ranges after filtering by max.dist from an ALT junction')

    if (length(query)==0 | length(subject)==0)
        return(gW(graph = gg))

  px.gg = gg$copy$disjoin(gr = sort(unique(grbind(query[, c()], subject[, c()]))), collapse = FALSE)

  if (verbose)
    message('disjoined graph has ', dim(px.gg)[1], ' nodes and ', dim(px.gg)[2], ' edges.')
  
  if (strict.ref) ## only use actual REF edges in graph
  {
    px.gg.ref = px.gg[, type == 'REF']
  }
  else ## create generic ref with current nodes following reference sequence
  {
    px.gg.ref = gG(breaks = px.gg$gr[, c()])
  }

  ## map query and subject ends onto disjoin graph
  query.alt = query %*% px.gg$nodes$gr[, 'node.id']
  subject.alt = subject %*% px.gg$nodes$gr[, 'node.id']

  query.ref = query %*% px.gg.ref$nodes$gr[, 'node.id']
  subject.ref = subject %*% px.gg.ref$nodes$gr[, 'node.id']

  if (verbose)
    message(sprintf('Running $dist on %s query and %s subject nodes lifted to proximity graph',
                    length(query.alt), length(subject.alt)))


  D.alt = px.gg$dist(query.alt$node.id, subject.alt$node.id)

  if (verbose)
    message(sprintf('Running $dist on %s query and %s subject nodes lifted to reference proximity graph',
                    length(query.ref), length(subject.ref)))

  D.ref = px.gg.ref$dist(query.ref$node.id, subject.ref$node.id)

  dt.alt = spmelt(D.alt, Inf)
  dt.alt[, qnid := px.gg$gr$snode.id[query.alt$node.id[i]]]
  dt.alt[, snid := px.gg$gr$snode.id[subject.alt$node.id[j]]]
  dt.alt[, qid := query.alt$id[i]][, sid := subject.alt$id[j]]
  
  dt.ref = spmelt(D.ref, Inf)
  dt.ref[, qid := query.ref$id[i]][, sid := subject.ref$id[j]]

  ## mark the node pair with minimum ALT distance for each qid and sid pair
  ## (there could be multiple chunks for each query and subjec and
  ## multiple haplotypes for complex graph with bubbles)  
  dt.alt[, is.min := 1:.N %in% which.min(val), by = .(qid, sid)]
  dt.alt = dt.alt[is.min == TRUE, ]
  dt.ref[, is.min := 1:.N %in% which.min(val), by = .(qid, sid)]
  dt.ref = dt.ref[is.min == TRUE, ]  

  setnames(dt.alt, 'val', 'altdist')
  setnames(dt.ref, 'val', 'refdist')

  ## merge ALT and REF distances using the original query and subject id as the key
  dt = merge(dt.alt, dt.ref[, .(qid, sid, refdist)], by = c('qid', 'sid'), allow.cartesian = TRUE, all = TRUE)
  dt[is.na(refdist), refdist := Inf]

  ## subset to qid / sid pairs that are closer in ALT genome than in reference
  dt = dt[altdist<pmin(max.dist, refdist), ]

  if (nrow(dt)==0) ## return blank graph if no proximal pairs
    return(gW(graph = gg))

  dt[, reldist := altdist/refdist]

  if (verbose)
    message(sprintf('Calculating %s ALT query to subject paths',
                    nrow(dt)))

  ## now let's gather the paths connecting these node pairs
  px = px.gg$paths(dt$qnid, dt$snid, cartesian = FALSE, mc.cores = mc.cores)

  if (verbose)
    message(sprintf('Populating metadata for %s paths', length(px)))

  ## reinstantiate px with this meta

  meta = cbind(px$dt[, .(source, sink, dist)],
               dt[ , .(qid, sid, reldist, altdist, refdist)])

  if (ncol(values(query.og))>0)
    meta = cbind(meta, values(query.og)[meta$qid, , drop = FALSE])

  if (ncol(values(subject.og))>0)
    meta = cbind(meta, values(subject.og)[meta$sid, , drop = FALSE])

  px = gW(snode.id = px$snode.id, graph = px$graph, meta = meta)

  px = px[order(dist), ]
  px$set(name = 1:length(px))
  return(px)
}



#' @name fusions
#' @description
#' fusions
#'
#' annotates all gene fusions in gGraph relative to cds definitions
#' 
#' cds = gencode cds GRanges gff3 / gtf or GRangesList the latter (converted via rtracklayer::import)
#' which has fields $exon_number
#'
#' "gene_name" GRangesList meta data field is used in annotation and in not creating "splice" fusions that arise from different transcripts of the same gene.
#'
#' @param graph input gGraph 
#' @param gencode  GFF containing gene boundaries and exons, in similar format to  https://www.gencodegenes.org/ 
#' @param query optional query limiting walks to specific regions of interest
#' @param prom.window window to use around each transcript to identify putative promoter if promoter is NULL
#' @return gWalks of gene fusions annotated with frame and gene(s)
#' @export
fusions = function(graph = NULL,
                   gencode = NULL,
                   genes = NULL,
                   annotate.graph = TRUE,  
                   mc.cores = 1,
                   verbose = FALSE)
{
  ## QC input graph or junctions
  if (!inherits(graph, "gGraph")){
    stop('Input must be be gGraph')
  }

  if (is.character(gencode) && file.url.exists(gencode))
  {
    gencode = rtracklayer::import(gencode)
  }

  GENCODE.FIELDS = c('type', 'transcript_id', 'gene_name', 'exon_number', 'exon_id')
  GENCODE.TYPES = c('CDS', 'UTR', 'exon', 'gene', 'start_codon', 'stop_codon', 'transcript')
  if (!inherits(gencode, 'GRanges') || !all(GENCODE.FIELDS %in% names(values(gencode))) || !all(GENCODE.TYPES %in% gencode$type))
    stop(sprintf('gencode argument must be either URL or path to a valid GENCODE gff or a GRanges object with the following metadata fields: %s\n and with field "type" containing elements of the following values: %s', paste(GENCODE.FIELDS, collapse = ', '), paste(GENCODE.TYPES, collapse = ', ')))

 
    tgg = make_txgraph(graph, gencode)
    if (is.null(tgg)){
        return(gWalk$new())
    }
  
  txp = get_txpaths(tgg, genes = genes, mc.cores = mc.cores, verbose = verbose)
  if (length(txp)>0)
    txp$set(fclass = 'path')
  
  txl = get_txloops(tgg, txp, genes = genes, mc.cores = mc.cores, verbose = verbose)

  if (length(txl)>0)
    txl$set(fclass = 'loop')

  ## concatenate and annotate loops and walks
  allp = annotate_walks(c(txp, txl))
  
  ## thread these annotated walks back on the original graph
  gw = gW(grl = allp$grl, meta = allp$meta, graph = graph)

  ## split graph further using gencode features --> mainly for cosmetics
  ## i.e. so that now walks will have embedded intersecting
  ## transcript / CDS features
  if (annotate.graph && length(gw)>0)
  {
    if (verbose)
      {
        message('Annotating gGraph with GENCODE elements')
      }
    ugene = unique(allp$nodes$dt$gene_name)
    gt = gt.gencode(gencode %Q% (gene_name %in% ugene))
    annotations = dat(gt)[[1]] ## stored in gTrack
    gw$disjoin(gr = annotations)
    gw$set(name = gw$dt$genes)
    gw$graph$set(colormap = colormap(gt)[1])

    ## make the intergenic nodes gray
    gw$graph$nodes[is.na(type)]$mark(col = 'gray')
  }

  return(gw)
}

#' @name make_txgraph
#' @name description
#'
#' Generates "transcript graph" from gGraph and GENCODE GRanges
#' The transcript graph is essentially a "natural join" of original gGraph 
#' with the set of all transcript that are broken by one or more ALT junctions.
#' Every node of the transcript graph is annotated by tx_strand and associated
#' with a 5p and 3p frame, cDNA coordinate, protein coordinate.  These nodes
#' are later composed into standard paths and paths with "loops" which are finally annotated
#' i.e. described with respect to their gene fusion patterns, frameness, and protein coordinate
#' fragments.
#'
#' "Intergenic nodes" are also incorporated, but these are strictly defined as those that
#' lack (protein coding) transcripts.
#' (Note: this may exclude fusion genes that "skip over" whole transcripts that themselves
#' are unrearranged, though such events would be require aberrant / alternate splicing
#' in addition to the somatic fusion).  
#' 
#' This is an internal function used in fusions upstream of get_txpaths and get_txloops.
#' 
#' @keywords internal
#' @param tgg gGraph output of get_txgraph
#' @param mc.cores number of cores across which to parallelize path search component of algorithm
#' @param verbose whether to provide verbose output
#' @return gWalk of paths representing derivative transcripts harboring "loops"
#' @author Marcin Imielinski
make_txgraph = function(gg, gencode)
  {
    tx = gencode %Q% (type == 'transcript')

    ## broken transcripts intersect at least one junction
    tx$in.break = tx %^% unlist(gg$junctions$grl)
    if (!any(tx$in.break)){
        warning("No breakpoint in any transcript.")
        return(NULL)
    }

    txb = tx[tx$in.break]

    ## we sort all cds associated with broken transcripts using their stranded coordinate value
    cds = gencode %Q% (type == 'CDS') %Q% (transcript_id %in% txb$transcript_id)

    ## remove any transcripts that lack a CDS (yes these exist)
    txb = txb %Q% (transcript_id %in% cds$transcript_id)

    ## now supplement cds with txends
#    cds = grbind(cds, gr.start(txb), gr.end(txb))

    cds = cds %Q% order(ifelse(strand == '+', 1, -1)*start)

    ## the "phase" of a CDS is a tricky concept, we want to convert this into the left / right and 5' / 3' "frame"
    ## of the two cds ends
    ## phase is described here http://gmod.org/wiki/GFF, quoting
    ## "For features of type "CDS", the phase indicates where the feature begins with reference to the reading frame.
    ## The phase is one of the integers 0, 1, or 2, indicating the <<number of bases that should be removed from the
    ## beginning of this feature to reach the first base of the next codon>>. In other words, a phase of "0" indicates
    ## that the next codon begins at the first base of the region described by the current line, a phase of "1"
    ## indicates that the next codon begins at the second base of this region, and a phase of "2" indicates
    ## that the codon begins at the third base of this region.
    ## This is NOT to be confused with the frame, which
    ## is simply start modulo 3. If there is no phase, put a "." (a period) in this field.
    ## For forward strand features, phase is counted from the start field. For reverse strand features,
    ## phase is counted from the end field. The phase is required for all CDS features."

    ## ok now we convert phase to frame 
    ## frame for cds coordinates is 0 for the first nucleotide in the codon, 1 for the second, 2 for the third
    ## (according to the above description frame) is just negative phase %% 3
    ## we convert "phase" (which refers to the phase of the leftmost base in the exon (with respect
    ## to the reference) into the first 5' base of the given exon with respect to the
    ## transcript into "left" and "right" frame with representing the frame of the bases on each side of the
    ## CDS exon with respect to the reference
    cds$phase = as.numeric(cds$phase)
    cds$fivep.frame = -cds$phase %% 3
    cds$threep.frame = (cds$fivep.frame + width(cds) -1 ) %% 3

    cds$left.frame = ifelse(strand(cds) == '+', cds$fivep.frame, cds$threep.frame)
    cds$right.frame = ifelse(strand(cds) == '+', cds$threep.frame, cds$fivep.frame)

    cds$exon_number = as.numeric(cds$exon_number)
    
    ## annotate protein and transcript coordinates of each spliced exon
    ## this is tricky since there are sometimes gaps in between exons .. eg in refseq
    ## i.e. the phase of the last base in the previous exon is not 1 + phase of the first 
    ## base in the last exon
    ## we quantify this in the "gap" i.e. distance between first and last exons
    cdsdt = gr2dt(cds)
    cdsdt[type == "CDS", gap := c(0, (fivep.frame[-1] - 1) %% 3 - threep.frame[-.N]), by = transcript_id]
    cdsdt[type == "CDS", fivep.cc := cumsum(c(1+fivep.frame[1], (width+gap)[-.N])), by = transcript_id]
    cdsdt[, threep.cc := fivep.cc + width-1]
    cdsdt[, fivep.pc := fivep.cc/3] ## 5' protein coordinate
    cdsdt[, threep.pc := threep.cc/3] ## 3' protein coordinate
    cdsdt[, left.cc := ifelse(strand == '+', fivep.cc, threep.cc)]  ## left cds coord 
    cdsdt[, right.cc := ifelse(strand == '+', threep.cc, fivep.cc)] ## right cds coord
    cdsdt[, left.pc := left.cc/3] ## left protein coordinate
    cdsdt[, right.pc := right.cc/3] ## right protein coordinate
    cdsdt[type != "CDS", exon_number := c(0, Inf), by = transcript_id] ## only two non CDS "ends" per transcript
    cdsdt[, is.start := exon_number == min(exon_number), by = transcript_id]
    cdsdt[, is.end := exon_number == max(exon_number), by = transcript_id]

    values(cds) = cbind(values(cds), cdsdt[, .(gap, fivep.cc, threep.cc, fivep.pc, threep.pc, left.cc, right.cc, left.pc, right.pc, is.start, is.end)])
   
    ## now do a left merge of gg nodes with broken transcripts ..
    ## where we keep our original gg nodes but just duplicate
    ## them based on transcript intersections
    nov = gg$nodes$gr %*% txb
    txnodes = unname(gg$nodes$gr[nov$query.id])
    txnodes$tx_strand = as.character(strand(txb)[nov$subject.id])
    values(txnodes) = cbind(values(txnodes), values(nov)[, c('transcript_id', 'gene_name', 'gene_id')])

    ## reorder the txnodes so they are in the direction of the given transcript
    tmpdt = gr2dt(txnodes[, c('transcript_id', 'tx_strand')])[, id := 1:.N][, start := ifelse(tx_strand == '+', start, -start)]
    setkeyv(tmpdt, c('transcript_id', 'seqnames', 'start'))

    txnodes = txnodes[tmpdt$id]

    ## compute start and end phase i.e. frame of txnodes
    ## by crossing with CDSs
    cdsov = gr2dt(gr.findoverlaps(txnodes, cds, by = 'transcript_id', qcol = 'tx_strand', scol = names(values(cds))))    
    
    ## for each query node we only keep the first or last cds exon
    cdsov[, is.min := exon_number == min(exon_number), by = .(query.id)]
    cdsov[, is.max := exon_number == max(exon_number), by = .(query.id)]
    cdsov = cdsov[is.min | is.max, ]

    ## we need to shift the left and right frame based on 
    ## the distance from the edge of that cds exon 
    ## (which will be zero unless the the right or left end of a node
    ## lies inside that cds exon)
    cdsov[, left.frame := (left.frame + start - start(cds)[subject.id]) %% 3]
    cdsov[, right.frame := (right.frame + end - end(cds)[subject.id]) %% 3]

    ## left cds exons are the min exon_number for + exons, and vice versa for negative exons
    cdsov[, is.left := ifelse(as.logical(strand(cds)[subject.id]=='+'), is.min, is.max)]
    cdsov[, is.right := ifelse(as.logical(strand(cds)[subject.id]=='+'), is.max, is.min)]
    cdsov[, is.threep := is.max]
    cdsov[, is.fivep := is.min]

    ## now merge this info back into txnodes
    txnode.ann = cdsov[, .(
      fivep.coord = ifelse(tx_strand[is.fivep]=='+', start[is.fivep], end[is.fivep]),
      threep.coord = ifelse(tx_strand[is.threep]=='+', end[is.threep], start[is.threep]),
      fivep.coord = start[is.fivep],
      fivep.frame = fivep.frame[is.fivep],
      threep.frame = threep.frame[is.threep],
      fivep.exon = exon_number[is.fivep],
      threep.exon = exon_number[is.threep],
      fivep.cc = fivep.cc[is.fivep],
      threep.cc = threep.cc[is.threep],
      fivep.pc = fivep.pc[is.fivep],
      threep.pc = threep.pc[is.threep],
      is.start = is.start[is.fivep],
      is.end = is.end[is.threep]
    ), keyby = query.id][.(seq_along(txnodes)), ]

    ## we need to fill in information for (internal) interexonic txnodes
    ## i.e. those don't contain exons .. which can be quite common
    ## we want to map them to their preceeding exon and then let them inherit
    ## the 3' properties of 

    ## here we number all runs of NA distinguish between transcripts that have NA runs at the very beginning
    ## (which we label 0) vs those that have their first NA run in the middle of the transcript
    txnode.ann[, transcript_id := txnodes$transcript_id[query.id]]
    txnode.ann[, na.run := ifelse(is.na(fivep.frame[1]), -1, 0) + label.runs(is.na(fivep.frame)),
               by = transcript_id]

    ## we also label all non NA runs
    txnode.ann[, nna.run := label.runs(!is.na(fivep.frame)), by = transcript_id]

    ## now label the fivep end of the transcript (different from is.start, which is CDS start)
    txnode.ann[, is.txstart := 1:.N %in% which.min(pmin(na.run, nna.run, na.rm = TRUE)), by = transcript_id]

    ## for each non NA run we collect the  <last> row and make a map
    ## we will key this map using na.runs, i.e. matching each na run to their last nna.run
    nna.map = txnode.ann[!is.na(nna.run), .(qid.last = query.id[.N]), keyby = .(k1 = transcript_id, k2 = nna.run)]

    ## now associate all NA runs with the features of their the last NNA run i.e. qid.last
    txnode.ann[!is.na(na.run), qid.last := nna.map[.(transcript_id, na.run), qid.last]]

    ## note that all the na.run == 0 will still have have an NA qid.last, which we want
    ## since they are 5' UTR should not be incorporated into any protein coding fusion
    ## but we need to also NA out qid.lasts for NA.runs that are 3' UTR, i.e. downstream of
    ## the translation end
    txnode.ann[!is.na(qid.last), qid.last := as.integer(ifelse(txnode.ann$is.end[qid.last], NA, qid.last))]

    ## now relabel both 5p and 3p internal txnode.ann features of internal NA runswith the <three prime>
    ## features of their qid.last

    txnode.ann[, is.cds := is.na(na.run)]
    txnode.ann[!is.na(na.run) & !is.na(qid.last),
               ":="(
                 fivep.frame = txnode.ann$threep.frame[qid.last],
                 threep.frame = txnode.ann$threep.frame[qid.last],
                 fivep.exon = txnode.ann$threep.exon[qid.last], 
                 threep.exon = txnode.ann$threep.exon[qid.last],
                 fivep.cc = txnode.ann$threep.cc[qid.last], 
                 threep.cc = txnode.ann$threep.cc[qid.last],
                 fivep.coord = txnode.ann$threep.coord[qid.last], 
                 threep.coord = txnode.ann$threep.coord[qid.last],
                 is.start = FALSE,
                 is.end = FALSE,
                 fivep.pc = txnode.ann$threep.pc[qid.last],
                 threep.pc = txnode.ann$threep.pc[qid.last])]
    txnode.ann[, twidth := ifelse(is.cds, (threep.cc - fivep.cc + 1)/3, 0)]

    values(txnodes) = cbind(values(txnodes),
                            txnode.ann[, .(is.cds, fivep.coord, threep.coord, fivep.frame,
                                           threep.frame, fivep.exon,
                                           threep.exon, fivep.cc, threep.cc, fivep.pc, threep.pc,
                                           is.txstart, is.start, is.end, twidth)])

    ## all the other nodes in the graph, which we include in case
    ## we have intergenic "bridging nodes" connecting different fusionsbu
    igrnodes = gg$nodes$gr[setdiff(1:length(gg), nov$query.id)]

    ## intergenic nodes should not intersect any transcript
    ## (may want to change in the future, since perhaps complex fusions could
    ## "skip over" transcripts ... though these would technically be splicing events)
    igrnodes = igrnodes[!(igrnodes %^% tx)] 

    newnodes = grbind(txnodes, igrnodes)
    newnodes$parent.node.id = newnodes$node.id
    newnodes$node.id = 1:length(newnodes)

    ## map of old nodes to new nodes
    nmap = data.table(new.node.id = newnodes$node.id, node.id = newnodes$parent.node.id, transcript.id = newnodes$transcript_id, tx_strand = newnodes$tx_strand)

    ## now merge edges using nmap
    edges = gg$edges$dt
    newedges = merge(merge(edges, nmap, by.x = 'n1', by.y = 'node.id', allow.cartesian = TRUE), nmap, by.x = 'n2', by.y = 'node.id', allow.cartesian = TRUE)

    ## only keep edges if
    ## (1) transcript_id.x or transcript_id.y are NA
    ## (2) transcript_id.x and transcript_id.y are both non-NA and equal
    ## (3) edge is ALT

    newedges = newedges[type == 'ALT' | is.na(transcript.id.x) | is.na(transcript.id.y) |
                        !is.na(transcript.id.x) & !is.na(transcript.id.y) & transcript.id.x == transcript.id.y, ]
    newedges[, n1 := new.node.id.x]
    newedges[, n2 := new.node.id.y]

    #' a transcript associated edge will NOT
    #' connect to a + is_start from the left or
    #' connect to a + is_end from the right or
    #' connect to a - is_start from the right or
    #' connect to a + is_end from the left
    #' ie we don't want to see an edge from the that is being "ended"

    newedges[, transcript_associated :=
                 !(is.na(tx_strand.x) | is.na(tx_strand.y))]

    ## make a note of splice variant edge
    ## i.e. those that connect two different transcripts of the same gene
    ## we will not be interested in these
    txdt = gr2dt(txb[, c('transcript_id', 'gene_name')])
    setkey(txdt, transcript_id)
    newedges$gene_name.x = txdt[.(newedges$transcript.id.x), gene_name]
    newedges$gene_name.y = txdt[.(newedges$transcript.id.y), gene_name]
    newedges[transcript_associated == TRUE,
             splicevar := gene_name.x == gene_name.y & transcript.id.x != transcript.id.y]


    ## first annotate whether transcript-associated edges emerge from 5' or 3' of node
    newedges[type == 'ALT' & transcript_associated,
             ":="(
               n1.threep = tx_strand.x == '+' & n1.side == 'right' |
                 tx_strand.x == '-' & n1.side == 'left',
               n2.threep = tx_strand.y == '+' & n2.side == 'right' |
                 tx_strand.y == '-' & n2.side == 'left')
             ]

    ## antisense = any ALT transcript associated edges that have 3' to 3' or 5' to 5' connections
    newedges[type == 'ALT' & transcript_associated, antisense := n1.threep == n2.threep]

    ## deadend = any ALT edges that emanate from the "stopped side"
    ## of a is.end == TRUE node ... any transcript that incorporates the 
    ## given node would not translate through the stop, i.e. would be
    ## truncated at that end, hence any path that proceeds beyond that 
    ## would have the remaining sequences truncated
    ## (*** unless we allow for splice variants in the rearranged allele
    ## ie those that skip over subsets of the included exons ... but this will
    ## greatly increase the combinatoric complexity of the allelic search ...
    ## ... hmm I suppose we could fix this by creating a non transcript_associated
    ## node for every unique transcript associated node ..  
    ## TBD FIXME)
    newnodes$end_dir = ifelse(newnodes$is.end, ifelse(newnodes$tx_strand == '+', 'right', 'left'), NA)
    newedges[,
             deadend := 
               newnodes$end_dir[n1] == n1.side |
               newnodes$end_dir[n2] == n2.side
             ]

    ## we also dead all starts, i.e. we don't allow edges that enter a node containing
    ## the 5p UTR of a transcript or a CDS start
    newnodes$start_dir = ifelse(newnodes$is.txstart | newnodes$is.start,
                         ifelse(newnodes$tx_strand == '+', 'left', 'right'), NA)
    newedges[,
             deadstart := 
               newnodes$start_dir[n1] == n1.side |
               newnodes$start_dir[n2] == n2.side
             ]

    ## in-frame = !antisense and the corresponding frame of the 5' node is 1 + the frame of the 3' node

    ## first annotate the frame of n1 and n2
    newedges[type == 'ALT' & transcript_associated & !antisense & !deadend,
             ":="(
               n1.frame = ifelse(n1.threep, newnodes$threep.frame[n1], newnodes$fivep.frame[n1]),
               n2.frame = ifelse(n2.threep, newnodes$threep.frame[n2], newnodes$fivep.frame[n2])
             )]

    newedges[type == 'ALT' & transcript_associated & !antisense & !deadend,
             in.frame := ifelse(n1.threep,
                                n2.frame == (n1.frame + 1) %% 3, ## n1 3'-> n2 5'
                                n1.frame == (n2.frame + 1) %% 3) ## n1 5'--> n2 3'
             ]

    ## set up edge weights on which we will run our shortest paths analysis

    ## just remove anti-sense and dead end edges
    ## (note: won't make our graph immune to these anti-sense paths via bridge nodes, see below)
    newedges = newedges[-which(deadend | deadstart | antisense | splicevar), ] 

    ## weigh edge
    ## (1) REF edges weigh = 1
    ## (2) ALT edges between non NA transcripts
    ##     (b) if sense and out of frame 100
    ##     (c) if sense and in frame 0
    ## (3) ALT edges between at least one NA transcript
    ##    ... just weigh 1

    newedges[, eweight := ifelse(type == 'ALT', 1, 10)] ## prefer ALT vs REF
    newedges[antisense == FALSE & in.frame == FALSE, eweight := 100] ## smaller out of frame penalty
    newedges[type == 'ALT' & in.frame == TRUE, eweight := 0] ## reward in frame ALT 

    ## FIXME: how do we enable non-transcript associated nodes i.e. "bridge nodes" in the graph to "transmit"
    ## the information about frame / strand from their transcript associated neighbors?
    ## might be an impossible problem given the current "dynamic programming formulation" involving
    ## "independent" edge weights, since edges entering and exiting non-transcript associated 
    ## nodes may have a different sense / frame status based on the path that they belong to ...
    ## i.e. that edge might be out of frame or antisense in one context but in frame or in-context in another ...
    ## we can of course verify frame-ness and sense-ness at the downstream "walk annotation" stage
    ## but there is still a risk is thus that we picked the wrong walk, i.e. we picked a "shortest path"
    ## with respect to the current formulation edge weight set that
    ## is hopelessly antisense or out-of-frame when for that given CDS start-end pair a perfectly
    ## good in-frame in-sense fusion exists
    ##
    ## currently we make the edge weights connecting to non-transcript associated nodes neutral,
    ## however one heuristic to minimize the occurrence of the above scenario would be to
    ## try to optimize the frame / tx_strand assignment of a non-transcript associated node
    ## i.e. propagate the "best" assignment according to some greedy heuristic

    tgg = gG(nodes = newnodes, edges = newedges)

    return(tgg)
  }

#' @name get_txpaths
#' @name description
#'
#' gets all transcript paths from a transcript graph (i.e. output of make_txgraph)
#' that connect CDS start and CDS end
#' 
#' 
#' @keywords internal
#' @param tgg gGraph output of get_txgraph
#' @param mc.cores number of cores across which to parallelize path search component of algorithm
#' @param verbose whether to provide verbose output
#' @return gWalk of txPaths
#' @author Marcin Imielinski
get_txpaths = function(tgg,
                       genes = NULL,
                       mc.cores = 1,
                       verbose = FALSE)
{
    starts = tgg$nodes$dt[is.start==TRUE, ifelse(tx_strand == '+', 1, -1)*node.id]
    ends = tgg$nodes$dt[is.end==TRUE, ifelse(tx_strand == '+', 1, -1)*node.id]

    if (!is.null(genes))
    {
      starts = starts[tgg$nodes[starts]$dt$gene_name %in% genes]
      ends = ends[tgg$nodes[ends]$dt$gene_name %in% genes]
    }    

    INF = 1e7
    ab.p = gW(graph = tgg)

    if (length(starts)>0 & length(ends)>0)
    {
      ## generate walks from aberrant shortest paths
      D = tgg$dist(starts, ends, weight = 'eweight', ignore.strand = FALSE)
        combos = as.data.table(which(D<INF & D>0, arr.ind = TRUE))
        
      combos[, start := starts[row]][, end := ends[col]][, dist := D[cbind(row, col)]]
      setkey(combos, 'start')
      
      ustart = unique(combos$start)
      if (nrow(combos)>0)
      {
        p = mclapply(1:length(ustart), function(k)
        {
          if (verbose)
            {
              message('traversing ', k, ' of ', length(ustart), ' candidate CDS start-end pairs')
            }

          p = tgg$paths(ustart[k],
                        combos[.(ustart[k]), end],
                        weight = 'eweight',
                        ignore.strand = FALSE)
          
          return(p)
        }, mc.cores = mc.cores,
        mc.preschedule = FALSE)
        
        paths = do.call('c', c(p, list(force = TRUE)))[dist<INF & length>1, ]
        
        ab.p = tryCatch(
        {
          paths$set(numchr = paths$eval(node = length(unique(seqnames))))
          paths$set(numab = paths$eval(edge = sum(type == 'ALT')))
          paths$set(numgenes = paths$eval(node = length(unique(gene_name[!is.na(gene_name)]))))
          paths$set(genes = paths$eval(node = paste(unique(gene_name[!is.na(gene_name)]), collapse = ',')))
          paths$set(maxcn = paths$eval(edge = min(cn)))
          ab.p = paths[numab>0] ## only include walks that contain one or more aberrant edges
          ## remove cryptic antisense paths
          ## (this can happen in highly rearranged genomes even in the absence of antisense
          ## edges, via ALT edges that go intergenic)
          ab.p = ab.p[!ab.p$eval(any(tx_strand != strand, na.rm = TRUE))]        
        }, error = function(e) paths)
      }
    }
    return(ab.p)
}
    
#' @name get_txloops
#' @name description
#'
#' gets internal "loops" that connect the right side of downstream nodes of
#' transcripts to left side of upstream nodes in the same transcript
#'
#' internal function used in fusions downstream of get_txpaths.
#' 
#' @keywords internal
#' @param tgg gGraph output of get_txgraph
#' @param mc.cores number of cores across which to parallelize path search component of algorithm
#' @param verbose whether to provide verbose output
#' @return gWalk of paths representing derivative transcripts harboring "loops"
#' @author Marcin Imielinski
get_txloops = function(tgg, 
                       other.p = NULL, genes = NULL,
                       mc.cores = 1, verbose = FALSE)
{       
    ## grab reference paths in tgg
    txgr = tgg$nodes$gr %Q% (!is.na(transcript_id)) %Q% (order(fivep.exon))

    if (!is.null(genes))
      {
        txgr = txgr %Q% (gene_name %in% genes)
      }

    snode.id = split(ifelse(txgr$tx_strand == '+', 1, -1) * txgr$node.id, txgr$transcript_id)
    ref.p  = gW(snode.id = snode.id, graph = tgg, drop = TRUE, meta = data.table(txid = names(snode.id)))

    if (!is.null(other.p))
      ref.p = c(ref.p, other.p)

    ## now identify reference path nodes that have an ALT edge
    ## leaving their right side 
    right = unique(ref.p$nodes$eright[type == 'ALT']$left$dt$snode.id)
    left = unique(ref.p$nodes$eleft[type == 'ALT']$right$dt$snode.id)
    
    ## exclude 3p edges associated with the last node in any transcript
    dt.right = ref.p$nodesdt[snode.id %in% right, ][walk.iid != lengths(ref.p$snode.id)[walk.id], ]
    dt.left = ref.p$nodesdt[snode.id %in% left, ]
    
    dtm = merge(dt.right, dt.left, by = "walk.id", allow.cartesian = TRUE)[walk.iid.x>=walk.iid.y, ]
    dtm[, dtmid := 1:.N][, found := FALSE]
    
    ## see if any of these are simple connections in the grapha
    ## using drop = TRUE functionality of gW
    loop.simple = gW(snode.id = split(rep(dtm$snode.id.x, 2), rep(dtm$dtmid, 2)), graph = tgg, drop = TRUE)
    if (length(loop.simple)>0)
    {
      simple.dtmid = as.numeric(loop.simple$meta$name)
      dtm[simple.dtmid, ":="(loop = loop.simple$snode.id)]
      dtm[simple.dtmid, found := TRUE]
    }
    
    setkeyv(dtm, c('snode.id.x', 'snode.id.y'))
    
    dtm.complex = unique(dtm[found == FALSE, ], by = c('snode.id.x', 'snode.id.y'))
    
    if (nrow(dtm.complex)>0)
    {
      starts = unique(dtm.complex$snode.id.x)
      ends = unique(dtm.complex$snode.id.y)
      D = tgg$dist(starts, ends, weight = 'eweight', ignore.strand = FALSE)
      dtm.complex$dist = D[cbind(match(dtm.complex$snode.id.x, starts), match(dtm.complex$snode.id.y, ends))]
      
      tmp = mclapply(1:nrow(dtm.complex), function(k)
      {
        if (verbose)
          message('traversing ', k, ' of ', nrow(dtm.complex), ' candidate cycles')

        so = dtm.complex$snode.id.x[k]
        si = dtm.complex$snode.id.y[k]
        if (so != si)
        {
          p = tgg$paths(so, si, weight = 'eweight', ignore.strand = FALSE)
        }
        else
        {
          ## not too edge case where we have a non-simple cycle that begins and ends in the same place
          ## to get the cyclic path we have to find a path between the ALT right neighbor of so

          edt = tgg$nodes[so]$eright[type == 'ALT']$dt

          ## we need the signed version of the source
          ## ?FIXME = develop better syntax eg include signed nodes in $dt
          soc = ifelse(edt$n2.side == 'left', edt$n2, -edt$n2)
          p = tgg$paths(soc, si, weight = 'eweight',
                        ignore.strand = FALSE)
          if (length(p)>0)
          {
            p = gW(split(c(so, p$snode.id[[1]]), 1), graph = tgg,
                   meta = data.table(source = so, sink = si))
          }
        }   
        return(p)
      }, mc.cores = mc.cores)

      loops.complex = do.call('c', c(tmp, list(force = TRUE)))

      ## remove any weird loops that cause parts of transcripts to be reversed (i.e. antisense)
      ## (this can happen in highly rearranged genomes even in the absence of antisense
      ## edges, via ALT edges that go intergenic)
      if (length(loops.complex)>0)
      {
        loops.complex =
          tryCatch(loops.complex[!loops.complex$eval(any(tx_strand != strand, na.rm = TRUE))],
                   error = function(e) loops.complex)
          
          ldt = cbind(loops.complex$dt[, .(snode.id.x = source, snode.id.y = sink)],
                      data.table(snode.id = loops.complex$snode.id, complex = TRUE))
          dtm = merge(dtm, ldt, all.x = TRUE, by = c("snode.id.x", "snode.id.y"), allow.cartesian = TRUE)
          dtm[!is.na(complex), ":="(loop = snode.id)]
          dtm[!is.na(complex), found := TRUE]
        }
    }
    
    loops = dtm[found == TRUE, ]
  ab.l = gW(graph = tgg)
    if (nrow(loops)>0)
    {
      ## now we need to create "prefix" and "suffix" subwalks from all of the
      ## walks using the beginning and end loop.ids
      wdt = ref.p$nodesdt
      wdt$txid = ref.p$meta[.(wdt$walk.id), txid]
      setkeyv(wdt, c('walk.id', 'snode.id'))
      loops[, begin := sapply(loop, "[", 1)]
      loops[, end := sapply(loop, function(x) x[length(x)])]
      loops$begini = wdt[.(loops$walk.id, loops$begin), walk.iid]
      loops$endi = wdt[.(loops$walk.id, loops$end), walk.iid]
      ## prefix is from beginning of walk to end of loop-1
      loops[, prefix := mapply(function(x,y)
        if (x==0) c() else y[1:x], begini-1, ref.p$snode.id[walk.id], SIMPLIFY = FALSE)]
      loops[, suffix := mapply(function(x,y)
        if (x>length(y)) c() else y[x:length(y)], endi+1, ref.p$snode.id[walk.id], SIMPLIFY = FALSE)]

      loops[, snode.id := mapply("c", prefix, loop, suffix, SIMPLIFY = FALSE)]
      ab.l = gW(snode.id = loops$snode.id, graph = tgg)

      ## dedup any loops with identical node strings
      ab.l = ab.l[!duplicated(sapply(ab.l$snode.id, paste, collapse = ', '))]

      tryCatch(
        {
          ab.l$set(numchr = ab.l$eval(node = length(unique(seqnames))))
          ab.l$set(numab = ab.l$eval(edge = sum(type == 'ALT')))
          ab.l$set(numgenes = ab.l$eval(node = length(unique(gene_name[!is.na(gene_name)]))))
          ab.l$set(genes = ab.l$eval(node = paste(unique(gene_name[!is.na(gene_name)]), collapse = ',')))
          ab.l$set(maxcn = ab.l$eval(edge = min(cn)))
        }, error = function(e) NULL)
    }
  return(ab.l)
}

annotate_walks = function(walks)
{
  if (length(walks)==0)
    return(walks)
  
  ## reinstantiate to "peel apart" walks .. i.e. reinstantiate with separate graph
  walks = gW(grl = walks$grl, meta = walks$meta, disjoin = FALSE)
  
  ## mark nodes with their walk.id
  walks$nodes$mark(wkid = walks$nodesdt$walk.id)
  
  ## annotate somatic coordinates and frames of fused transcripts
  gr = walks$nodes$gr
  grs = paste(gr.string(gr), gr$transcript_id)
  gr$uid = as.integer(factor(grs))
  ndt = gr2dt(gr)[is.cds == TRUE, ] ## only include cds chunks  
  ndt[is.na(twidth), twidth := 0]
  ndt[, threep.sc := cumsum(twidth), by = wkid]
  ndt[, fivep.sc := threep.sc-twidth+1/3]
  ndt[, threep.sc.frame := (round(threep.sc*3)-1) %% 3]
  ndt[, in.frame := threep.sc.frame == threep.frame]

  ## collect all edges to the left (i.e. 5p) of each node
  edt = walks$nodes$eleft$dt
  setkey(edt, n2) ## index on n2
  ndt[, fivep.alt := edt[.(node.id), class]!='REF'] ## determine if REF
  ndt[is.na(fivep.alt), fivep.alt := FALSE]

  ## use this to mark all distinct "chunks" within a transcript
  ## (two or more chunks will be created as a result of a deletion or a
  ## duplication)
  ndt[, chunk := cumsum(fivep.alt), by = .(wkid, transcript_id)]

  ## now want to label the protein coordinates and exons of each transcript "chunk"
  ## in "grl.string" format ... i.e. something parseable by parse.grl
  cdt = ndt[!is.na(fivep.frame) & is.cds,
            .(
              in.frame = in.frame[1],
              exon.start = fivep.exon[1],
              exon.end = threep.exon[.N],
              cc.start = ceiling(fivep.cc[1]),
              cc.end = ceiling(threep.cc[.N]),
              pc.start = ceiling(fivep.pc[1]),
              pc.end = floor(threep.pc[.N]),
              uids = paste(uid[1], uid[.N])
            ),
            by = .(wkid, gene_name, transcript_id, chunk)]

  cdt[, gene_name := paste0(ifelse(in.frame, '', '['), gene_name, ifelse(in.frame, '', ']fs'))]
  cdt[, transcript_id := paste0(ifelse(in.frame, '', '['), transcript_id, ifelse(in.frame, '', ']fs'))]
  
  
  .del = function(tx, start, end, label)
  {
    N = length(tx)
    ret = as.character(NA)
    if (tx[1] == tx[N])
    {
      ret = '' ## an empty deletion signals to us a "silent" deletion
      ## figure out what's missing
      ix = tx == tx
      del = setdiff(IRanges(start[1], end[N]),
                    IRanges(start[ix],
                            end[ix]))
      if (length(del)>0)
      {
          ret = paste0(label, ':',
                       ## round(start(del)/3,1), '-', 
                       ## round(end(del)/3,1), collapse = ';')          
                    ceiling(start(del)/3), '-', 
                    floor(end(del)/3), collapse = ';')      

      }
    }
    return(ret)
  }
  
  
  .amp = function(tx, start, end, label, uids)
  {
    N = length(tx)
    ret = as.character(NA)
    if (any(dup.ix <- duplicated(uids)))
    {
      ## figure out what's present in this transcript more than once
      ret = paste0(label[dup.ix], ':', 
                start[dup.ix], '-', 
                        end[dup.ix], collapse = ';')
    }
    return(ret)
  }

  cdt[, num.splice := length(unique(transcript_id)), by = .(wkid, gene_name)]

  adt = cdt[, .(
    gene.pc = paste0(gene_name, ':', pc.start, '-', pc.end, collapse = ';'),  
    del.pc = .del(transcript_id, cc.start, cc.end, gene_name[1]),
    amp.pc = .amp(transcript_id, pc.start, pc.end, gene_name, uids),
    splice.variant = any(num.splice>1),
    in.frame = all(in.frame, na.rm = TRUE),
    qin.frame = in.frame[1] & in.frame[.N],
    tx.cc = paste0(transcript_id, ':', cc.start, '-', cc.end, collapse = ';'),
    tx.ec = paste0(transcript_id, ':', exon.start, '-', exon.end, collapse = ';')
  ),
  keyby = wkid][.(1:length(walks)), ]

  adt[, silent := ifelse(is.na(del.pc), FALSE,
                  ifelse(nchar(del.pc)>0, FALSE, TRUE))]  

  adt[, frame.rescue := !in.frame & qin.frame]
  
  newmeta = cbind(walks$meta, adt[, .(gene.pc, amp.pc, del.pc, in.frame, frame.rescue, tx.cc, tx.ec, silent, splice.variant)])

  ## now we want to trim the first and last intervals in the walks according to the
  ## fivep.coord and threep.coord
  ngr = walks$nodes$gr
  ngrdt = gr2dt(ngr)[, ":="(is.first = 1:.N %in% 1, is.last = 1:.N %in% .N), by = wkid]
  ngrdt[is.first==TRUE, start := ifelse(tx_strand == '+', fivep.coord, start)]
  ngrdt[is.first==TRUE, end := ifelse(tx_strand == '+', end, fivep.coord)]
  ngrdt[is.last==TRUE, start := ifelse(tx_strand == '+', start, threep.coord)]
  ngrdt[is.last==TRUE, end := ifelse(tx_strand == '+', threep.coord, end)]
  
  end(ngr) = ngrdt$end
  start(ngr) = ngrdt$start

  walks = gW(grl = split(ngr, ngr$wkid), disjoin = FALSE, meta = newmeta)

  walks = walks[order(!in.frame, !is.na(amp.pc), !is.na(del.pc), !frame.rescue, splice.variant, silent)]

  return(walks)
}


#' @name chromoplexy
#' @title
#'
#' Finds "clean" choromoplexy chains eg paths and cycles of junctions with span > min.span
#' by first identifying jbp pairs within max.insert distance for which the "left" jbp is + and "right" jbp is -
#' (in reference coordinates) excluding any jbp that have more than one ALT jbp
#' (or optionally loose ends) with min.cushion distance.  These remaining jbp pairs
#' are then combined into a graph, and connected components in that graph are
#' scraped for paths and cycles, which are marked on the graph and added as metadata
#' to the outputted gGraph
#'
#' To differentiate from tic, chromoplexy cycles must have at least one
#' (quasi) reciprocal ie -+ junction pair with in the insert size. Other parameters
#' and algorithm almost identical to tic, except for allowing (and requiring at least one)
#' -+ bp pair to connect junctions
#'
#' @param gg gGraph
#' @param max.insert max insert to consider in a templated insertion 1e5
#' @param min.span min span for a TIC junction (1e6)
#' @param min.cushion minimum cushion between a TIC junction and any other nearby event (to ensure "clean" events), the bigger the cushion, the cleaner the calls
#' @param ignore.loose.ends logical flag (FALSE) determining whether we ignore loose ends when filtering on min.cushion
#' @param ignore.small.dups logical flag (FALSE) determining whether we ignore small dups when filtering on min.cushion
#' @param ignore.small.dels logical flag (FALSE) determining whether we ignore small dels when filtering on min.cushion
#' @param max.small threshold for calling a local dup or del "small" 1e4
#' @return gGraph with $meta annotated with gWalks corresponding to tic and tip and nodes and edges labeled with 'p1' through 'pk' for all k templated insertion paths and 'c1' through 'ck' for all k templated insertion cycles
chromoplexy = function(gg, max.insert = 1e5,
               min.cushion = 1e6,
               min.span = 1e6,
               ignore.loose.ends = FALSE,
               ignore.small.dups = FALSE,
               ignore.small.dels = FALSE,
               max.small = 1e4,
               mark = FALSE,
               mark.col = 'purple'
               )
{
  gg = refresh(gg)

  ## set empty output - in case we find no events can quit early
  gg.empty = gg$copy
  gg.empty$edges$mark(chromoplexy = as.integer(NA))
  gg.empty$set(chromoplexy = data.table())

  ed = gg$edges[type == 'ALT']
  
  if (length(ed) && ignore.small.dups)
    ed = ed[!(class == 'DUP-like' & ed$span<=max.small)]

  if (length(ed) && ignore.small.dels)
    ed = ed[!(class == 'DEL-like' & ed$span<=max.small)]

  if (length(ed)==0)
    return(gg.empty)

  jbp = grl.unlist(ed$grl)[, c('edge.id', 'grl.iix', 'snode.id')]
  jbp$span = gg$edges[jbp$edge.id]$span
  jbp$loose = FALSE

  if (!ignore.loose.ends)
  {
    ggl = gg$loose[, c()]
    ggl$loose = TRUE
    jbp = grbind(jbp, ggl)
  }

  jbpdt = as.data.table(jbp[order(gr.stripstrand(jbp))])
  jbpdt[, dist.to.next := c(start[-1]-end[-.N], NA), by = seqnames]

  ## mark those pairs that are within cushion distance of the next or previous and remove all groups of size != 2
  jbpdt[, in.cushion := dist.to.next<min.cushion, by = seqnames]
  jbpdt[, cushion.run := label.runs(in.cushion), by = seqnames]
  jbpdt[, last.cushion.run := c(NA, cushion.run[-.N]), by = seqnames] ## extend one node beyond
  jbpdt[is.na(cushion.run), cushion.run := last.cushion.run]
  jbpdt[, cushion.run.size := .N, by = .(seqnames, cushion.run)]

  ## identify bad cushions, and notify all junctions associated with bad cushions
  jbpdt[, bad := cushion.run.size>2 | any(span<min.span), by = .(seqnames, cushion.run)]

  if (any(jbpdt$bad)) ## propagate badness across junctions
    {
      badj = jbpdt[bad == TRUE, unique(edge.id)]
      jbpdt[, bad := bad | edge.id %in% badj]
    }

  ## keep only pairs or singletons
  jbpdt = jbpdt[cushion.run.size <= 2, ]

  if (nrow(jbpdt)==0)
    return(gg.empty)

  ## remaining bp pairs test to see if within max.insert of next and correct orientation
  jbpdt[, in.insert := dist.to.next<=min.insert, by = seqnames]
  jbpdt[, insert.run := label.runs(in.insert), by = seqnames] 
  jbpdt[, last.insert.run := c(NA, insert.run[-.N]), by = seqnames] ## extend one node beyond
  jbpdt[is.na(insert.run), insert.run := last.insert.run]
  jbpdt[!is.na(insert.run), insert.run.size := .N, by = .(seqnames, insert.run)]
  jbpdt[is.na(insert.run), insert.run.size := 0]

  jbpdt[, bad := bad | any(is.na(insert.run)), by = .(seqnames, cushion.run)] ## mark any cushion runs without an insert run

  ## mark any junction associated with a bad cushion run
  badj = unique(jbpdt[bad == TRUE, edge.id])

  ## keep paired breakpoints that are not associated with an uncushioned junctions
  ## though we keep track of "bad" junctions so they can propagate the badness
  jbpdt[, bad := bad | edge.id %in% badj]
  jbpdt = jbpdt[insert.run.size == 2, ]

  if (nrow(jbpdt)==0)
    return(gg.empty)

  ## now only pairs left
  ## any pairs with a loose end get NA sign
  jbpdt[, sign := ifelse(strand =='+', 1,-1)*ifelse(any(loose) | edge.id[1]==edge.id[2], NA, 1), by = .(seqnames, insert.run)]
  jbpdt[, ":="(pair.sign = sign[1]*sign[2],
               sign1 = sign[1],
               sign2 = sign[2])
      , by = .(seqnames, insert.run)]

  badj = unique(jbpdt[pair.sign>=0, edge.id])
#  jbpdt[, bad := bad | edge.id %in% badj] ## enough to disqualify a junction? maybe not 
  jbpdt = jbpdt[pair.sign<0, ]
  
  if (nrow(jbpdt)==0)
    return(gg.empty)

  ## since these are ordered by coordinate
  jbpdt[, edge.type := ifelse(sign1>0, 1, -1)]


  ## note here: absence of edge.type filter (i.e. different from tic)
  ## jbpdt = jbpdt[edge.type>0, ]

  ## also no need to compute actual paths on the original graph
  ## we only will keep track of paths on the junction graph 

  if (nrow(jbpdt)==0)
    return(gg.empty)

  ## now create a gGraph of a jgraph from the remaining breakpoints
  ## using jnodes that have coordinates edgeid:1-2
  ## and "sides" representing ends
  jnodes = unique(GRanges(jbpdt$edge.id, IRanges(1,2)))
  jnodes$edge.id = as.integer(as.character(seqnames(jnodes) ))
  jnodes$nodep = gg$edges[as.integer(as.character(seqnames(jnodes)))]$sdt$n1 ## left side of junction
  jnodes$rnodep = gg$edges[-as.integer(as.character(seqnames(jnodes)))]$sdt$n1 ## left side of flipped junction
  jnodes$bad = jnodes$edge.id %in% jbpdt[bad == TRUE, edge.id]

  jedges = jbpdt[
  , data.table(n1 = match(as.character(edge.id[1]), seqnames(jnodes)),
               n2 = match(as.character(edge.id[2]), seqnames(jnodes)),
               n1.side = sign(grl.iix[1] == 2),
               n2.side = sign(grl.iix[2] == 2),
               edge.type
               ),
  , by = .(seqnames, insert.run)][, edge.id := 1:.N]
  
  jgraph = gG(nodes = jnodes, edges = jedges)

  jgraph$clusters(mode = 'weak')

  ## weed out any bad clusters
  if (any(jgraph$nodes$dt$bad))
  {
    badcl = jgraph$nodes[bad == TRUE]$dt[, unique(c(cluster, rcluster))]
    jgraph = jgraph[!(cluster %in% badcl | rcluster %in% badcl), ]
  }

  ## graph may only be bad clusters
  if (length(jgraph)==0)
    return(gg.empty)

  ## map jgraph clusters back to gGraph
  wks = jgraph$walks()

  ## order walks by rev size
  wks = wks[rev(order(lengths(wks)))]

  ## for chromoplexy we need:
  ## 1) at least one edge type<0
  ## 2) at least two chromosomes
  ## 3) and if length = 2, then has to be a cycle
  ## check to make sure we have at least one edge.type<0 in each walk
  wks$set(num.neg.edge = wks$eval(edge = sum(edge.type<0)))
  wks$set(num.chrom = wks$eval(length(unique(seqnames(unlist(gg$edges$grl[edge.id]))))))
  wks = wks[num.neg.edge>0 & num.chrom>0 & (lengths(wks)>2 | wks$circular)]

  ## graph may only be bad clusters
  if (length(wks)==0)
    return(gg.empty)

  ## merge node and edge walk dts on the jgraph
  nedt = merge(wks$nodesdt, wks$edgesdt[!is.na(sedge.id), ], by = c('walk.id', 'walk.iid'), all = TRUE)

  ## note that unlike tic here we skip building strings of sedge.gnodes, since
  ## we won't be generating walks on the original graph ... 

  ## the "event" is the first junction, the nodes in the middle, and the last junction
  ## get the gedge corresponding to this jgraph snode
  nedt[, snode.gedges := sign(snode.id)*jgraph$nodes$gr$edge.id[abs(snode.id)]]

  wkstrings = nedt[, .(
    edgestr = paste(snode.gedges, collapse = ',')
  ), by = walk.id]

  ## only need to build edgedt 
  edgedt = dunlist(lapply(strsplit(wkstrings$edgestr, ','), as.integer))
  gg$edges[abs(edgedt$V1)]$mark(chromoplexy = edgedt$listid)
  gg$set(chromoplexy = wks$dts())

  if (mark)
  {
    gg$edges[!is.na(chromoplexy)]$mark(col = mark.col)
  }

  ## return walks with nodes 
  return(gg)
}




#' @name tic
#' @title
#'
#' Finds "clean" templated insertion chains eg paths and cycles of junctions with span > min.span
#' by first identifying jbp pairs within max.insert distance for which the "left" jbp is + and "right" jbp is -
#' (in reference coordinates) excluding any jbp that have more than one ALT jbp
#' (or optionally loose ends) with min.cushion distance.  These remaining jbp pairs
#' are then combined into a graph, and connected components in that graph are
#' scraped for paths and cycles, which are marked on the graph and added as metadata
#' to the outputted gGraph
#'
#' @param gg gGraph
#' @param max.insert max insert to consider in a templated insertion 1e5
#' @param min.span min span for a TIC junction (1e6)
#' @param min.cushion minimum cushion between a TIC junction and any other nearby event (to ensure "clean" events), the bigger the cushion, the cleaner the calls
#' @param ignore.loose.ends logical flag (FALSE) determining whether we ignore loose ends when filtering on min.cushion
#' @param ignore.small.dups logical flag (FALSE) determining whether we ignore small dups when filtering on min.cushion
#' @param ignore.small.dels logical flag (FALSE) determining whether we ignore small dels when filtering on min.cushion
#' @param max.small threshold for calling a local dup or del "small" 1e4
#' @return gGraph with $meta annotated with gWalks corresponding to tic and tip and nodes and edges labeled with 'p1' through 'pk' for all k templated insertion paths and 'c1' through 'ck' for all k templated insertion cycles
tic = function(gg, max.insert = 1e5,
               min.cushion = 1e6,
               min.span = 1e6,
               ignore.loose.ends = FALSE,
               ignore.small.dups = FALSE,
               ignore.small.dels = FALSE,
               max.small = 1e4,
               mark = FALSE,
               mark.col = 'purple'
               )
{
  gg = refresh(gg)

  ## set empty output - in case we find no events can quit early
  gg.empty = gg$copy
  gg.empty$nodes$mark(tic = as.integer(NA))
  gg.empty$edges$mark(tic = as.integer(NA))
  gg.empty$set(tic = data.table())

  ed = gg$edges[type == 'ALT']
  
  if (length(ed) && ignore.small.dups)
    ed = ed[!(class == 'DUP-like' & ed$span<=max.small)]

  if (length(ed) && ignore.small.dels)
    ed = ed[!(class == 'DEL-like' & ed$span<=max.small)]

  if (length(ed)==0)
    return(gg.empty)

  jbp = grl.unlist(ed$grl)[, c('edge.id', 'grl.iix', 'snode.id')]
  jbp$span = gg$edges[jbp$edge.id]$span
  jbp$loose = FALSE

  if (!ignore.loose.ends)
  {
    ggl = gg$loose[, c()]
    ggl$loose = TRUE
    jbp = grbind(jbp, ggl)
  }

  jbpdt = as.data.table(jbp[order(gr.stripstrand(jbp))])
  jbpdt[, dist.to.next := c(start[-1]-end[-.N], NA), by = seqnames]

  ## mark those pairs that are within cushion distance of the next or previous and remove all groups of size != 2
  jbpdt[, in.cushion := dist.to.next<min.cushion, by = seqnames]
  jbpdt[, cushion.run := label.runs(in.cushion), by = seqnames]
  jbpdt[, last.cushion.run := c(NA, cushion.run[-.N]), by = seqnames] ## extend one node beyond
  jbpdt[is.na(cushion.run), cushion.run := last.cushion.run]
  jbpdt[, cushion.run.size := .N, by = .(seqnames, cushion.run)]

  ## identify bad cushions, and notify all junctions associated with bad cushions
  jbpdt[, bad := cushion.run.size>2 | any(span<min.span), by = .(seqnames, cushion.run)]

  if (any(jbpdt$bad)) ## propagate badness across junctions
    {
      badj = jbpdt[bad == TRUE, unique(edge.id)]
      jbpdt[, bad := bad | edge.id %in% badj]
    }

  ## keep only pairs or singletons
  jbpdt = jbpdt[cushion.run.size <= 2, ]

  if (nrow(jbpdt)==0)
    return(gg.empty)

  ## remaining bp pairs test to see if within max.insert of next and correct orientation
  jbpdt[, in.insert := dist.to.next<=min.insert, by = seqnames]
  jbpdt[, insert.run := label.runs(in.insert), by = seqnames] 
  jbpdt[, last.insert.run := c(NA, insert.run[-.N]), by = seqnames] ## extend one node beyond
  jbpdt[is.na(insert.run), insert.run := last.insert.run]
  jbpdt[!is.na(insert.run), insert.run.size := .N, by = .(seqnames, insert.run)]
  jbpdt[is.na(insert.run), insert.run.size := 0]

  jbpdt[, bad := bad | any(is.na(insert.run)), by = .(seqnames, cushion.run)] ## mark any cushion runs without an insert run

  ## mark any junction associated with a bad cushion run
  badj = unique(jbpdt[bad == TRUE, edge.id])

  ## keep paired breakpoints that are not associated with an uncushioned junctions
  ## though we keep track of "bad" junctions so they can propagate the badness
  jbpdt[, bad := bad | edge.id %in% badj]
  jbpdt = jbpdt[insert.run.size == 2, ]

  if (nrow(jbpdt)==0)
    return(gg.empty)

  ## now only pairs left
  ## any pairs with a loose end get NA sign
  jbpdt[, sign := ifelse(strand =='+', 1,-1)*ifelse(any(loose) | edge.id[1]==edge.id[2], NA, 1), by = .(seqnames, insert.run)]
  jbpdt[, ":="(pair.sign = sign[1]*sign[2],
               sign1 = sign[1],
               sign2 = sign[2])
      , by = .(seqnames, insert.run)]

  badj = unique(jbpdt[pair.sign>=0, edge.id])
#  jbpdt[, bad := bad | edge.id %in% badj] ## enough to disqualify a junction? maybe not 
  jbpdt = jbpdt[pair.sign<0, ]
  
  if (nrow(jbpdt)==0)
    return(gg.empty)

  ## since these are ordered by coordinate
  ## the first bp in the pair has a strand
  jbpdt[, edge.type := ifelse(sign1>0, 1, -1)]
  badj = unique(jbpdt[edge.type<=0, edge.id])
#  jbpdt[, bad := bad | edge.id %in% badj] ## enough to disqualify a junction? maybe not
  jbpdt = jbpdt[edge.type>0, ]

  ## now let's make there is a path in the graph for every pair of surviving
  ## nodes

  ## orient edges correctly, so that the left bp (ie + and first in the pair) is n2 of the
  ## first edge and the second bp (ie second and - in the pair) is the $n1 of the second edge
  jbpdt[, seid := edge.id*
            ifelse(1:2 %in% 1,
            ifelse(grl.iix==2, 1, -1), ## flip if grl.iix == 2 so that $n2 of edge 1 enters the bridging node
            ifelse(grl.iix == 1, 1, -1)), by = .(insert.run, seqnames)] ## flip if grl.iix = 2 so that $n1 of edge 2 exits the bridging node
  jbpdt[, sid := gg$edges[seid]$sdt[, ifelse(1:2 %in% 1, n2, n1)],
        by = .(insert.run, seqnames)]
  jbpdt[, path := paste(unlist(gg$paths(sid[1], sid[2], ignore.strand = FALSE)$snode.id), collapse = ','), by = .(seqnames, insert.run)]

  badj = unique(jbpdt[nchar(path)==0, edge.id])
#  jbpdt[, bad := bad | edge.id %in% badj] ## enough to disqualify a junction? maybe not
  jbpdt = jbpdt[nchar(path)>0, ]

  if (nrow(jbpdt)==0)
    return(gg.empty)

  ## now create a gGraph of a jgraph from the remaining breakpoints
  ## using jnodes that have coordinates edgeid:1-2
  ## and "sides" representing ends
  jnodes = unique(GRanges(jbpdt$edge.id, IRanges(1,2)))
  jnodes$edge.id = as.integer(as.character(seqnames(jnodes) ))
  jnodes$nodep = gg$edges[as.integer(as.character(seqnames(jnodes)))]$sdt$n1 ## left side of junction
  jnodes$rnodep = gg$edges[-as.integer(as.character(seqnames(jnodes)))]$sdt$n1 ## left side of flipped junction
  jnodes$bad = jnodes$edge.id %in% jbpdt[bad == TRUE, edge.id]

  jedges = jbpdt[
  , data.table(n1 = match(as.character(edge.id[1]), seqnames(jnodes)),
               n2 = match(as.character(edge.id[2]), seqnames(jnodes)),
               n1.side = sign(grl.iix[1] == 2),
               n2.side = sign(grl.iix[2] == 2),
               path = path[1]
               ),
  , by = .(seqnames, insert.run)][, edge.id := 1:.N]

  jedges[, path :=  paste(as.character(strsplit(path, ',')[[1]]), collapse = ','), by = edge.id]
  jedges[, rpath := paste(-as.numeric(rev(strsplit(path, ',')[[1]])), collapse = ','), by = edge.id]
  
  jgraph = gG(nodes = jnodes, edges = jedges)

  jgraph$clusters(mode = 'weak')

  ## weed out any bad clusters
  if (any(jgraph$nodes$dt$bad))
  {
    badcl = jgraph$nodes[bad == TRUE]$dt[, unique(c(cluster, rcluster))]
    jgraph = jgraph[!(cluster %in% badcl | rcluster %in% badcl), ]
  }

  ## graph may only be bad clusters
  if (length(jgraph)==0)
    return(gg.empty)

  ## map jgraph clusters back to gGraph
  wks = jgraph$walks()

  ## order walks by rev size
  wks = wks[rev(order(lengths(wks)))]

  if (length(wks)==0)
    stop('something went wrong with jgraph walk calculation')

  ## merge node and edge walk dts on the jgraph
  nedt = merge(wks$nodesdt, wks$edgesdt[!is.na(sedge.id), ], by = c('walk.id', 'walk.iid'), all = TRUE)

  ## now use to build strings of nodes on the graph
  nedt[, sedge.gnodes := ifelse(sedge.id>0, jgraph$edges[abs(sedge.id)]$dt$path, jgraph$edges[abs(sedge.id)]$dt$rpath)]

  ## the "event" is the first junction, the nodes in the middle, and the last junction
  ## get the gedge corresponding to this jgraph snode
  nedt[, snode.gedges := sign(snode.id)*jgraph$nodes$gr$edge.id[abs(snode.id)]]

  wkstrings = nedt[, .(
    edgestr = paste(snode.gedges, collapse = ','),
    nodestr = paste(
#      ifelse(is.na(snode.gnodes), '', snode.gnodes),
      ifelse(is.na(sedge.gnodes), '', sedge.gnodes),
      sep = ' ', collapse = ',')
  ), by = walk.id]
  

  nodels = lapply(
    strsplit(gsub('^,', '',
                  gsub(',$', '',
                       gsub('\\,+', ',',
                            gsub('\\s+', ',', wkstrings$nodestr)))), ','), as.integer)

  edgedt = dunlist(lapply(strsplit(wkstrings$edgestr, ','), as.integer))
  nodedt = dunlist(nodels)

  gg$edges[abs(edgedt$V1)]$mark(tic = edgedt$listid)
  gg$nodes[abs(nodedt$V1)]$mark(tic = nodedt$listid)

  gwks = gW(nodels, graph = gg, circular = wks$dt$circular)
  gg$set(tic = gwks)

  if (mark)
  {
    gg$edges[!is.na(tic)]$mark(col = mark.col)
    gg$nodes[!is.na(tic)]$mark(col = mark.col)
  }

  ## return walks with nodes 
  return(gg)
}

               



#' @name chromothripsis
#' @title
#'
#' Finds chromothripsis as clusters of >= min.seg segments and >= min.jun
#' junctions, clusters defined as clusters of segs with  <= max.seg.width
#' spread across <= max.major footprints with cn <= max.cn and
#' cn amplitude cn <= max.cn.amplitude 
#' and also fractional width distribution consistent with dirichlet 
#' orientation distribution consistent with Uniform and
#' window distribution (if multiple windows) consistent with Uniform.
#'
#' "Major" footprints are defined as reduced footprints that
#' have width >= min.major.width.
#' 
#' We may also allow a number of (<= min.majo) "minor" footprints
#' (e.g. templated insertions) each <= min.major.width to prevent
#' clusters from being thrown out because certain TI junctions are mapped
#' in multiple parts.
#'
#' 
#' 
#' @param gg gGraph
#' @param width.thresh max width of a CT segment (in simplified gGraph) to consider (1e7)
#' @param min.seg minimum number of segments in a CT (8)
#' @param min.jun minimum number of juctnions in a CT (8)
#' @param max.cn max CN of a CT (4)
#' @param max.cn.ampltitude max difference between top and bottom CN in a CT (3)
#' @param max.major max number of "major" footprints of CT with width >= min.major.width
#' @param max.minor max number of "minor" footprints of CT with width < min.mmajor.width
#' @param min.major.width width threshold defining a major footprint
#' @param min.mean.stack  average number of "stacks" in event treating each junction span as a GRanges (2) 
#' @param mark logical flag whether to mark graph with CT events (FALSE)
#' @param mark.col logical flag of what to color events 
#' @param max.win max
#' @return gGraph with nodes and edges annotated with integer chromothripsis event or NA and metadata showing some statistics for the returns chromothripsis events 
chromothripsis = function(gg,
                          max.seg.width = 1e7,
                          min.seg = 8,
                          min.jun = 8,
                          max.cn = 4,
                          max.cn.amplitude = 3,
                          max.major = 3,
                          min.major.width = 1e5,
                          max.minor = 2,
                          min.width.p = 0.05,
                          min.orientation.p = 0.05,
                          min.mean.stack = 2, 
                          mark = FALSE,
                          mark.col = 'purple'
                          )
{
  ## hardcoded min CN amplitude
  min.cn.amplitude = 2

  ## make clusters
  ## note: all the code below will just be used to choose / filter among the clusters
  ## derived here
  oldclust = gg$nodes$dt$cluster
  oldrclust = gg$nodes$dt$rcluster
  gg.empty = refresh(gg)
  gg.empty$nodes$mark(chromothripsis = as.integer(NA))
  gg.empty$edges$mark(chromothripsis = as.integer(NA))
  gg = suppressWarnings(refresh(gg)$simplify()$clusters(width<max.seg.width))

  set.seed(10)
  ucl = sort(unique(gg$nodes$dt$cluster))
  if (length(ucl)==0)
    return(gg.empty)

  ## first round: get basic cluster stats to rule out obvious non CT based on CN
  ## and width distribution

  ## ks.test on beta distribution using quantile trimmed data
  .ksbeta = function(x, a=0.1, b = 0.9, alpha = 1) 
  {
    xtrim = x[x>=quantile(x, a) & x<=quantile(x, b)]
    if (length(xtrim)>5)
      suppressWarnings(ks.test(xtrim/sum(xtrim), pbeta, alpha*1, alpha*(length(xtrim)-1))$p)
    else
      1
  }

  cl.stats = gg$nodes$dt[!is.na(cluster), .(
                        nseg = .N,
                        cn.min = min(cn),
                        cn.max = max(cn),
                        median.wid = as.numeric(median(width)),
                        max.wid = as.numeric(quantile(width, 0.95)),
                        width.p = ## compare stick breaks to beta distribution
                          .ksbeta(width)
                      ), keyby = cluster]
  candidates = cl.stats[nseg>=min.seg & (cn.max-cn.min)<=max.cn.amplitude & cn.max <= max.cn & width.p>=0*min.width.p, ]

  if (nrow(candidates)==0)
      return(gg.empty)

  ## second round: edge stats, orientation distribution,
  cl.edges = gg$nodes[cluster %in% candidates$cluster]$edges[type == 'ALT']
  if (length(cl.edges)==0)
    return(gg.empty)

  grleft = cl.edges$junctions$left 
  grright = cl.edges$junctions$right
  
  ## orient junction ends from L to R 
  flip = gr.stripstrand(grleft)>=gr.stripstrand(grright)

  ## compiled edge stats
  ed.stats = data.table(
    cluster = pmin(gg$nodes[grleft$node.id]$dt$cluster, gg$nodes[grright$node.id]$dt$cluster, na.rm = TRUE),
    eid = cl.edges$dt$edge.id,
    orientation = ifelse(flip, paste(strand(grright), strand(grleft)),
                    paste(strand(grleft), strand(grright))))

  ed.stats[, njun := .N, by = cluster]

  orientations = data.table(prob = rep(0.25, 4), ori = as.data.table(expand.grid(c('-', '+'), c('-', '+')))[, paste(Var1, Var2)], key = 'ori')

  cluster.stats = ed.stats[,
                        .(
                          njun = .N,
                          orientation.p = {
                            data = merge(data.table(ori = orientation)[, .N, by = ori], orientations,by = 'ori', allow.cartesian = TRUE, all = TRUE)[is.na(N), N:= 0]
                            if (nrow(data)<=1 | sum(data$N)<=5) ## less than 5 junctions too little to call deviation
                              1
                            else
                              suppressWarnings(chisq.test(data$N, p = data$prob)$p.value)
                          }
                          ), by = cluster]

  cluster.stats = cluster.stats[njun>min.jun & orientation.p>=min.orientation.p, ]

  if (nrow(cluster.stats)==0)
      return(gg.empty)

  ## round 3:
  ## trim any surviving clusters by reducing width threshold
  ## maximizing "uniformity" across footprints
  ## while maintainin number of footprints and num.na below threshold
  ## num.na counts junctions that leave the cluster
  ## if any have too many footprints or not enough junction span overlap then
  ## remove

  ## max number of junctions we allow to "leave" the cluster e.g. to larger nodesxo
  max.na = 2

  cluster.stats[order(cluster), node.ids :=
                    {
                      cl.nodes = gg$nodes[gg$nodes$dt$cluster == cluster]$gr
                      cl.edges = gg$nodes[gg$nodes$dt$cluster == cluster]$edges[type == 'ALT']
                      grleft = cl.edges$junctions$left 
                      grright = cl.edges$junctions$right
                      wids = width(cl.nodes)
                      widths = rev(sort(wids[wids>=quantile(wids, 0.8)]))

                      ## p value for uniform distribution of junctions across
                      ## footprint pairs ... i.e. "loss" vs uniform model 
                      last.p = -1
                      last.na = length(setdiff(c(grleft$node.id, grright$node.id),
                                               cl.nodes$node.id))
                      p = 0
                      last.nodes = cl.nodes$node.id

                      ## check initial cluster footprints
                      ## if alreeady violating provided constraints then quit and return NA
                      cl.footprints = reduce(cl.nodes)
                                            
                      ## map junctions to cluster footprints
                      ## note some of these may be NA because some cluster junctions
                      ## may be connecting outside of the cluster
                      ## i.e. to larger nodes
                      cl.nodes$fpid = suppressWarnings(gr.match(cl.nodes, cl.footprints))
                      dt.nodes = as.data.table(cl.nodes)
                      setkey(dt.nodes, node.id)


                      grleft$fpid = dt.nodes[.(grleft$node.id), fpid]
                      grright$fpid = dt.nodes[.(grright$node.id), fpid]
                                            
                      span = ifelse(seqnames(grleft)==seqnames(grright),
                                    gr.string(
                                      GRanges(seqnames(grleft),
                                              IRanges(pmin(start(grleft), start(grright)),
                                                      pmax(start(grleft), start(grright))))), NA)

                      ## stack intrachromosomal junction spans
                      ## and for interchromosomal just use the cl.footprint span
                      stacks = grbind(
                        parse.gr(span[!is.na(span)]),
                        rep(cl.footprints, sum(is.na(span))))

                      stacksum = gr.sum(stacks) %&% cl.footprints                     
                      mean.stack = as.data.table(stacksum)[, sum(score*as.numeric(width))/sum(as.numeric(width))]
                                                       
                      num.major = sum(width(cl.footprints)>=min.major.width)
                      num.minor = sum(width(cl.footprints)<min.major.width)
                      
                      if (num.major>max.major | num.minor > max.minor | mean.stack<min.mean.stack)
                      {
                        ## will only get worse with decreasing footprints 
                        widths = c()
                        last.nodes = as.character(NA)
                      }

                      ## we go until 
                      while (length(widths)>0)
                      {
                        this.wid = widths[1]
                        these.nodes = cl.nodes %Q% (width<=this.wid)
                        cl.footprints = reduce(these.nodes)
                        
                        ## make table of expected fraction for all pairs of footprints inside cluster
                        ## using widths, will use this to determine whether distribution of
                        ## inter footprint junctions deviates from expectation
                        expected.frac = as.data.table(
                          expand.grid(
                            fp1 = 1:length(cl.footprints),
                            fp2 = 1:length(cl.footprints)))[fp1 <= fp2, ]
                        expected.frac[, area :=
                                          as.numeric(width(cl.footprints)[fp1])*
                                          as.numeric(width(cl.footprints)[fp2])]
                        expected.frac[, frac := area/sum(area)]
                        setkeyv(expected.frac, c('fp1', 'fp2'))
                      
                        ## map nodes back to footprints
                        cl.nodes$fpid = suppressWarnings(gr.match(cl.nodes, cl.footprints))
                        dt.nodes = as.data.table(cl.nodes)
                        setkey(dt.nodes, node.id)
                        
                        ## map junctions to cluster footprints
                        ## note some of these may be NA because some cluster junctions
                        ## may be connecting outside of the cluster
                        ## i.e. to larger nodes
                        grleft$fpid = dt.nodes[.(grleft$node.id), fpid]
                        grright$fpid = dt.nodes[.(grright$node.id), fpid]

                        efp1 = pmin(grright$fpid, grleft$fpid)
                        efp2 = pmax(grleft$fpid, grright$fpid)
                        ## compute goodness of fit of bp footprint
                        ## pair distribution vs expected / uniform
                        ## given footprint pair widths
                        observed.counts = data.table(
                          cluster,
                          fp1 = efp1,
                          fp2 = efp2)[, .N, keyby = .(cluster, fp1, fp2)][!is.na(fp1), ]
                        data = merge(observed.counts, expected.frac, allow.cartesian = TRUE, by =  c("fp1", "fp2"), all = TRUE)[is.na(N), N := 0]

                        last.p = p
                        if (nrow(data)<=1 | sum(data$N)<=5) ## not enough to rule out if <=5
                          p = 1
                        else
                          p = suppressWarnings(chisq.test(data$N, p = data$frac)$p.value)

                        ## now check whether we continue to improve p / loss vs model of
                        ## footprint uniformity, don't increase NA junctions,
                        ## and stay within our threshold of min vs max footprints
                        num.na = sum(is.na(efp1))
                        num.major = sum(width(cl.footprints)>=min.major.width)
                        num.minor = sum(width(cl.footprints)<min.major.width)
                        if (num.na<=last.na & p > last.p &
                           num.major<=max.major &
                           num.minor<=max.minor)                          
                        {
                          last.nodes = these.nodes$node.id
                          last.p = p
                          widths = widths[-1]
                        }
                        else ## finito
                          widths = c()
                      }
                      if (length(last.nodes)>1)
                        last.nodes = paste(last.nodes, collapse = ',')
                      last.nodes
                    }, by = cluster]

  cluster.stats = cluster.stats[!is.na(node.ids), ]

  gg$nodes$mark(chromothripsis = as.integer(NA))
  gg$edges$mark(chromothripsis = as.integer(NA))
  gg$set(chromothripsis = data.table())

  if (nrow(cluster.stats)>0)
  {
    cluster.stats$total.width = as.numeric(NA)
    for (i in 1:nrow(cluster.stats))
    {
      ## mark appropriate nodes and edges
      nids = as.numeric(unlist(strsplit(cluster.stats[i, node.ids], ',')))
      gg$nodes[nids]$mark(chromothripsis = i)
      gg$nodes[nids]$edges[type == 'ALT']$mark(chromothripsis = i)
      cluster.stats$total.width[i] = sum(as.numeric(width(reduce(gg$nodes[nids]$gr))))
    }
    if (mark)
    {
      gg$nodes[!is.na(chromothripsis)]$mark(col = alpha(mark.col, 0.3))
      gg$edges[!is.na(chromothripsis)]$mark(col = mark.col)
    }

    cluster.stats$node.ids = NULL
    ct = cbind(data.table(chromothripsis = 1:nrow(cluster.stats)),
               cluster.stats,
               cl.stats[.(cluster.stats$cluster), ]
               )

    gg$set(chromothripsis = ct)
  }

  ## restore old cluster (or NULL if not set)
  gg$nodes$mark(cluster = oldclust)
  gg$nodes$mark(rcluster = oldrclust)
  
  return(gg)
  
}


#' @name bfb
#' #' @export
#' #' @rdname internal
#' #' @description
#' #' Find the subgraph that is likely a BFB cycle event
#' #'
#' #' @param gg gGraph of the "alternate genome"
#' #'
#' #' @return gGraph object containing labeling the putative event
#' #' @export
bfb = function(gg){
    if (!inherits(gg, "gGraph")){
        stop("Input is not a gGraph object.")
    }
    
    if (length(gg$nodes)==0 |
        length(gg$edges)==0){
        return(gg)
    }
    if (!is.element("cn", colnames(gg$nodes$dt)) |
        !is.element("cn", colnames(gg$edges$dt)) |
        !any(gg$edges$dt[, type=="ALT"])){
        return(gg)
    }
    ## start from all FALSE, find one add one
    gg$annotate("bfb", data=0, id=gg$nodes$dt$node.id, class="node")
    gg$annotate("bfb", data=0, id=gg$edges$dt$edge.id, class="edge")
    ## label the fold-back junctions
    gg$annotate("fb",
                data=gg$junctions$span<=1e4 & gg$junctions$sign==1,
                id=gg$junctions$dt$edge.id,
                class="edge")
    if (!is.element("cn", colnames(gg$nodes$dt))){
        return(gg)
    }
    ## start with registering original node ix
    gg$nodes$mark(og.nid = seq_along(gg$nodes))
    gg$edges$mark(og.eid = seq_along(gg$edges))
    ## annotate the strongly connected components among amplicons
    amp.nid = gg$nodes$dt[cn>=5, node.id]
    if (length(amp.nid)==0){
        return(gg)
    }
    amp.cl = gg[amp.nid,] ## 5 is the baseline if a 2-round BFB event happened
    amp.cl$clusters(mode = "strong") ## FIXME: function does not handle empty grpah yet
    amp.cl.dt = copy(amp.cl$dt)
    cool.ix = amp.cl.dt[cluster==rcluster, unique(cluster)]
    if (length(cool.ix)==0){
        return(gg)
    }
    all.sg = list()
    out =
        lapply(seq_along(cool.ix),
               function(i){
                   x = cool.ix[[i]]
                   nix = which(amp.cl$nodes$dt$cluster==x)
                   this.sg = amp.cl[nix,]
                   all.sg[[as.character(i)]] <<- this.sg
                   palindromic.frac =
                       sum(width(this.sg$nodes[cluster==rcluster])) /
                       sum(width(this.sg$nodes))
                   this.juncs = this.sg$junctions
                   is.fb = this.sg$edges$dt[, fb]
                   n.fb = sum(is.fb, na.rm=T)
                   if (n.fb>0){
                       max.cn.fb = max(this.sg$edges$dt$cn[which(is.fb)], na.rm=T)
                   } else {
                       ## not zero because zero can ba an actual cn 0 fold back junction
                       max.cn.fb = -1
                   }
                   res = data.table(i = i,
                                    cix = x,
                                    palindromic.frac = palindromic.frac,
                                    n.fb = n.fb,
                                    max.cn.fb = max.cn.fb)
                   ## BFB criteria: at least two foldback juncs, max fb copy at least 2
                   ## palindromic fraction more than three quaters
                   res[, is.bfb := (n.fb>=2 &
                                    palindromic.frac>=0.75 &
                                    max.cn.fb>=2)]
                   return(res)
               })
    res = do.call(rbind, out)
    if (res[, any(is.bfb)]){
        ## annotate the nodes and edges
        for (i in res[is.bfb==TRUE, i]){
            this.sg = all.sg[[as.character(i)]]
            gg$annotate("bfb", i, this.sg$nodes$dt[, og.nid], "node")
            gg$annotate("bfb", i, this.sg$edges$dt[, og.eid], "edge")
        }
    }
    return(gg)
}


#' @name rigma
#' @description
#' "Ditches" of overlapping deletions
#' @param gg gGraph with $cn field annotated on nodes and edges
#' @param min.count minimum number of deletions to constitute a rigma (3)
#' @param min.frac minimum fraction of events in a cluster that must be deletions (0.9)
#' @param pad padding to allow between nearby deletions to constitute a cluster (1e5)
#' @param max.width max width of deletions to consider (1e6)
#' @param min.width min width of deletions to consider (1e3)
#' @export
rigma = function(gg,
                 min.count = 3,
                 min.frac = 0.9,
                 pad = 1e5,
                 mark = FALSE,
                 mark.col = 'purple', 
                 min.width = 1e3,
                 max.width = 1e7)
{
  if (!is(gg, 'gGraph'))
    stop('gg must be gGraph')

  if (length(gg)==0)
    return(gg)

  if (!is.element("cn", colnames(gg$nodes$dt)))
    {
      stop('nodes and edges must have $cn annotation for rigma function')
    }
  
  if (!any(gg$edges$dt[, type=="ALT"])){
    return(gg)
  }

  gg = refresh(gg)

  gg$nodes$mark(rigma = as.numeric(NA))
  gg$edges$mark(rigma = as.numeric(NA))


  ## lower "water level" from highest nonzero copy number to highest
  ## keep only nodes at or below that copy number that have not yet
  ## been marked as a rigma
  ucn = rev(sort(unique(gg$nodes$dt$cn)))[-1]

  ## keep track of junction breakpoints to score candidates
  juncs = gg$edges[type == 'ALT']$junctions
  grl = juncs$grl
  values(grl)$span = juncs$span
  jbp = grl.unlist(grl)[, c("edge.id", "class", "cn", "span")]
  
  ## define "rigma candidate" breakpoints
  jbp$rigmac = jbp$class == 'DEL-like' & jbp$span<=max.width & jbp$span>=min.width

  k = 1
  numevents = 0

  ## candidates as a function of "water level" and previous rigmas calls 
  .candidates = function(edges, cn.thresh, pad)
  {
    ## compute rigma candidates not yet labeled as rigma to identify remaining candidates regions
    if (length(edges)==0) ## return empty GRanges
      return(gg$nodes$gr[c()])

    ## compute td.spans
    gr1 = edges$junctions$left
    gr2 = edges$junctions$right
    td.dt = edges$dt
    td.dt$seqnames = as.character(seqnames(gr1))
    td.dt$start = pmin(start(gr1), start(gr2))
    td.dt$end = pmax(end(gr1), end(gr2))   
    td.dt[, cn := pmax(gg$nodes$dt$cn[n1], gg$nodes$dt$cn[n2])] ## min node cn associated with td
    td.dt = td.dt[cn<=cn.thresh, ]
    if (nrow(td.dt)==0) ## return empty GRanges
      return(gg$nodes$gr[c()])  
    
    td.spans = dt2gr(td.dt)
    
    ## reduce td.spans and keep only those reductions that intersect at least two td
    reduce(td.spans+pad)-pad
  }

  ## raise the "water level" until no more candidates
  while (length(candidates <- .candidates(
                  gg$edges[jbp$edge.id[jbp$rigmac]][is.na(rigma)],
                  ucn[k], pad)))
  {
    ## score candidate regions based on their junction mix
    ov = candidates %*% jbp
    values(candidates) = gr2dt(ov)[, .(np = sum(rigmac)/2, ntot = .N/2), keyby = query.id][.(1:length(candidates)), ]

    ix = which(candidates$np>=min.count &
               (candidates$np/candidates$ntot) >= min.frac) ## threshold for calling a rigma

    if (length(ix))
    {
      new.events = candidates[ix]
      new.events$pid = numevents + (1:length(new.events))
      numevents = numevents + length(new.events)
      ov = ov %Q% (rigmac & query.id %in% ix)
      ov$pid = new.events$pid[match(ov$query.id, ix)]
      gg$edges[ov$edge.id]$mark(rigma = ov$pid)
      grov = gg$nodes$gr[, 'node.id'] %*% new.events
      gg$nodes[grov$node.id]$mark(rigma = grov$pid)
    }
    k = k-1
  }

  urig = rev(sort(table(gg$edges$dt$rigma)))
  if (length(urig))
  {
    rnm = data.table(old = names(urig), new = as.numeric(1:length(urig)), key = 'old')
    
    ## rename so most popular (biggest) is first
    gg$edges$mark(rigma = rnm[gg$edges$dt$rigma, new])
    gg$nodes$mark(rigma = rnm[as.character(gg$nodes$dt$rigma), new])
    
    if (mark)
    {
      gg$nodes[!is.na(rigma)]$mark(col = alpha(mark.col, 0.3))
      gg$edges[!is.na(rigma)]$mark(col = mark.col)
    }
  }
  return(gg)
}


#' @name dm
#' @description
#' Discover the likely driver double minute contigs from the graph
#' they are defined as ultra high copy cycles with a mixture of ref and alt edges
#'
#' @param gg gGraph
#' @param min.amp the minimun amplitude to call a junction amplified
#' @return gg
#' @export
dm = function(gg,
              min.amp = 10){
    if (!is.element("cn", colnames(gg$nodes$dt)) |
        !is.element("cn", colnames(gg$edges$dt)) |
        !any(gg$edges$dt[, type=="ALT"])){
        return(gg)
    }

    ## mark the original node edge id
    gg$nodes$mark(og.nid = gg$nodes$dt$node.id)
    gg$edges$mark(og.eid = gg$edges$dt$edge.id)

    ## start from all FALSE, find one add one
    gg$annotate("dm", data=0, id=gg$nodes$dt$node.id, class="node")

    ## first get the high copy portion of the genome
    amp.gg = gg[, cn>=min.amp]
    if (length(amp.gg$nodes)==0 ||
        length(amp.gg$edges)==0 ||
        !amp.gg$edges$dt[, any(type=="ALT")]){
        return(gg)
    }

    ## find out the cycles from amp.gg
    amp.gg$clusters(mode = "strong")
    cool.cl = as.numeric(names(amp.gg$nodes$dt[cluster!=rcluster, which(table(cluster)>1)]))
    if (length(cool.cl)==0){
        return(gg)
    }
 
    tmp = amp.gg$nodes$dt[cluster %in% cool.cl, .(og.nid, cluster)]
    gg$annotate("dm", tmp$cluster, tmp$og.nid, "node")
    return(gg)
}


#' @name edgedist
#' @description
#' Computes distance matrix between junction endpoints, note that 0 = infinity,
#' and actual zero distance is represented as eps (e.g. 1e-9)
#'
#' @param edges gEdge object
#' @param mc.cores number of cores (1)
#' @param chunksize integer number of edges to consider when computing sparse distance matrix
#' @param ignore.strand whether to ignore.strand in distance computation
#' @param verbose logical flag (FALSE)
#' @param eps definition of epsilon in output matrix
#' @return sparse distance matrix where 0 represents infinity and actual 0 is eps
.edgedist = function(edges, mc.cores = 1, chunksize = 100, ignore.strand = FALSE, verbose = FALSE, eps = 1e-9)
{
  
  if (length(edges)==0){
    if (verbose){
      gmessage("No junction in this graph")                                     
    }
    return(NULL)
  }

  bp = grl.unlist(edges$grl)[, c("grl.ix", "grl.iix")]
  bp.dt = gr2dt(bp)

  adj = Matrix::sparseMatrix(1, 1, x = 0, dims = rep(length(bp), 2))
  ix = split(1:length(bp), ceiling(runif(length(bp))*ceiling(length(bp)/chunksize)))
  ixu = unlist(ix)

  if (verbose){
    message(sprintf('Computing junction graph across %s edges', length(edges)))
  }

  adj[ixu, ] = do.call(rbind,
                       mclapply(ix,
                                function(iix)
                                {
                                  if (verbose>1)
                                    cat('.')
                                  tmpm =
                                    gr.dist(bp[iix],
                                            gr.flipstrand(bp),
                                            ignore.strand = ignore.strand)
                                  return(as(tmpm, "Matrix"))
                                },
                                mc.cores = mc.cores))+eps
  adj[is.na(adj)] = 0
  ## get back to marcin's version
  colnames(adj) = rownames(adj) = paste(bp$grl.ix, bp$grl.iix, sep = '.')
  return(adj)
}



#' @name tornado
#' @title
#'
#' Identifies clusters of complex high copy amplicons with tens or hundreds
#' of high copy junctions.  Returns gGraph with nodes and edges annotated for
#' $tornado field and gGraph metadata $tornado data.table with stats for the
#' calls (if any)
#' 
#' @param gg gGraph
#' @param high.cn integer threshold for high copy number (20)
#' @param min.span minimum span for a contributing high copy rearrangement (1e6)
#' @param min.high min number of high copy junctions to qualify for a tornado (10)
#' @param min.low min number of low copy junctions to qualify for a tornado (10)
#' @param min.footprint min bp footprint of tornado to qualify (1e6)
#' @param pad pad to calculate number of footprints of tornado prior to reducing
#' @param mark logical flag whether to mark graph with CT events (FALSE)
#' @param mark.col logical flag of what to color events 
#' @param max.win max
#' @return gGraph with nodes and edges annotated with integer chromothripsis event or NA and metadata showing some statistics for the returns chromothripsis events
tornado = function(gg,
                   high.cn = 5,
                   min.span = 1e5,
                   min.high = 10,
                   min.low = 10,
                   min.footprint = 0.5e6,
                   min.nfootprint = 1,
                   mark = FALSE,
                   mark.col = 'purple',
                   mc.cores = 1, 
                   pad = 1e6)
{
  gg = refresh(gg)
  
  gg$nodes$mark(og.id = 1:length(gg$nodes))
  gg$nodes$mark(tornado = as.integer(NA))
  gg$set(tornado = data.table())
  gg.empty = gg$copy
  
  if (length(gg)==0)
    return(gg.empty)

  gg.high = gg[, cn>=high.cn]

  if (length(gg.high)==0)
    return(gg.empty)

  gg.high$clusters(mode = 'weak')

  ## don't test anything
  uclust = as.numeric(names(which(table(gg.high$nodes$dt$cluster)>=min.high)))

  if (length(uclust)==0)
    return(gg.empty)

  res = rbindlist(mclapply(uclust, function(clust)
  {
    seed = gg.high$nodes[cluster == clust]
    highcopy = seed$edges[type == 'ALT' & seed$edges$span>min.span]
    lowcopy = gg$nodes[seed$dt$og.id]$edges[type == 'ALT' & cn<high.cn]
    lowcopy = lowcopy[lowcopy$span>min.span]  
    ## check to see how many extend beyond the seed
    return(data.table(
      clust = clust,
      footprint = sum(width(seed$gr)),
      nfootprint = length(reduce(seed$gr+pad)),
      nchrom = length(unique(seqnames(seed$gr))),
      max.cn = max(seed$dt$cn),
      median.cn = median(seed$dt$cn),
      highcopy.edges = list(highcopy$dt$edge.id),
      lowcopy.edges = list(lowcopy$dt$edge.id),
      nhighcopy.edges = length(seed$edges[type == 'ALT']),
      nlowcopy.edges = length(lowcopy)))
  }, mc.cores = mc.cores))

  res = res[nhighcopy.edges>=min.high & nlowcopy.edges>=min.low & footprint>=min.footprint &
            nfootprint>=min.nfootprint, ]

  if (nrow(res)==0)
    return(gg.empty)

  res$tornado = 1:nrow(res)
  for (i in 1:nrow(res))
  {
    ## mark "seed" nodes in cluster
    ids = gg.high[cluster == res$clust[i]]$nodes$dt$og.id
    gg$nodes[ids]$mark(tornado = res$tornado[i])
    ## mark high and low copy edges
    eids = res[i, c(highcopy.edges[[1]], lowcopy.edges[[1]])]
    gg$edges[eids]$mark(tornado = res$tornado[i])
  }

  res$clust = NULL
  gg$set(tornado = res)

  if (mark)
  {
    gg$nodes[!is.na(tornado)]$mark(col = mark.col)
    gg$edges[!is.na(tornado)]$mark(col = mark.col)
  }

  return(gg)
}


#' @name pyrgos
#' @description
#'
#' Genomic "towers" formed from (mostly) DUP-like junctions
#' Find the subgraph with more than two overlaping, low cn, tDup junctions
#' and the middle part is highest copy number
#' @param gg gGraph
#' @param thresh Thresh max integer distance threshold to call pyrgos
#' @param cn.thresh maximum CN threshold for a junction to be called "pyrgos"
#' @param minj minimum number of junctions to qualify as a pyrgos (2)
#' @return gGraph with $pyrgo marking on nodes and edges labeling unique "events"
#' @export
pyrgos = function(gg,
                  min.count = 3,
                  mark = FALSE,
                  mark.col = 'purple',
                  min.frac = 0.9,
                  thresh = 1e7, cn.thresh = 2){
  if (!is(gg, 'gGraph'))
    stop('gg must be gGraph')

  if (length(gg)==0)
    return(gg)

  if (!is.element("cn", colnames(gg$nodes$dt)))
    {
      stop('nodes and edges must have $cn annotation for pyrgos function')
    }

  if (!any(gg$edges$dt[, type=="ALT"])){
    return(gg)
  }

  gg = refresh(gg)

  gg$nodes$mark(pyrgo = as.numeric(NA))
  gg$edges$mark(pyrgo = as.numeric(NA))

  ## raise the "water level" from lowest nonzero copy number to highest
  ## keep only nodes at or below that copy number that have not yet
  ## been marked as a pyrgo
  ucn = setdiff(sort(unique(gg$nodes$dt$cn)), 0)

  ## keep track of junction breakpoints to score candidates
  jbp = grl.unlist(gg$edges[type == 'ALT']$junctions$grl)[, c("edge.id", "class", "cn")]

  ## define "pyrgos candidate" breakpoints
  jbp$pyrgoc = jbp$class == 'DUP-like' & jbp$cn<=cn.thresh

  k = 1
  numevents = 0

  ## candidates as a function of water level and previous pyrgos calls 
  .candidates = function(edges, cn.thresh, thresh)
  {
    ## compute pyrgo candidates not yet labeled as pyrgo to identify remaining candidates regions    
    if (length(edges)==0) ## return empty GRanges
      return(gg$nodes$gr[c()])

    ## compute td.spans
    gr1 = edges$junctions$left
    gr2 = edges$junctions$right
    td.dt = edges$dt
    td.dt$seqnames = as.character(seqnames(gr1))
    td.dt$start = pmin(start(gr1), start(gr2))
    td.dt$end = pmax(end(gr1), end(gr2))   
    td.dt[, cn := pmin(gg$nodes$dt$cn[n1], gg$nodes$dt$cn[n2])] ## min node cn associated with td
    td.dt = td.dt[cn>=cn.thresh, ]
    if (nrow(td.dt)==0) ## return empty GRanges
      return(gg$nodes$gr[c()])  
    
    td.spans = dt2gr(td.dt)
    
    ## reduce td.spans and keep only those reductions that intersect at least two td
    reduce(td.spans) %&% (gr.sum(td.spans) %Q% (score>1)) %Q% (width<=thresh)
  }

  ## raise the "water level 
  while (length(candidates <- .candidates(gg$edges[jbp$edge.id[jbp$pyrgoc]][is.na(pyrgo)],
                                          ucn[k],
                                          thresh)))
  {
    ## score candidate regions based on their junction mix
    ov = candidates %*% jbp
    values(candidates) = gr2dt(ov)[, .(np = sum(pyrgoc)/2, ntot = .N/2), keyby = query.id][.(1:length(candidates)), ]
    
    ix = which(candidates$np>=min.count &
               (candidates$np/candidates$ntot) >= min.frac) ## threshold for calling a pyrgo
    
    if (length(ix))
    {
      new.events = candidates[ix]
      new.events$pid = numevents + (1:length(new.events))
      numevents = numevents + length(new.events)
      ov = ov %Q% (pyrgoc & query.id %in% ix)
      ov$pid = new.events$pid[match(ov$query.id, ix)]
      gg$edges[ov$edge.id]$mark(pyrgo = ov$pid)
      grov = gg$nodes$gr[, 'node.id'] %*% new.events
      gg$nodes[grov$node.id]$mark(pyrgo = grov$pid)
    }
    k = k+1
  }
  
  upyr = rev(sort(table(gg$edges$dt$pyrgo)))

  if (length(upyr))
    {
      rnm = data.table(old = names(upyr), new = as.numeric(1:length(upyr)), key = 'old')
      
      ## rename so most popular (biggest) is first
      gg$edges$mark(pyrgo = rnm[gg$edges$dt$pyrgo, new])
      gg$nodes$mark(pyrgo = rnm[as.character(gg$nodes$dt$pyrgo), new])
      
      if (mark)
      {
        gg$nodes[!is.na(pyrgo)]$mark(col = alpha(mark.col, 0.3))
        gg$edges[!is.na(pyrgo)]$mark(col = mark.col)
      }
    }

  return(gg)
}
