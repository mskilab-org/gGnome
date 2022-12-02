#' Proximity analysis of two genomic regions and rearrangment movement
#' @name proximity
<<<<<<< HEAD
#'
#' @description Takes a set of n "query" elements (GRanges object, e.g. genes) and determines 
#' their proximity to m "subject" elements (GRanges object, e.g. regulatory 
#' elements) subject to set of rearrangement adjacencies (GRangesList with width 1 range pairs)
#'
#' @param gg gGraph of the "alternate genome"
#' @param query GRanges of "intervals of interest" eg regulatory elements
#' @param subject GRanges of "intervals of interest" eg genes
#' @param strict.ref Boolean, only use actual reference edges in graph. 
#' Default: False
#' @param ignore.strand whether to ignore strand of input GRanges. Default: True
#' @param verbose logical flag, verbose output. Default: False
#' @param mc.cores how many cores to use for the path exploration step or if 
#' chunksize is provided, across chunks. Default: 1
#' @param chunksize chunks to split subject and query into to minimize memory 
#' usage, if mc.cores>1 then each chunk will be allotted a core. Default: NULL
#' @param max.dist maximum genomic distance to store and compute (1MB by default) 
#' should the maximum distance at which biological interactions may occur.
#' Default: 1e6
#' 
#' @details 
#' This analysis makes the (pretty liberal) assumption that all pairs of adjacencies 
#' that can be linked on a gGraph path are in cis (i.e. share a chromosome) in 
#' the tumor genome.
=======
#' @title proximity
#'
#' @description
#' Takes a set of n "query" elements (GRangse object, e.g. genes) and determines their proximity to m "subject" elements
#' (GRanges object, e.g. regulatory elements) subject to set of rearrangement adjacencies (GRangesList with width 1 range pairs)
#'
#' @details
#' This analysis makes the (pretty liberal) assumption that all pairs of adjacencies that can be linked on a gGraph path are in
#' cis (i.e. share a chromosome) in the tumor genome.
>>>>>>> upstream/master
#'
#' Each output proximity is a gWalk that connects query-subject on the genome
#' described by gGraph gg.  Each gWalk is  annotated by the metadata of the
#' corresponding query-subject GRanges pair as well as fields "altdist" and "refdist"
#' specifying the "alternate and "reference" gGraph distance of the
#' query-subject pair.  The gWalk metadata field "reldist" specifies the relative
#' distance (i.e. ratio of altdist to refdist) for that walk. 
#' 
#' For more details follow the Proximity Analysis in the gGnome Tutorial:
#' 
#' \href{http://mskilab.com/gGnome/tutorial.html#Proximity_analysis}{Proximity Analysis}
#' 
#' @return gWalk object each representing a proximity
#' @md
#' @export
proximity = function(gg,
                     query,
                     subject,
                     ignore.strand = TRUE,
                     verbose = F,
                     mc.cores = 1,
                     strict.ref = FALSE,
                     chunksize = NULL,
                     max.dist = 1e6)
{
    if (is.null(chunksize) && mc.cores >1){
      stop("chunksize must be specified if mc.cores > 1")
    }
    if (length(gg)==0)
        return(gW(graph = gg))
    ### check statement to see if query and subject are of the same grangelist seqnames
    if (!any(subject@seqnames@values %in% query@seqnames@values)){
      warning(" no matching seqnames between query and subject. 
              Check if annotation for both query and subject is correct.")
    }
    
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
            px = proximity(gg, query[qix], subject[six], ignore.strand = ignore.strand,
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
    ## dt.alt[, is.min := 1:.N %in% which.min(val), by = .(qid, sid)]
    dt.alt[, is.min := seq_len(.N) %in% which.min(val), by = .(qid, sid)]
  dt.alt = dt.alt[is.min == TRUE, ]
    ## dt.ref[, is.min := 1:.N %in% which.min(val), by = .(qid, sid)]
    dt.ref[, is.min := seq_len(.N) %in% which.min(val), by = .(qid, sid)]
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


#' find fusions
#' @name fusions
#' @description annotates all gene fusions in gGraph relative to cds definitions
#' 
#' cds = gencode cds GRanges gff3 / gtf or GRangesList the latter (converted via rtracklayer::import)
#' which has fields $exon_number
#'
#' "gene_name" GRangesList meta data field is used in annotation and in not 
#' creating "splice" fusionsthat arise from different transcripts of the same gene.
#'
#' @param graph input gGraph 
#' @param gencode  GFF containing gene boundaries and exons, in similar format to  
#' https://www.gencodegenes.org/ 
#' @param genes set of genes to pass for fusions.
#' @param mc.cores number of cores to run. Default: 1
#' @param annotate.graph Annotate the graph generated Default: True
#' @param verbose output verbose argument to function Default: False
#' 
#' @details 
#' For more info please follow the protein fusions analysis in the gGnome Tutorial:
#' 
#' \href{http://mskilab.com/gGnome/tutorial.html#Protein_fusions}{Protein fusions}
#' 
#' @return gWalks of gene fusions annotated with frame and gene(s)
#' @md
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
    gw = refresh(gw)
  }

  return(gw)
}

#' Make a transcript graph
#' @name make_txgraph
#'
#' @description  Generates "transcript graph" from gGraph and GENCODE GRanges
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
#' @param gg gGraph
#' @param verbose whether to provide verbose output
#' @return gWalk of paths representing derivative transcripts harboring "loops"
#' @author Marcin Imielinski
#' @noRd
#' @keywords internal
make_txgraph = function(gg, gencode)
  {
    tx = gencode %Q% (type == 'transcript')

    ## broken transcripts intersect at least one junction
    tx$in.break = tx %^% grbind(gg$loose, unlist(gg$edges[type == 'ALT']$junctions$grl))
    
    if (!any(tx$in.break)){
        warning("No breakpoint in any transcript.")
        return(NULL)
    }

    txb = tx[tx$in.break]

    ## we sort all cds associated with broken transcripts using their stranded coordinate value
    cds = gencode %Q% (type == 'CDS') %Q% (transcript_id %in% txb$transcript_id)

    ## remove any transcripts that lack a CDS (yes these exist)
    txb = txb %Q% (transcript_id %in% cds$transcript_id)

    ## now do a left merge of gg nodes with broken transcripts ..
    ## where we keep our original gg nodes but just duplicate
    ## them based on transcript intersections
    nov = gg$nodes$gr %*% txb
    txnodes = unname(gg$nodes$gr[nov$query.id])

    if (!length(txnodes))
      return(NULL)
   
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

#' get transcript paths from transcript graph
#' @name get_txpaths
#' 
#' @description
#' gets all transcript paths from a transcript graph (i.e. output of make_txgraph)
#' that connect CDS start and CDS end
#' 
#' 
#' @param tgg gGraph output of get_txgraph
#' @param mc.cores number of cores across which to parallelize path 
#' search component of algorithm
#' @param verbose whether to provide verbose output
#' @return gWalk of txPaths
#' @author Marcin Imielinski
#' @noRd
#' @keywords internal
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
    
#' get transcript loops
#' @name get_txloops
#' @description
#' gets internal "loops" that connect the right side of downstream nodes of
#' transcripts to left side of upstream nodes in the same transcript
#'
#' internal function used in fusions downstream of get_txpaths.
#' 
#' @param tgg gGraph output of get_txgraph
#' @param mc.cores number of cores across which to parallelize path search component of algorithm
#' @param verbose whether to provide verbose output
#' @return gWalk of paths representing derivative transcripts harboring "loops"
#' @author Marcin Imielinski
#' @noRd
#' @keywords internal
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
    if (!length(ref.p))
    {
      return(gW(graph = tgg))
    }

    right = unique(ref.p$nodes$eright[type == 'ALT']$left$dt$snode.id)
    left = unique(ref.p$nodes$eleft[type == 'ALT']$right$dt$snode.id)
    
    ## exclude 3p edges associated with the last node in any transcript
    dt.right = ref.p$nodesdt[snode.id %in% right, ][walk.iid != base::lengths(ref.p$snode.id)[walk.id], ]
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

        if (length(loops.complex))
          {
            
            ldt = cbind(loops.complex$dt[, .(snode.id.x = source, snode.id.y = sink)],
                        data.table(snode.id = loops.complex$snode.id, complex = TRUE))
            dtm = merge(dtm, ldt, all.x = TRUE, by = c("snode.id.x", "snode.id.y"), allow.cartesian = TRUE)
            dtm[!is.na(complex), ":="(loop = snode.id)]
            dtm[!is.na(complex), found := TRUE]
          }
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
        if (x==0) c() else list(y[1:x]), begini-1, ref.p$snode.id[walk.id], SIMPLIFY = FALSE)]
      loops[, suffix := mapply(function(x,y)
        if (x>length(y)) c() else list(y[x:length(y)]), endi+1, ref.p$snode.id[walk.id], SIMPLIFY = FALSE)]

#      loops[, snode.id := list(mapply("c", prefix, loop, suffix, SIMPLIFY = FALSE))]
      loops[, snode.id := list(mapply(function(x, y, z) c(unlist(x), unlist(y), unlist(z)), prefix, loop, suffix, SIMPLIFY = FALSE))]
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

#' annotate walks
#' @name annnotate_walks 
#' @description internal function to annotate the walks
#' @param walks walks input
#' @noRd
#' @keywords internal
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


#' Events calling
#' @name events
#'
#' @description Shortcut to call all simple and complex event types on JaBbA graph using 
#' standard settings on all event callers. 
#' 
#' @param gg gGraph
#' @param verbose verbose output Default: TRUE
#' @param mark Default: FALSE
#' @param QRP qrp events Default: FALSE
#' @return gGraph with nodes and edges annotated with complex events in their 
#' node and edge metadata and in the graph meta data field $events 
#' 
#' @return returns a gg with events found. A complete list of events include:
#' amplifcation events (tyfonas, dm, cpxdm, bfb), chromothripsis, deletions, duplications,
#' chromoplexies, tic's, and if QRP is True, qrp, qrpmix, qrppos, and qrpmin. 
#' 
#' For more information on event calling please follow the tutorial in gGnome:
#' 
#' \href{http://mskilab.com/gGnome/tutorial.html#Classifying_SV_events}{Classifying SV Events}
#' 
#' @seealso \code{\link{simple}}, \code{\link{amp}}, \code{\link{chromothripsis}},
#' \code{\link{del}}, \code{\link{dup}}, \code{\link{chromoplexy}}, 
#' \code{\link{tic}}, \code{\link{qrp}}
#' 
#' @md
#' @export
events = function(gg, verbose = TRUE, mark = FALSE, QRP = FALSE)
{
  gg = gg %>% simple(mark = TRUE)
  if (verbose)
    message('Finished simple')

  gg = gg %>% amp(mark = TRUE)
  if (verbose)
    message('Finished amp (tyfonas, dm, cpxdm, bfb)')
 
  gg = gg %>% chromothripsis(mark = TRUE)
  if (verbose)
      message('Finished chromothripsis')

  gg = gg %>% del(mark = TRUE)
  if (verbose)
    message('Finished del')

  gg = gg %>% dup(mark = TRUE)
  if (verbose)
    message('Finished dup')

  gg = gg %>% chromoplexy(mark = TRUE)
  if (verbose)
    message('Finished chromoplexy')

  gg = gg %>% tic(mark = TRUE)
  if (verbose)
    message('Finished tic')

  if (QRP){
    gg = gg %>% qrp(mark = TRUE)
    if (verbose)
      message('Finished qrp')
  }
  
    ev = rbind(
      gg$meta$simple,
      gg$meta$chromothripsis,
      gg$meta$chromoplexy,
      gg$meta$rigma,
      gg$meta$pyrgo,
      gg$meta$tic,
      gg$meta$amp,
      gg$meta$del,
      gg$meta$dup, fill = TRUE)[, ev.id := seq_len(.N)]
  if (QRP){
  ev = rbind(
      ev,
      gg$meta$qrp,
      ## gg$meta$qrpmix, #' ARRRGGGHHHH keh2019 Wednesday, Jan 26, 2022, Week 04, 07:53:16 PM
      gg$meta$qrppos,
      gg$meta$qrpmin,
      gg$meta$qrpmix, fill = TRUE)[, ev.id := seq_len(.N)]
  }

  gg$set(events = ev)
  return(gg)
}

#' Find Chromoplexy chains
#' @name chromoplexy
#'
#' @description Finds chromoplexy chains as clusters of "long distance" junctions 
#' that each span at least min.span (i.e. distant regions on the reference) have 
#' junctions nearby ie within max.dist. We filter to chains that have at least 
#' min.num footprints on the genome and involve at least min.num long distance 
#' junctions.  We keep track of how many "other" (non small dup and non small 
#' del) junctions there are in the vicinity for downstream filtering. 
#' 
#' @param gg gGraph
#' @param min.span minimimum span to define a "long distance" junction and also 
#' the span by which major footprints of the event must be separated 
#' Default: 1e7
#' @param max.dist maximum distance allowed in edge clusters
#' Default: 1e4
#' @param min.num minimum number of junctions and major footprints that define a 
#' chromoplexy
#' Default: 3
#' @param max.cn max copy number before filter
#' Default: 3
#' @param footprint.width padding around which to define the footprint of an 
#' event, note that the outputted footprint only includes the chromoplexy junction 
#' breakpoints.
#' Default: 1e6
#' @param ignore.small.dups logical flag determining whether we ignore 
#' small dups when filtering on min.cushion.
#' Default: True
#' @param ignore.small.dels logical flag determining whether we ignore small dels 
#' when filtering on min.cushion.
#' Default: True
#' @param max.small threshold for calling a local dup or del "small". 
#' Default: 5e4
#' @param mark mark chromoplexies. Default:FALSE
#' @param mark.col color to mark chromoplexies. Default: purple
#' 
#' @details For more details on what chromoplexies are and visualized:
#' \href{http://mskilab.com/gGnome/tutorial.html#Chromoplexy_and_TICs}{Chromoplexy and TICs} 
#' 
#' @return gGraph with $meta$chromoplexy annotated with chromoplexy event metadata
#' and edges labeled with $chromoplexy id or NA if the edge does not belong to a 
#' chromoplexy
#' @md 
#' @export
chromoplexy = function(gg,
                       min.span = 1e7,
                       max.dist = 1e4,
                       min.num = 3,
                       max.cn = 3,
                       footprint.width = 1e6,
                       ignore.small.dups = TRUE,
                       ignore.small.dels = TRUE,
                       max.small = 5e4,
                       mark = FALSE,
                       mark.col = 'purple'
               )
{
  gg = gGnome::refresh(gg)

  ## set empty output - in case we find no events can quit early
  gg.empty = gg$copy
  gg.empty$nodes$mark(chromoplexy = as.integer(NA))
  gg.empty$edges$mark(chromoplexy = as.integer(NA))
  gg.empty$set(chromoplexy = data.table())

  ggcopy = gg$copy
  ggcopy$edges$mark(og.id = ggcopy$edges$dt$edge.id)

  ed = ggcopy$edges[type == 'ALT']
  
  if (length(ed) && ignore.small.dups)
    ed = ed[!(class == 'DUP-like' & ed$span<=max.small)]
  
  if (length(ed) && ignore.small.dels)
    ed = ed[!(class == 'DEL-like' & ed$span<=max.small)]
     
  candidates = ed[ed$dt$cn<=max.cn]

  if (!length(candidates))
    return(gg.empty)

  gg.tmp = gG(si2gr(gg), junctions = candidates$junctions[, 'og.id'])$eclusters(thresh = max.dist)

  gg$edges$mark(ecluster = as.integer(NA))
  gg$edges[gg.tmp$edges$dt$og.id]$mark(ecluster = gg.tmp$edges$dt$ecluster)

  cl = gg$edges$dt[!is.na(ecluster), .N, by = ecluster][order(N), ][N>=min.num, ]

  cp = data.table()
  gg$edges$mark(chromoplexy = as.integer(NA))
  for (i in seq_along(cl$ecluster))
  {
    this.ed = gg$edges[gg$edges$span>=min.span & ecluster == cl$ecluster[i]]
    fp.minor = streduce(this.ed$footprint+footprint.width)
    fp.major = streduce(this.ed$footprint+min.span)
    if (length(fp.major)>=min.num)
    {
      num.other = ed[!(edge.id %in% this.ed$dt$edge.id)]$grl %&% fp.minor %>% length
      this.ed$mark(chromoplexy = i)
      this.ed$nodes$mark(chromoplexy = i)
      cp = rbind(
        cp,
        data.table(type = 'chromoplexy',
                   chromoplexy = i, 
                   njun = length(this.ed),
                   nfp.major = length(fp.major),
                   nfp.minor = length(fp.minor),
                   max.cn = max(this.ed$nodes$dt$cn),
                   min.cn = min(this.ed$nodes$dt$cn),
                   min.jcn = min(this.ed$dt$cn),
                   max.jcn = max(this.ed$dt$cn),
                   num.other = num.other,
                   footprint = paste(gr.string(sort(gr.stripstrand(this.ed$footprint))), collapse = ';'))[, frac.cp := njun/(njun+num.other)]
      )    
    }
  }

  gg$set(chromoplexy = cp)

  if (mark)
  {
    gg$edges[!is.na(chromoplexy)]$mark(col = mark.col)
  }

  ## return walks with nodes 
  return(gg)
}

## #' @name qrp
## #'
## #' Finds (quasi) reciprocal pairs of junctions.  Very related to chromoplexy or tic "cycles", but with
## #' exactly two junctions.
## #' 
## #' @param gg gGraph
## #' @param max.insert max insert to consider in a templated insertion (5e4)
## #' @param min.span min span for a TIC junction (1e6)
## #' @param min.cushion minimum cushion between a TIC junction and any other nearby event (to ensure "clean" events), the bigger the cushion, the cleaner the calls
## #' @param ignore.loose.ends logical flag (FALSE) determining whether we ignore loose ends when filtering on min.cushion
## #' @param ignore.small.dups logical flag (FALSE) determining whether we ignore small dups when filtering on min.cushion
## #' @param ignore.small.dels logical flag (FALSE) determining whether we ignore small dels when filtering on min.cushion
## #' @param max.small threshold for calling a local dup or del "small" 1e4
## #' @return gGraph with $meta annotated with gWalks corresponding to tic and tip and nodes and edges labeled with 'p1' through 'pk' for all k templated insertion paths and 'c1' through 'ck' for all k templated insertion cycles
## 
## qrp = function(gg, max.insert = 5e4,
##                    min.cushion = 1e6,
##                    min.span = 1e6,
##                    ignore.loose.ends = TRUE,
##                    ignore.small.dups = TRUE,
##                    ignore.small.dels = TRUE,
##                    max.small = 5e4,
##                mark = FALSE,
##                mark.col = 'purple'
##                )
## {
##   gg = gGnome::refresh(gg)

##   ## set empty output - in case we find no events can quit early
##   gg.empty = gg$copy
##   gg.empty$edges$mark(qrp = as.integer(NA))
##   gg.empty$nodes$mark(qrp = as.integer(NA))
##   gg.empty$set(qrp = data.table())

##   ed = gg$edges[type == 'ALT']
  
##   if (length(ed) && ignore.small.dups)
##     ed = ed[!(class == 'DUP-like' & ed$span<=max.small)]

##   if (length(ed) && ignore.small.dels)
##     ed = ed[!(class == 'DEL-like' & ed$span<=max.small)]

##   if (length(ed)==0)
##     return(gg.empty)

##   jbp = grl.unlist(ed$grl)[, c('edge.id', 'grl.iix', 'snode.id')]
##   jbp$span = gg$edges[jbp$edge.id]$span
##   jbp$loose = FALSE

##   if (!ignore.loose.ends)
##   {
##     ggl = gg$loose[, c()]
##     ggl$loose = TRUE
##     jbp = grbind(jbp, ggl)
##   }

##   jbpdt = as.data.table(jbp[order(gr.stripstrand(jbp))])
##   jbpdt[, dist.to.next := c(start[-1]-end[-.N], NA), by = seqnames]

##   ## mark those pairs that are within cushion distance of the next or previous and remove all groups of size != 2
##   jbpdt[, in.cushion := dist.to.next<min.cushion, by = seqnames]
##   jbpdt[, cushion.run := label.runs(in.cushion), by = seqnames]
##   jbpdt[, last.cushion.run := c(NA, cushion.run[-.N]), by = seqnames] ## extend one node beyond
##   jbpdt[is.na(cushion.run), cushion.run := last.cushion.run]
##   jbpdt[, cushion.run.size := .N, by = .(seqnames, cushion.run)]

##   ## identify bad cushions, and notify all junctions associated with bad cushions
##   jbpdt[, bad := cushion.run.size>2 | any(span<min.span), by = .(seqnames, cushion.run)]

##   if (any(jbpdt$bad)) ## propagate badness across junctions
##     {
##       badj = jbpdt[bad == TRUE, unique(edge.id)]
##       jbpdt[, bad := bad | edge.id %in% badj]
##     }

##   ## keep only pairs or singletons
##   jbpdt = jbpdt[cushion.run.size <= 2, ]

##   if (nrow(jbpdt)==0)
##     return(gg.empty)

##   ## remaining bp pairs test to see if within max.insert of next and correct orientation
##   jbpdt[, in.insert := dist.to.next<=max.insert, by = seqnames]
##   jbpdt[, insert.run := label.runs(in.insert), by = seqnames] 
##   jbpdt[, last.insert.run := c(NA, insert.run[-.N]), by = seqnames] ## extend one node beyond
##   jbpdt[is.na(insert.run), insert.run := last.insert.run]
##   jbpdt[!is.na(insert.run), insert.run.size := .N, by = .(seqnames, insert.run)]
##   jbpdt[is.na(insert.run), insert.run.size := 0]

##   jbpdt[, bad := bad | any(is.na(insert.run)), by = .(seqnames, cushion.run)] ## mark any cushion runs without an insert run

##   ## mark any junction associated with a bad cushion run
##   badj = unique(jbpdt[bad == TRUE, edge.id])

##   ## keep paired breakpoints that are not associated with an uncushioned junctions
##   ## though we keep track of "bad" junctions so they can propagate the badness
##   jbpdt[, bad := bad | edge.id %in% badj]
##   jbpdt = jbpdt[insert.run.size == 2, ]

##   if (nrow(jbpdt)==0)
##     return(gg.empty)

##   ## now only pairs left
##   ## any pairs with a loose end get NA sign
##   jbpdt[, sign := ifelse(strand =='+', 1,-1)*ifelse(any(loose) | edge.id[1]==edge.id[2], NA, 1), by = .(seqnames, insert.run)]
##   jbpdt[, ":="(pair.sign = sign[1]*sign[2],
##                sign1 = sign[1],
##                sign2 = sign[2])
##       , by = .(seqnames, insert.run)]

##   badj = unique(jbpdt[pair.sign>=0, edge.id])
##   ##  jbpdt[, bad := bad | edge.id %in% badj] ## enough to disqualify a junction? maybe not 
##   jbpdt = jbpdt[pair.sign<0, ]
  
##   if (nrow(jbpdt)==0)
##     return(gg.empty)

##   ## since these are ordered by coordinate
##   jbpdt[, edge.type := ifelse(sign1>0, 1, -1)]

##   ## note here: absence of edge.type filter (i.e. different from tic)
##   ## jbpdt = jbpdt[edge.type>0, ]

##   ## also no need to compute actual paths on the original graph
##   ## we only will keep track of paths on the junction graph 

##   if (nrow(jbpdt)==0)
##     return(gg.empty)

##   ## now create a gGraph of a jgraph from the remaining breakpoints
##   ## using jnodes that have coordinates edgeid:1-2
##   ## and "sides" representing ends
##   jnodes = unique(GRanges(jbpdt$edge.id, IRanges(1,2)))
##   jnodes$edge.id = as.integer(as.character(seqnames(jnodes) ))
##   jnodes$nodep = gg$edges[as.integer(as.character(seqnames(jnodes)))]$sdt$n1 ## left side of junction
##   jnodes$rnodep = gg$edges[-as.integer(as.character(seqnames(jnodes)))]$sdt$n1 ## left side of flipped junction
##   jnodes$bad = jnodes$edge.id %in% jbpdt[bad == TRUE, edge.id]

##   jedges = jbpdt[
##   , data.table(n1 = match(as.character(edge.id[1]), seqnames(jnodes)),
##                n2 = match(as.character(edge.id[2]), seqnames(jnodes)),
##                sn = seqnames,
##                n1.side = sign(grl.iix[1] == 2),
##                n2.side = sign(grl.iix[2] == 2),
##                edge.type
##                ),
##   , by = .(seqnames, insert.run)][, edge.id := 1:.N]
  
##   jgraph = gG(nodes = jnodes, edges = jedges)

##   jgraph$clusters(mode = 'weak')

##   ## weed out any bad clusters
##   if (any(jgraph$nodes$dt$bad))
##   {
##     badcl = jgraph$nodes[bad == TRUE]$dt[, unique(c(pcluster, ncluster))]
##     jgraph = jgraph[!(pcluster %in% badcl | cluster %in% badcl), ]
##   }

##   ## graph may only be bad clusters
##   if (length(jgraph)==0)
##     return(gg.empty)

##   ## map jgraph clusters back to gGraph
##   wks = jgraph$walks()

##   ## order walks by rev size
##   wks = wks[lengths(wks)==2 & wks$circular]

##   ## for reciprocal pair, unlike chromoplexy we don't need
##   ## 1) at least one edge type<0
##   ## 2) at least two chromosomes
##   ## 3) and if length = 2, then has to be a cycle
##   wks$set(num.neg.edge = wks$eval(edge = sum(edge.type<0)))
##   wks$set(num.chrom = wks$eval(edge = length(unique(sn))))
##   ## wks = wks[num.neg.edge>0 & num.chrom>0 & (lengths(wks)>2 | wks$circular)]

##   ## graph may only be bad clusters
##   if (length(wks)==0)
##     return(gg.empty)

##   ## merge node and edge walk dts on the jgraph
##   nedt = merge(wks$nodesdt, wks$edgesdt[!is.na(sedge.id), ], by = c('walk.id', 'walk.iid'), all = TRUE)

##   ## note that unlike tic here we skip building strings of sedge.gnodes, since
##   ## we won't be generating walks on the original graph ... 

##   ## the "event" is the first junction, the nodes in the middle, and the last junction
##   ## get the gedge corresponding to this jgraph snode
##   nedt[, snode.gedges := sign(snode.id)*jgraph$nodes$gr$edge.id[abs(snode.id)]]

##   wkstrings = nedt[, .(
##     edgestr = paste(snode.gedges, collapse = ',')
##   ), by = walk.id]

##   ## only need to build edgedt 
##   edgedt = dunlist(lapply(strsplit(wkstrings$edgestr, ','), as.integer))
##   gg$edges[abs(edgedt$V1)]$mark(qrp = edgedt$listid)

##   ## mark nodes with qrp
##   uid = unique(edgedt$listid)
##   for (i in uid)
##   {
##     gg$edges[which(qrp == i)]$nodes$mark(qrp = i)
##   }

##   summ = wks$dts()[, .(qrp = walk.id, length, circular, num.neg.edge, num.chrom)]
##   bps = sort(grl.unlist(gg$edges[!is.na(qrp)]$grl))
##   summ$footprint =  sapply(split(gr.string(bps), bps$qrp), paste, collapse = ';')[as.character(summ$qrp)]
##   summ$type = 'qrp'
##   gg$set(qrp = summ)

##   if (mark)
##   {
##     gg$edges[!is.na(chromoplexy)]$mark(col = mark.col)
##   }

##   ## return walks with nodes 
##   return(gg)
## }


#' Find templated insertion chains (tics)
#' @name tic
#' 
#' @description 
#' Finds "clean" templated insertion chains eg paths and cycles of junctions with span > min.span
#' by first identifying jbp pairs within max.insert distance for which the 
#' "left" jbp is + and "right" jbp is - (in reference coordinates) excluding any 
#' jbp that have more than one ALT jbp (or optionally loose ends) with min.cushion 
#' distance.  These remaining jbp pairs are then combined into a graph, and 
#' connected components in that graph are scraped for paths and cycles, which 
#' are marked on the graph and added as metadata to the outputted gGraph.
#'
#' @param gg gGraph
#' @param max.insert max insert to consider in a templated insertion 
#' Default: 5e5
#' @param min.span min span for a TIC junction 
#' Default: 1e6
#' @param min.cushion minimum cushion between a TIC junction and any other nearby 
#' event (to ensure "clean" events), the bigger the cushion, the cleaner the calls
#' Default: 5e5
#' @param ignore.loose.ends logical flag determining whether we ignore 
#' loose ends when filtering on min.cushion. Default: TRUE
#' @param ignore.small.dups logical flag determining whether we ignore 
#' small dups when filtering on min.cushion. Default: True
#' @param ignore.small.dels logical flag determining whether we ignore 
#' small dels when filtering on min.cushion. Default: True
#' @param max.small threshold for calling a local dup or del "small" 
#' Default: 5e4
#' @param mark Default: FALSE
#' @param mark.col Default: purple
#' 
#' @details For more detail on running tic and examples see the tutorial:
#' \href{http://mskilab.com/gGnome/tutorial.html#Chromoplexy_and_TICs}{tic}
#' 
#' @return gGraph with $meta annotated with gWalks corresponding to tic and tip 
#' and nodes and edges labeled with 'p1' through 'pk' for all k templated 
#' insertion paths and 'c1' through 'ck' for all k templated insertion cycles
#' @md
#' @export
tic = function(gg, max.insert = 5e4,
               min.cushion = 5e5,
               min.span = 1e6,
               min.length = 2,
               ignore.loose.ends = TRUE,
               ignore.small.dups = TRUE,
               ignore.small.dels = TRUE,
               max.small = 5e4,
               mark = FALSE,
               mark.col = 'purple'
               )
{
  gg = gGnome::refresh(gg)

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
  jbpdt[!is.na(cushion.run), cushion.run.size := .N, by = .(seqnames, cushion.run)]

  ## identify bad cushions, and notify all junctions associated with bad cushions
  jbpdt[, bad := FALSE]
  jbpdt[!is.na(cushion.run), bad := cushion.run.size>2 | any(span<min.span), by = .(seqnames, cushion.run)]

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
  jbpdt[, in.insert := dist.to.next<=max.insert, by = seqnames]
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

  if (nrow(jbpdt)==0)
    return(gg.empty)

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
    badcl = jgraph$nodes[bad == TRUE]$dt[, unique(c(pcluster, ncluster))]
    jgraph = jgraph[!(pcluster %in% badcl | ncluster %in% badcl), ]
  }

  ## graph may only be bad clusters
  if (length(jgraph)==0)
    return(gg.empty)

  ## map jgraph clusters back to gGraph
  wks = jgraph$walks()

  ## order walks by rev size
  wks = wks[rev(order(lengths(wks)))]

  ## filter walks on size adding +1 for length when not circular
  ## i.e. we allow templated insertion cycles but not paths with one junction
  ## pair if min.length = 2
  wks = wks[(sign(!wks$circular) + lengths(wks))>=min.length]
  
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
  meta = gwks$dt[, .(tic = walk.id, length, circular)]
  meta$footprint = sapply(1:length(gwks), function(x) grbind(gwks[x]$nodes$edges$footprint, gwks[x]$nodes$footprint) %>% streduce %>% gr.string %>% paste(collapse = ','))
  meta$num.chrom = gwks$eval(length(unique(seqnames)))
  meta$type = 'tic'
  gg$set(tic = meta)

  if (mark)
  {
    gg$edges[!is.na(tic)]$mark(col = mark.col)
    gg$nodes[!is.na(tic)]$mark(col = mark.col)
  }

  ## return walks with nodes 
  return(gg)
}
              
#' Find Chromothripsis
#' @name chromothripsis
#' 
#' @description Finds chromothripsis as clusters of >= min.seg segments and >= min.jun
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
#' @param min.seg minimum number of segments in a CT. 
#' Default: 8
#' @param min.jun minimum number of junctions in a CT. 
#' Default: 7
#' @param max.cn max CN of a CT. 
#' Default: 4
#' @param max.cn.ampltitude max difference between top and bottom CN in a CT. 
#' Default: 3
#' @param max.major max number of "major" footprints of CT with width >= min.major.width.
#' Default: 4
#' @param max.minor max number of "minor" footprints of CT with width < min.major.width.
#' Default: 2
#' @param min.major.width width threshold defining a major footprint.
#' Default: 1e5
#' @param min.mean.stack  average number of "stacks" in event treating each junction 
#' span as a GRanges.
#' Default: 2 
#' @param min.stack minimum stack threshold to make make clusters around ALT junction 
#' shadows.
#' Default: 3
#' @param min.p.orientation minimum pvalue orientation needed in cluster to keep
#' in cluster stats.
#' Default: 0.001
#' @param mark logical flag whether to mark graph with CT events 
#' Default: TRUE
#' @param mark.col logical flag of what to color events 
#' Default: purple
#' @param fbi.thresh Default: 5e4
#' @param remove.small.junctions removes small junctions found Default: TRUE
#' @param small.junction.thresh threshold for a junction to be considered small
#' Default: 1e4
#' @param scale.to.ploidy Scale to ploidy. if TRUE will double thresholds for 
#' amplitude and CN when ploidy is > 3.5
#' Default: True
#' @param mark logical flag to color eventsDefault: True
#' @param mark.col color of event Default: purple
#' 
#' @details for more details on chromothripsis:
#' \href{http://mskilab.com/gGnome/tutorial.html#Chromothripsis}{Chromothripsis}
#' 
#' @return gGraph with nodes and edges annotated with integer chromothripsis event 
#' or NA and metadata showing some statistics for the returns chromothripsis events
#' 
#' @md
#' @export
chromothripsis = function(gg,
              min.seg = 8,
              min.jun = 7,
              max.cn = 4,
              max.cn.amplitude = 3,
              max.major = 4,
              min.major.width = 1e5,
              max.minor = 2,
              min.p.orientation = 0.001,
              min.stack = 3,
              min.mean.stack = 3, 
              fbi.thresh = 5e4,
              remove.small.junctions = TRUE,
              small.junction.thresh = 1e4,
              scale.to.ploidy = TRUE, ## if TRUE will double thresholds for amplitude and CN when ploidy is > 3.5
              mark = TRUE,
              mark.col = 'purple'
              )
{
  ## initializations
  gg.empty = gGnome::refresh(gg)
  gg.empty$set(chromothripsis = data.table())
  gg.empty$nodes$mark(chromothripsis = as.integer(NA))
  gg.empty$edges$mark(chromothripsis = as.integer(NA))

  ## save current graph
  gg.og = gGnome::refresh(gg)

  ## copy and mark with original ids
  gg = gg.og$copy
  gg$nodes$mark(og.nodeid = gg$nodes$dt$node.id)
  gg$edges$mark(og.edgeid = gg$edges$dt$edge.id)

  ## simplify graph just in case it wasn't already
  gg = gg$simplify()

  ## scale thresholds to ploidy
  if (scale.to.ploidy)
  {
    ploidy = gg$nodes$dt[!is.na(cn), sum(cn*as.numeric(width))/sum(as.numeric(width))]
    if (ploidy>2)
      {
        max.cn = max.cn*ploidy/2
        max.cn.amplitude = max.cn.amplitude*ploidy/2
      }
  }

  gg$edges$mark(keep = TRUE)
  if (remove.small.junctions)
  {
    rem = which(gg$edges$dt$span<=small.junction.thresh &
                gg$edges$dt$class %in% c('DEL-like', 'DUP-like'))
    if (length(rem))
      gg$edges[rem]$mark(keep = FALSE)
  }

  if (!length(gg$edges[type == 'ALT' & keep == TRUE]))
    return(gg.empty)

  ## make clusters around ALT junction shadows that stack greater than min.stack
  cx.shadow = reduce((gr.sum(gg$edges[type == 'ALT' & keep == TRUE]$shadow) %Q% (score>=min.stack)))
  gg$nodes$mark(keep = (gg$nodes$gr %^% cx.shadow))
  gg = suppressWarnings(gg$clusters(keep == TRUE))

  set.seed(10)
  ucl = sort(unique(gg$nodes$dt$cluster))
  if (length(ucl)==0)
    return(gg.empty)

  ## all the remaining code is for filtering / annotating the above clusters

  ## first round: get basic cluster stats to rule out obvious non CT based on
  ## num segs and cn distribution
  maincol = c("seqnames", "start", "end")
  clgr = gg$nodes$gr %Q% (!is.na(cluster))
  cl.footprints = gr.reduce(clgr, by = 'cluster') %>% sort 
  cl.footprints$c.win = gr2dt(cl.footprints)[, c.win := 1:.N+0, by = cluster]$c.win
  nd = gr2dt(gr.val(clgr, cl.footprints, 'c.win'))

  nd = copy(gg$nodes$dt)[!is.na(cluster)][order(cluster), c("c.foot", "c.win") := {
    cluster.win = sort(gr.stripstrand(reduce(dt2gr(.SD))))
    mcols(cluster.win)$c.win = seq_along(cluster.win)
    gr2dt(gr.val(dt2gr(.SD[ , maincol, with = FALSE]), cluster.win, val = "c.win"))[, .(c.foot = paste0(seqnames, ":", start, "-", end), c.win)]
  }, by = cluster]

  ## mark nodes with their c.win in their respective cluster
  gg$nodes[match(nd$node.id, node.id)]$mark(c.win = nd$c.win)

  ## quantiles to use for computing "robust" cn.min and cn.max (and cn.amplitude)
  fp.stats = nd[, .(
    cn.min = gr.quantile(dt2gr(.SD), 0.01, 'cn') %>% as.integer,
    cn.max = gr.quantile(dt2gr(.SD), 0.99, 'cn') %>% as.integer,
    nseg = .N,
    median.wid = as.numeric(median(width)),
    max.wid = as.numeric(quantile(width, 0.95))
  ), keyby = .(cluster, c.win)]

  ## aggregate footprint stats to clusters
  cl.stats = fp.stats[, .(cn.amplitude  = max(cn.max-cn.min),    
                          nseg = sum(nseg),
                          cn.min = min(cn.min),
                          cn.max = max(cn.max)
                 ), keyby = cluster]

  candidates = cl.stats[nseg >= min.seg &
                        cn.amplitude <=
                        max.cn.amplitude & cn.max <= max.cn, ]

  if (nrow(candidates)==0)
      return(gg.empty)

  ## second round: edge stats, orientation distribution,
  cl.edges = gg$nodes[cluster %in% candidates$cluster]$edges[type == 'ALT' & keep == TRUE]

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
                             p.orientation = {
                               data = merge(data.table(ori = orientation)[, .N, by = ori], orientations,by = 'ori', allow.cartesian = TRUE, all = TRUE)[is.na(N), N:= 0]
                               if (nrow(data)<=1 | sum(data$N)<=5) ## less than 5 junctions too little to call deviation
                                 1
                               else
                                 suppressWarnings(chisq.test(data$N, p = data$prob, simulate = TRUE)$p.value)
                             }
                           ), by = cluster]

  cluster.stats = cluster.stats[p.orientation>=min.p.orientation  , ]

  if (nrow(cluster.stats)==0)
      return(gg.empty)

  ## third round:
  ## trim any surviving clusters by reducing width threshold
  ## maximizing "uniformity" across footprints
  ## while maintainin number of footprints and num.na below threshold
  ## num.na counts junctions that leave the cluster
  ## if any have too many footprints or not enough junction span overlap then
  ## remove

  ## p value via chisq based on widths distirbution of nodes and junctions joining them
  .pzigzag = function(cl.nodes, cl.edges)
  {
    ## declarations
    cl.footprints = reduce(cl.nodes+min.major.width)-min.major.width
    grleft = cl.edges$junctions$left 
    grright = cl.edges$junctions$right
    cl.nodes$fpid = suppressWarnings(gr.match(cl.nodes, cl.footprints))   
    dt.nodes = as.data.table(cl.nodes)
    setkey(dt.nodes, node.id)    
    grleft$fpid = dt.nodes[.(grleft$node.id), fpid] ## map junctions to cluster footprints
    grright$fpid = dt.nodes[.(grright$node.id), fpid]
    
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

    grleft = cl.edges$junctions$left 
    grright = cl.edges$junctions$right
    cl.nodes$fpid = suppressWarnings(gr.match(cl.nodes, cl.footprints))   
    dt.nodes = as.data.table(cl.nodes)
    setkey(dt.nodes, node.id)    
    grleft$fpid = dt.nodes[.(grleft$node.id), fpid] ## map junctions to cluster footprints
    grright$fpid = dt.nodes[.(grright$node.id), fpid]
    
    efp1 = pmin(grright$fpid, grleft$fpid)
    efp2 = pmax(grleft$fpid, grright$fpid)

    ## compute goodness of fit of bp footprint
    ## pair distribution vs expected / uniform
    ## given footprint pair widths
    observed.counts = data.table(
      fp1 = efp1,
      fp2 = efp2)[, .N, keyby = .(fp1, fp2)][!is.na(fp1), ]
    data = merge(observed.counts, expected.frac, allow.cartesian = TRUE, by =  c("fp1", "fp2"), all = TRUE)[is.na(N), N := 0][, N.expected := sum(N)*frac]

    nna = sum(is.na(efp1))
    if (nrow(data)<=1 | sum(data$N)<=5) ## not enough to rule out if <=5
      1
    else
      suppressWarnings(chisq.test(data$N, p = data$frac, simulate = TRUE)$p.value)
  }

  ## third round: filter on major, minor, njun, and annotate p.zigzag and mean.stack

  ## annotate surviing clusters with
  ## pzigzag, meanstack, num major and num minor 
  ## also determine whether mean.stack is sufficient
  cluster.stats = merge(cluster.stats,
                       cluster.stats[,
                    {
                      cl.nodes = gg$nodes[gg$nodes$dt$cluster == cluster]$gr
                      cl.edges = gg$nodes[gg$nodes$dt$cluster == cluster]$edges[type == 'ALT' & keep == TRUE]
                      is.internal = cl.edges$dt[, n1 %in% cl.nodes$node.id & n2 %in% cl.nodes$node.id]

                      is.fbi = cl.edges$class == 'INV-like' & cl.edges$span < fbi.thresh
                      ## check initial cluster footprints
                      ## if already violating provided constraints then quit and return NA
                      cl.footprints = reduce(cl.nodes+min.major.width)-min.major.width

                      stacksum = gr.sum(cl.edges$shadow %Q% (width>10)) %*% cl.nodes
                      mean.stack = gr2dt(stacksum)[, sum(as.numeric(width)*score)/sum(width)]
                      min.stack = gr.quantile(stacksum, 0.05, 'score') %>% as.integer
                      max.stack = gr.quantile(stacksum, 0.95, 'score') %>% as.integer

                      num.major = sum(width(cl.footprints)>=min.major.width)
                      num.minor = sum(width(cl.footprints)<min.major.width)

                      data.table(node.ids = paste(cl.nodes$node.id, collapse = ','),
                                 njun = sum(is.internal),
                                 njun.boundary = sum(!is.internal),
                                 njun.tot = length(is.internal),
                                 njun.fbi = sum(is.fbi),
                                 njun.frac = sum(is.internal)/length(is.internal),
                                 mean.stack = mean.stack,
                                 min.stack = min.stack,
                                 max.stack = max.stack,
                                 num.major = sum(width(cl.footprints)>=min.major.width),
                                 num.minor = sum(width(cl.footprints)<min.major.width),
                                 p.zigzag = .pzigzag(cl.nodes, cl.edges))
                    }, by = cluster], by = 'cluster')

  ## final filter
  cluster.stats = cluster.stats[njun >= min.jun &
                                num.major <= max.major &
                                mean.stack >= min.mean.stack & 
                                num.minor <= max.minor, ]

  ### wrap up 
  ## do final marking of nodes and event summarization, annotate gg object
  gg.og$nodes$mark(chromothripsis = as.integer(NA))
  gg.og$edges$mark(chromothripsis = as.integer(NA))
  gg.og$set(chromothripsis = data.table())

  ## label only the nodes that survived (ie were stored in cluster.stats$node.ids)
  if (nrow(cluster.stats)>0)
  {
    ## fill in stack and stats
    cluster.stats$total.width = as.numeric(NA)
    cluster.stats$footprint = as.character(NA)

    for (i in 1:nrow(cluster.stats))
    {
      ## mark appropriate nodes and edges
      nids = as.numeric(unlist(strsplit(cluster.stats[i, node.ids], ',')))
      footprint = gr.stripstrand(gg$nodes[nids]$footprint)
      cluster.stats$footprint[i] = paste(gr.string(footprint), collapse = ';')
      ## now reach back to mark original unsimplified graph
      nodes.og = gg.og$nodes %&% footprint
      nodes.og$mark(chromothripsis = i)
      nodes.og$edges[type == 'ALT']$mark(chromothripsis = i)
      cluster.stats$total.width[i] = sum(as.numeric(width(reduce(gg$nodes[nids]$gr))))
    }

    if (mark)
    {
      gg.og$nodes[!is.na(chromothripsis)]$mark(col = alpha(mark.col, 0.3))
      gg.og$edges[!is.na(chromothripsis)]$mark(col = mark.col)
    }

    cluster.stats$node.ids = NULL
    ct = cbind(data.table(chromothripsis = 1:nrow(cluster.stats)),
               cluster.stats,
               cl.stats[.(cluster.stats$cluster), ]
               )

    ct$type = 'chromothripsis'
    gg.og$set(chromothripsis = ct)
  }

  return(gg.og)  
}


#' Simple Event Calling
#' @name simple
#' @export
#' @description Call simple events in gGraph
#'
#' @param gg gGraph, must have 'cn' node and edge annotation
#' @param reciprocal.thresh reciprocal threshold added to alt.shadows to find
#' translocations
#' Default: 1e4
#' @param tra.pad translocation padding to add on when finding locations 
#' Default: 1e6
#' @param mark logical flag specifying whether to mark fold backs 
#' Default: TRUE
#' @param mark.col character specifying colors to mark graph with 
#' Default: 'purple'
#' 
#' @details simple event calling includes 3 events: inversion, inverted 
#' duplication, and translocation.
#' 
#' For how to run this function:
#' \href{http://mskilab.com/gGnome/tutorial.html#Simple}{Simple}
#' 
#' @return gGraph object containing labeling the putative event
#' 
#' 
#' @md
#' @export
simple = function(gg,
                  reciprocal.thresh = 1e4,
                  tra.pad = 1e6,
                  mark = TRUE,
                  mark.col = 'purple'
               )
{
  if (!is(gg, 'gGraph'))
    stop('gg must be gGraph')
  
  gg = gGnome::refresh(gg)
  gg$nodes$mark(simple = NULL) ## 
  gg$edges$mark(simple = NULL) ##
  gg.empty = gg$copy
  gg.empty$set(simple = data.table())
  gg.empty$nodes$mark(simple = NA_character_) ##
  gg.empty$edges$mark(simple = NA_character_) ##
  gg.empty$nodes$mark(simple = as.integer(NA))
  gg.empty$edges$mark(simple = as.integer(NA))
  
  if (length(gg)==0)
    return(gg)
  
  if (!is.element("cn", colnames(gg$nodes$dt)))
  {
    stop('nodes and edges must have $cn annotation for bfb function')
  }
  
  if (!any(gg$edges$dt[, type=="ALT"])){
    return(gg.empty)
  }

  gg$nodes$mark(simple = NA_character_) ##
  gg$edges$mark(simple = NA_character_) ##

  alt = gg$edges[type == 'ALT']
  alt.shadows = alt$shadow
  alt.shadows$class = alt$dt$class[alt.shadows$id]
  alt.shadows$jcn = alt$dt$cn[alt.shadows$id]
  alt.shadows$ncn.left = alt[alt.shadows$id]$left$dt$cn
  alt.shadows$ncn.right = alt[alt.shadows$id]$right$dt$cn
  sum.shadows = gr.sum(alt.shadows) %Q% (score>0)
  simple.inv = GRanges(seqinfo = seqinfo(gg))
  inv.shadows = alt.shadows %Q%
      (!is.na(ncn.left) & !is.na(ncn.right)) %Q%
      (class == 'INV-like') %Q% (ncn.left == ncn.right)
  inv.shadows$str2 = inv.shadows$str = strand(alt[inv.shadows$id]$junctions$left)
  inv.shadows$id2 = inv.shadows$id
  if (length(inv.shadows))
    {
      ov = inv.shadows[, c('id', 'str')] %*% inv.shadows[, c('id2', 'str2')] %Q%
        (id != id2 & str != str2 &
         start(inv.shadows)[query.id]<=start(inv.shadows)[subject.id] & 
         inv.shadows$jcn[query.id]==inv.shadows$jcn[subject.id] &
         inv.shadows$ncn.left[query.id]==inv.shadows$ncn.left[subject.id] &
     inv.shadows$ncn.right[query.id]==inv.shadows$ncn.right[subject.id]
        )
      
      simple.inv = suppressWarnings(
        dt2gr(gr2dt(gr2dt(ov)[, .(
                       seqnames = as.character(seqnames(inv.shadows)[query.id]),
                       start = pmin(start(inv.shadows)[query.id], start(inv.shadows)[subject.id]),
                       end = pmax(end(inv.shadows)[query.id], end(inv.shadows)[subject.id]),
                       inv.width = width,
                       type = ifelse(str == '+', 'INVDUP', 'INV'),
                       id1 = inv.shadows$id[query.id],
                       id2 = inv.shadows$id[subject.id],
                       jcn = inv.shadows$jcn[query.id],
                       ncn.left = inv.shadows$ncn.left[query.id],
                       ncn.right = inv.shadows$ncn.right[query.id]
                       )]), seqlengths = seqlengths(gg)))
      
      if (length(simple.inv))
      {
        simple.inv$inv.frac = ifelse(simple.inv$type == 'INV', simple.inv$inv.width/width(simple.inv), 1-simple.inv$inv.width/width(simple.inv))
        simple.inv$score = gr2dt(simple.inv %*% alt.shadows)[, .(ov = length(setdiff(id, c(id1, id2)))), keyby = 'query.id'][.(1:length(simple.inv)),ov ]
        simple.inv = simple.inv %Q% (score==0)
      }
    }

  simple = grbind(simple.inv)
  if (!is.null(simple))
    {
      simple$footprint = gr.string(simple)
      simple = as.data.table(values(simple))
    }
  else
    simple = data.table()

  ## find pure translocations
  ## i.e. reciprocal or unbalanced between two chromosomes that don't have other translocations
  ## and loci that are not in the shadow of another event 
  tra.shadows = alt.shadows %Q% (class == 'TRA-like')
  tra.sum = gr.sum(tra.shadows+reciprocal.thresh) %Q% (score>0)
  tra.keep = gr2dt(tra.sum)[, .(count = sum(score!=0)), keyby = seqnames][count == 1, ]$seqnames
  tra.candidates = tra.sum %Q% (score <= 2 & seqnames %in% tra.keep) %*% tra.shadows 
  tra.candidates = tra.candidates[!(tra.candidates %^% ((alt.shadows + tra.pad) %Q% (class != 'TRA-like')))]

  if (length(tra.candidates))
    {
      tra.id = gr2dt(tra.candidates)[, .N, keyby = id][N==2, id]
      tra.pairs = gr2dt(tra.candidates %Q% (id %in% tra.id))[, .(id1 = id[1], id2 = id[.N]), by = subject.id]
      tra.cl = clusters(graph.edgelist(cbind(as.character(tra.pairs$id1), as.character(tra.pairs$id2))))$membership
      simple.tra = suppressWarnings(merge(data.table(id = as.numeric(names(tra.cl)), simple = unname(tra.cl)),
                                          gr2dt(tra.candidates), by = 'id')[, .(
                                           jcn.min = max(jcn),
                                           type = 'TRA',
                                           ncn.left = max(ncn.left),
                                           ncn.right = max(ncn.right),
                                           footprint = paste0(seqnames, ':', start, '-', end, collapse = ',')
                                         ), by = simple])
      simple = rbind(simple, simple.tra, fill = TRUE)
    }

  if (nrow(simple))
  {
    ## concatenate simple translocations with rest
    simple = simple[, simple := paste0(type, 1:.N)]
    simple = simple[, intersect(names(simple), c("simple", "type", "jcn", "ncn.left", "ncn.right", "footprint")), with = FALSE]    
    simple$type = tolower(simple$type)
    gru = grl.unlist(parse.grl(simple$footprint, seqlengths = seqlengths(gg)))
    grun = gru %*% (gg$nodes$gr[, 'node.id']) %Q% (width>1)
    gg$nodes[grun$node.id]$mark(simple = simple$simple[grun$grl.ix])
    grue = gru %*% grl.unlist(gg$edges[type == 'ALT']$grl)[, 'edge.id']
    gg$edges[grue$edge.id]$mark(simple = simple$simple[grue$grl.ix])      
  }
  
  if (mark)
  {
    gg$nodes[!is.na(gg$nodes$dt$simple)]$mark(col = alpha(mark.col, 0.3))
    gg$edges[!is.na(gg$edges$dt$simple)]$mark(col = mark.col)
  }
  gg$set(simple = simple)

  return(gg)
}


#' Find simple deletions and rigmas
#' @name del
#' @description Calls simple deletions (del) and rigma, which are "rifts" or 
#' clusters of overlapping deletions
#' 
#' 
#' @param gg gGraph with $cn field annotated on nodes and edges
#' @param fdr.thresh False Discovery Rate threshold. 
#' Default: 0.5
#' @param tile.width bin width to use when computing deletion clustering.
#' Default: 1e6
#' @param cn.thresh node copy number threshold to call a deletion in ploidy units. 
#' Default: 2
#' @param jcn.thresh edge copy number threshold to call a deletion in ploidy units. 
#' Default: 1
#' @param min.count minimum number of overlapping deletions to constitute a rigma, i
#' ncludes the tile.width. 
#' Default: 2
#' @param max.width max width of deletions to consider for rigma. 
#' Default: 1e7
#' @param min.width min width of deletions to consider for rigma. 
#' Default: 1e4
#' @param max.width.flank max width flank each side of class dup-like to consider. 
#' Default: 1e4
#' @param mark color deletion events. Default: False
#' @param mark.col color of event. Default: purple
#' @param return.fish Parameter to return fishook::Fish output. Default: False
#' 
#' @details Deletions are defined as having low junction copy number and connect two nodes of low junction copy number
#' (with cn and jcn thresholds provided as parameters).  Simple deletions have no other overlapping junctions
#' in their shadow.  Rigma have a min.count and are also outliers (<fdr.thresh) in a negative binomial model
#' that incorporates the regional (non-DEL) junction count in tile.width genomic bins (set to 1 Mbp by default). 
#'
#' Note: Not all DEL-like junctions will be called a del or rigma.  
#' 
#' More details on how to run this function and examples:
#' \href{http://mskilab.com/gGnome/tutorial.html#Deletions_and_Rigma}{Deletions & Rigma}
#' 
#' @return gGraph with nodes and edges annotated with $del and $rigma metadata 
#' field, and data.tables $meta$rigma and $meta$set with event level statistics.
#' 
#' If return.fish = TRUE, returns FishHook::Fish() output.
#' @md
#' @export
del = function(gg,
               fdr.thresh = 0.5,
               tile.width = 1e6, ## pad around which we don't want to see any other junctions
               cn.thresh = 2,
               jcn.thresh = 1,
               min.count = 2, 
               mark = FALSE,
               return.fish = FALSE,
               mark.col = 'purple', 
               min.width = 1e4,
               max.width.flank = 1e4,
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

  gg = gGnome::refresh(gg) 
  gg$nodes$mark(rigma = as.integer(NA))
  gg$edges$mark(rigma = as.integer(NA))
  gg$nodes$mark(del = as.integer(NA))
  gg$edges$mark(del = as.integer(NA))
  gg$set(del = data.table())
  ploidy = ceiling(gg$nodes$dt[!is.na(cn), sum(cn*as.numeric(width))/sum(as.numeric(width))])

  all = gg$edges[type == 'ALT']
  all$mark(cn.max = pmax(all$left$dt$cn, all$right$dt$cn))

  ## define all dels 
  all.dels = all[class == 'DEL-like' &
                 all$dt$cn.max <= cn.thresh*ploidy &
                 all$dt$cn<=jcn.thresh*ploidy]

  ## candidate rigma dels be within certain width threshod
  dels = all.dels[ all.dels$span>=min.width  & 
                   all.dels$span<=max.width]

  if (length(dels)==0)
    return(gg)
  
  ## "other" are non-small dups and inversions / translocations 
  other = all[(class %in% c('DUP-like') & all$span>max.width.flank) | class %in% c('INV-like', 'TRA-like')]

  ## only consider dels that are not int he shadow of an other
  shadows = dels$shadow
  shadows$edge.id = dels$dt$edge.id
  shadows = shadows[!(shadows %^% other$shadow)]

  if (length(shadows)==0)
    return(gg)

  foci = gr.sum(shadows) %Q% (score>= min.count) 
  candidates = reduce((shadows+tile.width) %&% foci)-tile.width

  ## collect some stats simply for record keeping 
  candidates$depth = gr.val(candidates, foci, val = 'score', FUN = max, weighted = FALSE)$score
  candidates$min.cn = gr.val(candidates, gg$nodes$gr, val = 'cn', FUN = min, weighted = FALSE)$cn
  candidates$max.cn = gr.val(candidates, gg$nodes$gr, val = 'cn', FUN = max, weighted = FALSE)$cn  
  candidates$del.count = (candidates %N% unlist(dels$grl))/2
  candidates$other.count = (candidates %N% unlist(other$grl))/2
  candidates$del.frac = candidates$del.count/(candidates$del.count + candidates$other.count)
  candidates$flank.count = flank(candidates, tile.width) %N% unlist(other$grl) +
    flank(candidates, tile.width, start = FALSE) %N% unlist(all$grl)

  otherc = other$shadow %>% gr.sum
  otherc$score = log(otherc$score+1)
  cov = fishHook::Cov(otherc, field = 'score', type = 'numeric')
  tiles = gr.tile(otherc %>% si2gr, tile.width)
  tiles = c(tiles, tiles %+% (tile.width/2))
  fish = fishHook::Fish(hypotheses = tiles, events = shadows, cov = cov, verbose = FALSE)
  fish$score(nb = FALSE)
  if (return.fish)
    return(fish)
  fish$res %Q% (fdr<fdr.thresh) -> sig

  ## rigma are FDR significant outliers with >= min.count
  sig = gr.fix(sig, candidates, drop = TRUE)

  rigmas = candidates %&% sig
  ## mark nodes and edges associated with rigmas
  if (length(rigmas))
  {
    rigmas = rigmas %Q% rev(order(depth))
    rigmas$footprint = gr.string(rigmas)
    values(rigmas) = cbind(data.frame(rigma = 1:length(rigmas)), values(rigmas))
    rigma.trimmed = rigmas[, 'rigma'] - pmin(1, floor(width(rigmas[, 'rigma'])))
    gg$nodes$mark(rigma = as.integer((gg$nodes$gr %$% rigma.trimmed)$rigma))
    shadows = shadows %$% rigmas[, 'rigma']
    gg$edges[shadows$edge.id]$mark(rigma = shadows$rigma)
    rigmas$type = 'rigma'
    gg$set(rigma = as.data.table(values(rigmas)))
  }

  ## simple dels have no other junctions in vicinity
  ## and are also not sig
  simple.dels = all.dels[!(all.dels %^% (sig %Q% (count>1))) & !(all.dels %^% (other$shadow+tile.width))]
  if (length(simple.dels))
  {
    simple.dels$mark(del = 1:length(simple.dels))
    final.del.gr = simple.dels$shadow
    final.del.gr$del = 1:length(final.del.gr)
    final.del.gr$type = 'del'
    final.del.gr$footprint = final.del.gr %>% gr.stripstrand %>% gr.string
    gg$set(del = final.del.gr %>% values %>% as.data.table)
    final.del.gr.trimmed = (final.del.gr-pmin(floor(width(final.del.gr)/2), 1))
    gg$nodes$mark(del = as.integer((gg$nodes$gr %$% final.del.gr.trimmed)$del))
  }

  if (mark)
  {
    gg$nodes[!is.na(rigma)]$mark(col = mark.col)
    gg$edges[!is.na(rigma)]$mark(col = mark.col)
    gg$nodes[!is.na(del)]$mark(col = mark.col)
    gg$edges[!is.na(del)]$mark(col = mark.col)
  }

  return(gg)
}

#' Find duplications and pyrgos
#' @name dup
#' @description
#' Calls simple duplications (dup) and pyrgo, which are clusters or "towers" of overlapping duplications
#' 
#' Duplications are defined as having low junction copy number and connect two nodes of low junction copy number
#' (with cn and jcn thresholds provided as parameters).  Simple duplications have no other overlapping junctions
#' in their shadow.  Pyrgo have a min.count and are also outliers (<fdr.thresh) in a negative binomial modup
#' that incorporates the regional (non-DUP) junction count in tile.width genomic bins (set to 1 Mbp by default). 
#'
#' Note: Not all DUP-like junctions will be called a dup or pyrgo.  
#' 
#' @param gg gGraph with $cn field annotated on nodes and edges
#' @param fdr.thresh False discovery rate threshold for fishook calculated events.
#' Default: 0.5
#' @param tile.width bin width to use when computing duplication clustering 
#' Default: 1e6
#' @param jcn.thresh edge copy number threshold to call a duplication in ploidy units.
#' Default: 1
#' @param min.count minimum number of overlapping duplications to constitute a pyrgo, 
#' includes the tile.width.
#' Default 2
#' @param return.fish parameter to return FishHook::Fish() output. Default: False
#' @param max.width max width of duplications to consider for pyrgo. 
#' Default: 1e7
#' @param min.width min width of duplications to consider for pyrgo.
#' Default: 1e4
#' @param max.width.flank  max width flank each side of class DEL-like to consider.
#' Default: 1e4
#' @param mark color duplication events. 
#' Default: False
#' @param mark.col color of duplication events. 
#' Default: purple
#' 
#' @details For more details on how to run the function and examples:
#' \href{http://mskilab.com/gGnome/tutorial.html#Tandem_duplications_and_pyrgo}{Duplications & Pyrgo}
#' 
#' @return gGraph with nodes and edges annotated with $dup and $pyrgo metadata field, 
#' and data.tables $meta$pyrgo and $meta$set with event level statistics.
#' 
#' If return.fish = TRUE, returns FishHook::Fish() output.
#' @md
#' @export
dup = function(gg,
               fdr.thresh = 0.5,
               tile.width = 1e6, ## pad around which we don't want to see any other junctions
               jcn.thresh = 1,
               min.count = 2,
               return.fish = FALSE,
               mark = FALSE,
               mark.col = 'purple', 
               min.width = 1e4,
               max.width.flank = 1e4,
               max.width = 1e7)
{
  if (!is(gg, 'gGraph'))
    stop('gg must be gGraph')

  if (length(gg)==0)
    return(gg)

  if (!is.element("cn", colnames(gg$nodes$dt)))
  {
    stop('nodes and edges must have $cn annotation for pyrgo function')
  }

  if (!any(gg$edges$dt[, type=="ALT"])){
    return(gg)
  }

  gg = gGnome::refresh(gg) 
  gg$nodes$mark(pyrgo = as.integer(NA))
  gg$edges$mark(pyrgo = as.integer(NA))
  gg$nodes$mark(dup = as.integer(NA))
  gg$edges$mark(dup = as.integer(NA))
  gg$set(dup = data.table())

  ploidy = ceiling(gg$nodes$dt[!is.na(cn), sum(cn*as.numeric(width))/sum(as.numeric(width))])

  all = gg$edges[type == 'ALT']
  all$mark(cn.max = pmax(all$left$dt$cn, all$right$dt$cn))

  ## define all dups 
  all.dups = all[class == 'DUP-like' &
                 all$dt$cn<=jcn.thresh*ploidy]

  ## candidate pyrgo dups be within certain width threshod
  dups = all.dups[ all.dups$span>=min.width  & 
                   all.dups$span<=max.width]

  if (length(dups)==0)
    return(gg)
  
  ## "other" are non-small dels and inversions / translocations 
  other = all[(class %in% c('DEL-like') & all$span>max.width.flank) | class %in% c('INV-like', 'TRA-like')]

  shadows = dups$shadow
  shadows$edge.id = dups$dt$edge.id

  ## only consider shadows that don't intersect an "other" in model
  shadows = shadows[!(shadows %^% other$shadow)]

  if (length(shadows)==0)
    return(gg)

  foci = gr.sum(shadows) %Q% (score>= min.count) 
  candidates = reduce((shadows+tile.width) %&% foci)-tile.width

  ## collect some stats simply for record keeping 
  candidates$depth = gr.val(candidates, foci, val = 'score', FUN = max, weighted = FALSE)$score
  candidates$min.cn = gr.val(candidates, gg$nodes$gr, val = 'cn', FUN = min, weighted = FALSE)$cn
  candidates$max.cn = gr.val(candidates, gg$nodes$gr, val = 'cn', FUN = max, weighted = FALSE)$cn  
  candidates$dup.count = (candidates %N% unlist(dups$grl))/2
  candidates$other.count = (candidates %N% unlist(other$grl))/2
  candidates$dup.frac = candidates$dup.count/(candidates$dup.count + candidates$other.count)
  candidates$flank.count = flank(candidates, tile.width) %N% unlist(other$grl) +
    flank(candidates, tile.width, start = FALSE) %N% unlist(all$grl)

  otherc = other$shadow %>% gr.sum
  otherc$score = log(otherc$score+1)
  cov = fishHook::Cov(otherc, field = 'score', type = 'numeric')
  tiles = gr.tile(otherc %>% si2gr, tile.width)
  tiles = c(tiles, tiles %+% (tile.width/2))
  fish = fishHook::Fish(hypotheses = tiles, events = shadows, cov = cov, verbose = FALSE)
  fish$score(nb = FALSE)
  if (return.fish)
    return(fish)
  fish$res %Q% (fdr<fdr.thresh) -> sig

  ## pyrgo are FDR significant outliers with >= min.count
  sig = gr.fix(sig, candidates, drop = TRUE)

  pyrgos = candidates %&% sig
  ## mark nodes and edges associated with pyrgos
  if (length(pyrgos))
  {
    pyrgos = pyrgos %Q% rev(order(depth))
    pyrgos$footprint = gr.string(pyrgos)
    values(pyrgos) = cbind(data.frame(pyrgo = 1:length(pyrgos)), values(pyrgos))
    pyrgo.trimmed = pyrgos[, 'pyrgo'] - pmin(1, floor(width(pyrgos[, 'pyrgo'])))
    gg$nodes$mark(pyrgo = as.integer((gg$nodes$gr %$% pyrgo.trimmed)$pyrgo))
    shadows = shadows %$% pyrgos[, 'pyrgo']
    gg$edges[shadows$edge.id]$mark(pyrgo = shadows$pyrgo)
    pyrgos$type = 'pyrgo'
    gg$set(pyrgo = as.data.table(values(pyrgos)))
  }

  ## simple dups have no other junctions in vicinity
  ## and are also not sig
  simple.dups = all.dups[!(all.dups %^% (sig %Q% (count>1))) & !(all.dups %^% (other$shadow+tile.width))]
#  simple.dups = all.dups[!(all.dups %^% sig) & !(all.dups %^% (other$shadow+tile.width))]
  if (length(simple.dups))
  {
    simple.dups$mark(dup = 1:length(simple.dups))
    final.dup.gr = simple.dups$shadow
    final.dup.gr$dup = 1:length(final.dup.gr)
    final.dup.gr$type = 'dup'
    final.dup.gr$footprint = final.dup.gr %>% gr.stripstrand %>% gr.string
    gg$set(dup = final.dup.gr %>% values %>% as.data.table)
    final.dup.gr.trimmed = (final.dup.gr-pmin(floor(width(final.dup.gr)/2), 1))
    gg$nodes$mark(dup = as.integer((gg$nodes$gr %$% final.dup.gr.trimmed)$dup))
  }
  
  if (mark)
  {
    gg$nodes[!is.na(pyrgo)]$mark(col = mark.col)
    gg$edges[!is.na(pyrgo)]$mark(col = mark.col)
    gg$nodes[!is.na(dup)]$mark(col = mark.col)
    gg$edges[!is.na(dup)]$mark(col = mark.col)
  }

  return(gg)
}

#' Amplifications
#' @name amp
#' @description
#' Classifies high-level amplifications in copy number annotated gGraph.  Default parameters
#' should be used. 
#' 
#' @param gg gGraph
#' @param jcn.thresh minimal ALT edge junction threshold to classify a high copy cluster. 
#' Default: 9
#' @param cn.thresh minimal node copy number in ploidy units to classify a high copy cluster.
#' Default: 2
#' @param fbi.cn.thresh fraction of total CN in cluster that is contributed to by 
#' fold back inversions, if higher than this will call a BFB.
#' Default: 0.5
#' @param n.jun.high.bfb.thresh max number of high copy junctions in a fbi.cn 
#' high cluster, if fbi.cn is high and high copy junctions exceed this, then will call a tyfonas.
#' Default: 26
#' @param n.jun.high.dm.thresh double minute threshold set for junctions. If fbi.cn 
#' is low, and high copy junctions exceed this thresh, it will call it cpxdm, else
#' dm.
#' Default: 31
#' @param width.thresh minimum width to consider for an amplification event. 
#' Default: 1e5
#' @param mark.nos (logical) Default: FALSE
#' @param min.nodes (numeric) minimum number of nodes for a cluster to be 
#' designated amp-NOS. Default: 3
#' @param min.jun (numeric) minimum number of aberrant junctions for a cluster 
#' to be designated amp-NOS. Default: 2
#' 
#' @details Amplification events are defined into 3 groups, bfb, tyfonas, and dm events.
#' More details can be found in the tutorial:
#' \href{http://mskilab.com/gGnome/tutorial.html#Complex_amplicons_(bfb,_dm,_tyfonas)}{Complex Amplicons}
#' 
#' @return gg of amplification events found.
#' @md
#' @export
amp = function(gg, jcn.thresh = 8, cn.thresh = 2, fbi.cn.thresh = 0.5,  
               n.jun.high.bfb.thresh = 26, n.jun.high.dm.thresh = 31, width.thresh = 1e5, 
               fbi.width.thresh = 1e5, mc.cores = 1, mark = TRUE, mark.col = 'purple', 
               mark.nos = FALSE, min.nodes = 3, min.jun = 2)
{
    if (mark.nos) {
        gg$nodes$mark(nos = as.integer(NA))
        gg$edges$mark(nos = as.integer(NA))
    }

    gg$nodes$mark(cpxdm = as.integer(NA))
    gg$edges$mark(cpxdm = as.integer(NA))
    gg$nodes$mark(tyfonas = as.integer(NA))
    gg$edges$mark(tyfonas = as.integer(NA))
    gg$nodes$mark(dm = as.integer(NA))
    gg$edges$mark(dm = as.integer(NA))
    gg$nodes$mark(bfb = as.integer(NA))
    gg$edges$mark(bfb = as.integer(NA))
    gg$edges$mark(fbi = gg$edges$class == 'INV-like' & gg$edges$span < fbi.width.thresh)
    gg$set(amp = data.table())
    ploidy = gg$nodes$dt[!is.na(cn), sum(cn*as.numeric(width))/sum(as.numeric(width))]
    keep = (gg$nodes$dt$cn/ploidy) > cn.thresh
    gg$clusters(keep)
    if (!any(!is.na(gg$nodes$dt$cluster)))
        return(gg)

  tiny = gg$edges$mark(tiny = gg$edges$dt$class %in% c('DEL-like', 'DUP-like') & gg$edges$span <1e4)
  ucl = gg$nodes$dt[!is.na(cluster), .(wid = sum(width)), by = cluster][wid > width.thresh, cluster] %>% sort

  amps = mclapply(ucl, function(cl, ploidy) {
    cl.nodes = gg$nodes[cluster == cl]
    cl.edges = cl.nodes$edges[type == "ALT" & tiny == FALSE]
    if (!length(cl.edges)) 
      return(NULL)
    if (!length(cl.edges)) 
      return(NULL)
    data.table(cluster = cl,
               nodes = paste(cl.nodes$dt$node.id,
                             collapse = ","),
               edges = paste(cl.edges$dt$edge.id, 
                             collapse = ","),
               fbi.cn = 2 * sum(cl.edges$dt[fbi == TRUE, sum(cn)]),
               n.jun = length(cl.edges),
               n.jun.high = sum(cl.edges$dt[, sum(cn > 3)]), 
               max.jcn = max(c(0, cl.edges$dt$cn)),
               max.cn = max(cl.nodes$dt$cn),
               footprint = paste(gr.string(cl.nodes$footprint),
                                 collapse = ","))
  }, ploidy, mc.cores = mc.cores) %>% rbindlist

  if (nrow(amps))
  {
      if (!mark.nos) {
          amps = amps[max.jcn >= jcn.thresh,]
          ##amps[max.jcn >= jcn.thresh, ]
      } else {
          ## keep only clusters with a sufficient number of nodes but don't filter by jcn
          amps = amps[max.jcn >= jcn.thresh | 
                      (n.jun >= min.jun &
                       (strsplit(nodes, ",") %>% lapply(length) %>% unlist) >= min.nodes),]
      }
  }


  ## implementing decision tree in https://tinyurl.com/srlbkh2
  if (nrow(amps))
  {
      ## order / rename and mark
      gg$set(amp = amps)

      ## call and mark event types
      amps[, type := ifelse(
                 max.jcn < jcn.thresh,
                 "nos",
                     ifelse(
                         ## few high copy junctions, high fbi cn -> BFB, otherwise tyfonas
                         fbi.cn / max.cn >= fbi.cn.thresh,
                     ifelse(n.jun.high < n.jun.high.bfb.thresh, 
                            'bfb',
                            'tyfonas'),
                     ## few high copy junctions, low fbi cn -> DM, otherwise CPXDM
                     ifelse(n.jun.high >= n.jun.high.dm.thresh, 
                            'cpxdm',     
                            'dm')
                     )
             )]
      
      amps[, ev.id := 1:.N, by = type]

      ## unlist node and edge ids and map back to type and ev label
      nodelist = strsplit(amps$nodes, ',') %>% lapply(as.integer) %>% dunlist
      edgelist = strsplit(amps$edges, ',') %>% lapply(as.integer) %>% dunlist
      nodelist = cbind(nodelist, amps[nodelist$listid, .(type, ev.id)])
      edgelist = cbind(edgelist, amps[edgelist$listid, .(type, ev.id)])
    
      nodelist[, {
          if (type == 'dm')
          {
              gg$nodes[V1]$mark(dm = ev.id)
          }
          else if (type == 'tyfonas')
          {
              gg$nodes[V1]$mark(tyfonas = ev.id)
          }
          else if (type == 'cpxdm')
          {
              gg$nodes[V1]$mark(cpxdm = ev.id)
          }
          else if (type == "bfb")
          {
              gg$nodes[V1]$mark(bfb = ev.id)
          }
          else
          {
              gg$nodes[V1]$mark(nos = ev.id)
          }
      }, by = type]
      
      edgelist[, {
          if (type == 'dm')
          {
              gg$edges[V1]$mark(dm = ev.id)
          }
          else if (type == 'tyfonas')
          {
              gg$edges[V1]$mark(tyfonas = ev.id)
          }
          else if (type == 'cpxdm')
          {
              gg$edges[V1]$mark(cpxdm = ev.id)
          }
          else if (type == 'bfb')
          {
              gg$edges[V1]$mark(bfb = ev.id)
          }
          else
          {
              gg$edges[V1]$mark(nos = ev.id)
          }
      }, by = type]
          
      if (mark)
      {
          gg$nodes[!is.na(tyfonas)]$mark(col = mark.col)
          gg$edges[!is.na(tyfonas)]$mark(col = mark.col)
          
          gg$nodes[!is.na(dm)]$mark(col = mark.col)
          gg$edges[!is.na(dm)]$mark(col = mark.col)

          gg$nodes[!is.na(bfb)]$mark(col = mark.col)
          gg$edges[!is.na(bfb)]$mark(col = mark.col)

          gg$nodes[!is.na(cpxdm)]$mark(col = mark.col)
          gg$edges[!is.na(cpxdm)]$mark(col = mark.col)

          if (mark.nos) {
              gg$nodes[!is.na(nos)]$mark(col = mark.col)
              gg$edges[!is.na(nos)]$mark(col = mark.col)
          }
      }
  }

  return(gg)
}

#' get microhomology
#' @name microhomology
<<<<<<< HEAD
#' @description Computes microhomology at 5bp, 10bp, 50bp, and 100bp windows
#'  around ALT junctions of input gGraph (or Junction object)
#' gg and adds these as an edge annotation to the appropriate edges.
=======
#' @title microhomology
#' 
#' @description
#' Computes microhomology at junction breakends
#'
#' @details
#' Computes microhomology at 5bp, 10bp, 50bp, and 100bp windows around ALT junctions of input gGraph (or Junction object) gg and adds these as an edge annotation to the appropriate edges.
#'
#' The default behavior is to compute the maximum microhomology using local alignment across the entire window. However, the longest common prefix within each window can be specified by setting the argument prefix_only to TRUE.
#'
#' Care should be taken that the sequence names of junctions are consistent with those provided in the reference. There will be an error if the sequence names of the junction are not a subset of those of the reference, if ignore_missing is FALSE (default). If ignore_missing is TRUE, then those junctions with missing seqnames will be assigned score -1.
>>>>>>> upstream/master
#'
#' Requires Biostrings.
#' 
#' @param gg gGraph or Junctions
#' @param hg DNAStringSet or path to reference fasta
<<<<<<< HEAD
#' @return gGraph with $pyrgo marking on nodes and edges labeling unique "events"
#' 
=======
#' @param prefix_only (logical) default FALSE. if TRUE, considers only the longest common prefix. if FALSE, considers the longest local alignment.
#' @param pad (numeric) default NA (use the default window lengths of 5, 10, 50, and 100). otherwise, an integer specifying window length.
#' @param ignore_missing (logical) ignore junctions where at least one breakend is not found on the reference, and return -1 for microhomology. default FALSE, which will cause an error.
#' 
#' @return gGraph with edges augmented with metadata mh labeling unique "events"
>>>>>>> upstream/master
#' @export
microhomology = function(gg, hg, prefix_only = FALSE, pad = c(5, 10, 50, 100), ignore_missing = FALSE)
{
  if (!requireNamespace("Biostrings", quietly = TRUE)) {
      stop('You must have the package "Biostrings" installed in order for this function to work. Please install it.')
  }
  if (inherits(gg, 'gGraph'))
  {
    gg = gg$clone()
    
    if ((!length(gg$edges)) || (!length(gg$edges[type == 'ALT'])))
    {
      gg$edges$mark(mh5 = NA_integer_)
      gg$edges$mark(mh10 = NA_integer_)
      gg$edges$mark(mh50 = NA_integer_)
      gg$edges$mark(mh100 = NA_integer_)
      return(gg)
    }

    ed = gg$edges[type == 'ALT']
    
    ## note: the region why this should be flipped
    ## is that the homology should be on the SAME strand
    ## if the junction joins "opposite" strands (e.g. deletions, duplications)
    ## and on the OPPOSITE strand
    ## if the junction joins the "same" strand (e.g. inversions)
    bp1 = ed$junctions$left %>% gr.flipstrand
    bp2 = ed$junctions$right
  }
  else if (inherits(gg, 'Junction'))
  {
    if (!length(gg))
      return(gg)
    bp1 = gg$left %>% gr.flipstrand
    bp2 = gg$right
  }
  else
    stop('Input must be either gGraph or Junction object')

  if (is.character(hg))
    hg = Biostrings::readDNAStringSet(hg)

  ## use seqlengths from supplied reference
  bp1 = dt2gr(gr2dt(bp1), seqlengths = seqlengths(hg), seqinfo = seqinfo(hg))
  bp2 = dt2gr(gr2dt(bp2), seqlengths = seqlengths(hg), seqinfo = seqinfo(hg))

  ## check whether all seqnames in bp1 and bp2 are present in the supplied reference
  if (ignore_missing) {
      ## if we are simply ignoring missing seqnames
      ## subset bp1 and bp2 to include only breakends on valid contigs
      bp1.og = bp1[, c()]
      bp2.og = bp2[, c()]
      keep = (as.character(seqnames(bp1)) %in% seqlevels(hg)) & (as.character(seqnames(bp2)) %in% seqlevels(hg))
      new.index = match(ed$dt$edge.id, ed$dt$edge.id[keep])
      bp1 = bp1[keep]
      bp2 = bp2[keep]
  } else {

      ## otherwise error out if there are discrepancies
      if (length(setdiff(c(seqnames(bp1), seqnames(bp2)), seqlevels(hg)))) {
          stop('seqnames in breakpoints missing from the provided reference, plesae check and fix the seqlevels of the provided graph / junctions / and/or reference')
      }
      new.index = 1:length(bp1)##match(ed$dt$edge.id, ed$dt$edge.id)
  }
      

  ## define some internal functions
  ## workaround weird sudden biostrings xstring 2^31 unlist problem ?!?!!?
  dodo.call = function (FUN, args) 
  {
    if (!is.character(FUN)) 
      FUN = substitute(FUN)
    cmd = paste(FUN, "(", paste("args[[", 1:length(args), "]]", 
                                collapse = ","), ")", sep = "")
    return(eval(parse(text = cmd)))
  }

  ## grab sequence associated with certain genomic range
  ## reverse complementing if the strand of the range is negative
  .getseq = function(hg, gr)
    {
      res = dodo.call('c', mapply(function(c,s,e) Biostrings::subseq(hg[c], start = s, end = e), seqnames(gr) %>% as.character, start(gr), end(gr)))
      res = ifelse(strand(gr)=='+', res, Biostrings::reverseComplement(res))
      res = Biostrings::DNAStringSet(res)
      return(res)
    }

  ## create base substitution penalty matrix for local alignment
  ## (kind of overkill - basically allow exact matches, and give anything else a score of zero, lol)
  .mat = function(match = 1, mismatch = 0, baseOnly = FALSE, type = "DNA", letters = NULL) 
  {
      "%safemult%" <- function(x, y) ifelse(is.infinite(x) & y == 
                                            0, 0, x * y)
      nLetters <- length(letters)
      splitLetters <- strsplit(letters, split = "")
      submat <- matrix(0, nrow = nLetters, ncol = nLetters, dimnames = list(names(letters), 
                                                                            names(letters)))
      for (i in 1:nLetters) for (j in i:nLetters) submat[i, j] <- submat[j, 
                                                                         i] <- mean(outer(splitLetters[[i]], splitLetters[[j]], 
                                                                                          "=="))
      abs(match) * submat - abs(mismatch) %safemult% (1 - submat)
  }

  ## create substitution matrix
  letters = Biostrings::alphabet(hg)
  names(letters) = letters
  mat = .mat(match = 1, mismatch = -1000, baseOnly = TRUE, letters = letters)

  ## loop over desired alignment length
  for (pad.length in pad) {

      ## grab the sequence at each breakend
      ## but specifically only the fused side
      ## since the strand was flipped for bp1, we want to check the end
      bp1.gr = trim(resize(bp1, width = pad.length, fix = "end")[, c()])
      ## for bp2 we can keep the original strand
      bp2.gr = trim(resize(bp2, width = pad.length, fix = "start")[, c()])
      ## get sequence as character vector (needed for lcprefix and lcsubstr)
      seq1 = .getseq(hg, bp1.gr)
      seq2 = .getseq(hg, bp2.gr)

      if (!is.na(prefix_only) && (prefix_only)) {
          ## need character vectors
          seq1 = as.character(seq1)
          seq2 = as.character(seq2)
          mh.score = sapply(1:length(seq1),
                            function(ix) {
                                ## what is the longest suffix of the left breakend
                                ## that is the prefix of the right breakend?
                                ## O(N^2) but whatever hehe :)
                                ## hehe RIP algorithms
                                mh.match = sapply(1:pad.length,
                                                  function(suffix_length) {
                                                      return(base::substring(seq2[ix], 1, suffix_length) == base::substring(seq1[ix], pad.length - suffix_length + 1, pad.length))
                                                  })
                                return(sum(mh.match))
                                ## because the strand was flipped
                                ## Biostrings::lcprefix(paste(rev(strsplit(seq1[ix], "")[[1]]),
                                   ##                         collapse = ""),
                                      ##                seq2[ix])
                            })
      } else {
          mh.score = Biostrings::pairwiseAlignment(seq1,
                                                   seq2,
                                                   substitutionMatrix = mat,
                                                   gapOpening = 1000, ## note: high penalty effectively forces contiguity
                                                   gapExtension = 1000,
                                                   type = 'local',
                                                   scoreOnly = TRUE)
      }

      ## map back to original index
      mh.score = mh.score[new.index]
      ## if there were some invalid breakends, set the microhomology score of those to -1
      mh.score[is.na(mh.score)] = -1

      ## set score in original gGraph
      if (inherits(gg, "gGraph")) {
          cmd = sprintf("ed$mark(mh%s = mh.score)", pad.length)
      } else {
          cmd = sprintf("gg$set(mh%s = mh.score)", pad.length)
      }
      eval(parse(text = cmd))
  }

  return(gg)
}

#' Get reciprical connected junctions
#' @name reciprocal
<<<<<<< HEAD
#' @description Identifies reciprocally connected junctions,
#' i.e. breakends from non-identical junctions that are "linked"
#' by an inter-breakpoint distance less than a given threshold.
#' Edges and nodes are marked by the "ecluster" metadata field
#'
=======
#' @title
#' @description
#'
#' Identifies reciprocally connected junctions,
#'
#' @details
#' Reciprocal junctions are junctions with breakends that are mutually adjacent and opposite.
#' 
#' i.e. breakends from non-identical junctions that are "linked"
#' by an inter-breakpoint distance less than a given threshold.
#' Edges and nodes are marked by the "ecluster" metadata field
#' 
>>>>>>> upstream/master
#' @param gg gGraph
#' @param thresh threshold for edge clusters. Default: 5e5
#' @param max.small max small threshold for edge clusters. Default: 1e4 
#' 
#' @return gGraph with $ecluster marking on nodes and edges labeling unique reciprocal events
#' @export
reciprocal = function(gg, thresh = 5e5, max.small = 1e4) {
  gg = gGnome::refresh(gg)
  gg$eclusters(thresh = thresh, max.small = max.small, only_chains = TRUE)
  gg$nodes$mark(ecluster = NA_integer_)
  gg$set(recip_event = data.table())
  eclust.edge = gg$edges[!is.na(ecluster)]
  if (length(eclust.edge)) {
    edt = eclust.edge$dt
    ndt = melt(edt, id.vars = c("ecluster"),
               measure.vars = c("n1", "n2"))[!duplicated(cbind(ecluster, value))]
    ndt$footprint = gr.string(gg$nodes[ndt$value]$gr)
    ndt = ndt[, .(footprint = paste(unique(footprint), collapse = ","),
                  nnodes = length(unique(value))), by = ecluster]
    recip_bp = merge(gg$meta$recip_bp[!is.na(ecluster)], ndt, by = "ecluster")
    recip_event = recip_bp[
     ,.(njuncs = nclust[1], nnodes = nclust[1],
        all_positive = unique(all_positive),
        all_negative = unique(all_negative),
        mixed = unique(mixed),
        bridge = unique(bridge),
        footprint = footprint[1]), by = ecluster]
    gg$set(recip_event = recip_event)
    gg$nodes[ndt$value]$mark(ecluster = ndt$ecluster)
  }
  return(gg)
}

# Finds Quasi-Reciprocal Pairs of Junctions
#' @name qrp
#' 
#' @description Finds (quasi) reciprocal pairs of junctions.  Very related to 
#' chromoplexy or tic "cycles", but with exactly two junctions.
#'
#' @param gg gGraph
#' @param thresh integer - maximal size of bridge.
#' Default: 1e6
#' @param max.small integer indicating maximum size of candidate SVs for clustering.
#' Default: 1e5
#' @param breakend_pairing "strict": can only be "monogamously" matched to one 
#' breakend and if the nearest breakend is of the wrong orientation, 
#' it is thrown out; "one_to_one:" breakends can only be coupled to a single 
#' "monogamous" match but without considering breakends of the wrong orientation. 
#' If breakend B is nearest to C, B will only be matched to C. If A's nearest breakend 
#' is B but is further away than C, A will not be matched to B. "loose": the nearest 
#' breakend in the correct orientation under the threshold is considered. The same 
#' as one_to_one except A will be matched to B, while B will be matched to C. 
#' Default: c("strict", "one_to_one", "loose")
#' @param mark logical, mark the edges and nodes in color specified by mark.col.
#' Default: True
#' @param mark.col character, color to mark nodes/edges involved in any QRP.
#' Default: purple
#' @return gGraph with $ecluster marking on nodes and edges labeling unique reciprocal events
#' @author Kevin Hadi
#' @export
qrp = function(gg, thresh = 1e6, max.small = 1e5,
               breakend_pairing = c("strict", "one_to_one", "loose"),
               mark = TRUE, mark.col = "purple") {
    if (identical(breakend_pairing, c("strict", "one_to_one", "loose"))){
        breakend_pairing = "strict"
    } else if (length(breakend_pairing) > 1) {
        breakend_pairing = intersect(breakend_pairing[1], c("strict", "one_to_one", "loose"))
        stopifnot(length(breakend_pairing) > 0)
    }
    
    gg$eclusters2(thresh = thresh,
                  max.small = max.small,
                  only_chains = TRUE,
                  ignore.small = TRUE,
                  strict = breakend_pairing,
                  ignore.isolated = TRUE)

    recip_event = copy3(gg$meta$recip_event)

    gg$edges$mark(qrpmix = NA_integer_)
    gg$edges$mark(qrpmin = NA_integer_)
    gg$edges$mark(qrppos = NA_integer_)

    gg$nodes$mark(qrpmix = NA_integer_)
    gg$nodes$mark(qrpmin = NA_integer_)
    gg$nodes$mark(qrppos = NA_integer_)
    
    if (NROW(recip_event)) {

        qrppos = recip_event[bridge == FALSE][njuncs == 2][all_positive == TRUE][
           ,.(ecluster = ecluster, qrppos = rleseq(ecluster,clump=T)$idx, footprint = footprint)] %>% unique
        qrppos[["type"]] = rep_len2("qrppos", qrppos)
        gg$set(qrppos = qrppos[seq_along2(qrppos)])

        qrpmix = recip_event[mixed == TRUE][num_negative == 1][bridge == FALSE][njuncs == 2][
           ,.(ecluster = ecluster, qrpmix = rleseq(ecluster,clump=T)$idx, footprint = footprint)] %>% unique
        qrpmix[["type"]] = rep_len2("qrpmix", qrpmix)
        gg$set(qrpmix = qrpmix)
        
        qrpmins = recip_event[bridge == FALSE][njuncs == 2][all_negative == TRUE][
           ,.(ecluster = ecluster, qrpmin = rleseq(ecluster,clump=T)$idx, footprint = footprint)] %>% unique
        qrpmins[["type"]] = rep_len2("qrpmin", qrpmins)
        gg$set(qrpmin = qrpmins)

        qrpmixix = gg$meta$qrpmix[match3(gg$edges$dt$ecluster, gg$meta$qrpmix$ecluster)]$qrpmix
        qrpminix = gg$meta$qrpmin[match3(gg$edges$dt$ecluster, gg$meta$qrpmin$ecluster)]$qrpmin
        qrpposix = gg$meta$qrppos[match3(gg$edges$dt$ecluster, gg$meta$qrppos$ecluster)]$qrppos

        gg$edges$mark(qrpmix = qrpmixix)
        gg$edges$mark(qrpmin = qrpminix)
        gg$edges$mark(qrppos = qrpposix)

        if (isTRUE(mark)) {
            eid = which(gg$edges$dt[, !is.na(qrpmix) | !is.na(qrpmin) | !is.na(qrppos)]) ## can't seem to index with gg$edges[] directly with certain cases?
            gg$edges[eid]$mark(col = mark.col)
        }

        qrpmixix = gg$meta$qrpmix[match3(gg$nodes$dt$ecluster, gg$meta$qrpmix$ecluster)]$qrpmix
        qrpminix = gg$meta$qrpmin[match3(gg$nodes$dt$ecluster, gg$meta$qrpmin$ecluster)]$qrpmin
        qrpposix = gg$meta$qrppos[match3(gg$nodes$dt$ecluster, gg$meta$qrppos$ecluster)]$qrppos

        gg$nodes$mark(qrpmix = qrpmixix)
        gg$nodes$mark(qrpmin = qrpminix)
        gg$nodes$mark(qrppos = qrpposix)

        if (isTRUE(mark)) {
            nid = which(gg$nodes$dt[, !is.na(qrpmix) | !is.na(qrpmin) | !is.na(qrppos)]) ## can't seem to index with gg$nodes[] directly with certain cases?
            gg$nodes[nid]$mark(col = mark.col)
        }
    }

    return(gg)
    
}

#' Events to GRanges
#' @name events.to.gr
#' @description Extract event annotation as a GRanges
#'
#' @param gg gGraph
#' 
#' @return GRanges containing ranges of annotated events along with all metadata 
#' from gg$meta$events
#' @author Alon Shaiber
#' @export
events.to.gr = function(gg){
    if (!inherits(gg, 'gGraph')){
        stop('Expected gGraph, but got: ', class(gg))
    }
    if (!('events' %in% names(gg$meta))){
        stop('Missing events field in gGraph meta. Are you sure you ran gGraph::events()?')
    }
    if (gg$meta$events[,.N] == 0){
        warning('No events annotated in this gGraph')
        return(GRanges())
    }
    ggrl = parse.grl(gg$meta$events$footprint)
    mcols(ggrl) = gg$meta$events
    ggr = grl.unlist(ggrl)
    return(ggr)
}
