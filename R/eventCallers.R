#' @name fusions
#' @export
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
fusions = function(graph = NULL,
                   gencode = NULL,
                   prom.window = 1e3,
                   exhaustive = FALSE,
                   annotate.graph = TRUE, 
                   max.chunk = 1e10, ## parameters for gr.findoverlaps
                   mc.cores = 1,
                   verbose = FALSE){
    ## QC input graph or junctions

  if (!inherits(graph, "gGraph")){
    stop('Input must be be gGraph')
  }

  seg = graph$gr ## added
  A = graph$adj ## added       


  if (is.null(A) | is.null(seg))
    stop('Some essential args are NULL')
  
  ## loading default cds
  if (is.null(gencode))
  {
    default.path = Sys.getenv("GENCODE_DIR")
    
    if (file.url.exists(default.path))
    {
      warning(paste('GENCODE object missing, attempting to read default from file path / URL of gtf, gff3, or .rds provided in', default.path))
      if (grepl('.rds', default.path))
      {
        read.rds.url(default.path)
      } else
      {
        gencode = rtracklayer::import()
      }
    } else
    {
      
          stop('GENCODE object not provided (as .gtf, gtf.gz, .gff, .gff.gz, or .rds file) and GENCODE_DIR is not set or pointing to a non-existent path')
    }
    }
  
  ## parse GENCODE into CDS
  seg = gr.fix(seg, gencode)
  gencode = gr.fix(gencode, seg)

  cds = gencode[gencode$type == 'CDS']
  genes = gencode[gencode$type == 'gene']
  transcripts = as.data.table(gencode[gencode$type == 'transcript'])
  setkey(transcripts, transcript_id)

  cds = split(cds[, c('exon_number')], cds$transcript_id)
  tx.span = dt2gr(transcripts[names(cds), ], seqlengths = seqlengths(cds))

  if (verbose){
    message('got transcript boundaries')
  }

  ## determine set of transcript fragments
  ## these correspond to transcripts that intersect a segment boundary
  cds.frag.left = gUtils::gr.findoverlaps(tx.span, gr.start(seg),
                                          qcol = c('gene_name', 'transcript_id'),
                                          ignore.strand = F, max.chunk = max.chunk)
  strand(cds.frag.left) = strand(tx.span)[cds.frag.left$query.id]
  cds.frag.right = gUtils::gr.findoverlaps(tx.span, gr.end(seg),
                                           qcol = c('gene_name', 'transcript_id'),
                                           ignore.strand = F, max.chunk = max.chunk)
  strand(cds.frag.right) = strand(tx.span)[cds.frag.right$query.id]

  ## I want to find all unique walks that involve tx fragments
  if (length(cds.frag.left)>0 & length(cds.frag.right) > 0 )
    tmp = merge(data.frame(i = 1:length(cds.frag.left), key1 = cds.frag.left$query.id, key2 = cds.frag.left$subject.id),
                data.frame(j = 1:length(cds.frag.right), key1 = cds.frag.right$query.id, key2 = cds.frag.right$subject.id), all = T)
  else
    return(gWalk$new(graph = graph))

  pos.right = which(as.logical( strand(cds.frag.right)=='+'))
  pos.left = which(as.logical(strand(cds.frag.left)=='+'))
  neg.right = which(as.logical(strand(cds.frag.right)=='-'))
  neg.left = which(as.logical( strand(cds.frag.left)=='-'))

  ## positive start fragments will be "right" fragments
  cds.start.frag.pos = cds.frag.right[tmp[is.na(tmp$i) & tmp$j %in% pos.right, ]$j]
  start(cds.start.frag.pos) = start(tx.span)[cds.start.frag.pos$query.id]
  if (length(cds.start.frag.pos)>0)
    cds.start.frag.pos$type = 'start'

  ## positive end fragments will be "left" fragments
  cds.end.frag.pos = cds.frag.left[tmp[is.na(tmp$j) & tmp$i %in% pos.left, ]$i]
  end(cds.end.frag.pos) = end(tx.span)[cds.end.frag.pos$query.id]
  if (length(cds.end.frag.pos)>0)
    cds.end.frag.pos$type = 'end'

  ## negative start fragments will be "right" fragments
  cds.start.frag.neg = cds.frag.left[tmp[is.na(tmp$j) & tmp$i %in% neg.left, ]$i]
  end(cds.start.frag.neg) = end(tx.span)[cds.start.frag.neg$query.id]
  if (length(cds.start.frag.neg)>0)
    cds.start.frag.neg$type = 'start'

  ## negative end fragments will be "left" fragments
  cds.end.frag.neg = cds.frag.right[tmp[is.na(tmp$i) & tmp$j %in% neg.right, ]$j]
  start(cds.end.frag.neg) = start(tx.span)[cds.end.frag.neg$query.id]
  if (length(cds.end.frag.neg)>0)
    cds.end.frag.neg$type = 'end'

  ## remaining will be "middle" fragments
  middle.frag = cds.frag.left[tmp[!is.na(tmp$i) & !is.na(tmp$j),]$i]
  end(middle.frag) = end(cds.frag.right[tmp[!is.na(tmp$i) & !is.na(tmp$j),]$j])
  if (length(middle.frag)>0)
    middle.frag$type = 'middle'

  ## concatenate fragments
  ## subject.id of frags is the id of the node on the graph

  all.frags = c(cds.start.frag.pos, cds.end.frag.pos, cds.start.frag.neg, cds.end.frag.neg, middle.frag)


  ## whatever segments are left out... include as extra "bridging" nodes
  extra_seg_ids = setdiff(seq_along(seg), all.frags$subject.id)
  extra_ids = integer(0)
  if (length(extra_seg_ids) > 0) {
    tmp_seg = seg[extra_seg_ids][,c()]
    tmp_seg$subject.id = extra_seg_ids
    tmp_seg$type = "middle"
    extra_ids = seq(length(all.frags) + 1, by = 1, length.out = length(extra_seg_ids))
    all.frags = grbind(all.frags, tmp_seg)
  }


  ## now connect all.frags according to A
  ## i.e. apply A connections to our fragments, so draw an edge between fragments
  ## if
  ## (1) there exists an edge connecting segment and
  ## (2) only allowable connections are 'start' --> 'middle' --> 'middle' --> 'end'
  seg.edges = as.data.frame(Matrix::which(A!=0, arr.ind = T))
  colnames(seg.edges) = c('from.seg', 'to.seg')
  edges = merge(merge(data.frame(i = 1:length(all.frags), from.seg = all.frags$subject.id),
                      seg.edges), data.frame(j = 1:length(all.frags), to.seg = all.frags$subject.id))

  edges = edges[all.frags$type[edges$i] == 'start' & all.frags$type[edges$j] == 'middle' |
                all.frags$type[edges$i] == 'start' & all.frags$type[edges$j] == 'end' |
                all.frags$type[edges$i] == 'middle' & all.frags$type[edges$j] == 'middle' |
                all.frags$type[edges$i] == 'middle' & all.frags$type[edges$j] == 'end', ]

  ## this removes splice variants .. keeping only links that fuse different genes or same transcripts
  ## edges = edges[tx.span$gene_name[all.frags$query.id[edges$i]] != tx.span$gene_name[all.frags$query.id[edges$j]] |
  ##               all.frags$query.id[edges$i] == all.frags$query.id[edges$j],]
  na2true = function(v) {
    v[is.na(v)] = TRUE
    as.logical(v)
  }
  edges = edges[na2true(tx.span$gene_name[all.frags$query.id[edges$i]] != tx.span$gene_name[all.frags$query.id[edges$j]] |
                        all.frags$query.id[edges$i] == all.frags$query.id[edges$j]),]

  if (nrow(edges)==0){
    return(GRangesList())
  }

  if (verbose){
    message('computed subgraph')
  }

  A.frag = sparseMatrix(edges$i, edges$j, x = 1, dims = rep(length(all.frags),2))
#  keep.nodes = which(Matrix::rowSums(A.frag)>0 | Matrix::colSums(A.frag)>0)
#  A.frag = A.frag[keep.nodes, keep.nodes]
#  all.frags = all.frags[keep.nodes]

  sources = which(all.frags$type == 'start')
  sinks = which(all.frags$type == 'end')

  G = graph.adjacency(A.frag)
  C = igraph::clusters(G, 'weak')
  vL = split(1:nrow(A.frag), C$membership)
  vL = vL[elementNROWS(vL) > 1] ## new
  paths = do.call('c', mclapply(1:length(vL), function(i) {
    if (verbose & (i %% 10)==0){
      message(i, ' of ', length(vL))
    }
    x = vL[[i]]

    tmp.source = setdiff(match(sources, x), NA)
    tmp.sink = setdiff(match(sinks, x), NA)
    tmp.mat = A.frag[x, x, drop = FALSE]!=0

    if (length(x)<=1){
      return(NULL)
    }

    if (length(x)==2){
      list(x[c(tmp.source, tmp.sink)])
    }else if (all(Matrix::rowSums(tmp.mat)<=1) & all(Matrix::colSums(tmp.mat)<=1)){
      get.shortest.paths(G, from = intersect(x, sources), to = intersect(x, sinks))$vpath
    } else {
      if (exhaustive){
        lapply(JaBbA:::all.paths(A.frag[x,x, drop = FALSE],
                                 source.vertices = tmp.source,
                                 sink.vertices = tmp.sink,
                                 verbose = verbose)$paths,
               function(y) x[y])
      } else {
        ## ALERT: possible wrong syntax!!!!!!!
        out = do.call('c',
                      lapply(intersect(x, sources),
                             function(x, sinks) suppressWarnings(get.shortest.paths(G, from = x, to = sinks)$vpath), sinks = intersect(x, sinks))
                      )
        out = out[sapply(out, length)!=0]
        if (length(out)>0)
          out = out[!duplicated(sapply(out, paste, collapse = ','))]
        return(out)
      }
    }
  }, mc.cores = mc.cores))

  if (verbose){
    message('computed paths')
  }

  paths.u = unlist(paths)
  ## paths.i = unlist(lapply(1:length(paths), function(x) rep(x, length(paths[[x]]))))
  paths.i = rep(seq_along(paths), times = elementNROWS(paths))
  walks = split(seg[all.frags$subject.id[paths.u]], paths.i)
  values(walks)$seg.id = split(all.frags$subject.id[paths.u], paths.i)

  ## for now just want to pick the non duplicated paths on the original graph and send these to the walk annotation module
  walks = walks[!duplicated(sapply(values(walks)$seg.id, function(x) paste(x, collapse = ',')))]

  if (verbose){
    message(sprintf('Annotating %s walks', length(walks)))
  }

  if (length(walks)==0){
    return(gWalk$new(graph = graph)) ## return empty gWalks
  }

  names(walks) = 1:length(walks)
  awalks = annotate.walks.with.cds(walks,
                                   cds,
                                   tx.span,
                                   prom.window = prom.window,
                                   verbose = verbose,
                                   exhaustive = FALSE,
                                   mc.cores = mc.cores)

  if (verbose)
    message('Finished annotation walks, creating gWalk')

  values(awalks)$seg.id = NULL
  values(awalks)$coords = NULL

  gw = gW(grl = awalks, graph = graph, meta = as.data.table(values(awalks)))

  if (annotate.graph)
  {  
    if (verbose)
      message('Annotating gGraph with GENCODE elements')    
    gt = gt.gencode(gencode[gr.in(gencode, f$footprint)])
    annotations = dat(gt)[[1]] ## stored in gTrack
    gw$disjoin(gr = annotations)
    gw$set(name = gw$dt$genes)
    gw$graph$set(colormap = colormap(gt)[1])
  }

  return(gw)
}


#' @name annotate.walks.with.cds
#' @rdname internal
#' @description
#' 
#'
#' Internal function that annotates grl with cds info from gencode
#' 
#' @param walks GRangesList of walks
#' @param cds  GRangesList of cds
#' @param transcript GRanges of transcripts, same length as cds
#' @param filter.splice logical flag (default = TRUE) specifying whether to remove splice isoforms of each transcript (ie keep one arbitrary)
#' @param prom.window window to use around each transcript to identify putative promoter if promoter is NULL
#' @param verbose logical flag whether to print verbose
#' @param mc.cores integer number of cores to parallelize over
#' @param exhaustive logical flag whether to report all possible transcript variants 
#' @return GRangesList of walks annotated with cds info
#' @keywords internal
#' @noRd
annotate.walks.with.cds = function(walks, cds, transcripts, filter.splice = T, verbose = F, prom.window = 1e3, max.chunk = 1e9, mc.cores = 1, exhaustive = exhaustive)
{
  if (inherits(walks, 'GRanges'))
    walks = GRangesList(walks)

  if (is(walks, 'list'))
    walks = do.call(GRangesList, walks)

  tx.span = transcripts

  cdsu = gr2dt(grl.unlist(cds)[, c('grl.ix')])
  setkey(cdsu, grl.ix)

  ## There are negative width CDS in the new GENCODE v27!!!!
  ## KH: how was this determined?
  ## any(mcols(cds)$End - mcols(cds)$Start <= 0) == FALSE


  ## cds.span comprises the start of the first and coordinates of the last cds in each transcript 
  cds.span = cdsu[, list(start = min(start), end = max(end)), keyby = grl.ix][list(1:length(tx.span)), ]
  setkey(cds.span, grl.ix)

  ## figuring out the ranges of the left and right UTRs (if they exist)
  ## basically doing this by taking the interval between the left and right ends of the tx.span
  ## and the beginning of cds.span
  utr.left.dt = gr2dt(tx.span)[
  , list(seqnames = seqnames,
         start = start,
         strand = strand,
         end = cds.span[list(1:length(start)), start],
         transcript_id,
         transcript_name,
         gene_name)]

  utr.right.dt = gr2dt(tx.span)[
  , list(seqnames = seqnames,
         start = cds.span[list(1:length(start)), end],
         strand = strand,
         end = end,
         transcript_id,
         transcript_name,
         gene_name)]

  ## DO a check of eligibility before converting to GRanges
  ## MOMENT
  trash.ix = which(utr.left.dt[, start>end] |
                   utr.right.dt[, start>end])

  if (length(trash.ix)>0){
    message(sprintf("some of the provided CDS annotations are outside the bounds of the provided transcripts, removing all %s of such transcripts", length(trash.ix)))
    utr.left = dt2gr(utr.left.dt[-trash.ix], seqlengths = seqlengths(cds))
    utr.right = dt2gr(utr.right.dt[-trash.ix], seqlengths = seqlengths(cds))
  }  else
  {
    utr.left = dt2gr(utr.left.dt, seqlengths = seqlengths(cds))
    utr.right = dt2gr(utr.right.dt, seqlengths = seqlengths(cds))    
  }


  utr = c(utr.left, utr.right)
  
  promoters = flank(tx.span, prom.window)
  values(promoters) = values(tx.span)

  tx.span$type = 'cds'
  promoters$type = 'gene'
  utr.left$type = 'gene'
  utr.right$type = 'gene'
  tx.span$cds.id = 1:length(tx.span)
  promoters$cds.id = 1:length(promoters) ## KMH added
  tx.span = seg2gr(as.data.table(rrbind(as.data.frame(tx.span), as.data.frame(promoters)))[, list(seqnames = seqnames[1], start = min(start), end = max(end), strand = strand[1], gene_name = gene_name[1], transcript_id = transcript_id[1], transcript_name = transcript_name[1], cds.id = cds.id), keyby = cds.id][!is.na(cds.id), ], seqlengths = seqlengths(tx.span)) ## promoters get thrown out?? what is the point of the is.na(cds.id) filter when promoters will not have cds.id??
  

                                        # match up tx.span to walks
  walks.u = grl.unlist(walks)

  ## these are fragments of transcripts that overlap walks
  this.tx.span = gr.findoverlaps(tx.span, walks.u, qcol = c('transcript_id', 'transcript_name', 'gene_name', 'cds.id'), verbose = verbose, max.chunk = max.chunk)
  this.tx.span$tx.id = this.tx.span$query.id

  strand(this.tx.span) = strand(tx.span)[this.tx.span$query.id]
  this.tx.span$left.broken = start(this.tx.span) != start(tx.span)[this.tx.span$query.id]
  this.tx.span$right.broken = end(this.tx.span) != end(tx.span)[this.tx.span$query.id]

                                        # remove elements that are "unbroken" by the window boundaries
  this.tx.span = this.tx.span[this.tx.span$left.broken | this.tx.span$right.broken]
  this.tx.span$cds.sign = c('-'= -1, '+' = 1)[as.character(strand(this.tx.span))]
  this.tx.span$window.sign = c('-'= -1, '+' = 1)[as.character(strand(walks.u)[this.tx.span$subject.id])]

                                        # annotate left and right ends (if information available)
                                        # i.e. UTR, CDS, Promoter

  ## we trim cds by 1 nucleotide at both ends so that 'cds' annotation only refers to a cds fragment (not including start and stop)
  ## annotate ends with feature types so that positions internal to cds bases will be 'cds',
  ## external to cds but in cds will 'utr'
  ## and 5' flank (with respect to cds orientation) will be promoter
  this.tx.span$left.feat = NA
  this.tx.span$left.feat[gr.findoverlaps(gr.start(this.tx.span), tx.span, by = 'transcript_id', max.chunk = max.chunk)$query.id] = 'cds'
  this.tx.span$left.feat[gr.findoverlaps(gr.start(this.tx.span), utr, by = 'transcript_id', max.chunk = max.chunk)$query.id] = 'utr'
  this.tx.span$left.feat[gr.findoverlaps(gr.start(this.tx.span), promoters, by = 'transcript_id', max.chunk = max.chunk)$query.id] = 'promoter'

  this.tx.span$right.feat = NA
  this.tx.span$right.feat[gr.findoverlaps(gr.end(this.tx.span), tx.span, by = 'transcript_id', max.chunk = max.chunk)$query.id] = 'cds'
  this.tx.span$right.feat[gr.findoverlaps(gr.end(this.tx.span), utr, by = 'transcript_id', max.chunk = max.chunk)$query.id] = 'utr'
  this.tx.span$right.feat[gr.findoverlaps(gr.end(this.tx.span), promoters, by = 'transcript_id', max.chunk = max.chunk)$query.id] = 'promoter'

  ## if lands in CDS annotate this.tx.span ends with right and/or left exon frame
  ## (if lands in intron then annotate left end with frame of next exon on right and right end
  ## with frame of next exon on left)
  ## otherwise annotate as NA
  ## we will eventually integrate frames across walks and call a transition "in frame" if the
  ## frame of the right (left) end of the previous + (-) interval

  ## now we want to find the first exon to the right of the left boundary and
  ## the first exon to the left of the right boundary for each fragment
  ##tix = match(this.tx.span$transcript_id, tx.span$transcript_id)

  ## cds.u = grl.unlist(cds[this.tx.span$tx.id])
  ## cds.u = grl.unlist(cds[this.tx.span$cds.id])
  cds.u = grl.unlist(cds)
  cds.u$transcript_id = names(cds)[cds.u$grl.ix]
  cds.u$start.local = gr2dt(cds.u)[, id := 1:length(cds.u)][, tmp.st := 1+c(0, cumsum(width)[-length(width)]), by = transcript_id][, keyby = id][, tmp.st] 
  cds.u$end.local = cds.u$start.local + width(cds.u) -1


                                        #        ranges(cds.u) =  ranges(pintersect(cds.u, tx.span[this.tx.span$tx.id[cds.u$grl.ix]], resolve.empty = 'start.x'))
  ## ranges(cds.u) =  ranges(pintersect(cds.u, tx.span[this.tx.span$tx.id[cds.u$grl.ix]]))
  ## ranges(cds.u) =  ranges(pintersect(cds.u, tx.span[this.tx.span$cds.id[cds.u$grl.ix]]))
  ranges(cds.u) =  ranges(pintersect(cds.u, tx.span[cds.u$grl.ix]))

  tmp = gr.findoverlaps(this.tx.span, cds.u, scol = c('start.local', 'end.local', 'exon_number'), by = 'transcript_id', verbose = verbose, max.chunk = max.chunk)

  leftmost.cds.exon = data.table(id = 1:length(tmp), qid = tmp$query.id, start = start(tmp))[, id[which.min(start)], by = qid][, V1]
  rightmost.cds.exon = data.table(id = 1:length(tmp), qid = tmp$query.id, end = end(tmp))[, id[which.max(end)], by = qid][, V1]

  ## now we want to get the frame of the left and right base
  ## of each leftmost and rightmost exon (etc.
  ## this will depend on orientation of the exon and side that we are querying
  ## for left side of - exon, (exonFrame + width) %% 3
  ## for right side of - exon, (exonFrame + width(og.exon) - width %% 3
  ## for left side of + exon, (exonFrame + width(og.exon) - width) %%3
  ## for right side of - exon, (exonFrame + int.exon) %% 3

  leftmost.coord = ifelse(as.logical(strand(cds.u[tmp$subject.id[leftmost.cds.exon]])=='+'),
  (start(tmp)[leftmost.cds.exon] - start(cds.u)[tmp$subject.id[leftmost.cds.exon]] + cds.u$start.local[tmp$subject.id[leftmost.cds.exon]]),
  (end(cds.u)[tmp$subject.id[leftmost.cds.exon]] - start(tmp)[leftmost.cds.exon] + cds.u$start.local[tmp$subject.id[leftmost.cds.exon]]))

  rightmost.coord = ifelse(as.logical(strand(cds.u[tmp$subject.id[rightmost.cds.exon]])=='+'),
  (end(tmp)[rightmost.cds.exon] - start(cds.u)[tmp$subject.id[rightmost.cds.exon]] + cds.u$start.local[tmp$subject.id[rightmost.cds.exon]]),
  (end(cds.u)[tmp$subject.id[rightmost.cds.exon]] - end(tmp)[rightmost.cds.exon] + cds.u$start.local[tmp$subject.id[rightmost.cds.exon]]))

  leftmost.frame = ifelse(as.logical(strand(cds.u[tmp$subject.id[leftmost.cds.exon]])=='+'),
  (start(tmp)[leftmost.cds.exon] - start(cds.u)[tmp$subject.id[leftmost.cds.exon]] + cds.u$phase[tmp$subject.id[leftmost.cds.exon]]) %% 3,
  (end(cds.u)[tmp$subject.id[leftmost.cds.exon]] - start(tmp)[leftmost.cds.exon] + cds.u$phase[tmp$subject.id[leftmost.cds.exon]]) %% 3)

  rightmost.frame = ifelse(as.logical(strand(cds.u[tmp$subject.id[rightmost.cds.exon]])=='+'),
  (end(tmp)[rightmost.cds.exon] - start(cds.u)[tmp$subject.id[rightmost.cds.exon]] + cds.u$phase[tmp$subject.id[rightmost.cds.exon]]) %% 3,
  (end(cds.u)[tmp$subject.id[rightmost.cds.exon]] - end(tmp)[rightmost.cds.exon] + cds.u$phase[tmp$subject.id[rightmost.cds.exon]]) %% 3)

  leftmost.frame = leftmost.coord %% 3
  rightmost.frame = rightmost.coord %% 3

  this.tx.span$left.coord = this.tx.span$left.boundary = this.tx.span$right.coord = this.tx.span$right.boundary = NA

  this.tx.span$left.coord[tmp$query.id[leftmost.cds.exon]] = leftmost.coord
  this.tx.span$right.coord[tmp$query.id[rightmost.cds.exon]] = rightmost.coord

  this.tx.span$left.boundary[tmp$query.id[leftmost.cds.exon]] =
    ifelse(strand(this.tx.span)[tmp$query.id[leftmost.cds.exon]]=="+",
           tmp$start.local[leftmost.cds.exon], tmp$end.local[leftmost.cds.exon])

  this.tx.span$right.boundary[tmp$query.id[rightmost.cds.exon]] =
    ifelse(strand(this.tx.span)[tmp$query.id[rightmost.cds.exon]]=="+",
           tmp$end.local[rightmost.cds.exon], tmp$start.local[rightmost.cds.exon])

  this.tx.span$right.exon_del = this.tx.span$left.exon_del = NA;
  this.tx.span$right.frame = this.tx.span$left.frame = NA;
  this.tx.span$right.exon_id= this.tx.span$left.exon_id = NA;

  ## keep track of exon frames to left and right
  this.tx.span$left.frame[tmp$query.id[leftmost.cds.exon]] = leftmost.frame
  this.tx.span$right.frame[tmp$query.id[rightmost.cds.exon]] = rightmost.frame

  this.tx.span$left.exon_id[tmp$query.id[leftmost.cds.exon]] = cds.u[tmp$subject.id[leftmost.cds.exon]]$exon_number
  this.tx.span$right.exon_id[tmp$query.id[rightmost.cds.exon]] = cds.u[tmp$subject.id[rightmost.cds.exon]]$exon_number

  ## keep track of which exons have any sort of sequence deletion
  this.tx.span$left.exon_del[tmp$query.id[leftmost.cds.exon]] =
    start(this.tx.span)[tmp$query.id[leftmost.cds.exon]] > start(cds.u)[tmp$subject.id[leftmost.cds.exon]]
  this.tx.span$right.exon_del[tmp$query.id[rightmost.cds.exon]] =
    end(this.tx.span)[tmp$query.id[rightmost.cds.exon]] <  end(cds.u)[tmp$subject.id[rightmost.cds.exon]]


  ## now traverse each walk and connect "broken transcripts"
  ## each walk is a list of windows, connections will be made with respect to the window walk
  ## paying attention to (1) side of breakage, (2) orientation of adjacent windows in walk (3) orientation of transcript
  ##
  ## right broken + cds upstream of ++ junction attaches to left  broken + cds in next window
  ##              - cds upstream of ++ junction attached to left  broken - cds in next window
  ##              + cds upstream of +- junction attaches to right broken - cds in next window
  ##              - cds upstream of +- junction attaches to right broken + cds in next window
  ##
  ## left  broken + cds upstream of -- junction attaches to right broken + cds in next window
  ##              - cds upstream of -- junction attached to right broken - cds in next window
  ##              + cds upstream of -+ junction attaches to left  broken - cds in next window
  ##              - cds upstream of -+ junction attaches to left  broken + cds in next window
  ##
  ##
  ##
  
  ## to achieve this create a graphs on elements of this.tx.span
  ## connecting elements i and j if a window pair k k+1 in (the corresponding) input grl
  ## produces a connection with the correct orientation
  
  ## match up on the basis of the above factors to determine edges in graph
  
  dt_a = data.table(
    i = 1:length(this.tx.span),
    key1 = walks.u$grl.ix[this.tx.span$subject.id],
    key2 = walks.u$grl.iix[this.tx.span$subject.id],
    key3 = this.tx.span$cds.sign*this.tx.span$window.sign,
    key4 = ifelse(this.tx.span$window.sign>0, this.tx.span$right.broken, this.tx.span$left.broken))


  dt_b = data.table(
    j = 1:length(this.tx.span),
    key1 = walks.u$grl.ix[this.tx.span$subject.id],
    key2 = walks.u$grl.iix[this.tx.span$subject.id]-1,
    ## key2 = walks.u$grl.iix[this.tx.span$subject.id]+1,
    key3 = this.tx.span$cds.sign*this.tx.span$window.sign,
    key4 = ifelse(this.tx.span$window.sign>0, this.tx.span$left.broken, this.tx.span$right.broken))

  extra_walks_ids = setdiff(seq_along(walks.u), this.tx.span$subject.id) ## probably more robust alternative to negative indexing as below
  
  if (length(extra_walks_ids) > 0) {
    dt_c = data.table(
      ## key1 = walks.u$grl.ix[-this.tx.span$subject.id],
      ## key2 = walks.u$grl.iix[-this.tx.span$subject.id],
      key1 = walks.u$grl.ix[extra_walks_ids],
      key2 = walks.u$grl.iix[extra_walks_ids],
      key3 = 1,
      key4 = TRUE)
    dt_c[, i := seq(length(this.tx.span) + 1, length(this.tx.span) + dt_c[, .N])]
    dt_c = rbind(dt_c, copy(dt_c)[, key3 := -1])


    dt_d = data.table(
      ## key1 = walks.u$grl.ix[-this.tx.span$subject.id],
      ## key2 = walks.u$grl.iix[-this.tx.span$subject.id]-1,
      key1 = walks.u$grl.ix[extra_walks_ids],
      key2 = walks.u$grl.iix[extra_walks_ids]-1,
      key3 = 1,
      key4 = TRUE)
    dt_d[, j := seq(length(this.tx.span) + 1, length(this.tx.span) + dt_d[, .N])]
    dt_d = rbind(dt_d, copy(dt_d)[, key3 := -1])
    
  } else {
    dt_c = NULL
    dt_d = NULL
  }

  dt_1 = rbind(dt_a, dt_c)
  dt_2 = rbind(dt_b, dt_d)

  edges = merge(dt_1, dt_2, by = c('key1', 'key2', 'key3', 'key4'), allow.cartesian = TRUE)[key4 == TRUE,]  

  ## tmp_tx_span = grbind(this.tx.span, walks.u[-this.tx.span$subject.id][,c()])
  tmp_tx_span = grbind(this.tx.span, walks.u[extra_walks_ids][,c()])
  ## extra_ids = seq(length(this.tx.span) + 1, length(this.tx.span) + length(walks.u[-this.tx.span$subject.id][,c()]))
  extra_ids = seq(length(this.tx.span) + 1, by = 1, length.out = length(extra_walks_ids))
  if (length(extra_ids) > 0) {
    mcols(tmp_tx_span)[extra_ids,][["transcript_id"]] = NA ## these are redundant lines, but for insurance... keep
    mcols(tmp_tx_span)[extra_ids,][["gene_name"]] = NA
  }
  
  ## remove edges that link different transcripts of same gene
  if (filter.splice) {
    na2true = function(v) {
      v[is.na(v)] = TRUE
      as.logical(v)
    }
    ## edges = edges[!(this.tx.span$gene_name[edges$i] == this.tx.span$gene_name[edges$j] & this.tx.span$transcript_id[edges$i] != this.tx.span$transcript_id[edges$j]), ]
    edges = edges[na2true(!(tmp_tx_span$gene_name[edges$i] == tmp_tx_span$gene_name[edges$j] & tmp_tx_span$transcript_id[edges$i] != tmp_tx_span$transcript_id[edges$j])), ]
    ## NA will be any gene_name or transcript ids which are NA, i.e. segments outside of coding region
  }

  if (nrow(edges)==0)
    return(GRangesList())

  ## dim_to_rep =  length(this.tx.span) + length(walks.u[-this.tx.span$subject.id])
  dim_to_rep = length(tmp_tx_span)
  A = sparseMatrix(edges$i, edges$j, x = 1, dims = rep(dim_to_rep,2))
  sources = Matrix::which(Matrix::colSums(A!=0)==0)
  sinks = Matrix::which(Matrix::rowSums(A!=0)==0)

  G = graph.adjacency(A)
  C = igraph::clusters(G, 'weak')
  vL = split(1:nrow(A), C$membership)
  vL = vL[elementNROWS(vL) > 1]

  ## collate all paths through this graph
  paths = do.call('c', mclapply(1:length(vL), function(i) {
    if (verbose & (i %% 10)==0)
      message(i, ' of ', length(vL))
    x = vL[[i]]
    tmp.source = setdiff(match(sources, x), NA)
    tmp.sink = setdiff(match(sinks, x), NA)
    tmp.mat = A[x, x, drop = FALSE]!=0
    if (length(x)<=1)
      return(NULL)
    if (length(x)==2)
      list(x[c(tmp.source, tmp.sink)])
    else if (all(Matrix::rowSums(tmp.mat)<=1) & all(Matrix::colSums(tmp.mat)<=1))
      get.shortest.paths(G, from = intersect(x, sources), intersect(x, sinks))$vpath
    else
    {
      if (exhaustive)
        lapply(JaBbA:::all.paths(A[x,x, drop = FALSE], source.vertices = tmp.source, sink.vertices = tmp.sink, verbose = FALSE)$paths, function(y) x[y])
      else
      {
        out = do.call('c', lapply(intersect(x, sources),
                                  function(x, sinks) suppressWarnings(get.shortest.paths(G, from = x, to = sinks)$vpath), sinks = intersect(x, sinks)))
        out = out[sapply(out, length)!=0]
        if (length(out)>0)
          out = out[!duplicated(sapply(out, paste, collapse = ','))]
        return(out)
      }
    }
  }, mc.cores = mc.cores))

  ## fus.sign = this.tx.span$cds.sign * this.tx.span$window.sign
  ## paths = lapply(paths, function(x) if (fus.sign[x][1]<0) rev(x) else x) ## reverse "backward paths" (i.e. those producing fusions in backward order)
  ## paths.first = sapply(paths, function(x) x[1])
  ## paths.last = sapply(paths, function(x) x[length(x)])
  ## paths.broken.start = ifelse(as.logical(strand(this.tx.span)[paths.first] == '+'), this.tx.span$left.broken[paths.first], this.tx.span$right.broken[paths.first])
  ## paths.broken.end = ifelse(as.logical(strand(this.tx.span)[paths.last] == '+'), this.tx.span$right.broken[paths.last], this.tx.span$left.broken[paths.last])
  paths.u = unlist(paths)
  paths.i = unlist(lapply(1:length(paths), function(x) rep(x, length(paths[[x]]))))
  ## tmp.gr = this.tx.span[paths.u]
  tmp.gr = tmp_tx_span[paths.u]


  ## annotate steps of walk with out of frame vs in frame (if cds is a grl)

  ## left and right exon frame
  paths.u.lec = tmp.gr$left.coord
  paths.u.rec = tmp.gr$right.coord
  paths.u.lef = tmp.gr$left.frame
  paths.u.ref = tmp.gr$right.frame
  paths.u.str = as.character(strand(tmp.gr))
  paths.u.lcds = tmp.gr$left.feat == 'cds'
  paths.u.rcds = tmp.gr$right.feat == 'cds'
  paths.u.lout = tmp.gr$left.feat %in% c('promoter', 'utr')
  paths.u.rout = tmp.gr$right.feat %in% c('promoter', 'utr')

  ## a fragment is in frame either if (1) it begins a walk at frame 0 outside of the cds
  ## or if (2) its frame is concordant with previous cds
  paths.u.inframe = rep(NA, length(paths.u))
  paths.u.cdsend = paths.u.cdsstart = rep(FALSE, length(paths.u)) # this keeps track of cds starts and cds ends in the fusion

  outside = TRUE
  for (i in 1:length(paths.u.inframe))
  {
    if (i == 1)
      outside = TRUE
    else if (paths.i[i] != paths.i[i-1] | (paths.u.str[i] == '+' & paths.u.lout[i]) | (paths.u.str[i] == '-' & paths.u.rout[i]))
      outside = TRUE

    if (outside)
    {
      if (paths.u.str[i] == '+')
        paths.u.inframe[i] = paths.u.lout[i] & paths.u.lef[i] == 0
      else
        paths.u.inframe[i] = paths.u.rout[i] & paths.u.ref[i] == 0

      paths.u.cdsstart[i] = paths.u.inframe[i]
      outside = F
    }
    else
    {
      if (paths.u.str[i] == '+' & paths.u.str[i-1] == '+')
        paths.u.inframe[i] = paths.u.lec[i] != 1 & paths.u.lef[i] == ((paths.u.ref[i-1]+1) %% 3) & paths.u.lcds[i] & paths.u.rcds[i-1]
      else if (paths.u.str[i] == '+' & paths.u.str[i-1] == '-')
        paths.u.inframe[i] = paths.u.lec[i] != 1 & paths.u.ref[i]  == ((paths.u.ref[i-1]+1) %% 3) & paths.u.rcds[i] & paths.u.rcds[i-1]
      else if (paths.u.str[i] == '-' & paths.u.str[i-1] == '-')
        paths.u.inframe[i] = paths.u.rec[i] != 1 & paths.u.ref[i] == ((paths.u.lef[i-1]+1) %% 3) & paths.u.rcds[i] & paths.u.lcds[i-1]
      else if (paths.u.str[i] == '-' & paths.u.str[i-1] == '+')
        paths.u.inframe[i] = paths.u.rec[i] != 1 & paths.u.lef[i] == ((paths.u.lef[i-1]+1) %% 3) & paths.u.lcds[i] & paths.u.lcds[i-1]
    }

    if ((paths.u.str[i] == '+' & paths.u.rout[i]) | (paths.u.str[i] == '-' & paths.u.lout[i]))
    {
      paths.u.cdsend[i] = paths.u.inframe[i]
      outside = T
    }
  }

  tmp.gr$in.frame = paths.u.inframe;
  tmp.gr$cds.start = paths.u.cdsstart
  tmp.gr$cds.end = paths.u.cdsend
  tmp.gr$del5 = ifelse(paths.u.str == '+', tmp.gr$left.exon_del, tmp.gr$right.exon_del)
  tmp.gr$del3 = ifelse(paths.u.str == '+', tmp.gr$right.exon_del, tmp.gr$left.exon_del)
  tmp.gr$first.coord = ifelse(paths.u.str == '+', tmp.gr$left.coord, tmp.gr$right.coord)
  tmp.gr$last.coord  = ifelse(paths.u.str == '+', tmp.gr$right.coord, tmp.gr$left.coord)
  tmp.gr$first.boundary = ifelse(paths.u.str == '+', tmp.gr$left.boundary, tmp.gr$right.boundary)
  tmp.gr$last.boundary  = ifelse(paths.u.str == '+', tmp.gr$right.boundary, tmp.gr$left.boundary)
  tmp.gr$first.exon = ifelse(paths.u.str == '+', tmp.gr$left.exon_id, tmp.gr$right.exon_id)
  tmp.gr$last.exon = ifelse(paths.u.str == '+', tmp.gr$right.exon_id, tmp.gr$left.exon_id)
  fusions = split(tmp.gr, paths.i)

                                        # now annotate fusions

  values(fusions)$walk.id = data.table(wid = walks.u$grl.ix[tmp.gr$subject.id], fid = paths.i)[, wid[1], keyby = fid][, V1]
                                        #    values(fusions)$walk.id = vaggregate(walks.u$grl.ix[tmp.gr$subject.id], by = list(paths.i), FUN = function(x) x[1])

  tmp.g = tmp.gr$gene_name
  tmp.cds = tmp.gr$transcript_id
  tmp.fe = as.numeric(tmp.gr$first.exon)
  tmp.le = as.numeric(tmp.gr$last.exon)
  tmp.fb = tmp.gr$first.boundary
  tmp.lb = tmp.gr$last.boundary
  tmp.fc = tmp.gr$first.coord
  tmp.lc = tmp.gr$last.coord
  tmp.5d = tmp.gr$del5
  tmp.3d = tmp.gr$del3
  tmp.5d[is.na(tmp.5d)] = FALSE
  tmp.3d[is.na(tmp.3d)] = FALSE
  paths.u.cdsend[is.na(paths.u.cdsend)] = FALSE
  paths.u.cdsstart[is.na(paths.u.cdsstart)] = FALSE
  paths.u.inframe[is.na(paths.u.inframe)] = FALSE

  totpaths = max(paths.i)

  if (verbose)
    message('Populating coordinates')

  values(fusions)[, 'coords'] = mcmapply(function(x) paste(unique(x), collapse = '; '),
                                         split(gr.string(tmp_tx_span[paths.u], mb = TRUE, round = 1), paths.i), mc.cores = mc.cores)

  if (verbose)
    message('Populating transcript names')
  values(fusions)[, 'transcript_names'] = mcmapply(function(x, y) paste(x, ' (', y, ')', sep = '', collapse = '; '),
                                                   split(values(tx.span)[, 'gene_name'][this.tx.span$query.id[paths.u]], paths.i),
                                                   split(values(tx.span)[, 'transcript_name'][this.tx.span$query.id[paths.u]], paths.i), mc.cores = mc.cores)

  if (verbose)
    message('Populating transcript ids')
  values(fusions)[, 'transcript_ids'] = mcmapply(function(x, y) paste(x, ' (', y, ')', sep = '', collapse = '; '),
                                                 split(values(tx.span)[, 'gene_name'][this.tx.span$query.id[paths.u]], paths.i),
                                                 split(values(tx.span)[, 'transcript_id'][this.tx.span$query.id[paths.u]], paths.i), mc.cores = mc.cores)

  if (verbose)
    message('Populating gene names')
  values(fusions)[, 'genes'] = mcmapply(function(x) paste(unique(x), collapse = '; '),
                                        split(values(tx.span)[, 'gene_name'][this.tx.span$query.id[paths.u]], paths.i), mc.cores = mc.cores)

  if (verbose)
    message('Populating alteration')
  values(fusions)$alteration = vaggregate(1:length(paths.i), by = list(paths.i),
                                          FUN = function(x)
                                          {
                                            if (verbose & (x[1] %% 10)==0)
                                              message('Path ', unique(paths.i[x]), ' of ', totpaths)
                                            if (length(unique((tmp.cds[x])))==1) ## single transcript event
                                            {
                                              out = NULL
                                              x = x[!is.na(tmp.fe[x]) & !is.na(tmp.le[x])]
                                              if (length(x)>0)
                                              {
                                        #                                                   browser()
                                                ir = IRanges(pmin(tmp.le[x], tmp.fe[x]), pmax(tmp.fe[x], tmp.le[x]))
                                                if (length(del <- setdiff(IRanges(min(tmp.fe[x]), max(tmp.le[x])), ir))>0)
                                                {
                                                  del.fc = pmax(tmp.lc[x[match(start(del)-1, tmp.le[x])]]+1, 1, na.rm = TRUE)
                                                  del.lc = pmin(tmp.fc[x[match(end(del)+1, tmp.fe[x])]]-1, max(tmp.lc[x]), na.rm = TRUE)
                                                  out = c(out,## some portion deleted
                                                          ifelse(start(del)==end(del),
                                                                 paste('deletion of exon ', start(del),
                                                                       ' [', del.fc, '-', del.lc, 'bp]',
                                                                       sep = '', collapse = ', '),
                                                                 paste('deletion of exons ', start(del), '-', end(del),
                                                                       ' [', del.fc, '-', del.lc, 'bp]',
                                                                       sep = '', collapse = ', ')))
                                                }

                                                if (length(amp <- IRanges(coverage(ir)>1))>0)
                                                {
                                                  amp.fc = tmp.lc[x[match(start(amp), tmp.le[x])]]
                                                  amp.lc = tmp.fc[x[match(end(amp), tmp.fe[x])]]
                                                  out = c(out,   ## some portion duplicated
                                                          ifelse(start(amp)==end(amp),
                                                                 paste('duplication of exon ', end(amp),
                                                                       '[', amp.fc, '-', amp.lc, 'bp]',
                                                                       sep = '', collapse = ', '),
                                                                 paste('duplication of exons ', start(amp), '-', end(amp),
                                                                       ' [', amp.fc, '-', amp.lc, 'bp]',
                                                                       sep = '', collapse = ', ')))
                                                }

                                                if (any(ix <- tmp.5d[x]))
                                                {
                                                  out = c(out, paste("partial 5' deletion of exon ", tmp.fe[x[ix]],
                                                                     ' [', tmp.fb[x[ix]], '-', tmp.fc[x[ix]], 'bp]',
                                                                     sep = '', collapse = ', '))  ## some portion duplicated
                                                }
                                                if (any(ix <- tmp.3d[x]))
                                                {
                                                  del.fc = pmax(tmp.lc[x[ix-1]] + 1, 1, na.rm = TRUE)
                                                  del.lc = pmin(tmp.fc[x[ix+1]]-1, max(tmp.lc[x]), na.rm = TRUE)
                                                  out = c(out, paste("partial 3' deletion of exon ", tmp.fe[x[ix]],
                                                                     ' [', tmp.lc[x[ix]], '-', tmp.lb[x[ix]], 'bp]',
                                                                     sep = '', collapse = ', '))  ## some portion duplicated
                                                }
                                              }

                                              if (length(out)>0)
                                                paste(out, collapse = '; ')
                                              else
                                                ''
                                            }
                                            else
                                            {
                                              return(paste(tmp.g[x], ' ', ifelse(paths.u.cdsstart[x], 'S', ''), ifelse(tmp.5d[x],  'tr', ''),
                                                           ifelse(is.na(tmp.fe[x]), 'UTR',
                                                           ifelse(tmp.le[x]==tmp.fe[x],
                                                                  paste('exon ', tmp.le[x], sep = ''),
                                                                  paste('exons ', tmp.fe[x], '-', tmp.le[x], sep = ''))),
                                                           ' [', tmp.fc[x], '-', tmp.lc[x], 'bp]',
                                                           ifelse(tmp.3d[x],  'tr', ''), ifelse(paths.u.cdsend[x], 'E', ''), ' ',
                                                           ifelse(c(paths.u.inframe[x[-1]], FALSE), '-',
                                                           ifelse((1:length(x))!=length(x), '-X', '')), sep = '', collapse = '-> '))
                                            }
                                          })

  values(fusions)$max.inframe = vaggregate(paths.u.inframe, by = list(paths.i),
                                           FUN = function(x) return(max(c(0, rle(x)$lengths[which(rle(x)$values == T)]))))
  values(fusions)$num.win = vaggregate(paths.u.inframe, by = list(paths.i), length)
  values(fusions) = cbind(values(walks)[values(fusions)$walk.id, , drop = FALSE], values(fusions))

  fusions = fusions[nchar(values(fusions)$alteration)>0, ]
  fusions = fusions[!gUtils::grl.eval(fusions, length(na.omit(gene_name)) == 1)]
  return(fusions)
}


#' @name gt.gencode
#' @description
#'
#' internal function to format transcript annotations for fusion output
#'
#' @param gencode GRanges output of rtracklayer::import of GENCODE gtf
#' @param bg.col character representing color to put in colormap
#' @param cds.col color for CDS
#' @param utr.col color for UTR
#' @param st.col scalar character representing color of CDS start
#' @param en.col scalar character representing color of CDS end
#' @keywords internal
#' @noRd
#' @return 
gt.gencode = function(gencode, bg.col = alpha('blue', 0.1), cds.col = alpha('blue', 0.6), utr.col = alpha('purple', 0.4), st.col = 'green',
  en.col = 'red')  
{
  tx = gencode[gencode$type =='transcript']
  genes = gencode[gencode$type =='gene']
  exons = gencode[gencode$type == 'exon']
  utr = gencode[gencode$type == 'UTR']
  ## ut = unlist(utr$tag)
  ## utix = rep(1:length(utr), sapply(utr$tag, length))
  ## utr5 = utr[unique(utix[grep('5_UTR',ut)])]
  ## utr3 = utr[unique(utix[grep('3_UTR',ut)])]
  ## utr5$type = 'UTR5'
  ## utr3$type = 'UTR3'
  startcodon = gencode[gencode$type == 'start_codon']
  stopcodon = gencode[gencode$type == 'stop_codon']
  OUT.COLS = c('gene_name', 'transcript_name', 'transcript_id', 'type', 'exon_number', 'type')
  tmp = c(genes, tx, exons, utr, startcodon, stopcodon)[, OUT.COLS]
  
  ## compute tx ord of intervals
  ord.ix = order(tmp$transcript_id, match(tmp$type, c('gene', 'transcript', 'exon', 'UTR', 'start_codon','stop_codon')))
  tmp.rle = rle(tmp$transcript_id[ord.ix])
  tmp$tx.ord[ord.ix] = unlist(lapply(tmp.rle$lengths, function(x) 1:x))
  tmp = tmp[rev(order(match(tmp$type, c('gene', 'transcript', 'exon', 'UTR', 'start_codon','stop_codon'))))] 
  tmp.g = tmp[tmp$type != 'transcript']
  cmap = list(type = c(gene = bg.col, transcript = bg.col, exon = cds.col, start_codon = st.col, stop_codon = en.col, UTR = utr.col))
  tmp.g = gr.disjoin(gr.stripstrand(tmp.g))
  return(gTrack(tmp.g[, c('type', 'gene_name')], colormap = cmap))
}
