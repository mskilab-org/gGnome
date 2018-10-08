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
                   verbose = FALSE){
    ## QC input graph or junctions

  if (!inherits(graph, "gGraph")){
    stop('Input must be be gGraph')
  }


  if (is.character(gencode) && file.url.exists(gencode))
  {
    gencode = rtracklayer::import(gencode)
  }
  
  gencode = skidb::read_gencode()
  
  tgg = make_txgraph(graph, gencode)
  
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
  if (inherits(walks, 'GRanges')){
    walks = GRangesList(walks)
  }

  if (is(walks, 'list')){
    walks = do.call(GRangesList, walks)
  }

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
  tx.span = seg2gr(
    as.data.table(rrbind(as.data.frame(tx.span), as.data.frame(promoters)))[, list(seqnames = seqnames[1], start = min(start), end = max(end), strand = strand[1], gene_name = gene_name[1], transcript_id = transcript_id[1], transcript_name = transcript_name[1], cds.id = cds.id), keyby = cds.id][!is.na(cds.id), ], seqlengths = seqlengths(tx.span)) ## promoters get thrown out?? what is the point of the is.na(cds.id) filter when promoters will not have cds.id??

  ## create stranded intergenic regions from in between transcripts
  genomep = genomen = si2gr(tx.span)
  strand(genomep) = '+'
  strand(genomen) = '-'
  igr = setdiff(grbind(genomep, genomen), tx.span)
  tx.span = grbind(tx.span, igr)

                                        # match up tx.span to walks
  walks.u = grl.unlist(walks)

  ## these are fragments of transcripts that overlap walks
  ## grl.ix and grl.iix will keep track of original walk id 
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


  ##        ranges(cds.u) =  ranges(pintersect(cds.u, tx.span[this.tx.span$tx.id[cds.u$grl.ix]], resolve.empty = 'start.x'))
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

   edges = merge(dt_a, dt_b, by = c('key1', 'key2', 'key3', 'key4'), allow.cartesian = TRUE)[key4 == TRUE,]  

  ## mimielinski Saturday, Sep 01, 2018 09:58:47 AM
  ## commenting out because creates fake edges
  ## instead added intergenic regions to tx.span above so they
  ## just get merged like any other features above

  ## extra_walks_ids = setdiff(seq_along(walks.u), this.tx.span$subject.id) ## probably more robust alternative to negative indexing as below
  
  ## if (length(extra_walks_ids) > 0) {
  ##   dt_c = data.table(
  ##     ## key1 = walks.u$grl.ix[-this.tx.span$subject.id],
  ##     ## key2 = walks.u$grl.iix[-this.tx.span$subject.id],
  ##     key1 = walks.u$grl.ix[extra_walks_ids],
  ##     key2 = walks.u$grl.iix[extra_walks_ids],
  ##     key3 = 1,
  ##     key4 = TRUE)
  ##   dt_c[, i := seq(length(this.tx.span) + 1, length(this.tx.span) + dt_c[, .N])]
  ##   dt_c = rbind(dt_c, copy(dt_c)[, key3 := -1])


  ##   dt_d = data.table(
  ##     ## key1 = walks.u$grl.ix[-this.tx.span$subject.id],
  ##     ## key2 = walks.u$grl.iix[-this.tx.span$subject.id]-1,
  ##     key1 = walks.u$grl.ix[extra_walks_ids],
  ##     key2 = walks.u$grl.iix[extra_walks_ids]-1,
  ##     key3 = 1,
  ##     key4 = TRUE)
  ##   dt_d[, j := seq(length(this.tx.span) + 1, length(this.tx.span) + dt_d[, .N])]
  ##   dt_d = rbind(dt_d, copy(dt_d)[, key3 := -1])
    
  ## } else {
  ##   dt_c = NULL
  ##   dt_d = NULL
  ## }

  ## dt_1 = rbind(dt_a, dt_c)
  ## dt_2 = rbind(dt_b, dt_d)

 ## edges = merge(dt_1, dt_2, by = c('key1', 'key2', 'key3', 'key4'), allow.cartesian = TRUE)[key4 == TRUE,]  

  ## tmp_tx_span = grbind(this.tx.span, walks.u[-this.tx.span$subject.id][,c()])
  ## tmp_tx_span = grbind(this.tx.span, walks.u[extra_walks_ids][,c()])
  ## ## extra_ids = seq(length(this.tx.span) + 1, length(this.tx.span) + length(walks.u[-this.tx.span$subject.id][,c()]))
  ## extra_ids = seq(length(this.tx.span) + 1, by = 1, length.out = length(extra_walks_ids))
  ## if (length(extra_ids) > 0) {
  ##   mcols(tmp_tx_span)[extra_ids,][["transcript_id"]] = NA ## these are redundant lines, but for insurance... keep
  ##   mcols(tmp_tx_span)[extra_ids,][["gene_name"]] = NA
  ## }
  
  ## remove edges that link different transcripts of same gene
  if (filter.splice) {
    na2true = function(v) {
      v[is.na(v)] = TRUE
      as.logical(v)
    }
    ##edges = edges[!(this.tx.span$gene_name[edges$i] == this.tx.span$gene_name[edges$j] & this.tx.span$transcript_id[edges$i] != this.tx.span$transcript_id[edges$j]), ]
    edges = edges[na2true(!(this.tx.span$gene_name[edges$i] == this.tx.span$gene_name[edges$j] & this.tx.span$transcript_id[edges$i] != this.tx.span$transcript_id[edges$j])), ]
#    edges = edges[na2true(!(tmp_tx_span$gene_name[edges$i] == tmp_tx_span$gene_name[edges$j] & tmp_tx_span$transcript_id[edges$i] != tmp_tx_span$transcript_id[edges$j])), ]
    ## NA will be any gene_name or transcript ids which are NA, i.e. segments outside of coding region
  }

  if (nrow(edges)==0){
      return(GRangesList())
  }

  ## dim_to_rep =  length(this.tx.span) + length(walks.u[-this.tx.span$subject.id])
  ##  dim_to_rep = length(tmp_tx_span)
  dim_to_rep = length(this.tx.span)
  A = sparseMatrix(edges$i, edges$j, x = 1, dims = rep(dim_to_rep,2))
  sources = Matrix::which(Matrix::colSums(A!=0)==0)
  sinks = Matrix::which(Matrix::rowSums(A!=0)==0)

  G = graph.adjacency(A)
  C = igraph::clusters(G, 'weak')
  vL = split(1:nrow(A), C$membership)
  vL = vL[elementNROWS(vL) > 1]

  ## collate all paths through this graph
  paths = do.call('c', mclapply(1:length(vL), function(i) {
    if (verbose & (i %% 10)==0){
      message(i, ' of ', length(vL))
    }
    x = vL[[i]]
    tmp.source = setdiff(match(sources, x), NA)
    tmp.sink = setdiff(match(sinks, x), NA)
    tmp.mat = A[x, x, drop = FALSE]!=0
    if (length(x)<=1){
      return(NULL)
    }
    if (length(x)==2){
      list(x[c(tmp.source, tmp.sink)])
    }
    else if (all(Matrix::rowSums(tmp.mat)<=1) & all(Matrix::colSums(tmp.mat)<=1)){
      get.shortest.paths(G, from = intersect(x, sources), intersect(x, sinks))$vpath
    }
    else
    {
      if (exhaustive){
        lapply(all.paths(A[x,x, drop = FALSE], source.vertices = tmp.source, sink.vertices = tmp.sink, verbose = FALSE)$paths, function(y) x[y])}
      else
      {
        out = do.call('c', lapply(intersect(x, sources),
                                  function(x, sinks) suppressWarnings(get.shortest.paths(G, from = x, to = sinks)$vpath), sinks = intersect(x, sinks)))
        out = out[sapply(out, length)!=0]
        if (length(out)>0){
          out = out[!duplicated(sapply(out, paste, collapse = ','))]
        }
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
  tmp.gr = this.tx.span[paths.u]
##  tmp.gr = tmp_tx_span[paths.u]


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
    if (i == 1){
      outside = TRUE
    }
    else if (paths.i[i] != paths.i[i-1] | (paths.u.str[i] == '+' & paths.u.lout[i]) | (paths.u.str[i] == '-' & paths.u.rout[i])){
      outside = TRUE
    }
    if (outside)
    {
      if (paths.u.str[i] == '+'){
          paths.u.inframe[i] = paths.u.lout[i] & paths.u.lef[i] == 0
      }
      else{
          paths.u.inframe[i] = paths.u.rout[i] & paths.u.ref[i] == 0
      }
      paths.u.cdsstart[i] = paths.u.inframe[i]
      outside = F
    }
    else
    {
      if (paths.u.str[i] == '+' & paths.u.str[i-1] == '+'){
        paths.u.inframe[i] = paths.u.lec[i] != 1 & paths.u.lef[i] == ((paths.u.ref[i-1]+1) %% 3) & paths.u.lcds[i] & paths.u.rcds[i-1]}
      else if (paths.u.str[i] == '+' & paths.u.str[i-1] == '-'){
          paths.u.inframe[i] = paths.u.lec[i] != 1 & paths.u.ref[i]  == ((paths.u.ref[i-1]+1) %% 3) & paths.u.rcds[i] & paths.u.rcds[i-1]
          }
      else if (paths.u.str[i] == '-' & paths.u.str[i-1] == '-'){
          paths.u.inframe[i] = paths.u.rec[i] != 1 & paths.u.ref[i] == ((paths.u.lef[i-1]+1) %% 3) & paths.u.rcds[i] & paths.u.lcds[i-1]
          }
      else if (paths.u.str[i] == '-' & paths.u.str[i-1] == '+'){
          paths.u.inframe[i] = paths.u.rec[i] != 1 & paths.u.lef[i] == ((paths.u.lef[i-1]+1) %% 3) & paths.u.lcds[i] & paths.u.lcds[i-1]
          }
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

  ## now annotate fusions

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

  if (verbose){
    message('Populating coordinates')
    }

  values(fusions)[, 'coords'] = mcmapply(function(x) paste(unique(x), collapse = '; '),
                                         split(gr.string(this.tx.span[paths.u], mb = TRUE, round = 1), paths.i), mc.cores = mc.cores)

  if (verbose){
    message('Populating transcript names')
  }
  values(fusions)[, 'transcript_names'] = mcmapply(function(x, y) paste(x, ' (', y, ')', sep = '', collapse = '; '),
                                                   split(values(tx.span)[, 'gene_name'][this.tx.span$query.id[paths.u]], paths.i),
                                                   split(values(tx.span)[, 'transcript_name'][this.tx.span$query.id[paths.u]], paths.i), mc.cores = mc.cores)

  if (verbose){
    message('Populating transcript ids')}
  values(fusions)[, 'transcript_ids'] = mcmapply(function(x, y) paste(x, ' (', y, ')', sep = '', collapse = '; '),
                                                 split(values(tx.span)[, 'gene_name'][this.tx.span$query.id[paths.u]], paths.i),
                                                 split(values(tx.span)[, 'transcript_id'][this.tx.span$query.id[paths.u]], paths.i), mc.cores = mc.cores)

  if (verbose){
    message('Populating gene names')}
  values(fusions)[, 'genes'] = mcmapply(function(x) paste(unique(x), collapse = '; '),
                                        split(values(tx.span)[, 'gene_name'][this.tx.span$query.id[paths.u]], paths.i), mc.cores = mc.cores)

  if (verbose){
    message('Populating alteration')}
  values(fusions)$alteration = vaggregate(1:length(paths.i), by = list(paths.i),
                                          FUN = function(x)
                                          {
                                            if (verbose & (x[1] %% 10)==0){
                                              message('Path ', unique(paths.i[x]), ' of ', totpaths)}
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

                                              if (length(out)>0){
                                                paste(out, collapse = '; ')}
                                              else{
                                                  ''
                                                  }
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
                                        #  fusions = fusions[!gUtils::grl.eval(fusions, length(na.omit(gene_name)) == 1)]
  browser()
  return(fusions)
}


#' @name proximity
#' @export
#' @rdname internal
#' proximity
#'
#' Takes a set of n "query" elements (GRanges object, e.g. genes) and determines their proximity to m "subject" elements
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
#' @param ref gGraph of the "reference genome", by default is the reference genome but can be any gGraph
#' @param ignore.strand whether to ignore strand of input GRanges
#' @param verbose logical flag
#' @param mc.cores how many cores (default 1)
#' @param max.dist maximum genomic distance to store and compute (1MB by default) should the maximum distance at which biological interactions may occur
#' @return gWalk object each representing a proximity
proximity = function(gg, query, subject, ref = NULL, reduce = TRUE, ignore.strand = TRUE,
                     verbose = F, mc.cores = 1,
  max.dist = 1e6 ## max distance to store / compute in the output matrix.cores
  )
{
  if (!ignore.strand)
    stop('strand-aware proximity is TBD')

  if (!is.null(ref))
    stop('proximity with arbitrary reference is TBD')
  
  if (length(query)==0 | length(subject)==0)
    return(gW(graph = gg))

  if (is.null(names(query)))
    names(query) = 1:length(query)

  if (is.null(names(subject)))
    names(subject) = 1:length(subject)

  ra = gg$edges[type == 'ALT', ]$junctions$grl

  query.og = query
  subject.og = subject

  query.nm = names(query);
  subject.nm = names(subject);

  query = query
  subject = subject

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

  if (length(query)==0 | length(subject)==0)
    return(gW(graph = gg))

  query$type = 'query'
  subject$type = 'subject'

  gr = gr.fix(grbind(query, subject), gg)

  ## node.start and node.end delinate the nodes corresponding to the interval start and end
  ## on both positive and negative tiles of the karyograph
  gr$node.start = gr$node.end = gr$node.start.n = gr$node.end.n = NA;

  if (reduce) ## 
  {
    ## we do (super cheap) "reduce" on graph to only include ALT junctions and our breaks
    values(ra) = NULL
    px.gg = gG(breaks = gr[, c()], junc = ra)
  } else ## we do a simplify to reduce complexity of path search
  {
    px.gg = gg$copy$simplify(FUN = NULL)

    ## break px.gg into nodes defined by our new tiling
    ## without collapsing the input graphs
    px.gg$disjoin(gr = gr[,c()], collapse = FALSE)
  }
  
  values(ra) = NULL

  ## start and end indices of nodes
  tip = which(as.character(strand(px.gg$gr))=='+')
  tin = which(as.character(strand(px.gg$gr))=='-')

  ## to run proximity we
  ## just to identify indices of "node.start" and "node.end"
  ## for each query / subject gr
  ## i.e. the nodes in the graph that correspond to the
  ## ends of nodes in our query and subject
  ## we do this separately for negative and positive side

  ## lifting gr onto px.gg
  ovs = gr.start(gr,1)[, c()] %*% gr.start(px.gg$nodes$gr[, c('snode.id', 'index')])
  ove = gr.end(gr,1)[, c()] %*% gr.end(px.gg$nodes$gr[, c('snode.id', 'index')])
  ove$rindex = px.gg$queryLookup(ove$snode.id)$rindex
  ovs$rindex = px.gg$queryLookup(ovs$snode.id)$rindex

  ov = merge(gr2dt(ovs)[, .(node.start = index, node.end.n = rindex, query.id)],
             gr2dt(ove)[, .(node.end = index, node.start.n = rindex, query.id)],
             by = "query.id", allow.cartesian = TRUE)

  ## temp check: every input gr should be represented in every
  ## ovs and ove above ... 
  if (length(setdiff(1:length(gr), ovs$query.id))!=0 |
      length(setdiff(1:length(gr), ove$query.id))!=0)
    stop('error lifting query or subject nodes onto gGraph')

  ## sync gr to ovs ... i.e. replicate for
  ## every overlap
  ##
  ## every gr here now represents a lifted chunk of query or subject
  gr = gr[ov$query.id, ]

  ## we annotate the start and ends of the lifted intervals in "node coordinates"
  gr$node.start = ov$node.start
  gr$node.end = ov$node.end
  gr$node.start.n = ov$node.start.n
  gr$node.end.n = ov$node.end.n

  ## so for each query end we will find the shortest path to all subject starts
  ## and for each query start we will find the shortest.path from all subject ends
  ix.query = which(gr$type == 'query')
  ix.subj = which(gr$type == 'subject')
  
  ## so now we build distance matrices from query ends to subject starts
  ## and subject ends to query starts
  node.start = gr$node.start
  node.end = gr$node.end
  node.start.n = gr$node.start.n
  node.end.n = gr$node.end.n

  w = width(px.gg$gr)

  G = px.gg$igraph
  ed = px.gg$sedgesdt
  to = ed[.(igraph::E(G)$sedge.id), to]
  igraph::E(G)$weight = width(px.gg$gr)[to]

  ## ix.query and ix.subj give the indices of query / subject in gr
  ## node.start, node.end map gr to graph node ids
  ##
  ## these matrices are in dimensions of query and subject, and will hold the pairwise distances between
  ##
  D.rel = D.ra = D.ref = D.which = Matrix::Matrix(data = 0, nrow = length(ix.query), ncol = length(ix.subj))

  ## "REF" graph (missing ALT edges)

  G.ref = subgraph.edges(G, Matrix::which(igraph::E(G)$type == 'REF'), delete.vertices = F)

  EPS = 1e-9

  tmp = mclapply(ix.query, function(i)
  {
    if (verbose)
      cat('starting interval', i, 'of', length(ix.query), '\n')

    ## D1 = shortest query to subject path, D2 = shortest subject to query path, then take shortest of D1 and D2
    ## for each path, the edge weights correspond to the interval width of the target node, and to compute the path
    ## length we remove the final node since we are measuring the distance from the end of the first vertex in the path
    ## to the beginning of the final vertex

    u.node.start = unique(node.start[ix.subj]) ## gets around annoying igraph::shortest.path issue (no dups allowed)
    u.node.end = unique(node.end[ix.subj])

    uix.start = match(node.start[ix.subj], u.node.start)
    uix.end = match(node.end[ix.subj], u.node.end)

    tmp.D1 = (shortest.paths(G, node.end[i], u.node.start, weights = igraph::E(G)$weight, mode = 'out') - w[u.node.start])[uix.start]
    tmp.D2 = (shortest.paths(G, node.start[i], u.node.end, weights = igraph::E(G)$weight, mode = 'in') - w[node.start[i]])[uix.end]
    tmp.D3 = (shortest.paths(G, node.end.n[i], u.node.start, weights = igraph::E(G)$weight, mode = 'out') - w[u.node.start])[uix.start]
    tmp.D4 = (shortest.paths(G, node.start.n[i], u.node.end, weights = igraph::E(G)$weight, mode = 'in') - w[node.start.n[i]])[uix.end]
    tmp.D = pmin(tmp.D1, tmp.D2, tmp.D3, tmp.D4)
    ix = Matrix::which(tmp.D<max.dist)
    D.ra[i, ix] = tmp.D[ix]+EPS
    D.which[i, ix] = apply(cbind(tmp.D1[ix], tmp.D2[ix], tmp.D3[ix], tmp.D4[ix]), 1, which.min)

    u.node.start = unique(node.start[ix.subj][ix]) ## gets around annoying igraph::shortest.path issue (no dups allowed)
    u.node.end = unique(node.end[ix.subj][ix])

    uix.start = match(node.start[ix.subj][ix], u.node.start)
    uix.end = match(node.end[ix.subj][ix], u.node.end)

    tmp.D1 = (shortest.paths(G.ref, node.end[i], u.node.start, weights = E(G.ref)$weight, mode = 'out') - w[u.node.start])[uix.start]
    tmp.D2 = (shortest.paths(G.ref, node.start[i], u.node.end, weights = E(G.ref)$weight, mode = 'in') - w[node.start[i]])[uix.end]
    tmp.D3 = (shortest.paths(G.ref, node.end.n[i], u.node.start, weights = E(G.ref)$weight, mode = 'out') - w[u.node.start])[uix.start]
    tmp.D4 = (shortest.paths(G.ref, node.start.n[i], u.node.end, weights = E(G.ref)$weight, mode = 'in') - w[node.start.n[i]])[uix.end]
    tmp.D = pmin(tmp.D1, tmp.D2, tmp.D3, tmp.D4)
    D.ref[i, ix] = tmp.D+EPS

    ## if subject and query intersect (on the reference) then we count both RA and Ref distance as 0
    ## (easier to do a simple range query here)
    ix.zero = gr.in(gr[ix.subj[ix]], gr[ix.query[i]])
    if (any(ix.zero))
    {
      D.ra[i, ix[ix.zero]] = 0
      D.ref[i, ix[ix.zero]] = 0
    }
    D.rel[i, ix] = ((D.ra[i, ix]-EPS) / (D.ref[i, ix]-EPS)) + EPS

    if (verbose)
      cat('finishing interval', i, 'of', length(ix.query), ':', paste(round(D.rel[i, ix],2), collapse = ', '), '\n')

    return(list(D.rel = D.rel, D.ref = D.ref, D.ra = D.ra, D.which = D.which))
  }, mc.cores = mc.cores)

  for (i in 1:length(tmp))
  {
    if (class(tmp[[i]]) != 'list')
    {
      warning(sprintf('Query %s failed', ix.query[i]))
    }
    else
    {
      D.rel = D.rel + tmp[[i]]$D.rel
      D.ra = D.ra + tmp[[i]]$D.ra
      D.ref = D.ref + tmp[[i]]$D.ref
      D.which = D.which + tmp[[i]]$D.which
    }
  }

  ## sparse melt yum
  .spmelt = function(A) {
    ij = Matrix::which(A!=0, arr.ind = TRUE);
    dt = data.table(i = ij[,1], j = ij[,2], val = A[ij])
  }

   ## "full" size matrix
  ## rel = ra = ref = ra.which =
  ##   Matrix::Matrix(data = 0, nrow = length(qix.filt), ncol = length(six.filt), dimnames = list(dedup(query.nm), dedup(names(subject.nm))))

  
  Dt = .spmelt(D.rel)[, .(i, j, reldist = val)]
  Dt = merge(Dt, .spmelt(D.ra)[, .(i, j, altdist = val)], by = c("i", "j"))
  Dt = merge(Dt, .spmelt(D.ref)[, .(i, j, refdist = val)], by = c("i", "j"))
  Dt = merge(Dt, .spmelt(D.which)[, .(i, j, which = val)], by = c("i", "j"))
  Dt[, ":="( query.id = gr$id[ix.query[i]],
            subject.id = gr$id[ix.subj[j]])]
  Dt = Dt[reldist<1, ] ## only keep those with reldist<1
  setkey(Dt, altdist) ## sort by ALT distance
  Dt = unique(Dt, by = c('query.id', 'subject.id')) ## keep only unique (ie nearest) i j pairs

  setkey(Dt, reldist)

  if (nrow(Dt)==0)
    return(gW(graph = gg)) ## return empty gWalk

  paths = mcmapply(function(x, y, rw, i)
  {
    if (verbose)
      {
        message('path ', i, ' of ', nrow(Dt), '\n')
      }
    if (rw == 1)
      {
        out = get.shortest.paths(G, values(gr)[x, 'node.end'], values(gr)[y, 'node.start'], weights = igraph::E(G)$weight, mode = 'out')$vpath[[1]]
      }
    else if (rw == 2)
      {
        out = rev(get.shortest.paths(G, values(gr)[x, 'node.start'], values(gr)[y, 'node.end'], weights = igraph::E(G)$weight, mode = 'in')$vpath[[1]])
      }
    else if (rw == 3)
      {
        out = get.shortest.paths(G, values(gr)[x, 'node.end.n'], values(gr)[y, 'node.start'], weights = igraph::E(G)$weight, mode = 'out')$vpath[[1]]
      }
    else if (rw == 4)
      {
        out = rev(get.shortest.paths(G, values(gr)[x, 'node.start.n'], values(gr)[y, 'node.end'], weights = igraph::E(G)$weight, mode = 'in')$vpath[[1]])
      }
    return(px.gg$gr$snode.id[as.integer(out)])
  }, Dt$i, Dt$j, Dt$which, 1:nrow(Dt), SIMPLIFY = F, mc.cores = mc.cores)

  Dt = Dt[, .(reldist, altdist, refdist, query.id, subject.id)]

  if (ncol(values(query))>0)
  {
    Dt = cbind(Dt, as.data.table(values(query.og)[Dt$query.id, , drop = FALSE]))
  }
  
  if (ncol(values(subject))>0)
  {
    Dt = cbind(Dt, as.data.table(values(subject.og)[Dt$subject.id, , drop = FALSE]))
  }

  ## create gWalk object
  gw = gW(snode.id = paths, graph = px.gg, meta = Dt)

  ## merge gw back into original gGraph to retrieve metadata
  ## and validate that these are indeed valid walks
  gw$disjoin(graph = gg)

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
#' @keyword internal
#' @param tgg gGraph output of get_txgraph
#' @param mc.cores number of cores across which to parallelize path search component of algorithm
#' @param verbose whether to provide verbose output
#' @return gWalk of paths representing derivative transcripts harboring "loops"
#' @author Marcin Imielinski
#' @noRd
make_txgraph = function(gg, gencode)
  {
    tx = gencode %Q% (type == 'transcript')

    ## broken transcripts intersect at least one junction
    tx$in.break = tx %^% unlist(gg$junctions$grl)

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
    cds$fivep.frame = -cds$phase %% 3
    cds$threep.frame = (cds$fivep.frame + width(cds) -1 ) %% 3

    cds$left.frame = ifelse(strand(cds) == '+', cds$fivep.frame, cds$threep.frame)
    cds$right.frame = ifelse(strand(cds) == '+', cds$threep.frame, cds$fivep.frame)

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

    ## now relabel internal txnode.ann features of internal NA runswith the three prime features of their qid.last
    txnode.ann[!is.na(na.run) & !is.na(qid.last),
               ":="(
                 fivep.frame = txnode.ann$threep.frame[qid.last],
                 threep.frame = txnode.ann$threep.frame[qid.last],
                 fivep.exon = txnode.ann$threep.frame[qid.last],
                 threep.exon = txnode.ann$threep.frame[qid.last],
                 fivep.cc = txnode.ann$threep.cc[qid.last],
                 threep.cc = txnode.ann$threep.cc[qid.last],
                 is.start = FALSE,
                 is.end = FALSE,
                 fivep.pc = txnode.ann$threep.pc[qid.last],
                 threep.pc = txnode.ann$threep.pc[qid.last])]
    txnode.ann[, twidth := (threep.cc - fivep.cc + 1)/3]

    values(txnodes) = cbind(values(txnodes),
                            txnode.ann[, .(fivep.coord, threep.coord, fivep.frame,
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
    newedges = newedges[-which(deadend | deadstart | antisense), ] 

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
#' @noRd
get_txpaths = function(tgg,
                       genes = NULL, mc.cores = 1, verbose = FALSE)
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

          p = tgg$paths(ustart[k], combos[.(ustart[k]), end], weight = 'eweight', ignore.strand = FALSE)
          
          return(p)
        }, mc.cores = mc.cores)
        
        paths = do.call('c', c(p, list(force = TRUE)))[dist<INF & length>1, ]
        
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
#' @noRd
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
    dt.right = ref.p$nodesdt[snode.id %in% right, ][walk.iid != lengths(ref.p)[walk.id], ]
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
          loops.complex = loops.complex[!loops.complex$eval(any(tx_strand != strand, na.rm = TRUE))]
          
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
      
      loops[, snode.id := mapply("c", prefix, loop, suffix)]
      ab.l = gW(snode.id = loops$snode.id, graph = tgg)

      ## dedup any loops with identical node strings
      ab.l = ab.l[!duplicated(sapply(ab.l$snode.id, paste, collapse = ', '))]

      ab.l$set(numchr = ab.l$eval(node = length(unique(seqnames))))
      ab.l$set(numab = ab.l$eval(edge = sum(type == 'ALT')))
      ab.l$set(numgenes = ab.l$eval(node = length(unique(gene_name[!is.na(gene_name)]))))
      ab.l$set(genes = ab.l$eval(node = paste(unique(gene_name[!is.na(gene_name)]), collapse = ',')))
      ab.l$set(maxcn = ab.l$eval(edge = min(cn)))
    }
    return(ab.l)
  }

annotate_walks = function(walks)
{
  ## reinstantiate to "peel apart" walks .. i.e. reinstantiate with separate graph
  walks = gW(grl = walks$grl, meta = walks$meta, disjoin = FALSE)

  ## mark nodes with their walk.id
  walks$nodes$mark(wkid = walks$nodesdt$walk.id)

  ## annotate somatic coordinates and frames of fused transcripts
  ndt = gr2dt(walks$nodes$gr)
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
  cdt = ndt[!is.na(fivep.frame),
            .(
              in.frame = in.frame[1],
              exon.start = fivep.exon[1],
              exon.end = threep.exon[.N],
              pc.start = ceiling(fivep.pc[1]),
              pc.end = ceiling(threep.pc[.N])
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
      del = setdiff(IRanges(start[1], end[N]),
                    IRanges(start[tx == tx[1]],
                            end[tx == tx[1]]))
      if (length(del)>0)
      {
        ret = grl.string(split(GRanges(label, del), 1))
      }
    }
    return(ret)
  }

  .amp = function(tx, start, end, label)
  {
    N = length(tx)
    ret = as.character(NA)
    if (tx[1] == tx[N])
    {
      amp = as(coverage(
        GRanges(label,
                IRanges(start[tx == tx[1]],
                        end[tx == tx[1]])))>1, 'GRanges')
      amp = split(amp[amp$score], 1)
      if (length(amp)>0)
      {
        ret = grl.string(amp)
      }
    }
    return(ret)
  }

  adt = cdt[, .(
    gene.pc = paste0(gene_name, ':', pc.start, '-', pc.end, collapse = ';'),  
    del.pc = .del(transcript_id, pc.start, pc.end, gene_name),
    amp.pc = .amp(transcript_id, pc.start, pc.end, gene_name),
    in.frame = all(in.frame, na.rm = TRUE),
    qin.frame = in.frame[1] & in.frame[.N],
    tx.pc = paste0(transcript_id, ':', pc.start, '-', pc.end, collapse = ';'),
    gene.ec = paste0(gene_name, ':', exon.start, '-', exon.end, collapse = ';'),
    tx.ed = paste0(transcript_id, ':', exon.start, '-', exon.end, collapse = ';')
  ),
  keyby = wkid][.(1:length(walks)), ]

  adt[, frame.rescue := !in.frame & qin.frame]

  newmeta = cbind(walks$meta, adt[, .(gene.pc, amp.pc, del.pc, in.frame, frame.rescue, gene.ec, tx.pc)])

  ## now we want to trim the first and last intervals in the walks according to the
  ## fivep.coord and threep.coord
  ngr = walks$nodes$gr
  ngrdt = gr2dt(ngr)[, ":="(is.first = 1:.N %in% 1, is.last = 1:.N %in% .N), by = wkid]
  ngrdt[is.first==TRUE, start := ifelse(tx_strand == '+', fivep.coord, start)]
  ngrdt[is.first==TRUE, end := ifelse(tx_strand == '+', end, fivep.coord)]
  ngrdt[is.last==TRUE, start := ifelse(tx_strand == '+', start, threep.coord)]
  ngrdt[is.last==TRUE, end := ifelse(tx_strand == '+', threep.coord, end)]

  start(ngr) = ngrdt$start
  end(ngr) = ngrdt$end

  walks = gW(grl = split(ngr, ngr$wkid), disjoin = FALSE, meta = newmeta)

  return(walks)
}


