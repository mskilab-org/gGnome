##########
########## KH additional functions
##########

## #' @name eclusters2
## #' @title gGraph R6 public method eclusters2
## #' 
## #' @description
## #' Marks ALT edges belonging (quasi) reciprocal cycles
## #'
## #' 
## #' @param juncs GRangesList of junctions
## #' @param mc.cores parallel
## #' @param only_chains TRUE will only pair breakend to its nearest nearest neighbor IFF the nearest neighbor is reciprocal, see arguments to "strict" for 3 different matching heuristics
## #' @param max.small size below which simple dups and dels are excluded
## #' @param weak logical flag if TRUE will not differentiate between cycles and paths and will return all weakly connected clusters in the junction graph [FALSE]
## #' @param strict Only active if only_chains = TRUE. Can be one of "strict", "one_to_one", or "loose". \cr
## #' "strict": each breakend can only be "monogamously" matched to one other breakend and if the nearest breakend is of the wrong orientation, it is thrown out. \cr
## #' "one_to_one" breakends can only be coupled to a single "monogamous" match but without considering breakends of the wrong orientation. If breakend B is nearest to C, B will only be matched to C. If A's nearest breakend is B but is further away than C, A will not be matched to B. \cr
## #' "loose": the nearest breakend in the correct orientation under the threshold is considered. The same as one_to_one except A will be matched to B, while B will be matched to C. \cr
## #' @param ignore.isolated If TRUE, all simple duplications, duplications nested with only other duplications, and simple deletions without any breakends will be thrown out.
## #' @return gGraph object with edges marked with ecluster id and metadata in gEdge, gNode, and $meta
## #' @author Marcin Imielinski
## eclusters2 = function (thresh = 1000, weak = TRUE, paths = !weak,
##                            mc.cores = 1, verbose = FALSE, chunksize = 1e+30, method = "single",
##                            return_pairs = FALSE, ignore.small = TRUE,
##                            max.small = 1e4, ignore.isolated = TRUE,
##                            strict = c("strict", "one_to_one", "loose"),
##                            min.isolated = max.small,
##                            only_chains = FALSE) {


##   if (!is.character(strict) || any(!strict %in% c("strict", "one_to_one", "loose"))) {
##     stop("strict must be one of 'strict', 'one_to_one', or 'loose'")
##   } else if (length(strict) > 1) {
##     strict = "one_to_one"
##   }
##   self$edges$mark(ecluster = as.integer(NA))
##   altedges = self$edges[type == "ALT", ]
##   if (length(altedges) == 0) {
##     if (verbose) {
##       gmessage("No junction in this graph")
##     }
##     return(NULL)
##   }
##   if (ignore.small) {
##       altedges = altedges[!((class == "DUP-like" | class == "DEL-like") & altedges$span <= max.small)]
##   }
##   if (length(altedges) == 0) {
##     if (verbose) {
##       gmessage("No junction in this graph")
##     }
##     return(NULL)
##   }
##   deldup = altedges[class %in% c("DUP-like", "DEL-like")]
##   ## below removes non-nested events below size threshold
##   ## and simple or nested duplications with no subsumed breakends
##   if (length(deldup) > 0 && ignore.isolated) {
##     altes = deldup$shadow
##     ## altes$sedge.id = altedges[class %in% c("DUP-like", "DEL-like")]$dt[altes$id]$sedge.id
##     bp = grl.unlist(altedges$grl)[, c("grl.ix", "grl.iix", "class", "sedge.id")]
##     bp$sedge.id.y = bp$sedge.id; bp$sedge.id = NULL
##     addon = deldup$dt[altes$id][, .(sedge.id, class)]
##     altes$sedge.id = addon$sedge.id
##     altes$class = addon$class
##     altes$nbp = altes %N% bp # number of breakpoints of any SV that fall within segment
##     numsum = altedges$shadow %>% gr.sum # using the shadows of all of the SVs not just dels and dups
##     altes = altes %$% numsum
##     iso = ((altes) %Q% (score == 1.0))$id
##     ## rm.edges = unique(altes[iso] %Q% (width < thresh))$sedge.id ## old
##     rm.edges = unique(altes[iso] %Q% (width < min.isolated))$sedge.id
##     rm.dups = S4Vectors::with(altes, sedge.id[class == "DUP-like" & nbp <= 2])
##     rm.dups = c(rm.dups, dedup.cols(gr2dt(altes %*% bp))[sedge.id != sedge.id.y][class == "DUP-like"][, .(all(1:2 %in% grl.iix), class.1 = class.1[1]), by = .(sedge.id, sedge.id.y)][, all(V1 == TRUE) & all(class.1 == "DUP-like"), by = sedge.id][V1 == TRUE]$sedge.id) # removing dups that have only other nested dups 
##     rm.edges = union(rm.edges, rm.dups)
##     keepeid = setdiff(altedges$dt$sedge.id, rm.edges)
##     altedges = altedges[as.character(keepeid)]
##   } # ignoring isolated dup and del edges that are smaller than threshold
##   if (verbose & weak)
##     message("Computing weak eclusters")
##   if (length(altedges) == 0) {
##     if (verbose) {
##       gmessage("No junction in this graph")
##     }
##     return(NULL)
##   }
##   bp = grl.unlist(altedges$grl)[, c("grl.ix", "grl.iix", "edge.id")]
##   bp$m.ix = seq_along(bp)
##   bp.dt = gr2dt(bp)
##   ix = split(1:length(bp), ceiling(runif(length(bp)) * ceiling(length(bp)/chunksize)))
##   ixu = unlist(ix)
##   eps = 1e-09
##   ## ij = do.call(rbind, split(1:length(bp), bp$grl.ix))
##   xt.adj = xt.adj0 = Matrix::sparseMatrix(1, 1, x = 0, dims = rep(length(bp),
##                                                         2))
##   if (verbose) {
##     message(sprintf("Computing junction graph across %s ALT edges with distance threshold %s",
##                     length(altedges), thresh))
##   }
##   if (!exists(".INF")) {
##     .INF = pmax(sum(seqlengths(self)), 1e+09)
##   }
##   bp.pair = as.matrix(dcast(data.table(ix = bp$m.ix,
##                              grl.ix = bp$grl.ix,
##                              grl.iix = bp$grl.iix),
##                   grl.ix ~ grl.iix, value.var = "ix")[,2:3])
##   bp.pair = rbind(bp.pair, cbind(bp.pair[,2], bp.pair[,1]))
##   ifun = function(iix, ignore.strand = FALSE,
##                   verbose = FALSE, eps = 1e-9) {
##     if (verbose > 1)
##       cat(".")
##     tmpm = gr.dist(bp[iix], gr.flipstrand(bp), ignore.strand = ignore.strand) +
##       eps
##     return(as(tmpm, "Matrix"))
##   }
##   xt.adj[ixu, ] = do.call(rbind,
##                           mclapply(ix, ifun, ignore.strand = FALSE,
##                                    mc.cores = mc.cores))
##   diag(xt.adj) = NA_real_
##   ## only_chains = TRUE, enforcing that nearest breakpoints are only considered
##   if (only_chains) {
    
##     ## enforcing that no distances between breakends from the same junction  are considered
##     xt.adj[bp.pair] = NA_real_
##     xt.adj0[ixu, ] = do.call(rbind,
##                              mclapply(ix, ifun, ignore.strand = TRUE,
##                                       mc.cores = mc.cores)) # this is to find nearest breakends, regardless of orientation
##     ## enforcing that self-to-self breakend distances are not considered
##     diag(xt.adj0) = NA_real_
##     ## enforcing that no distances between breakends from the same junction  are considered
##     xt.adj0[bp.pair] = NA_real_
##     suppressWarnings({
##       nearest_ix = dunlist(lapply(
##         seq_len(nrow(xt.adj)),
##         function(x) which(xt.adj[x,] == min(xt.adj[x,], na.rm = T)))) %>%
##         as.matrix
##       nearest0_ix = dunlist(lapply(
##         seq_len(nrow(xt.adj0)),
##         function(x) which(xt.adj0[x,] == min(xt.adj0[x,], na.rm = T)))) %>%
##         as.matrix
##     })
##     if (nrow(nearest_ix) & nrow(nearest0_ix)) {
##       if (strict == "strict") {
##         nearest_ix = cbind(rowMins(nearest_ix), rowMaxs(nearest_ix))
##         nearest0_ix = cbind(rowMins(nearest0_ix), rowMaxs(nearest0_ix))
##         nearest_ix = nearest_ix[duplicated(nearest_ix),,drop = FALSE]
##         nearest0_ix = nearest0_ix[duplicated(nearest0_ix),,drop = FALSE]
##         nearest_ix = rbind(nearest_ix, cbind(nearest_ix[,2], nearest_ix[,1]))
##         nearest0_ix = rbind(nearest0_ix, cbind(nearest0_ix[,2], nearest0_ix[,1]))
##         nearest_ix = as.matrix(merge(nearest_ix, nearest0_ix)) # if there is a breakend closer but in the wrong orientation, that distance will be thrown out downstream, i.e. there is a one to one match of breakend and any nearest breakend that is in the wrong orientation disqualifies the clustser
##       } else if (strict == "one_to_one") {
##         nearest_ix = cbind(rowMins(nearest_ix), rowMaxs(nearest_ix))
##         nearest_ix = nearest_ix[duplicated(nearest_ix),,drop = FALSE]
##         nearest_ix = rbind(nearest_ix, cbind(nearest_ix[,2], nearest_ix[,1]))
##         nearest_ix = as.matrix(nearest_ix) # only one-to-one breakends are considered
##       } else if (strict == "loose") {
##         nearest_ix = rbind(nearest_ix, cbind(nearest_ix[,2], nearest_ix[,1]))
##         nearest_ix = nearest_ix[!duplicated(nearest_ix),,drop = FALSE]
##       }
##     }
##   } else if (!only_chains) {
##     xt.adj[bp.pair] = 1
##   }
##   ## rm(xt.adj0)
##   adj = xt.adj
##   xt.adj[which(is.na(as.matrix(xt.adj)))] = .INF + 1
##   adj[which(is.na(as.matrix(adj)))] = 0
##   adj[which(as.matrix(adj) > thresh)] = 0
##   if (only_chains && nrow(nearest_ix)) {
##     tmp = xt.adj[nearest_ix]
##     xt.adj[] = .INF + 1
##     adj[] = 0
##     xt.adj[nearest_ix] = tmp
##     adj[nearest_ix] = tmp
##   }
##   dt = Matrix::which(xt.adj < thresh, arr.ind = T)
##   dt = unique(data.table(cbind(rowMins(dt), rowMaxs(dt))))
##   dt[, bp.dist := xt.adj[dt[, cbind(V1, V2)]]]
##   dt$sign = ifelse(strand(bp[dt$V1]) == "+",
##             ifelse(gr.flipstrand(bp[dt$V1]) > bp[dt$V2], -1, 1),
##             ifelse(gr.flipstrand(bp[dt$V1]) < bp[dt$V2], -1, 1))
##   dt2 = dt[,idx := seq_len(.N)] %>% melt(measure.vars = c("V1", "V2"))
##   meta = cbind(gr2dt(bp)[dt2$value], dt2)[order(idx)]
##   meta = rbind(meta,
##               gr2dt(bp)[paste(grl.ix, grl.iix) %nin% meta[, paste(grl.ix, grl.iix)]],
##               fill = T)
##   if (return_pairs) {
##     return(meta)
##   }
##   hcl = stats::hclust(as.dist(xt.adj), method = "single")
##   hcl.lbl = cutree(hcl, h = thresh)
##   bp.dt$hcl = hcl.lbl
##   bp.hcl = bp.dt[, .(hcl.1 = .SD[grl.iix == 1, hcl], hcl.2 = .SD[grl.iix ==
##                                                                  2, hcl]), keyby = grl.ix]
##   altedges$mark(hcl.1 = bp.hcl[.(seq_along(altedges)), hcl.1])
##   altedges$mark(hcl.2 = bp.hcl[.(seq_along(altedges)), hcl.2])
##   hcl.ig = igraph::graph_from_edgelist(bp.hcl[, unique(cbind(hcl.1,
##                                                              hcl.2))], directed = FALSE)
##   hcl.comp = components(hcl.ig)
##   altedges$mark(ehcl = as.integer(hcl.comp$membership)[bp.hcl[,
##                                                               hcl.1]])
##   adj[adj > thresh] = 0
##   refg = self[, type == "REF"]
##   bpp = Matrix::which(adj != 0, arr.ind = TRUE)
##   dref = pdist(bp[bpp[, 1]], bp[bpp[, 2]])
##   drefg = diag(refg$dist(bp[bpp[, 1]], bp[bpp[, 2]]))
##   ix = which(drefg > dref)
##   if (length(ix))
##     adj[bpp[ix, , drop = FALSE]] = FALSE
##   if (verbose > 1)
##     cat("\n")
##   adj = adj | t(adj)
##   junpos = bp1 = bp$grl.iix == 1
##   junneg = bp2 = bp$grl.iix == 2
##   adj2 = adj & FALSE
##   adj2[junpos, junpos] = adj[bp2, bp1]
##   adj2[junpos, junneg] = adj[bp2, bp2]
##   adj2[junneg, junpos] = adj[bp1, bp1]
##   adj2[junneg, junneg] = adj[bp1, bp2]
##   if (verbose)
##     message(sprintf("Created basic junction graph using distance threshold of %s",
##                     thresh))
##   cl = split(1:length(bp), igraph::clusters(graph.adjacency(adj2),
##                                             ifelse(weak, "weak", "strong"))$membership)
##   cl = cl[S4Vectors::elementNROWS(cl) > 1]
##   cl = cl[order(S4Vectors::elementNROWS(cl))]
##   ## browser()
##   ## jcl = lapply(cl, function(x) unique(sort(bp$grl.ix[x])))
##   jcl = lapply(cl, function(x) unique(sort(bp$edge.id[x])))
##   jcls = sapply(jcl, paste, collapse = " ")
##   jcl = jcl[!duplicated(jcls)]
##   adj3 = adj2
##   altedges$mark(ecycle = as.character(NA))
##   if (length(jcl) > 0) {
##     dcl = dunlist(unname(jcl))[, `:=`(listid, paste0(ifelse(weak,
##                                                             "", "c"), listid))]
##     if (!weak)
##         ## altedges[dcl$V1]$mark(ecycle = dcl$listid)
##         altedges[as.character(dcl$V1)]$mark(ecycle = dcl$listid)
##     ## altedges[dcl$V1]$mark(ecluster = dcl$listid)
##     altedges[as.character(dcl$V1)]$mark(ecluster = dcl$listid)
##     meta = merge(meta, altedges$dt[, .(edge.id, ecluster, hcl.1, hcl.2, ehcl, ecycle)], by = "edge.id")
##     meta[!is.na(ecluster),
##          `:=`(
##            nclust = length(unique(edge.id)),
##            all_positive = all(replace(sign, is.na(sign),  3e9) > 0),
##            all_negative = all(replace(sign, is.na(sign), -3e9) < 0),
##            mixed = {naom = na.omit(sign); any(naom > 0) & any(naom < 0)},
##            bridge = anyNA(bp.dist)
##            ),
##          by = ecluster]
##     meta[, `:=`(
##         num_positive = sum(sign[!duplicated(idx)] > 0, na.rm = T),
##         num_negative = sum(sign[!duplicated(idx)] < 0, na.rm = T)
##     ),
##     by = ecluster]
##     self$set(recip_bp = meta)
##     if (verbose) {
##         message("Annotated weakly connected junction clusters and added ecluster pair metadata")
##     }
##   }
##   if (verbose)
##     message(sprintf("Annotated %s junction cycles in edge field $ecycle",
##                     length(jcl)))
##   if (paths & !weak) {
##     if (verbose)
##       message("Analyzing paths")
##     if (length(jcl) > 0) {
##       adj3[unlist(jcl), unlist(jcl)] = FALSE
##     }
##     sinks = Matrix::which(Matrix::rowSums(adj3) == 0)
##     sources = Matrix::which(Matrix::colSums(adj3) == 0)
##     cl2 = split(1:length(bp), igraph::clusters(graph.adjacency(adj3),
##                                                "weak")$membership)
##     cl2 = cl2[S4Vectors::elementNROWS(cl2) > 1]
##     if (any(ix <- S4Vectors::elementNROWS(cl2) > 2)) {
##       cl3 = do.call(c, mclapply(cl2[ix], function(x) {
##         tmp.adj = adj3[x, x]
##         lapply(all.paths(tmp.adj, sources = sources,
##                          sinks = sinks, verbose = verbose)$paths, function(i) x[i])
##       }, mc.cores = mc.cores))
##       cl2 = c(cl2[!ix], cl3)
##     }
##     jcl2 = lapply(cl2, function(x) unique(sort(bp$grl.ix[x])))
##     jcls2 = sapply(jcl2, paste, collapse = " ")
##     jcl2 = jcl2[!duplicated(jcls2)]
##     altedges$mark(epath = as.character(NA))
##     if (length(jcl2) > 0) {
##       dcl2 = dunlist(unname(jcl2))[, `:=`(listid, paste0("p",
##                                                          listid))]
##       altedges[dcl2$V1]$mark(epath = dcl2$listid)
##       self$edges$mark(ecluster =
##                         ifelse(is.na(self$edges$dt$ecycle) &
##                                is.na(self$edges$dt$epath),
##                                as.character(NA),
##                                paste0(ifelse(is.na(self$edges$dt$ecycle), "",
##                                              self$edges$dt$ecycle),
##                                       ifelse(is.na(self$edges$dt$epath),
##                                              "", self$edges$dt$epath))))
##     }
##     if (verbose)
##       message(sprintf("Annotated %s paths in edge field $epath",
##                       length(jcl2)))
##   }
##   return(invisible(self))
## }
## ## gGraph$public_methods$eclusters2 = tmpeclustpairs # if replacing a binding

## #' @description
## #' make eclusters
## gGraph$public_methods$eclusters2 = NULL; gGraph$set("public", "eclusters2", eclusters2)



#' @name copy3
#' @title make deep copy, recursively
#'
#' useful for dev
#' makes deep copy of R6 object, S4 object, or anything else really
#'
copy3 = function (x, recurse_list = TRUE) {
    if (inherits(x, "R6")) {
        x2 = rlang::duplicate(x$clone(deep = T))
        for (name in intersect(names(x2$.__enclos_env__), c("private", 
            "public"))) for (nname in names(x2$.__enclos_env__[[name]])) tryCatch({
            x2$.__enclos_env__[[name]][[nname]] = copy3(x2$.__enclos_env__[[name]][[nname]])
        }, error = function(e) NULL)
        return(x2)
    } else if (isS4(x)) {
        x2 = rlang::duplicate(x)
        slns = slotNames(x2)
        for (sln in slns) {
            tryCatch({
                slot(x2, sln) = copy3(slot(x2, sln))
            }, error = function(e) NULL)
        }
        return(x2)
    } else if (inherits(x, c("list"))) {
        x2 = rlang::duplicate(x)
        x2 = rapply(x2, copy3, how = "replace")
        return(x2)
    } else {
        x2 = rlang::duplicate(x)
        return(x2)
    }
}


#' @name rep_len2
#' @title recycle vector along length OR nrow of object
#'
#'
#' @author Kevin Hadi
#' @param x data
#' @param objalong any object to recycle x along if uselen = TRUE, or an actual integer value if uselen = FALSE
#' @return vector
rep_len2 = function(x, objalong, uselen = TRUE) {
    if (uselen)
        rep(x, length.out = NROW(objalong))
    else
        rep(x, length.out = objalong)
}

#' @name seq_along2
#' @title seq along either row of table or length of vector
#'
#'
#' @author Kevin Hadi
#' @param x data
#' @return vector
seq_along2 = function(x)  {
  seq_len(NROW(x))
}


#' @name match3
#' @title similar to setkey except a general use utility
#'
#' very slow version of keying a la data.table
#' but for general/interactive use
#' @author Kevin Hadi
match3 = function(x, table, nomatch = NA_integer_, old = TRUE, use.data.table = TRUE) {
  out = if (use.data.table) {
    tryCatch({
      dx = data.table(x = x)[, id.x := seq_len(.N)]
      dtb = data.table(table = table)[, id.tb := seq_len(.N)]
      ## setkey(dx, x)[list(dtb$table)]$id.x
      setkey(dtb, table)[list(dx$x)]$id.tb
    }, error = function(e) structure("err", class = "err"))
  }
  if (!is.null(out) && !inherits(out, "err")) return(out)
  if (old) {
    dx = within(data.frame(x = x), {id.x = seq_along(x)})
    dtb = within(data.frame(table = table), {id.tb = seq_along(table)})
    res = merge(dx, dtb, by.x = "x", by.y = "table", all.x = TRUE,
      allow.cartesian = TRUE)
    return(res$id.tb[order(res$id.x)])
  } else  {
    m = match(table,x)
    mat = cbind(m, seq_along(m))
    mat = mat[!is.na(mat[,1]),,drop=FALSE]
    mat = mat[order(mat[,1], na.last = FALSE),,drop = FALSE]
    mat = cbind(mat, seq_len(dim(mat)[1]))
    m2 = match(x,table)
    ix = which(!duplicated(m2) & !is.na(m2))
    mat_rix = unlist(rep(split(mat[,3], mat[,1]), base::tabulate(m2)[m2][ix]))
    ## mat_rix = unlist(rep(split(mat[,3], mat[,1]), base::tabulate(m2)[m2][ix]))
    ix = rep(1, length.out = length(m2))
    ## original line
    ## ix[!is.na(m2)] = base::tabulate(m)[!is.na(m2)]
    ix[!is.na(m2)] = base::tabulate(m)[m][m2][!is.na(m2)]
    out = rep(m2, ix)
    out[!is.na(out)] = mat[mat_rix,,drop=F][,2]
    return(out)
    ## m = match(table, x)
    ## mat = cbind(m, seq_along(m))
    ## mat = mat[!is.na(mat[, 1]), , drop = FALSE]
    ## mat = mat[order(mat[, 1]), , drop = FALSE]
    ## mat = cbind(mat, seq_len(dim(mat)[1]))
    ## m2 = match(x, table)
    ## ix = which(!duplicated(m2))
    ## mat_rix = unlist(rep(split(mat[, 3], mat[, 1]), base::tabulate(m2)[m2][ix]))
    ## mat[mat_rix, , drop = F][, 2]
  }
}


#' @name dedup.cols
#' @title applies dedup to colnames
#'
#' dedup the column names of a data.frame/data.table
#'
#' @return A data.table or data.frame
dedup.cols = function(tbl, remove = FALSE) {
    if (remove) {
        if (!inherits(tbl, "data.table"))
            return(tbl[, match(unique(colnames(tbl)), colnames(tbl))])
        else
            return(tbl[, match(unique(colnames(tbl)), colnames(tbl)), with = FALSE])
    } else {
            colnames(tbl) = base::make.unique(colnames(tbl))
            return(tbl)
    }
}


#' @name rowMins
#' @title rowMins hack
#'
#'
#' @return A vector
rowMins = function(x) {
  do.call(pmin, as.data.frame(x))
}

#' @name rowMaxs
#' @title rowMaxs hack
#'
#'
#' @return A vector
rowMaxs = function(x) {
  do.call(pmax, as.data.frame(x))
}

#' @name not.in
#' @title not.in
#'
#' @description
#' not in
#'
#' @param x value to test
#' @param table table to test x against
#'
#' @noRd
`%nin%` = function (x, table)
{
    match(x, table, nomatch = 0L) == 0L
}

#' @name gr.noval
#' @title get rid of mcols on GRanges/GRangesLists
#' @description
#'
#' remove all metadata from GRanges or GRangesList
#'
#' @return GRanges or GRangesList
#' @author Kevin Hadi
gr.noval = function(gr, keep.col = NULL, drop.col = NULL) {
    if (is.null(keep.col) & is.null(drop.col)) {
        select_col = NULL
    } else {
        all_col = colnames(gr@elementMetadata)
        if (inherits(gr, "GRangesList")) {
            all_col = c(all_col, colnames(gr@unlistData@elementMetadata))
        }

        if (!is.null(keep.col) & is.null(drop.col)) {
            select_col = intersect(all_col, keep.col)
        } else if (is.null(keep.col) & !is.null(drop.col)) {
            select_col = setdiff(all_col, drop.col)
        } else if (!is.null(keep.col) && !is.null(drop.col)) {
            if (intersect(keep.col, drop.col) > 0) {
                warning("drop.col and keep.col args have overlapping elements\nkeeping the columns that overlap")
                select_col = intersect(setdiff(all_col, setdiff(drop.col, keep.col)), keep.col)
            }
        }
    }
    if (inherits(gr, "GRangesList")) {
        tmp_query = intersect(select_col, colnames(gr@unlistData@elementMetadata))
        gr@unlistData@elementMetadata = gr@unlistData@elementMetadata[,c(tmp_query), drop = FALSE]
    }
    tmp_query = intersect(select_col, colnames(gr@elementMetadata))
    gr@elementMetadata = gr@elementMetadata[,c(tmp_query),drop = FALSE]
    return(gr)
}


#' @name setcols
#' @title convenience function to set columns
#'
#' sets columns of an object
#'
#' @param dt data frame/table or matrix
#' @param old integer or character or logical vector corresponding to current colnames in dt
#' @param new character vector for new column names
#' @return colnamed object
#' @author Kevin Hadi
setcols = function(dt, old, new) {
  if (inherits(dt, c("GRanges", "GRangesList"))) {
    mcols(dt) = setcols(mcols(dt), old, new)
    return(dt)
  }
  cnames = colnames2(dt)
  if (missing(new) || missing(old)) {
    if (missing(old)) {
      old = new
    }
    if (is.character(old) && length(old) == length(cnames)) {
      colnames(dt) = old
      return(dt)
    } else {
      stop("names provided must be same length as ncol(dt)")
    }
  }
  if (is.character(old)) {
    out = merge(data.frame(cnames, seq_along(cnames)), data.frame(cnames = old, new = new),
      allow.cartesian = T)
    cnames[out[[2]]] = out[[3]]
    colnames(dt) = cnames
    return(dt)
  }
  if (is.logical(old)) {
    if (! length(old) == length(cnames)) stop("logical vector must be same length as ncol(dt)")
    old = which(old)
  }
  cnames[old] = new
  colnames(dt) = cnames
  return(dt)
}

#' @name colnames2
#' @title robust colnames
#'
#' gives back character vector same number of columns of input regardless whether named or not
#'
#' @param str a path string
#' @return a string with multiple parentheses replaced with a single parenthesis
#' @author Kevin Hadi
colnames2 = function(x) {
    nm = colnames(x)
    if (is.null(nm))
        return(rep("", length.out = NCOL(x)))
    else
        return(nm)
}




#' @name %K%
#' @title similar to setkey except a general use utility
#'
#' slower version of setkey, but for interactive use
#'
#' @author Kevin Hadi
`%K%` = function(thisx,thisy, old = TRUE) {
    if (old)
        return(match3(table = thisx, x = thisy, old = TRUE))
    else
        return(match3(table = thisx, x = thisy, old = FALSE))
}


##' @name enframe_list
##' @title data table-ize list elements and add name column
##'
##' @description
##'
##' @author Kevin Hadi
#enframe_list = function(lst, name = "name", value = "value", as.data.table = TRUE, rbind = TRUE, mc.cores = 1) {
#    nms = names(lst)
#    expr = parse(text = sprintf("cbind(%s = nm, df)", name, value))
#    out = mcmapply(function(el, nm) {
#        if (NROW(el)) {
#            if (is.null(dim(el)))
#                df = setColnames(as.data.frame(el), value)
#            else
#                df = as.data.frame(el)
#            if (as.data.table)
#                setDT(df)
#            eval(expr)
#        }
#    }, lst, nms, SIMPLIFY = FALSE, mc.cores = mc.cores)
#    if (rbind) {
#        if (as.data.table)
#            return(rbindlist(out))
#        else
#            return(do.call(rbind, out))
#    }
#}

##' @name setColnames
##' @title convenience function to set column names
##'
##' @param object tabled object
##' @param nm names of the new columns
##' @return colnamed object
##' @author Kevin Hadi
#setColnames = function(object = nm, nm = NULL, pattern = NULL, replacement = "") {
#    if (!is.null(nm)) {
#        if (is.null(names(nm)))
#            colnames2(object)  = nm
#        else {
#            ix = match3(names(nm), colnames(object))
#            colnames2(object)[ix] = nm
#        }
#    } else if (!is.null(pattern)) {
#        colnames2(object) = gsub(pattern, replacement, colnames2(object))
#    }
#    return(object)
#}

#' @name isNA
#' @title is.na but also tests for "NA" character
#'
#' @description
#'
#' @author Kevin Hadi
isNA = function(x, na.char = c("NA", "NULL", "na", "null")) {
    if (is.character(x)) {
        return(is.na(x) | x %in% na.char)
    } else {
        return(is.na(x))
    }
}

#' @name na2false
#' @title replace logical vector with NA to FALSE
#'
#' @description
#' A convenience function to set a logical vector with NAs to false
#'
#' @return A logical vector with NAs set to FALSE
#' @author Kevin Hadi
na2false = function(v)
{
    v[isNA(v)] = FALSE
    ## mode(v) = "logical"
    v
}

#' @name rand.string
#' @title make a random string
#'
#' @return random string
#' @author Someone from Stackoverflow
rand.string <- function(n=1, length=12)
{
    randomString <- c(1:n)                  # initialize vector
    for (i in 1:n)
    {
        randomString[i] <- paste(sample(c(0:9, letters, LETTERS),
                                        length, replace=TRUE),
                                 collapse="")
    }
    return(randomString)
}

#' @name gr.spreduce
#' @title reduce based on a field(s) to split by in elementMetadata of GRanges, or given vector
#' @description
#'
#' split and reduce GRanges by field(s)
#' if providing a variable not already within the GRanges,
#' may need to use dynget(variable_name)
#'
#' @return GRanges
#' @author Kevin Hadi
gr.spreduce = function(gr,  ..., ignore.strand = FALSE, pad = 0, return.grl = FALSE, sep = paste0(" ", rand.string(length = 8), " ")) {
  lst = as.list(match.call())[-1]
  ix = which(!names(lst) %in% c("gr", "sep", "pad", "ignore.strand", "return.grl"))
  vars = unlist(sapply(lst[ix], function(x) unlist(sapply(x, toString))))
  if (length(vars) == 1) {
    if (!vars %in% colnames(mcols(gr)))
      vars = tryCatch(unlist(list(...)), error = function(e) vars)
  }
  if (!all(vars %in% colnames(mcols(gr))))
    stop("Must specify valid metadata columns in gr")
  tmpix = do.call(
    function(...) paste(..., sep = sep),
    as.list(mcols(gr)[,vars, drop = F]))
  unix = which(!duplicated(tmpix))
  tmpix = factor(tmpix, levels = tmpix[unix])
  grl = unname(gr.noval(gr) %>% GenomicRanges::split(tmpix))
  grl = GenomicRanges::reduce(grl + pad, ignore.strand = ignore.strand)
  if (return.grl) {
    mcols(grl) = mcols(gr)[unix,vars,drop = F]
    return(grl)
  } else {
    out = unlist(grl)
    mcols(out) = mcols(gr)[rep(unix, times = IRanges::width(grl@partitioning)),
      vars,drop = F]
    return(out)
  }
}

#' @name duped
#' @title duped
#'
#'
#'
#' @param ... vectors to paste by
#' @author Kevin Hadi
duped = function(..., binder = "data.table") {
    duplicated(tryCatch(et(sprintf("%s(...)", binder)),
                        ## error = function(e) do.call(cbind, list(...))))
                        error = function(e) paste(...)))
}

#' @name merge.repl
#' @title merging data tables with collapsing columns with the same name
#'
#' Merge two data tables with various replacing strategies
#' for columns common between x and y that are not used to merge
#' (i.e. not specified in the "by" argument)
#'
#' @param replace_NA logical, only use values in dt.y, any dt.x not in dt.y is clobbered (NA)
#' @param force_y logical, should x and y common columns be merged?
#' @param overwrite_x logical, if force_y = TRUE, should NA values in y replace x?
#' @return A data.table
#' @author Kevin Hadi
merge.repl = function(dt.x,
                      dt.y,
                      sep = "_",
                      replace_NA = TRUE,
                      force_y = TRUE,
                      overwrite_x = FALSE,
                      keep_order = FALSE,
                      keep_colorder = TRUE,
                      keep_factor = TRUE,
                      ...) {
    arg_lst = as.list(match.call())
    by.y = eval(arg_lst[['by.y']], parent.frame())
    by.x = eval(arg_lst[['by.x']], parent.frame())
    by = eval(arg_lst[['by']], parent.frame())
    all.x = eval(arg_lst[['all.x']], parent.frame())
    all.y = eval(arg_lst[['all.y']], parent.frame())
    all = eval(arg_lst[['all']], parent.frame())
    allow.cartesian = eval(arg_lst[['allow.cartesian']])
    key_x = key(dt.x)
    if (is.null(all.x)) {
        all.x = TRUE
    }
    if (is.null(all.y)) {
        all.y = FALSE
    }
    if (!is.null(all) && all) {
        all.y = TRUE
        all.x = TRUE
    }
    if (is.null(allow.cartesian)) {
        allow.cartesian = FALSE
    }
    if (!inherits(dt.x, "data.table")) {
        dt.x = as.data.table(dt.x)
    }
    if (!inherits(dt.y, "data.table")) {
        dt.y = as.data.table(dt.y)
    }
    if (keep_order == TRUE) {
        dt.x[['tmp.2345098712340987']] = seq_len(nrow(dt.x))
    }

    dt.x[['in.x.2345098712340987']] = rep(TRUE, length.out = nrow(dt.x))
    dt.y[['in.y.2345098712340987']] = rep(TRUE, length.out = nrow(dt.y))

    new_ddd_args = list(by = by, by.x = by.x, by.y = by.y, all.x = all.x, all.y = all.y, allow.cartesian = allow.cartesian)

    if (is.null(by.y) & is.null(by.x) & is.null(by)) {

        if (length(attributes(dt.x)[['sorted']]) > 0 &&
            length(attributes(dt.y)[['sorted']]) > 0) {
            k.x = key(dt.x)
            k.y = key(dt.y)
        } else {
            k.y = k.x = intersect(names2(dt.x), names2(dt.y))
            if (length(k.x) == 0)
                stop("no common columns to merge by!")
            message("intersecting by: ", paste(k.x, collapse = ", "))
            new_ddd_args[['by']] = k.x
        }
        if (is.null(k.x) | is.null(k.y) || (k.x != k.y)) {
            stop("neither by.x/by.y  nor by are supplied, keys of dt.x and dt.y must be identical and non NULL")
        }
        x.cols = setdiff(names(dt.x), k.x)
        y.cols = setdiff(names(dt.y), k.y)

    } else if (!is.null(by.x) & !is.null(by.y)) {

        x.cols = setdiff(names(dt.x), by.x)
        y.cols = setdiff(names(dt.y), by.y)
        new_ddd_args = new_ddd_args[setdiff(names(new_ddd_args), c("by"))]

    } else if (!is.null(by)) {

        x.cols = setdiff(names(dt.x), by)
        y.cols = setdiff(names(dt.y), by)
        if (! all(by %in% colnames(dt.x)) | ! all(by %in% colnames(dt.y))) {
            stop("column ", by, " does not exist in one of the tables supplied \nCheck the column names")
        }
        new_ddd_args = new_ddd_args[setdiff(names(new_ddd_args), c("by.y", "by.x"))]

    }
    these_cols = intersect(x.cols, y.cols)
    ## if (replace_in_x) {
    if (!replace_NA) {
        dt.x.tmp = copy(dt.x)
        for (this_col in these_cols) {
            data.table::set(dt.x.tmp, i = NULL, j = this_col, value = NULL)
        }
        dt.repl = suppressWarnings(do.call("merge", args = c(list(x = dt.x.tmp, y = dt.y), new_ddd_args)))
        ## dt_na2false(dt.repl, c("in.x.2345098712340987", "in.y.2345098712340987"))
    } else {
        dt.repl = suppressWarnings(do.call("merge", args = c(list(x = dt.x, y = dt.y), new_ddd_args)))
        dt_na2false(dt.repl, c("in.x.2345098712340987", "in.y.2345098712340987"))
        in.x = which(dt.repl[["in.x.2345098712340987"]])
        in.y = which(dt.repl[["in.y.2345098712340987"]])
        this_env = environment()
        for (this_col in these_cols) {
            x_cname = paste0(this_col, ".x")
            y_cname = paste0(this_col, ".y")
            x_col = dt.repl[[x_cname]]
            y_col = dt.repl[[y_cname]]
            xf = inherits(x_col, "factor")
            yf = inherits(y_col, "factor")
            if ( {xf || yf} && keep_factor) {
                if (!xf) { x_col = factor(x_col); xf = TRUE }
                if (!yf) { y_col = factor(y_col); yf = TRUE }
            }
            if (xf && !keep_factor) { x_col = as.character(x_col); xf = FALSE } 
            if (yf && !keep_factor) { y_col = as.character(y_col); yf = FALSE }
            if (force_y) {
                if (!overwrite_x) {
                    ## if (inherits(x_col, "factor") & inherits(y_col, "factor")) {
                    ##     new_col = factor(y_col, forcats::lvls_union(list(y_col, x_col)))
                    ##     new_col[is.na(new_col)] = x_col[is.na(new_col)]
                    ## } else {
                    ##     new_col = ifelse(!is.na(y_col), y_col, x_col)
                    ## }                    
                    if (xf || yf) {
                        new_col = factor(y_col, forcats::lvls_union(list(y_col, x_col)))
                        new_col[is.na(new_col)] = x_col[is.na(new_col)]
                    } else {
                        new_col = ifelse(!is.na(y_col), y_col, x_col)
                    }
                } else {
                    ## if (inherits(x_col, "factor") & inherits(y_col, "factor")) {
                    ##     new_col = factor(x_col, forcats::lvls_union(list(y_col, x_col)))
                    ## } else {
                    ##     new_col = x_col
                    ## }
                    ## new_col[dt.repl[['in.y.2345098712340987']]] = y_col[dt.repl[['in.y.2345098712340987']]]
                    if (xf || yf) {
                        new_col = factor(x_col, forcats::lvls_union(list(y_col, x_col)))
                    } else {
                        new_col = x_col
                    }
                    new_col[in.y] = y_col[in.y]
                }
            } else {
                ## if (inherits(x_col, "factor") & inherits(y_col, "factor")) {
                ##     new_col = factor(x_col, forcats::lvls_union(list(x_col, y_col)))
                ##     new_col[is.na(new_col) & !is.na(y_col)] = y_col[is.na(new_col) & !is.na(y_col)]
                ## } else {
                ##     new_col = ifelse(is.na(x_col) & !is.na(y_col), y_col, x_col)
                ## }
                if (xf | yf) {
                    new_col = factor(x_col, forcats::lvls_union(list(x_col, y_col)))
                    new_col[is.na(new_col) & !is.na(y_col)] = y_col[is.na(new_col) & !is.na(y_col)]
                } else {
                    new_col = ifelse(is.na(x_col) & !is.na(y_col), y_col, x_col)
                }
            }
            data.table::set(dt.repl, j = c(x_cname, y_cname, this_col), value = list(NULL, NULL, this_env[["new_col"]]))
        }
    }
    ## } else if (!replace_in_x & !is.null(suffix)) {
    ##     y.suff.cols = paste0(y.cols, sep, suffix)
    ##     ## dt.y.tmp = copy(dt.y)[, eval(dc(y.suff.cols)) := eval(dl(y.cols))][, eval(dc(y.cols)) := NULL]
    ##     dt.y.tmp = copy(dt.y)
    ##     data.table::set(dt.y, j = y.suff.cols, value = dt.y[, y.cols, with = FALSE])
    ##     data.table::set(dt.y, j = y.cols, value = NULL)
    ##     ## dt.repl = merge(dt.x, dt.y.tmp, all.x = TRUE, ...)
    ##     dt.repl = do.call("merge", args = c(list(x = dt.x, y = dt.y.tmp), new_ddd_args))
    ## }
    if (keep_order == TRUE) {
        data.table::setorderv(dt.repl, "tmp.2345098712340987")
        dt.repl[['tmp.2345098712340987']] = NULL
    }
    data.table::set(dt.repl, j = c("in.y.2345098712340987", "in.x.2345098712340987"),
                    value = list(NULL, NULL))
    if (keep_colorder) {
        x_cols = colnames(dt.x)
        ## get the order of columns in dt.repl in order of X with
        ## additional columns tacked on end
        setcolorder(dt.repl,
                    intersect(union(colnames(dt.x), colnames(dt.repl)),
                              colnames(dt.repl)))
    }
    return(dt.repl)
}


#' @name dt_na2false
#' @title convert columns with NA to false
#'
#' coerce NA in columns of class "logical" to FALSE
#'
#' @param dt data.table
#' @param these_cols NULL by default, will select columns of class logical, otherwise will be specified
#' @return A data.table
#' @author Kevin Hadi
dt_na2false = function(dt, these_cols = NULL) {
    na2false = function(v)
    {
        ## v = ifelse(is.na(v), v, FALSE)
        v[is.na(v)] = FALSE
        as.logical(v)
    }
    if (is.null(these_cols)) {
        these_cols = which(sapply(dt, class) == "logical")
    }
    for (this_col in these_cols) {
        ## this_val = as.data.frame(dt[, this_col, with = FALSE])[,1]
        this_val = dt[[this_col]]
        data.table::set(dt, j = this_col, value = na2false(this_val))
    }
    return(dt)
}
