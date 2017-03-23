#' Break GRanges at given breakpoints into disjoint gr
#'
#' @import gTrack
#' @import GenomicRanges
#' @import Matrix
#' @import parallel
#' @import data.table
#' @import gUtils
#'
#' @param query a disjoint \code{GRanges} object to be broken
#' @param bps \code{GRanges} of width 1, locations of the bp
#'
#' @return \code{GRanges} disjoint object at least the same length as query,
#' with a metadata column \code{qid} indicating input index where new segment is from
#'
#' @examples
#'
#'
#' @export
gr.breaks = function(query, bps=NULL){
    ## preprocess query
    if (!isDisjoint(query)){
        warning("Query GRanges not disjoint.")
        queryDj = disjoin(query)
        queryDj$qid = queryDj %N% query ## only retain the first occurence
        values(queryDj) = cbind(values(queryDj),
                                as.data.table(values(query))[queryDj$qid])
        query = queryDj
    } else {
        if ("qid" %in% colnames(values(query))){
            warning("'qid' col in query overwritten.")
        }
        query$qid = seq_along(query)
    }

    ## if bps not provided, return back-traced disjoin wrapper
    if (is.null(bps)) {
        return(query)
    } else {
        ## preprocess bps
        ## having strand info? remove it!
        if (any(strand(bps)!="*")){
            warning("Some breakpoints have strand info. Force to '*'.")
            bps = gr.stripstrand(bps)
        }
        ## some not a point? turn it into a point
        if (any(width(bps)!=1)){
            warning("Some breakpoint width>1.")
            bps = reduce(c(gr.start(bps), gr.end(bps)))
        }

        bps$inQuery = bps %^% query
        if (any(bps$inQuery==F)){
            warning("Some breakpoint not within query ranges.")
        }

        ## label and only consider breakpoints not already at the boundary of query
        bps$inner = bps$inQuery
        bps$inner[which(bps %^% gr.start(query) | bps %^% gr.end(query))]=F
        ## maybe no inner bp at all, then no need to proceed
        if (!any(bps$inner)){
            return(query)
        }
        bpsInner = bps %Q% (inner==T)
        ## map query and inner breakpoints
        qbMap = gr.findoverlaps(query, bpsInner)
        mappedQ = seq_along(query) %in% qbMap$query.id
        ## raw coors to construct ranges from
        tmpRange = data.table(qid2 = qbMap$query.id,
                              startFrom = start(query[qbMap$query.id]),
                              breakAt = start(bpsInner[qbMap$subject.id]),
                              upTo = end(query[qbMap$query.id]))
        tmpCoor = tmpRange[, .(pos=sort(unique(c(startFrom, breakAt, upTo)))), by=qid2]

        ## construct new ranges
        newRange = tmpCoor[, .(start=pos[-which.max(pos)],
                               end=pos[-which.min(pos)]), by=qid2]
        newRange[, ":="(chr = as.vector(seqnames(query)[qid2]),
                        strand = as.vector(strand(query)[qid2]))]
        newRange$start = newRange[, ifelse(start==min(start), start, start+1)]

        ## put together the mapped and broken
        newGr = GRanges(newRange)
        values(newGr) = values(query)[newGr$qid2, , drop=F] ## preserve the input metacol
        ## with the intact not mapped part of query
        output = c(newGr, query[!mappedQ]) %Q% (order(strand, seqnames, start))
        return(output)
    }

}
