#' applications of gGraph

#' @name balance
#' @title balance gGnome graphs 
#' @description
#'
#' Here we analyze gGraphs with "cn" (copy number) field to enforce integer
#' cn and junction balance, ie sum of incoming (or outgoing) edge
#' cn should be equal to node copy cn.
#'
#' The goal is to find a balaned assignment of "cn" to the nodes and edges of the gGraph
#' that maximally resemble the input weights while minimizing the loose end penalty.
#' The similarity / distance function can be weighted by optional node / edge
#' metadata field $weight (when weighted = TRUE). 
#'
#' To output this gGraph, we design a MIP with
#' (1) objective function that minimizes (weighted sum of square) distance of fit node and junction copy number to provided values in
#'     $cn field 
#' (2) and lambda* the sum of copy number at non terminal loose ends subject to 
#' (3) junction balance constraint
#' (4) fixing copy number of a subset of nodes and junctions
#'
#' Objective weight can be modulated at nodes and edges with $weight metadata
#' field (default node weight is node width, and edge weight is 1).
#' These fields will then set the penalty incurred to a fit of x to that node / edge
#' with copy number c and weight w as (x-c)^2/w.
#' 
#' Lambda can be modulated at nodes with $lambda node metadata field (default 1)
#'
#' For "haplographs" ie graphs that have more than one node overlapping a given location, it may
#' be important to constrain total copy number using a haploid read depth signal.
#' The marginal parameter enables this through a GRanges annotated with $cn and optionally $weight
#' field that provides a target total copy number that the optimization will attempt to satisfy.
#' This provided copy number c and weight w (default 1) will be evaluated against the
#' sum s of the fit copy numbers of all nodes overlapping that location by adding a penalty
#' of (c-s)^2/w to the corresponding solution. marginal can also have an optional logical field
#' $fix that will actually constrain the marginal copy number to be equal to the provided value
#' (note: that the optimization may be infeasible, and function will error out)
#' 
#' Additional controls can be inputted by changing the graph metadata - e.g. adding fields
#' $lb and $ub to nodes and edges will constrain their fit copy number to those bounds.
#' Adding $reward field to edges will add a reward for each copy of that edge in the solution.
#' 
#' 
#' @param gg gGraph with field $cn, can be NA for some nodes and edges, optional field $weight which will adjust the quadratic penalty on the fit to x as (x-$cn)^2/weight
#' @param lambda positive number specifying loose end penalty, note if gg$node metadata contain $lambda field then this lambda will be multiplied by the node level lambda (default 10)
#' @param marginal GRanges with field $cn and optional $weight field will be used to fit the summed values at each base of the genome to optimally fit the marginal value, optional field $fix will actually constrain the marginal to be the provided value
#' @param emarginal Junctions object with marginal CN in the $cn field (and optionally $weight in the weight field). optional field $fix will actually constrain the marginal to be the provided value.
#' @param tight indices or epxression on node metadata specifying at which nodes to disallow loose ensd
#' @param nfix indices or expression on node metadata specifying which node cn to fix
#' @param efix indices or expression on edge metadata specifying which edge cn to fix
#' @param nrelax indices or expression on node metadata specifying which nodes cn to relax
#' @param erelax  indices or expression on edge metadata specifying which edges cn to relax
#' @param L0  flag whether to apply loose end penalty as L1 (TRUE)
#' @param loose.collapse (parameter only relevant if L0 = TRUE) will count all unique (by coordinate) instances of loose ends in the graph as the loose end penalty, rather than each instance alone ... useful for fitting a metagenome graph   (FALSE)
#' @param phased (logical) indicates whether to run phased/unphased. default = FALSE
#' @param ism  (logical) additional ISM constraints (FALSE)
#' @param lp (logical) solve as linear program using abs value (default TRUE)
#' @param M  (numeric) big M constraint for L0 norm loose end penalty (default 1e3)
#' @param verbose (integer)scalar specifying whether to do verbose output, value 2 will spit out MIP (1)
#' @param tilim (numeric) time limit on MIP in seconds (10)
#' @param epgap (numeric) relative optimality gap threshhold between 0 and 1 (default 1e-3)
#' @param nsol (integer) number of solutions (default 1)
#' @param debug (logical) returns list with names gg and sol. sol contains full RCPLEX solution. (default FALSE)
#' 
#' @return balanced gGraph maximally resembling input gg in CN while minimizing loose end penalty lambda.
#' @author Marcin Imielinski
#' 
#' @export 
balance = function(gg,
                   lambda = 10,
                   marginal = NULL,
                   emarginal = NULL,
                   tight = NULL,
                   nfix = NULL, efix = NULL, nrelax = NULL, erelax = NULL,
                   L0 = TRUE,
                   loose.collapse = FALSE,
                   M = 1e3,
                   phased = FALSE,
                   ism = FALSE,
                   lp = TRUE,
                   verbose = 1,
                   tilim = 10,
                   epgap = 1e-3,
                   nsol = 1,
                   debug = FALSE)
{
    if (verbose) {
        message("creating copy of input gGraph")
    }

    gg = gg$copy

    if (verbose) {
        message("Checking inputs")
    }

    if (ism) {
        if (!L0) {
            stop("ISM can only be set to true if using L0 penalty")
        }
    }
    
    if (!('cn' %in% names(gg$nodes$dt)))
    {
        warning('cn field not defined on nodes, setting to NA')    
        gg$nodes$mark(cn = NA_real_)
    }

    if (!('cn' %in% names(gg$edges$dt)))
    {
        warning('cn not defined on edges, providing NA')    
        gg$edges$mark(cn = NA_real_)
    }

    if (phased) {
        if (!("allele" %in% names(gg$nodes$dt))) {
            stop("cannot run phased balance without $allele field in nodes")
        }
    }
        
    if (!is.null(marginal)) {
        if (!inherits(marginal, 'GRanges') || is.null(marginal$cn)) {
            stop('marginal must be a GRanges with field $cn')
        }
        if (is.null(marginal$fix)) {
            if (verbose) {
                message("$fix not supplied. marginals not fixed by default.")
            }
            marginal$fix = 0
        }
        if (is.null(marginal$weight)) {
            if (verbose) {
                message("$weight not supplied. set to range width in Mbp by default.")
            }
            marginal$weight = width(marginal)
        }
    }

    if (!is.null(emarginal)) {
        if (!inherits(emarginal, 'Junction') || is.null(emarginal$dt$cn)) {
            stop('emarginal must be Junction with field $cn')
        }
        ## don't mutate?
        ## emarginal = emarginal$copy
        if (is.null(emarginal$dt$fix)) {
            if (verbose) {
                message('$fix not supplied in emarginal. not fixed by default')
            }
            emarginal$set(fix = 0)
        }
        if (is.null(emarginal$dt$weight)) {
            if (verbose) {
                message("$weight not supplied in emarginal. set to 1 by default")
            }
            emarginal$set(weight = 1)
        }
    }
             
    ## default local lambda: default local lambda is 1 for consistency with JaBbA
    if (!('lambda' %in% names(gg$nodes$dt)))
        gg$nodes$mark(lambda = 1)

    ## default node weight is its width
    if (!('weight' %in% names(gg$nodes$dt)))
    {
        gg$nodes$mark(weight = width(gg$nodes$gr))
    }

    ## default edge weight is its width
    if (!('weight' %in% names(gg$edges$dt)))
    {
        gg$edges$mark(weight = 1)
    }

    ## default reward is 0 
    if (!('reward' %in% names(gg$edges$dt)))
    {
        gg$edges$mark(reward = 0)
    }
  
    ## handle parsing of efix, nfix, nrelax, erelax
    if (!any(deparse(substitute(nfix)) == "NULL")) ## R voodo to allow "with" style evaluation 
        nfix = tryCatch(eval(eval(parse(text = substitute(deparse(substitute(nfix)))), parent.frame()), gg$nodes$dt, parent.frame(2)), error = function(e) NULL)

    if (!any(deparse(substitute(nrelax)) == "NULL")) ## R voodo to allow "with" style evaluation 
        nrelax = tryCatch(eval(eval(parse(text = substitute(deparse(substitute(nrelax)))), parent.frame()), gg$nodes$dt, parent.frame(2)), error = function(e) NULL)

    if (!any(deparse(substitute(efix)) == "NULL")) ## R voodo to allow "with" style evaluation 
        efix = tryCatch(eval(eval(parse(text = substitute(deparse(substitute(efix)))), parent.frame()), gg$edges$dt, parent.frame(2)), error = function(e) NULL)

    if (!any(deparse(substitute(erelax)) == "NULL")) ## R voodo to allow "with" style evaluation 
        erelax = tryCatch(eval(eval(parse(text = substitute(deparse(substitute(erelax)))), parent.frame()), gg$edges$dt, parent.frame(2)), error = function(e) NULL)

    if (!any(deparse(substitute(tight)) == "NULL")) ## R voodo to allow "with" style evaluation 
        tight = tryCatch(eval(eval(parse(text = substitute(deparse(substitute(tight)))), parent.frame()), gg$nodes$dt, parent.frame(2)), error = function(e) NULL)


    if (is.logical(nfix))
        nfix = which(nfix)

    if (is.logical(efix))
        efix = which(efix)
    
    if (is.logical(nrelax))
        nrelax = which(nrelax)

    if (is.logical(erelax))
        erelax = which(erelax)

    if (length(nfix) & verbose)
        message('Fixing ', length(nfix), ' nodes')

    if (length(efix) & verbose)
        message('Fixing ', length(efix), ' edges')

    if (length(nrelax) & verbose)
        message('Relaxing ', length(nrelax), ' nodes')

    gg$nodes[nrelax]$mark(weight = 0)
    
    if (length(erelax) & verbose)
        message('Relaxing ', length(erelax), ' edges')

    gg$nodes[erelax]$mark(weight = 0)

    if (!is.logical(tight))
        tight = 1:length(gg$nodes) %in% tight

    if (any(tight) & verbose)
        message('Leaving ', sum(tight), ' nodes tight')
    
    gg$nodes$mark(tight = tight)

    if (is.null(gg$nodes$dt$lb))
        gg$nodes$mark(lb = 0)

    if (is.null(gg$nodes$dt$ub))
        gg$nodes$mark(ub = Inf)

    if (is.null(gg$edges$dt$lb))
        gg$edges$mark(lb = 0)

    if (is.null(gg$edges$dt$ub))
        gg$edges$mark(ub = Inf)

  if (loose.collapse)
  {
    if (verbose)
      message('Collapsing loose ends')

    uleft = unique(gr.start(gg$nodes$gr))
    uright = unique(gr.end(gg$nodes$gr))
    
    gg$nodes$mark(loose.left.id = paste0(gr.match(gr.start(gg$nodes$gr), uleft), 'l'))
    gg$nodes$mark(loose.right.id = paste0(gr.match(gr.end(gg$nodes$gr), uright), 'r'))      
  }
  else
  {
    gg$nodes$mark(loose.left.id = paste0(1:length(gg$nodes), 'l'))
    gg$nodes$mark(loose.right.id = paste0(1:length(gg$nodes), 'r'))
  }

  ########
  ## VARIABLES
  ########

  ## create state space, keeping track of graph ids
  vars = rbind(
    gg$dt[, .(cn, snode.id, lb, ub, weight, gid = index, type = 'node', vtype = 'I')], ## signed nodes
    gg$sedgesdt[, .(from, to, lb, ub, sedge.id,  cn, reward, gid = sedge.id, type = 'edge', vtype = 'I')], ## signed edges

    ## for loose ends lid marks all "unique" loose ends (which if loose.collapse = TRUE
    ## will be defined on the basis of coordinate overlap)
    gg$dt[tight == FALSE, .(cn = NA, snode.id, lambda, gid = index,
                            ulid = paste0(index, 'i'),
                            lid = ifelse(strand == '+', loose.left.id, paste0('-', loose.right.id)),
                            type = 'loose.in', vtype = 'I')], ## incoming loose ends
    gg$dt[tight == FALSE, .(cn = NA, snode.id, lambda, gid = index,
                            ulid = paste0(index, 'o'),
                            lid = ifelse(strand == '+', loose.right.id, paste0('-', loose.left.id)),
                            type = 'loose.out', vtype = 'I')], ## outgoing loose ends
    gg$dt[tight == FALSE, .(gid = index, cn, weight, type = 'nresidual', vtype = 'C')], ## node residual 
    gg$sedgesdt[, .(gid = sedge.id, cn, weight, type = 'eresidual', vtype = 'C')], ## edge residual 
    fill = TRUE)

  if (L0)
  {
    ## loose ends are labeled with lid and ulid, lid is only relevant if loose.collapse is true
    ## (i.e. we need indicator.sum and indicator.sum.indicator
    vars = rbind(vars, 
                 rbind( 
                   vars[type == 'loose.in', ][ , type := 'loose.in.indicator'][, vtype := 'B'][, gid := lid],
                   vars[type == 'loose.out', ][ , type := 'loose.out.indicator'][, vtype := 'B'][, gid := lid]
                 ))

    if (loose.collapse)
    {
      ## sum will sum all the loose ends assocaited with the same lid
      vars = rbind(vars, 
                   unique(rbind( 
                     vars[type == 'loose.in', ][ , type := 'loose.in.indicator.sum'][, vtype := 'I'][, gid := lid],
                     vars[type == 'loose.out', ][ , type := 'loose.out.indicator.sum'][, vtype := 'I'][, gid := lid]
                   ), by = 'gid'))
      
      ## sum.indicator is an binary indicator on the sum
      vars = rbind(vars, 
                   rbind( 
                     vars[type == 'loose.in.indicator.sum', ][ , type := 'loose.in.indicator.sum.indicator'][, vtype := 'B'][, gid := lid],
                     vars[type == 'loose.out.indicator.sum', ][ , type := 'loose.out.indicator.sum.indicator'][, vtype := 'B'][, gid := lid]
                   ))        
    }        
  }
  
    if (!is.null(marginal)) {
        ## first disjoin marginal against the nodes
        ## ie wee ned to create a separate residual variable for every unique
        ## disjoint overlap of marginal with the nodes
        dmarginal = gg$nodes$gr %>% gr.stripstrand %*% grbind(marginal %>% gr.stripstrand) %>%
            disjoin %$% marginal[, c('cn', 'weight', 'fix')] %Q%
            (!is.na(cn)) %Q% (!is.na(weight)) %Q% (!is.infinite(weight))

        vars = rbind(vars,
                     gr2dt(dmarginal)[, .(cn, weight, mfix = fix>0,
                                          rid = 1:.N, type = 'mresidual', vtype = 'C')][, gid := rid],
                     fill = TRUE
                     )
    }

    if (!is.null(emarginal)) {
        ## we need to identify which junction in the marginal each junction in the phased graph corresponds to
        junction.map = merge.Junction(
            phased = gg$junctions[, c()],
            emarginal = emarginal[, c("cn", "weight", "fix")],
            cartesian = TRUE,
            all.x = TRUE)$dt
        ## match this back with edge id and add this to vars
        vars[type == "edge", emarginal.id := junction.map[abs(sedge.id), seen.by.emarginal]]
        ## add weight and target total CN
        emtch = match(emarginal.id, junction.map$seen.by.emarginal)
        emarginal = unique(
            vars[type == "edge",][, type := "emresidual"][, cn := junction.map$cn[emtch]][, weight := junction.map$weight[emtch]][, fix := junction.map$fix[emtch]], ## lol change to merge
            by = "emarginal.id")
        vars = rbind(vars, emarginal, emresidual, fill = TRUE)
    }

    if (lp) {
        ## need delta plus and delta minus for nodes and edges
        delta.node = gg$dt[tight == FALSE, .(gid = index, cn, weight, vtype = 'C')] ## node residual 
        delta.edge = gg$sedgesdt[, .(gid = sedge.id, cn, weight, vtype = 'C')] ## edge residual 

        deltas = rbind(
            delta.node[, .(gid, weight, vtype, type = "ndelta.plus")],
            delta.node[, .(gid, weight, vtype, type = "ndelta.minus")],
            delta.edge[, .(gid, weight, vtype, type = "edelta.plus")],
            delta.edge[, .(gid, weight, vtype, type = "edelta.minus")]
        )

        deltas[, lb := 0] ## must be greater than zero

        vars = rbind(
            vars,
            deltas,
            fill = TRUE
        )

        ## add deltas for marginals if marginals are supplied
        if (!is.null(marginal)) {
            mdeltas = rbind(
                vars[type == "mresidual", .(rid, weight, vtype, type = "mdelta.plus")][, gid := rid],
                vars[type == "mresidual", .(rid, weight, vtype, type = "mdelta.minus")][, gid := rid]
            )
            vars = rbind(vars, mdeltas, fill = TRUE)
        }

        ## add deltas for emresiduals if emarginals are supplied
        if (!is.null(emarginal)) {
            emdeltas = rbind(
                vars[type == "emresidual", .(emarginal.id, weight, type = "emdelta.plus")][, gid := emarginal.id],
                vars[type == "emresidual", .(emarginal.id, weight, type = "emdelta.minus")][, gid := emarginal.id]
            )
            vars = rbind(vars, emdeltas, fill = TRUE)
        }
    }

    if (phased) {
        ## add allele information and og.node.id
        node.match = match(vars[, snode.id], gg$dt$snode.id)
        vars[, ":="(allele = gg$dt$allele[node.match],
                    og.node.id = gg$dt$og.node.id[node.match])]

        ## add ref/alt information and og.edge.id
        edge.match = match(vars[, sedge.id], gg$sedgesdt$sedge.id)
        vars[, ":="(ref.or.alt = gg$sedgesdt$type[edge.match], ## need type info but rename column...
                    og.edge.id = gg$sedgesdt$og.edge.id[edge.match])]

        edge.indicator.vars = vars[type == "edge"][, type := "edge.indicator"][, vtype := "B"][, gid := sedge.id]
        vars = rbind(vars, edge.indicator.vars, fill = TRUE)
    }

    if (ism) {
        ## if not phased, must add edge indicators (for just the ALT edges)
        if (!phased) {
            edge.match = match(vars[, sedge.id], gg$sedgesdt$sedge.id)
            vars[, ":="(ref.or.alt = gg$sedgesdt$type[edge.match])] ## need ref.or.alt information
            edge.indicator.vars = vars[type == "edge" & ref.or.alt == "ALT"][, type := "edge.indicator"][, vtype := "B"][, gid := sedge.id]
            vars = rbind(vars, edge.indicator.vars, fill = TRUE)
        }

        vars[type == "loose.in.indicator" & sign(snode.id) == 1, ee.id := paste(snode.id, "left")]
        vars[type == "loose.out.indicator" & sign(snode.id) == 1, ee.id := paste(snode.id, "right")]

        vars[type == "edge.indicator" & sign(sedge.id) == 1 & ref.or.alt == "ALT",
             ":="(ee.id.n1 = paste(gg$edges$dt$n1[match(sedge.id, gg$edges$dt$sedge.id)],
                                   gg$edges$dt$n1.side[match(sedge.id, gg$edges$dt$sedge.id)]),
                  ee.id.n2 = paste(gg$edges$dt$n2[match(sedge.id, gg$edges$dt$sedge.id)],
                                   gg$edges$dt$n2.side[match(sedge.id, gg$edges$dt$sedge.id)]))]

        if (phased) {
            ## homologous extremity exclusivity (only for phased graphs)
            ## get stranded breakpoint ID's associated with the start and end of each node

            ## number of unique starts should be equal to number of snodes in the original unphased graph
            ## aka 2 * number of og edge ids
            
            vars[type == "loose.in.indicator", hee.id := paste(og.node.id, "in")]
            vars[type == "loose.out.indicator", hee.id := paste(og.node.id, "out")]

            vars[type == "edge.indicator" & ref.or.alt == "ALT" & sign(sedge.id) == 1,
                 ":="(og.n1 = gg$dt$og.node.id[from],
                      og.n1.side = gg$edges$dt$n1.side[match(abs(sedge.id), gg$edges$dt$edge.id)],
                      og.n2 = gg$dt$og.node.id[to],
                      og.n2.side = gg$edges$dt$n2.side[match(abs(sedge.id), gg$edges$dt$edge.id)])]

            vars[type == "edge.indicator" & ref.or.alt == "ALT" & sign(sedge.id) == 1,
                 ":="(hee.id.n1 = ifelse(og.n1.side == "left",
                                         paste(og.n1, "in"),
                                         paste(og.n1, "out")),
                      hee.id.n2 = ifelse(og.n2.side == "left",
                                         paste(og.n2, "in"),
                                         paste(og.n2, "out")))]

            ## reciprocal homologous extremity exclusivity
            ## implement config indicators. there is one per og.edge.id per configuration
            straight.config = unique(vars[type == "edge.indicator" & ref.or.alt == "REF" & sedge.id > 0, ][, type := "straight.config"][, config.id := paste("straight", og.edge.id)], by = "og.edge.id")
            cross.config = unique(vars[type == "edge.indicator" & ref.or.alt == "REF" & sedge.id > 0, ][, type := "cross.config"][, config.id := paste("cross", og.edge.id)], by = "og.edge.id")

            ## add straight/cross to REF edges
            vars[type == "edge.indicator" & ref.or.alt == "REF",
                 connection := gg$sedgesdt$connection[match(sedge.id, gg$sedgesdt$sedge.id)]]
            
            ## add config ID's to corresponding edge indicators
            vars[type == "edge.indicator" & ref.or.alt == "REF" & sedge.id > 0,
                 config.id := paste(connection, og.edge.id)]

            vars = rbind(vars, straight.config, cross.config, fill = TRUE)

            ## add straight edge id e.g. for each n1 and n2, add the sedge.id of the corresponding straight edge
            straight.sedges = gg$edges$dt[type == "REF" & connection == "straight" & sedge.id > 0,
                                          .(n1.full = paste(n1, n1.side), n2.full = paste(n2, n2.side), sedge.id)]
            cross.sedges = gg$edges$dt[type == "REF" & connection == "cross" & sedge.id > 0,
                                       .(n1.full = paste(n1, n1.side), n2.full = paste(n2, n2.side), sedge.id)]


            ## pull alt edges from sedgesdt
            alt.sedges = gg$edges$dt[type == "ALT" & sedge.id > 0,
                                     .(n1.full = paste(n1, n1.side), n2.full = paste(n2, n2.side), sedge.id)]

            alt.sedges[, ":="(s1 = straight.sedges$sedge.id[match(n1.full, straight.sedges$n1.full)],
                              s2 = straight.sedges$sedge.id[match(n2.full, straight.sedges$n1.full)],
                              s3 = straight.sedges$sedge.id[match(n1.full, straight.sedges$n2.full)],
                              s4 = straight.sedges$sedge.id[match(n2.full, straight.sedges$n2.full)])]

            alt.sedges[, ":="(c1 = cross.sedges$sedge.id[match(n1.full, cross.sedges$n1.full)],
                              c2 = cross.sedges$sedge.id[match(n2.full, cross.sedges$n1.full)],
                              c3 = cross.sedges$sedge.id[match(n1.full, cross.sedges$n2.full)],
                              c4 = cross.sedges$sedge.id[match(n2.full, cross.sedges$n2.full)])]

            ## pull loose ends
            vars[type == "loose.in.indicator" & snode.id > 0, n2.full := paste(snode.id, "left")]
            vars[type == "loose.out.indicator" & snode.id > 0, n1.full := paste(snode.id, "right")]

            ## merge sedge.id
            vars[type == "loose.in.indicator" & snode.id > 0, ":="(s = straight.sedges$sedge.id[match(n2.full, straight.sedges$n2.full)])]
            vars[type == "loose.out.indicator" & snode.id > 0, ":="(s = straight.sedges$sedge.id[match(n1.full, straight.sedges$n1.full)])]

            vars[type == "loose.in.indicator" & snode.id > 0, ":="(c = cross.sedges$sedge.id[match(n2.full, cross.sedges$n2.full)])]
            vars[type == "loose.out.indicator" & snode.id > 0, ":="(c = cross.sedges$sedge.id[match(n1.full, cross.sedges$n1.full)])]

            ## merge this info into vars
            vars = merge(vars,
                         alt.sedges[, .(sedge.id, s1, s2, s3, s4, c1, c2, c3, c4)],
                         by = "sedge.id",
                         all.x = TRUE,
                         all.y = FALSE)

        }
    }
    
    vars[, id := 1:.N] ## set id in the optimization
    vars[is.na(lb), lb := -Inf]
    vars[is.na(ub), ub := Inf]
    vars[, relax := FALSE][, fix := FALSE]
    if ("mresidual" %in% vars$type) {
        vars[type == 'mresidual' & mfix == TRUE, ":="(lb = 0, ub = 0)]
        message("Number of fixed marginals: ", nrow(vars[type == 'mresidual' & mfix == TRUE,]))
    }
    if ("emresidual" %in% vars$type) {
        vars[type == "emresidual" & fix == TRUE, ":="(lb = 0, ub = 0)]
    }
    vars[type %in% c('node', 'edge'), lb := pmax(lb, 0, na.rm = TRUE)]
    vars[type %in% c('node', 'edge'), ub := ifelse(is.na(ub), M, pmax(ub, M, na.rm = TRUE))]
    vars[type %in% c('loose.in', 'loose.out'), ":="(lb = 0, ub = Inf)]
  
    vars[type %in% c('edge'), reward := pmax(reward, 0, na.rm = TRUE)]


    ## figure out junctions and nodes to fix

    vars[!is.na(cn) & type == 'node' & abs(snode.id) %in% nfix, ":="(lb = cn, ub = cn, fix = TRUE)]
    vars[!is.na(cn) & type == 'edge' & abs(sedge.id) %in% efix, ":="(lb = cn, ub = cn, fix = TRUE)]

    ## figure out terminal node sides for in and out loose ends
    ## these will not have loose ends penalized
    qtips = gr.end(si2gr(seqlengths(gg$nodes))) ## location of q arm tips
    term.in = c(which(start(gg$nodes$gr) == 1), ## beginning of chromosome
                -which(gg$nodes$gr %^% qtips)) ## flip side of chromosome end
    term.out = -term.in
    vars$terminal = FALSE
    vars[(type %in% c('loose.in', 'loose.in.indicator')) & (snode.id %in% term.in), terminal := TRUE]
    vars[(type %in% c('loose.out', 'loose.out.indicator')) & (snode.id %in% term.out), terminal := TRUE]

  ########
  ## CONSTRAINTS
  ## the key principle behind this "melted" form of constraint building is the cid 
  ## (constraint id) which is the key that will group coefficients into constraints
  ## when we finally build the matrices.  So all we need to do is make sure that
  ## that value / cid pairs make sense and that every cid has an entry in b
  ########

  ## we need one junction balance constraint per loose end

  ## constraints indexed by cid
  constraints = rbind(
    vars[type == 'loose.in', .(value = 1, id, cid = paste('in', gid))],
    vars[type == 'edge', .(value = 1, id, cid = paste('in', to))],
    vars[type == 'node', .(value = -1, id, cid = paste('in', gid))],
    vars[type == 'loose.out', .(value = 1, id, cid = paste('out', gid))],
    vars[type == 'edge', .(value = 1, id, cid = paste('out', from))],
    vars[type == 'node', .(value = -1, id, cid = paste('out', gid))],
   fill = TRUE)

  b = rbind(
    vars[type == 'node', .(value = 0, sense = 'E', cid = paste('in', gid))],
    vars[type == 'node', .(value = 0, sense = 'E', cid = paste('out', gid))],
    fill = TRUE)

  ## add to the constraints the definitions of the node and edge
  ## residuals
  constraints = rbind(
    constraints,
    rbind(
      vars[type == 'node', .(value = 1, id, cid = paste('nresidual', gid))],
      vars[type == 'nresidual', .(value = -1, id, cid = paste('nresidual', gid))],
      vars[type == 'edge', .(value = 1, id, cid = paste('eresidual', gid))],
      vars[type == 'eresidual', .(value = -1, id, cid = paste('eresidual', gid))],
      fill = TRUE)
  )

  b = rbind(b,
            vars[type == 'node', .(value = cn, sense = 'E', cid = paste('nresidual', gid))],
            vars[type == 'edge', .(value = cn, sense = 'E', cid = paste('eresidual', gid))],
            fill = TRUE)
  
  ## add the reverse complement equality constraints on nodes and edges
  constraints = rbind(
    constraints,
    rbind( ## +1 coefficient for positive nodes, -1 for negative nodes, matched by abs (snode.id)
      vars[type == 'node', .(value = sign(snode.id), id, cid = paste('nrc', abs(snode.id)))],
      vars[type == 'edge', .(value = sign(sedge.id), id, cid = paste('erc', abs(sedge.id)))],
      fill = TRUE)
  )
  
  b = rbind(b,
            vars[type == 'node' & snode.id>0, .(value = 0, sense = 'E', cid = paste('nrc', abs(snode.id)))],
            vars[type == 'edge' & sedge.id>0, .(value = 0, sense = 'E', cid = paste('erc', abs(sedge.id)))],
            fill = TRUE)

  ## if solving as LP, add deltas constraints (absolute value trick)

  if (lp) {
    if (verbose) {
      message("adding delta constraints for LP")
    }

    ## ## constrain deltas to be at least zero
    ## delta.lbs = rbind(
    ##   vars[type == "ndelta.minus", .(value = 1, id, cid = paste("ndelta.minus.lb", gid))],
    ##   vars[type == "ndelta.plus", .(value = 1, id, cid = paste("ndelta.plus.lb", gid))],
    ##   vars[type == "edelta.minus", .(value = 1, id, cid = paste("edelta.minus.lb", gid))],
    ##   vars[type == "edelta.plus", .(value = 1, id, cid = paste("edelta.plus.lb", gid))],
    ##   vars[type == "mdelta.minus", .(value = 1, id, cid = paste("mdelta.minus.lb", gid))],
    ##   vars[type == "mdelta.plus", .(value = 1, id, cid = paste("mdelta.plus.lb", gid))]
    ## )

    ## delta.lbs.rhs = rbind(
    ##   vars[type == "ndelta.minus", .(value = 0, sense = "G", cid = paste("ndelta.minus.lb", gid))],
    ##   vars[type == "ndelta.plus", .(value = 0, sense = "G", cid = paste("ndelta.plus.lb", gid))],
    ##   vars[type == "edelta.minus", .(value = 0, sense = "G", cid = paste("edelta.minus.lb", gid))],
    ##   vars[type == "edelta.plus", .(value = 0, sense = "G", cid = paste("edelta.plus.lb", gid))],
    ##   vars[type == "mdelta.minus", .(value = 0, sense = "G", cid = paste("mdelta.minus.lb", gid))],
    ##   vars[type == "mdelta.plus", .(value = 0, sense = "G", cid = paste("mdelta.plus.lb", gid))]
    ## )


    ## constraints = rbind(constraints, delta.lbs, fill = TRUE)
    ## b = rbind(b, delta.lbs.rhs, fill = TRUE)

    ## ## add upper bound to prevent problem from becoming unbounded
    ## ## constrain deltas to be at least zero
    ## delta.ubs = rbind(
    ##   vars[type == "ndelta.minus", .(value = 1, id, cid = paste("ndelta.minus.ub", gid))],
    ##   vars[type == "ndelta.plus", .(value = 1, id, cid = paste("ndelta.plus.ub", gid))],
    ##   vars[type == "edelta.minus", .(value = 1, id, cid = paste("edelta.minus.ub", gid))],
    ##   vars[type == "edelta.plus", .(value = 1, id, cid = paste("edelta.plus.ub", gid))],
    ##   vars[type == "mdelta.minus", .(value = 1, id, cid = paste("mdelta.minus.ub", gid))],
    ##   vars[type == "mdelta.plus", .(value = 1, id, cid = paste("mdelta.plus.ub", gid))]
    ## )

    ## delta.ubs.rhs = rbind(
    ##   vars[type == "ndelta.minus", .(value = M, sense = "L", cid = paste("ndelta.minus.ub", gid))],
    ##   vars[type == "ndelta.plus", .(value = M, sense = "L", cid = paste("ndelta.plus.ub", gid))],
    ##   vars[type == "edelta.minus", .(value = M, sense = "L", cid = paste("edelta.minus.ub", gid))],
    ##   vars[type == "edelta.plus", .(value = M, sense = "L", cid = paste("edelta.plus.ub", gid))],
    ##   vars[type == "mdelta.minus", .(value = M, sense = "L", cid = paste("mdelta.minus.ub", gid))],
    ##   vars[type == "mdelta.plus", .(value = M, sense = "L", cid = paste("mdelta.plus.ub", gid))]
    ## )

    ## constraints = rbind(constraints, delta.ubs, fill = TRUE)
    ## b = rbind(b, delta.ubs.rhs, fill = TRUE)

    vars[type %like% "delta.plus" | type %like% "delta.minus", ":="(ub = M, lb = 0)]

    ## add the residual constraints
    ## kind of gross code, should just write a function for this
    ndelta.slack = rbind(
      vars[type == "nresidual", .(value = -1, id, cid = paste("ndelta.minus.slack", gid))],
      vars[type == "ndelta.minus", .(value = -1, id, cid = paste("ndelta.minus.slack", gid))],
      vars[type == "nresidual", .(value = 1, id, cid = paste("ndelta.plus.slack", gid))],
      vars[type == "ndelta.plus", .(value = -1, id, cid = paste("ndelta.plus.slack", gid))]
    )

    ndelta.slack.rhs = rbind(
      vars[type == "ndelta.minus", .(value = 0, sense = "L", cid = paste("ndelta.minus.slack", gid))],
      vars[type == "ndelta.plus", .(value = 0, sense = "L", cid = paste("ndelta.plus.slack", gid))]
    )

    edelta.slack = rbind(
      vars[type == "eresidual", .(value = -1, id, cid = paste("edelta.minus.slack", gid))],
      vars[type == "edelta.minus", .(value = -1, id, cid = paste("edelta.minus.slack", gid))],
      vars[type == "eresidual", .(value = 1, id, cid = paste("edelta.plus.slack", gid))],
      vars[type == "edelta.plus", .(value = -1, id, cid = paste("edelta.plus.slack", gid))]
    )

    edelta.slack.rhs = rbind(
      vars[type == "edelta.minus", .(value = 0, sense = "L", cid = paste("edelta.minus.slack", gid))],
      vars[type == "edelta.plus", .(value = 0, sense = "L", cid = paste("edelta.plus.slack", gid))]
    )

    mdelta.slack = rbind(
      vars[type == "mresidual", .(value = -1, id, cid = paste("mdelta.minus.slack", gid))],
      vars[type == "mdelta.minus", .(value = -1, id, cid = paste("mdelta.minus.slack", gid))],
      vars[type == "mresidual", .(value = 1, id, cid = paste("mdelta.plus.slack", gid))],
      vars[type == "mdelta.plus", .(value = -1, id, cid = paste("mdelta.plus.slack", gid))]
    )

    mdelta.slack.rhs = rbind(
      vars[type == "mdelta.minus", .(value = 0, sense = "L", cid = paste("mdelta.minus.slack", gid))],
      vars[type == "mdelta.plus", .(value = 0, sense = "L", cid = paste("mdelta.plus.slack", gid))]
    )

    constraints = rbind(constraints, ndelta.slack, edelta.slack, mdelta.slack, fill = TRUE)
    b = rbind(b, ndelta.slack.rhs, edelta.slack.rhs, mdelta.slack.rhs, fill = TRUE)

    ## browser()

  }

    if (phased) {
        
        ## add constraints that force indicators to be 1 if edge CN > 0
        
        ## add constraints for upper bound (same setup as L0 penalty) - one per edge
        iconstraints = vars[type == "edge", .(value = 1, id,
                                              sedge.id, 
                                              cid = paste("edge.indicator.ub", sedge.id))]

        ## add matching indicator variables, matching by cid
        iconstraints = rbind(
            iconstraints,
            vars[type == "edge.indicator", ][
                sedge.id %in% iconstraints$sedge.id, .(value = -M, id, cid = iconstraints$cid, sedge.id)],
            fill = TRUE)

        ## upper bound is M if indicator is positive, and zero otherwise
        constraints = rbind(
            constraints,
            iconstraints,
            fill = TRUE)

        ## add the RHS of this constraint (upper bound)
        b = rbind(
            b,
            vars[type == "edge", .(value = 0, sense = "L", cid = paste("edge.indicator.ub", sedge.id))],
            fill = TRUE
        )

        ## add constraints for the lower bound
        iconstraints = vars[type == "edge",
                            .(value = 1, id, sedge.id, cid = paste("edge.indicator.lb", sedge.id))]

        ## add matching indicator variables for LB
        iconstraints = rbind(
            iconstraints,
            vars[type == "edge.indicator", ][sedge.id %in% iconstraints$sedge.id,
                                             .(value = -0.1, id, cid = iconstraints$cid, sedge.id)],
            fill = TRUE)

        constraints = rbind(
            constraints,
            iconstraints,
            fill = TRUE)

        ## add the RHS of this constraint (upper bound)
        b = rbind(
            b,
            vars[type == "edge", .(value = 0, sense = "G", cid = paste("edge.indicator.lb", sedge.id))],
            fill = TRUE
        )
    }

    if (ism) {

        ## implement edge edge indicators if not already (e.g. if not doing phasing)
        if (!phased) {

            ## importantly, we only want to add these for ALT edges
            iconstraints = vars[type == "edge" & ref.or.alt == "ALT" & sign(sedge.id) == 1,
                                .(value = 1, id,
                                  sedge.id, 
                                  cid = paste("edge.indicator.ub", sedge.id))]

            ## add matching indicator variables, matching by cid
            iconstraints = rbind(
                iconstraints,
                vars[type == "edge.indicator" & ref.or.alt == "ALT" & sign(sedge.id) == 1, ][
                    sedge.id %in% iconstraints$sedge.id,
                    .(value = -M, id, cid = iconstraints$cid, sedge.id)],
                fill = TRUE)

            ## upper bound is M if indicator is positive, and zero otherwise
            constraints = rbind(
                constraints,
                iconstraints,
                fill = TRUE)

            ## add the RHS of this constraint (upper bound)
            b = rbind(
                b,
                vars[type == "edge" & ref.or.alt == "ALT" & sign(sedge.id) == 1,
                     .(value = 0, sense = "L", cid = paste("edge.indicator.ub", sedge.id))],
                fill = TRUE
            )

            ## add constraints for the lower bound
            iconstraints = vars[type == "edge" & ref.or.alt == "ALT" & sign(sedge.id) == 1,
                                .(value = 1, id, sedge.id, cid = paste("edge.indicator.lb", sedge.id))]

            ## add matching indicator variables for LB
            iconstraints = rbind(
                iconstraints,
                vars[type == "edge.indicator" & ref.or.alt == "ALT" & sign(sedge.id) == 1, ][
                    sedge.id %in% iconstraints$sedge.id,
                    .(value = -0.1, id, cid = iconstraints$cid, sedge.id)],
                fill = TRUE)

            constraints = rbind(
                constraints,
                iconstraints,
                fill = TRUE)

            ## add the RHS of this constraint (upper bound)
            b = rbind(
                b,
                vars[type == "edge" & ref.or.alt == "ALT" & sign(sedge.id) == 1,
                     .(value = 0, sense = "G", cid = paste("edge.indicator.lb", sedge.id))],
                fill = TRUE
            )
        }

        ## fix loose ends at zero if there's a junction there (only valid if not phasing)
        if (!phased) {
            ## extremity exclusivity (relevant for ALL graphs)
            loose.constraints = rbind(
                vars[type == "loose.in.indicator" & sign(snode.id) == 1,
                     .(value = 1, id, cid = paste("extremity.exclusivity", ee.id))],
                vars[type == "loose.out.indicator" & sign(snode.id) == 1,
                     .(value = 1, id, cid = paste("extremity.exclusivity", ee.id))]
            )

            edge.constraints = rbind(
                vars[type == "edge.indicator" & ref.or.alt == "ALT" & sign(sedge.id) == 1,
                     .(value = 1, id, cid = paste("extremity.exclusivity", ee.id.n1))],
                vars[type == "edge.indicator" & ref.or.alt == "ALT" & sign(sedge.id) == 1,
                     .(value = 1, id, cid = paste("extremity.exclusivity", ee.id.n2))]
            )

            constraints = rbind(constraints, loose.constraints, edge.constraints, fill = TRUE)

            loose.b = unique(loose.constraints[, .(cid, value = 1, sense = "L")], by = "cid")
            edge.b = unique(edge.constraints[, .(cid, value = 1, sense = "L")], by = "cid")

            b = rbind(b, edge.b, loose.b, fill = TRUE)

            edge.ee.ids = unique(c(vars[type == "edge.indicator", ee.id.n1], vars[type == "edge.indicator", ee.id.n2]))
            edge.ee.ids = edge.ee.ids[!is.na(edge.ee.ids)]

            loose.zeros = rbind(
                vars[type == "loose.in.indicator" & sign(snode.id) == 1 & ee.id %in% edge.ee.ids,
                     .(value = 1, id, cid = paste("extremity.exclusivity", ee.id))],
                vars[type == "loose.out.indicator" & sign(snode.id) == 1 & ee.id %in% edge.ee.ids,
                     .(value = 1, id, cid = paste("extremity.exclusivity", ee.id))]
            )

            loose.zeros.rhs = unique(loose.zeros[, .(cid, value = 0, sense = "E")], by = "cid")

            constraints = rbind(constraints, loose.zeros, fill = TRUE)
            b = rbind(b, loose.zeros.rhs, fill = TRUE)
        }

        if (phased) {
            ## homologous extremity exclusivity
            loose.constraints = rbind(
                vars[type == "loose.in.indicator" & sign(snode.id)==1,
                     .(value = 1, id, cid = paste("homol.extremity.exclusivity", hee.id))],
                vars[type == "loose.out.indicator" & sign(snode.id)==1,
                     .(value = 1, id, cid = paste("homol.extremity.exclusivity", hee.id))]
            )

            edge.constraints = rbind(
                vars[type == "edge.indicator" & ref.or.alt == "ALT" & sign(sedge.id)==1,
                     .(value = 1, id, cid = paste("homol.extremity.exclusivity", hee.id.n1))],
                vars[type == "edge.indicator" & ref.or.alt == "ALT" & sign(sedge.id)==1,
                     .(value = 1, id, cid = paste("homol.extremity.exclusivity", hee.id.n2))]
            )

            constraints = rbind(constraints, loose.constraints, edge.constraints, fill = TRUE)

            rhs = unique(rbind(
                vars[type == "loose.in.indicator" & sign(snode.id)==1,
                     .(value = 1, sense = "L", cid = paste("homol.extremity.exclusivity", hee.id))],
                vars[type == "loose.out.indicator" & sign(snode.id)==1,
                     .(value = 1, sense = "L", cid = paste("homol.extremity.exclusivity", hee.id))]
            ), by = "cid")
            
            b = rbind(b, rhs, fill = TRUE)
        
            ## reciprocal homologous extremity exclusivity
            ## implement configuration indicators (OR constraint)
            config.dt = vars[type == "straight.config" | type == "cross.config",]
            
            config.constraints.lt = rbind(
                vars[type == "edge.indicator" & sedge.id > 0 & ref.or.alt == "REF",
                     .(value = -1, id, cid = paste("config lt", sedge.id))],
                vars[type == "edge.indicator" & sedge.id > 0 & ref.or.alt == "REF",
                     .(value = 1, id = config.dt$id[match(config.id, config.dt$config.id)],
                       cid = paste("config lt", sedge.id))])

            config.constraints.gt = rbind(
                vars[type == "edge.indicator" & sedge.id > 0 & ref.or.alt == "REF",
                     .(value = -1, id, cid = config.id)],
                vars[type == "straight.config" & sedge.id > 0 & ref.or.alt == "REF",
                     .(value = 1, id, cid = config.id)],
                vars[type == "cross.config" & sedge.id > 0 & ref.or.alt == "REF",
                     .(value = 1, id, cid = config.id)])

            rhs = unique(rbind(
                config.constraints.lt[, .(cid, value = 0, sense = "G")],
                config.constraints.gt[, .(cid, value = 0, sense = "L")]),
                by = "cid")

            constraints = rbind(constraints, config.constraints.lt, config.constraints.gt, fill = TRUE)
            b = rbind(b, rhs, fill = TRUE)

            ## implement reciprocal homologous extremity exclusivity
            straight.config.dt = vars[type == "straight.config",]
            cross.config.dt = vars[type == "cross.config",]
            
            rhomol.constraints = rbind(
                ## corresponding cross indicator
                vars[type == "edge.indicator" & sedge.id > 0 & ref.or.alt == "REF" & connection == "straight",
                     .(value = 1, id = cross.config.dt$id[match(og.edge.id, cross.config.dt$og.edge.id)],
                       cid = paste("rhee", sedge.id))],
                
                ## corresponding cross indicator
                vars[type == "edge.indicator" & sedge.id > 0 & ref.or.alt == "REF" & connection == "cross",
                     .(value = 1, id = straight.config.dt$id[match(og.edge.id, straight.config.dt$og.edge.id)],
                       cid = paste("rhee", sedge.id))],

                ## actual ALT edges
                vars[type == "edge.indicator" & sedge.id > 0 & ref.or.alt == "ALT" & !is.na(s1),
                     .(value = 1, id, cid = paste("rhee", s1))],
                vars[type == "edge.indicator" & sedge.id > 0 & ref.or.alt == "ALT" & !is.na(s2),
                     .(value = 1, id, cid = paste("rhee", s2))],
                vars[type == "edge.indicator" & sedge.id > 0 & ref.or.alt == "ALT" & !is.na(s3),
                     .(value = 1, id, cid = paste("rhee", s3))],
                vars[type == "edge.indicator" & sedge.id > 0 & ref.or.alt == "ALT" & !is.na(s4),
                     .(value = 1, id, cid = paste("rhee", s4))],
                vars[type == "edge.indicator" & sedge.id > 0 & ref.or.alt == "ALT" & !is.na(c1),
                     .(value = 1, id, cid = paste("rhee", c1))],
                vars[type == "edge.indicator" & sedge.id > 0 & ref.or.alt == "ALT" & !is.na(c2),
                     .(value = 1, id, cid = paste("rhee", c2))],
                vars[type == "edge.indicator" & sedge.id > 0 & ref.or.alt == "ALT" & !is.na(c3),
                     .(value = 1, id, cid = paste("rhee", c3))],
                vars[type == "edge.indicator" & sedge.id > 0 & ref.or.alt == "ALT" & !is.na(c4),
                     .(value = 1, id, cid = paste("rhee", c4))],

                ## loose indicators
                vars[type == "loose.in.indicator" & snode.id > 0 & !is.na(s),
                     .(value = 1, id, cid = paste("rhee", s))],
                vars[type == "loose.in.indicator" & snode.id > 0 & !is.na(c),
                     .(value = 1, id, cid = paste("rhee", c))],
                vars[type == "loose.out.indicator" & snode.id > 0 & !is.na(s),
                     .(value = 1, id, cid = paste("rhee", s))],
                vars[type == "loose.out.indicator" & snode.id > 0 & !is.na(c),
                     .(value = 1, id, cid = paste("rhee", c))]
                )

            rhs = unique(rhomol.constraints[, .(value = 2, sense = "L", cid)], by = "cid")

            constraints = rbind(constraints, rhomol.constraints, fill = TRUE)
            b = rbind(b, rhs, fill = TRUE)

        }
    }

    if (phased) {

        ## add the edge indicator sum constraints (ISM consistency)
        iconstraints = unique(
            vars[type == "edge.indicator" & ref.or.alt == "ALT",
                 .(value = 1, id, og.edge.id,
                   edge.id = abs(sedge.id),
                   cid = paste("edge.indicator.sum.ub", og.edge.id))],
            by = "edge.id"
        )

        constraints = rbind(
            constraints,
            iconstraints[, .(value, id, cid)],
            fill = TRUE)

        edge.indicator.b = unique(
            vars[type == "edge.indicator" & ref.or.alt == "ALT",
                 .(value = 1, sense = "L", cid = paste("edge.indicator.sum.ub", og.edge.id))],
            by = "cid"
        )

        b = rbind(b, edge.indicator.b, fill = TRUE)

        ## force nonzero CN for ALT edges (because these have nonzero CN in original JaBbA output)
        ## can become infeasible ...
        iconstraints = unique(
            vars[type == "edge.indicator" & ref.or.alt == "ALT",
                 .(value = 1, id, og.edge.id,
                   edge.id = abs(sedge.id),
                   cid = paste("edge.indicator.sum.lb", og.edge.id))],
            by = "edge.id"
        )

        constraints = rbind(
            constraints,
            iconstraints[, .(value, id, cid)],
            fill = TRUE)

        edge.indicator.b = unique(
            vars[type == "edge.indicator" & ref.or.alt == "ALT",
                 .(value = 1, sense = "G", cid = paste("edge.indicator.sum.lb", og.edge.id))],
            by = "cid"
        )

        b = rbind(b, edge.indicator.b, fill = TRUE)

        ## REF edge configuration constraint (added by default basically)
        ## only add this if there are no unphased nodes
        iconstraints.from = unique(
            vars[type == "edge.indicator" & ref.or.alt == "REF",
                 .(value = 1, id,
                   edge.id = abs(sedge.id),
                   snode.id = from, ## this is actually a misleading name because from is the row in gg$dt
                   cid = paste("ref.configuration.constraint.from", from))],
            by = "edge.id"
        )

        iconstraints.to = unique(
            vars[type == "edge.indicator" & ref.or.alt == "REF",
                 .(value = 1, id,
                   edge.id = abs(sedge.id),
                   snode.id = to,
                   cid = paste("ref.configuration.constraint.to", to))],
            by = "edge.id"
        )

        iconstraints = rbind(iconstraints.from, iconstraints.to)

        ## sum to at most 1 if phased, unconstrained if unphased
        iconstraints[, ":="(allele = gg$dt$allele[iconstraints$snode.id])]
        
        edge.indicator.b = unique(iconstraints[allele %in% c("major", "minor"),
                                               .(value = 1, sense = "L", cid)],
                                  by = "cid")
## rbind(
##             unique(iconstraints[allele %in% c("major", "minor"),
##                                 .(value = 1, sense = "L", cid)], by = "cid"),
##             unique(iconstraints[!(allele %in% c("major", "minor")),
##                                 .(value = 2, sense = "L", cid)], by = "cid")
##         )

        constraints = rbind(
            constraints,
            iconstraints[allele %in% c("major", "minor"),
                         .(value, id, cid)],
            fill = TRUE)
        
        ## add to b
        b = rbind(b, edge.indicator.b, fill = TRUE)
    }
  
  if (L0) ## add "big M" constraints
  {
    ## indicator constraints ie on ulids 
    iconstraints = rbind(
      vars[type == 'loose.out', .(value = 1, id, ulid, cid = paste('loose.out.indicator.ub', ulid))],
      vars[type == 'loose.in', .(value = 1, id, ulid, cid = paste('loose.in.indicator.ub', ulid))],
      fill = TRUE)

    ## add the matching indicator variables, matching to the cid from above
    iconstraints = rbind(
      iconstraints,
      vars[type %in% c('loose.out.indicator', 'loose.in.indicator'), ][
        match(iconstraints$ulid, ulid), .(value = -M, id, cid = iconstraints$cid)],
      fill = TRUE)
                         
    ## upper bounds "infinity" ie M if indicator positive, 0 otherwise
    constraints = rbind(
      constraints,
      iconstraints,
      fill = TRUE)

    ## upper bound sense is 'L' i.e. less than because -M on left hand side
    b = rbind(b,
              vars[type == 'loose.in', .(value = 0, sense = 'L', cid = paste('loose.in.indicator.ub', ulid))],
              vars[type == 'loose.out', .(value = 0, sense = 'L', cid = paste('loose.out.indicator.ub', ulid))],
              fill = TRUE)

    ## lower bound 0.1 if indicator positive, 0 otherwise
    iconstraints = rbind(
      vars[type == 'loose.out', .(value = 1, id, ulid, cid = paste('loose.out.indicator.lb', ulid))],
      vars[type == 'loose.in', .(value = 1, id, ulid, cid = paste('loose.in.indicator.lb', ulid))],
      fill = TRUE)
    
    ## add the matching indicator variables, matching to the cid from above
    iconstraints = rbind(
      iconstraints,
      vars[type %in% c('loose.out.indicator', 'loose.in.indicator'), ][
        match(iconstraints$ulid, ulid), .(value = -.1, id, cid = iconstraints$cid)],
      fill = TRUE)

    ## upper bounds "infinity" ie M if indicator positive, 0 otherwise
    constraints = rbind(
      constraints,
      iconstraints,
      fill = TRUE)

    ## lower bound sense is 'G' i.e. greater than because -M on left hand side
    b = rbind(b,
              vars[type == 'loose.in', .(value = 0, sense = 'G', cid = paste('loose.in.indicator.lb', ulid))],
              vars[type == 'loose.out', .(value = 0, sense = 'G', cid = paste('loose.out.indicator.lb', ulid))],
              fill = TRUE)

    if (loose.collapse)
    {
      ##################
      ## loose indicator sum  = sum of indicators
      ##################
      iconstraints = rbind(
        vars[type == 'loose.out.indicator', .(value = 1, id, lid, cid = paste('loose.out.indicator.sum', lid))],
        vars[type == 'loose.in.indicator', .(value = 1, id, lid, cid = paste('loose.in.indicator.sum', lid))],
        fill = TRUE)
      
      ## indicator sum is the sum of all indicators mapping to that loose end
      iconstraints = rbind(
          iconstraints,
        unique(vars[type %in% c('loose.out.indicator.sum', 'loose.in.indicator.sum'), ][
          match(iconstraints$lid, lid), .(value = -1, id, lid, cid = iconstraints$cid)], by = 'lid'),
        fill = TRUE)
      
      constraints = rbind(
        constraints,
        iconstraints,
        fill = TRUE)
      
      b = rbind(b,
                vars[type == 'loose.in.indicator.sum', .(value = 0, sense = 'E', cid = paste('loose.in.indicator.sum', lid))],
                vars[type == 'loose.out.indicator.sum', .(value = 0, sense = 'E', cid = paste('loose.out.indicator.sum', lid))],
                fill = TRUE)
     
      ##################
      ## now we make new indicator variables on the sum of the individual loose end indicators      
      ## upper bound bound 0.1 if indicator positive, 0 otherwise
      ##################

      iconstraints = rbind(
        vars[type == 'loose.out.indicator.sum', .(value = 1, id, lid, cid = paste('loose.out.indicator.sum.indicator.ub', lid))],
        vars[type == 'loose.in.indicator.sum', .(value = 1, id, lid, cid = paste('loose.in.indicator.sum.indicator.ub', lid))],
        fill = TRUE)

      ## add the matching indicator variables, matching to the cid from above
      iconstraints = rbind(
        iconstraints,
        vars[type %in% c('loose.out.indicator.sum.indicator', 'loose.in.indicator.sum.indicator'), ][
          match(iconstraints$lid, lid), .(value = -M, id, lid, cid = iconstraints$cid)],
        fill = TRUE)
      
      ## upper bounds "infinity" ie M if indicator positive, 0 otherwise
      constraints = rbind(
        constraints,
        iconstraints,
        fill = TRUE)

      ## upper bound sense is 'L' i.e. less than because -M on left hand side
      b = rbind(b,
                vars[type == 'loose.in.indicator.sum', .(value = 0, sense = 'L', cid = paste('loose.in.indicator.sum.indicator.ub', lid))],
                vars[type == 'loose.out.indicator.sum', .(value = 0, sense = 'L', cid = paste('loose.out.indicator.sum.indicator.ub', lid))],
                fill = TRUE)

      ## lower bound 0.1 if indicator positive, 0 otherwise
      iconstraints = rbind(
        vars[type == 'loose.out.indicator.sum', .(value = 1, id, lid, cid = paste('loose.out.indicator.sum.indicator.lb', lid))],
        vars[type == 'loose.in.indicator.sum', .(value = 1, id, lid, cid = paste('loose.in.indicator.sum.indicator.lb', lid))],
        fill = TRUE)
      
      ## add the matching indicator variables, matching to the cid from above
      iconstraints = rbind(
        iconstraints,
        vars[type %in% c('loose.out.indicator.sum', 'loose.in.indicator.sum'), ][
          match(iconstraints$lid, lid), .(value = -.1, id, lid, cid = iconstraints$cid)],
        fill = TRUE)

      ## upper bounds "infinity" ie M if indicator positive, 0 otherwise
      constraints = rbind(
        constraints,
        iconstraints,
        fill = TRUE)

      ## lower bound sense is 'G' i.e. greater than because -M on left hand side
      b = rbind(b,
                vars[type == 'loose.in.indicator.sum', .(value = 0, sense = 'G', cid = paste('loose.in.indicator.sum.indicator.lb', lid))],
                vars[type == 'loose.out.indicator.sum', .(value = 0, sense = 'G', cid = paste('loose.out.indicator.sum.indicator.lb', lid))],
                fill = TRUE)

    }    
  }

    if (!is.null(marginal) && length(dmarginal)) 
    {
        ## match against nodes and store query.id as rid
        ## this will be the constraint id that will allow us
        ## to sum the appropriate nodes to constrain to the residual
        ov = dmarginal[, c('cn', 'weight')] %*% gg$nodes$gr %>% gr2dt

        ov[, rid := query.id]

        constraints = rbind(
            constraints,
            rbind(
                ## match up vars and marginal by snode.id and populate coefficients
                merge(vars[type == 'node', !"rid"], ov, by = 'snode.id')[, .(value = 1, id , cid = paste('mresidual', rid))],
                ## the residual is the difference between the sum and marginal cn
                vars[type == 'mresidual' & rid %in% ov$rid, .(value = -1, id, cid = paste('mresidual', rid))],        
                fill = TRUE),
            fill = TRUE
        )

        b = rbind(b,
                  vars[type == 'mresidual' & rid %in% ov$rid, .(value = cn, sense = 'E', cid = paste('mresidual', rid))],
                  fill = TRUE)
    }

    if (!is.null(emarginal)) {

        emconstraints = rbind(
            vars[type == "edge", .(value = 1, id, cid = paste("emresidual", emarginal.id))],
            vars[type == "emresidual", .(value = -1, id, cid = paste("emresidual", emarginal.id))]
        )

        constraints = rbind(constraints, emconstraints, fill = TRUE)

        emb = vars[type == "emresidual", .(value = cn, sense = "E", cid = paste("emresidual", emarginal.id))]

        b = rbind(emb, b, fill = TRUE)
    }

  ########
  ## MAKE MATRICES
  ########

  ## now Rcplex time
  ## remove any rows with b = NA

    ## get rid of any constraints with NA values
    keep.constraints = intersect(b[!is.na(value), cid], constraints[!is.na(value), cid])
    b = b[cid %in% keep.constraints,]
    constraints = constraints[cid %in% keep.constraints,]

    ## convert constraints to integers
    ucid = unique(b$cid)
    b[, cid.char := cid]
    b[, cid := cid %>% factor(ucid) %>% as.integer]
    constraints[, cid.char := cid]
    constraints[, cid := cid %>% factor(ucid) %>% as.integer]

    pmt = match(ucid, b$cid.char) ## get right permutation
    bvec = b[pmt, value]
    sense = b[pmt, sense]
    if (verbose) {
        message("Unique cids (A): ", length(unique(constraints$cid)))
        message("Unique cids (b): ", length(unique(b$cid)))
        message("Number of variables: ", length(unique(constraints$id)))
    }

    ## create constraint matrix, Qmat, and cobj, lb, ub from vars and constraints  lambda = 10
    Amat = sparseMatrix(constraints$cid, constraints$id, x = constraints$value, dims = c(length(ucid), nrow(vars)))
    vars[is.na(weight), weight := 0]

  if (verbose) {

    message("bvec length: ", length(bvec))
    message("Amat nrow: ", nrow(Amat))

  }
  if (any(ix <- is.infinite(vars$weight)))
  {
    warning('nodes with infinite weight, setting to 0, please check inputs')
    vars[ix, weight := 0]
  }
  Qmat = vars[, weight * (type %in% c('nresidual', 'eresidual', 'mresidual'))] %>% as.numeric %>% Diagonal(x = .) %>% as('CsparseMatrix')

  ## set lambda to 0 at terminal or other non NA nodes
  vars[is.na(lambda), lambda := 0]

 
  ## set cvec by multiplying global lambda by local lambda for non-terminal loose end
  ## vars (or their indicators if L0 is TRUE)
  if (L0)
    {
      if (loose.collapse)
      {
          cvec = lambda*(vars[, lambda*(type %in% c('loose.in.indicator.sum.indicator', 'loose.out.indicator.sum.indicator') & !terminal)] %>% as.numeric)
          ## cvec = lambda*(vars[, lambda*(type %in% c('loose.in.indicator.sum.indicator', 'loose.out.indicator.sum.indicator', 'loose.in.indicator', 'loose.out.indicator') & !terminal)] %>% as.numeric)
        }
      else
        {
          cvec = lambda*(vars[, lambda * (type %in% c('loose.in.indicator', 'loose.out.indicator') & !terminal)] %>% as.numeric)
        }
    } else {
    cvec = lambda*(vars[, lambda*(type %in% c('loose.in', 'loose.out') & !terminal)] %>% as.numeric)
  }

  ## implement reward if provided
  if (length(ix <- which(vars$reward!=0)))
  {
    if (verbose)
      message('Applying reward')
    cvec[ix] = -vars$reward[ix]

  }

    if (lp) {
        ## add weights of stuff
        indices = which(vars$type %in% c("mdelta.plus", "mdelta.minus",
                                         "ndelta.plus", "ndelta.minus",
                                         "edelta.plus", "edelta.minus"))
        wts = vars$weight[indices]
        cvec[indices] = wts
        Qmat = NULL ## no Q if solving LP
    }

  lb = vars$lb
  ub = vars$ub

  control = list(trace = ifelse(verbose>=2, 1, 0), tilim = tilim, epgap = epgap, round = 1)
  ## sol = Rcplex::Rcplex(cvec = cvec, Amat = Amat, bvec = bvec, Qmat = Qmat, lb = lb, ub = ub, sense = sense, vtype = vars$vtype, objsense = 'min', control = control)

    ## call our wrapper for CPLEX
    sol =  Rcplex2(cvec,
                   Amat,
                   bvec,
                   Qmat = Qmat,
                   lb = lb,
                   ub = ub,
                   sense = sense,
                   vtype = vars$vtype,
                   objsense = "min",
                   control = control,
                   tuning = FALSE)
    
  vars$cvec = cvec
  vars$x = sol$x

  ## for debugging
  ppc = function(x) (x %>% merge(vars, by = 'id') %>% merge(b, by = 'cid.char'))[, paste(paste(round(value.x, 1), '*', paste(type, gid, sep=  '_'), '(', signif(x, 2), ')', collapse = ' + '), ifelse(sense[1] == 'E', '=', ifelse(sense[1] == 'G', '>=', '<=')), round(value.y[1],2)), by = cid.char]
  
  ppv = function(x) {tmp = x %>% merge(constraints, by = 'id'); constraints[cid %in% tmp$cid, ] %>% ppc}

  .check = function(x) data.table(obs = sign(as.numeric(round(Amat %*% x - bvec))),
                                  sense)
  chk = .check(sol$x)

  if (any(is.na(sol$x)))
    stop('Rcplex did not converge or failed to find a solution, please run with verbose = 2 to get more detailed output')

  if (chk[sense == 'E', any(obs != 0, na.rm = TRUE)] |
      chk[sense == 'G', any(obs < 0, na.rm = TRUE)] |
      chk[sense == 'L', any(obs > 0, na.rm = TRUE)])
    stop('Constraint violation likely due to M parameter being too large for problem causing CPLEX numerical instability, consider lowering M parameter')

   ##.obj = function(x) 0.5 * rbind(x) %*% Qmat %*% cbind(x) + cvec %*% x

 
  ## update graph  
  nmark = vars[type == 'node', .(nid = abs(snode.id), cn = round(x))]
  emark = vars[type == 'edge', .(eid = abs(sedge.id), cn = round(x))]

  loosei = vars[type == 'loose.in' & snode.id>0, .(cn = round(x)), keyby = snode.id]
  looseo = vars[type == 'loose.out' & snode.id>0, .(cn = round(x)), keyby = snode.id]

  nodes = gg$nodes[loosei$snode.id] ## need to do this to use nodes active binding settings
  nodes$loose.left = loosei$cn>0

  nodes = gg$nodes[looseo$snode.id] ## need to do this to use nodes active binding settings
  nodes$loose.right = looseo$cn>0

  gg$nodes$mark(loose.cn.left = 0, loose.cn.right = 0)
  gg$nodes[loosei$snode.id]$mark(loose.cn.left = loosei$cn)
  gg$nodes[looseo$snode.id]$mark(loose.cn.right = looseo$cn)

  ## cache old cn values
  gg$nodes$mark(cn.old = gg$nodes$dt$cn)
  gg$edges$mark(cn.old = gg$edges$dt$cn)
  gg$nodes$mark(cn = NULL) ## reset to avoid weird type casting issue
  gg$edges$mark(cn = NULL) ## reset to avoid weird type casting issue
  gg$nodes[nmark$nid]$mark(cn = nmark$cn)
  gg$edges[emark$eid]$mark(cn = emark$cn)
  gg$set(y.field = 'cn')

  gg$set(obj = sol$obj)

##  fix loose ends
  nodes = gg$nodes 
  nodes$loose.left = nodes$dt$loose.cn.left>0
  nodes$loose.right = nodes$dt$loose.cn.right>0

  ## if phased, mark edges with different colors to make it easier to visualize
  if (phased) {
    if (verbose) {
      message("formatting phased graph...")
    }
    ## edge formatting
    ref.edge.col = alpha("blue", 0.5)
    alt.edge.col = alpha("red", 0.5)
    ref.edge.lwd = 1.0
    alt.edge.lwd = 1.0
    edge.col = ifelse(gg$edges$dt$type == "REF", ref.edge.col, alt.edge.col)
    edge.lwd = ifelse(gg$edges$dt$type == "REF", ref.edge.lwd, alt.edge.lwd)
    gg$edges$mark(col = edge.col, lwd = edge.lwd)

    ## mark zero cn edges
    zero.cn.col = alpha("gray", 0.1)
    zero.cn.lwd = 0.5
    zero.cn.edges = which(gg$edges$dt$cn == 0)
    gg$edges[zero.cn.edges]$mark(col = zero.cn.col, lwd = zero.cn.lwd)
  }
    if (debug) {
        return(list(gg = gg, sol = sol))
    }    
    return(gg)
}

#' @name jbaLP
#' @title jbaLP
#' @description jbaLP
#'
#' Simple (probably temporary) wrapper around balance for JaBbA LP
#' Takes karyograph as input and balances it.
#'
#' @param kag.file (character)
#' @param kag (karyograph object)
#' @param cn.field (character) column in karyograph with CN guess, default cnmle
#' @param var.field (character) column in karyograph with node variance estimate, default loess.var
#' @param bins.field (character) column in karyograph containing number of bins
#' @param min.var (numeric) min allowable variance default 1e-3
#' @param min.bins (numeric) min allowable bins default 5
#' @param lambda (numeric) slack penalty, default 10
#' @param L0 (logical) default TRUE
#' @param loose.collapse (logical) default FALSE
#' @param M (numeric) max CN
#' @param verbose (numeric) 0 (nothing) 1 (everything but MIP) 2 (print the MIP), default 1
#' @param tilim (numeric) default 1e3
#' @param ism (logical
#' @param epgap (numeric) default 1e-3
#'
#' @return
#' karyograph with modified segstats/adj. Adds fields epgap, cl, ecn.in, ecn.out, eslack.in, eslack.out to $segstats and edge CNs to $adj
#' 
#' @author Marcin Imielinski, Zi-Ning Choo
#' @export
jbaLP = function(kag.file = NULL,
                 kag = NULL,
                 cn.field = "cnmle",
                 var.field = "loess.var",
                 bins.field = "nbins",
                 min.var = 1e-3,
                 min.bins = 3,
                 lambda = 10,
                 L0 = TRUE,
                 loose.collapse = FALSE,
                 M = 1e3,
                 verbose = 1,
                 tilim = 1e3,
                 ism = FALSE,
                 epgap = 1e-3)
{
    if (is.null(kag.file) & is.null(kag)) {
        stop("one of kag or kag.file must be supplied")
    }
    if (!is.null(kag.file) & !is.null(kag)) {
        warning("both kag.file and kag supplied. using kag.")
    }
    if (!is.null(kag)) {
        if (verbose) {
            message("using supplied karyograph")
        }
    } else {
        if (file.exists(kag.file)) {
            if (verbose) {
                message("reading karyograph from file")
            }
            kag = readRDS(kag.file)
        } else {
            stop("kag.file does not exist and kag not supplied")
        }
    }
    kag.gg = gG(jabba = kag)

    if (verbose) {
        message("Marking nodes with cn contained in column: ", cn.field)
    }
    
    if (is.null(values(kag.gg$nodes$gr)[[cn.field]])) {
        stop("karyograph must have field specified in cn.field")
    }
    kag.gg$nodes$mark(cn  = values(kag.gg$nodes$gr)[[cn.field]])

    if (verbose) {
        message("Computing node weights using variance contained in column: ", var.field)
    }
    
    if (is.null(values(kag.gg$nodes$gr)[[var.field]]) | is.null(values(kag.gg$nodes$gr)[[bins.field]])) {
        warning("karyograph missing var.field. setting weights to node widths")
        wts = width(kag.gg$nodes$gr)
    } else {
        
        ## process variances
        vars = values(kag.gg$nodes$gr)[[var.field]]
        vars = ifelse(vars < min.var, NA, vars) ## filter negative variances
        sd = sqrt(vars) * kag$beta ## rel2abs the standard deviation

        ## process bins
        bins = values(kag.gg$nodes$gr)[[bins.field]]
        bins = ifelse(bins < min.bins, NA, bins)

        ## compute node weights
        wts = bins / (sd / sqrt(2)) ## for consistency with Laplace distribution
        wts = ifelse(is.infinite(wts) | is.na(wts) | wts < 0, NA, wts)
    }
    kag.gg$nodes$mark(weight = wts)
    
    ## no edge CNs
    kag.gg$edges$mark(cn = NULL)
    kag.gg$nodes[cn > M]$mark(cn = NA, weight = NA)

    if (verbose) {
        message("Starting LP balance on gGraph with...")
        message("Number of nodes: ", length(kag.gg$nodes))
        message("Number of edges: ", length(kag.gg$edges))
    }

    res = balance(kag.gg,
                  debug = TRUE,
                  lambda = lambda,
                  L0 = TRUE,
                  verbose = verbose,
                  tilim = tilim,
                  epgap = epgap,
                  lp = TRUE,
                  ism = ism)
    
    bal.gg = res$gg
    sol = res$sol

    ## just replace things in the outputs
    out = copy(kag)
    new.segstats = bal.gg$gr
    nnodes = length(new.segstats)

    new.segstats$cl = 1 ## everything same cluster
    new.segstats$epgap = sol$epgap ## add epgap from genome-side opt
    
    ## weighted adjacency
    adj = sparseMatrix(i = bal.gg$sedgesdt$from, j = bal.gg$sedgesdt$to,
                       x = bal.gg$sedgesdt$cn, dims = c(nnodes, nnodes))
    ## add the necessary columns
    new.segstats$ecn.in = Matrix::colSums(adj)
    new.segstats$ecn.out = Matrix::rowSums(adj)
    target.less = (Matrix::rowSums(adj, na.rm = T) == 0)
    source.less = (Matrix::colSums(adj, na.rm = T) == 0)
    new.segstats$eslack.out[!target.less] = new.segstats$cn[!target.less] - Matrix::rowSums(adj)[!target.less]
    new.segstats$eslack.in[!source.less] =  new.segstats$cn[!source.less] - Matrix::colSums(adj)[!source.less]
    out$adj = adj

    ## add metadata
    out$segstats = new.segstats
    out$status = sol$status
    out$epgap = sol$epgap
    return(out)
}

#' @name balance.alleles
#' @description balance.alleles
#'
#' Infers parental haplotype graph given simplified unphased junction-balanced graph and het SNP counts
#'
#' @param jab JaBbA object with fields $agtrack, $asegstats, $aadj
#' @param major.count.field character specifying alt.count meta data field in input het.sites (default $alt)
#' @param minor.count.field character specifying ref.count meta data field in input het.sites (default $ref)
#' @param fix.marginals (logical) whether to fix marginals (default TRUE)
#' @param fix.width.thres (logical) fix marginals with width above this threshold (default 1e6)
#' @param ism (logical) ism contraints (default TRUE)
#' @param lambda (numeric) slack penalty (default NULL, set to median node weight)
#' @param postprocess logical, collapse unphased segments (default TRUE)
#' @param epgap (numeric) default 1e-4
#' @param tilim (numeric) default 1e3
#' @param verbose logical (default 0)
#' 
#' @return
#' balanced parental allelic gGraph
#'
#' @export
balance.alleles = function(jab,
                           major.count.field = "high",
                           minor.count.field = "low",
                           fix.marginals = TRUE,
                           fix.width.thres = 1e6,
                           ism = TRUE,
                           lambda = NULL,
                           postprocess = TRUE,
                           epgap = 1e-3,
                           tilim = 1000,
                           verbose = 0) {
    if (!(all(c("agtrack", "purity", "ploidy") %in% names(jab)))) {
        stop("jab must contain agtrack")
    }
    ## extract het counts from agtrack
    hets.gr = jab$agtrack@data[[1]][, c("count", "type")]
    hets.gr$allele = ifelse(hets.gr$type == "high", "major", "minor")

    ## get unphased gGraph
    unphased.gg = gG(jabba = jab)

    marginals.gr = unphased.gg$nodes$gr[, "cn"] %>% gr.stripstrand %Q% (!is.na(cn) & cn > 0)

    if (fix.marginals) {
        marginals.gr$fix = ifelse(width(marginals.gr) > fix.width.thres, 1, 0)
    } else {
        marginals.gr$fix = 0
    }

    ## make unbalanced phased graph
    phased.gg = phased.binstats(unphased.gg, hets.gr, purity = jab$purity, ploidy = jab$ploidy)

    ## infer lambda and M
    if (is.null(lambda)) {
        lambda = median(phased.gg$nodes$dt$weight, na.rm = TRUE)
    }

    ## balance (gives back solution status, etc)
    res = balance(phased.gg,
                  phased = TRUE,
                  marginal = marginals.gr,
                  lambda = lambda,
                  verbose = verbose,
                  M = 1000,
                  ism = ism,
                  epgap = epgap,
                  tilim = tilim,
                  debug = TRUE)
    
    ## unpack
    balanced.gg = res$gg
    sol = res$sol
    
    ## postprocess
    postprocessed.gg = phased.postprocess(balanced.gg)

    ## store epgap and solution status as node metadata i guess
    postprocessed.gg$nodes$mark(epgap = sol$epgap, cl = 1)

    return(postprocessed.gg)
}


#' @name nodestats
#' @description nodestats
#'
#' Computes copy number (CN) "stats" using on graph gg using (binned) copy number data in GRanges data to prep a
#' graph for balance function (see above).
#'
#' The function outputs a graph whose nodes are populated with $cn and $weight fields.  The cn field is computed as
#' the centroid (measured by function FUN) and the weight is computed as number of bins / 2 * (sample variance) across
#' the bins overlapping that node.
#'
#' If loess = TRUE is specified, then a mean vs variance curve is computed across all nodes and the sample variance in each node
#' is replaced with the mapping assigned to the centroid
#' 
#' The node to bin mapping can be additionally restricted by the columns in "by", which must exist in both the node metadata
#' of gg and data
#' 
#' @param gg gGraph to compute node stats across
#' @param data GRanges of cn data to compute node stats against
#' @param cn.field numeric field corresponding to cn data which is used to aggregate (default is first column in the data)
#' @param loess flag whether to use LOESS to compute variance as mean to variance 
#' @param FUN function to compute centroids with (mean)
#' @param by by field to match data against node data, the column specified in this field must be present in both the gGraph (NULL)
#' @return gGraph with fields $cn
#' @export
#' @author Marcin Imielinski
nodestats = function(gg,
                   data,
                   cn.field = names(values(data))[[1]],
                   FUN = mean,
                   loess = FALSE, 
                   by = NULL)
{
  gg = gg$copy
  data$cn = values(data)[cn.field]
  ov = gr.findoverlaps(gg$nodes$gr, data[, "cn"], scol = 'cn', by = by)  %>% gr2dt
  dt = ov[, .(mean = FUN(cn, na.rm = TRUE), var = var(cn, na.rm = TRUE), nbins = .N), keyby = query.id][.(1:length(gg$nodes)), ]
  dt$weight = dt$nbins/(2*dt$var)
  dt[is.infinite(weight), weight := NA]
  gg$nodes$mark(cn = dt$mean, weight = dt$weight)
  return(gg)
}



#' @name transplant
#' @title transplant donor subggraph into recipient
#' @description
#'
#' Here we transplant a subgraph (eg representing an event)
#' into a recipient graph.  The transplant operation will
#' (1) fix the donor ALT junction copy numbers
#' (2) leave the subgraph node copy numbers unconstrained
#' (3) minimize the loose ends of the resulting balanced graph
#'
#' Donor and recipient graph should both have $cn field defined on
#' nodes and edges. Donor can also just be junctions, in which case
#' the junction shadow will be taken as the footprint, or a
#' footprint can be manually provided.  The footprint
#' will then provide a region where the copy number constraints
#' are relaxed while loose ends strictly enforced. 
#' 
#' Currently transplant is only well defined for haploid / unphased graphs ie those
#' where there is at most a single interval representing a given
#' reference genomic regions. 
#' 
#' @param gg "haploid" gGraph with field $cn
#' @param donor donor (haploid) gGraph or Junction object with field $cn 
#' @param footprint if gGraph not provided then footprint can be directly provided 
#' @param lambda loose end penalty to apply to balance operation (1000) to nodes within the donor region
#' @param L0 flag whether to use L0 penalty (FALSE)
#' @return balanced gGraph with donor junctions incorporated at their donor cn
#' @author Marcin Imielinski 
#' @export
transplant = function(gg, donor, footprint = NULL, lambda = 1000, L0 = FALSE)
{
  if (!('cn' %in% names(gg$nodes$dt) & 'cn' %in% names(gg$edges$dt)))
    stop('cn field must be defined on the nodes and edges recipient graph')

  tmp = (gg$nodes$gr %>% gr.sum)
  is.haploid = sum(tmp$score>1)==0

  if (!is.haploid)
    stop('transplant is only defined on haploid graphs, please check your graph for overlapping nodes')

  gg = gg$copy

  if (inherits(donor, 'gGraph'))
  {
    footprint = donor$footprint
    junctions = donor$junctions[type == 'ALT']
    gg$disjoin(donor$gr)
  } else
  {
    junctions = donor
    if (is.null(footprint)) ## in this case gg is a junction
    {
      footprint = junctions$shadow
    }
    gg$disjoin(grbind(junctions$grl %>% unlist, footprint))
  }
  
  if (length(junctions)>0 && !('cn' %in% names(junctions$dt)))
    stop('cn field must be defined on the edges of the donor graph')

  combined.junctions = gg$junctions[type == 'ALT']
  if (length(junctions)) ## mark as donor and combine with gg $junctions
    {
      junctions$set(donor = TRUE)
      combined.junctions = c(combined.junctions, junctions[, c('cn', 'donor')])
    }
  
  ## mark edges of the donor subgraph as donor = TRUE
  gg$edges$mark(donor = FALSE)

  ## make new gGraph combining junctions from recipient and donor
  ggn = suppressWarnings(gG(breaks = gg$nodes$gr,
                            junctions = combined.junctions))


  ## fix cn on the edges of the donor subgraph
  efix = ggn$edges[type == 'ALT']$dt$edge.id

  fp.nodes = ggn$nodes %&% footprint

  ## relax loose ends everywhere but in the subgraph
  this.lambda = lambda
  ggn$nodes$mark(lambda = 1)
  fp.nodes$mark(lambda = this.lambda)

  ## balance combined gGraph
  ggn = balance(ggn, efix = efix, nrelax = fp.nodes$dt$node.id, L0 = L0)
  
  return(ggn)
}

#' @name bcheck
#' @title bcheck
#' @description
#'
#' Checks if genome graph with node and edge cn's annotated and loose cn's is balanced.
#' 
#' @param gg gGraph with node and edge fields 'cn', node field loose.cn.left and loose.cn.right
bcheck = function(gg)
{
  lmerge = merge(gg$nodes$dt, gg$nodes$eleft$dt, by.x = 'node.id', by.y = 'n2', all.x = TRUE)
  rmerge = merge(gg$nodes$dt, gg$nodes$eright$dt, by.x = 'node.id', by.y = 'n1', all.x = TRUE)

  if (!('loose.cn.left' %in% names(lmerge)))
    {
      lmerge[, loose.cn.left := 0]
    }

  if (!('loose.cn.right' %in% names(rmerge)))
    {
      rmerge[, loose.cn.right := 0]
    }

  out = merge(lmerge[, .(left = cn.x[1] - sum(cn.y, na.rm = TRUE) - loose.cn.left[1]), keyby = node.id],
              rmerge[, .(right = cn.x[1] - sum(cn.y, na.rm = TRUE) - loose.cn.right[1]), keyby = node.id])[.(gg$nodes$dt$node.id), ]
  out[, balanced := left == 0 & right == 0]

  return(out)
}

#' @name loosefix
#' @title loosefix
#' @description
#'
#' Updates loose.left, loose.right, loose.cn.left, loose.cn.right based on value of cn field
#' 
#' @param gg gGraph with node and edge fields 'cn', node field loose.cn.left and loose.cn.right
#' @author Marcin Imielinski
loosefix = function(gg)
{
  lmerge = merge(gg$nodes$dt, gg$nodes$eleft$dt, by.x = 'node.id', by.y = 'n2', all.x = TRUE)
  rmerge = merge(gg$nodes$dt, gg$nodes$eright$dt, by.x = 'node.id', by.y = 'n1', all.x = TRUE)

  ## mark any cn NA as 0
  gg$nodes[is.na(cn)]$mark(cn = 0)
  gg$edges[is.na(cn)]$mark(cn = 0)

  out = merge(lmerge[, .(loose.cn.left = cn.x[1] - sum(cn.y, na.rm = TRUE)), keyby = node.id],
              rmerge[, .(loose.cn.right = cn.x[1] - sum(cn.y, na.rm = TRUE)), keyby = node.id])[.(gg$nodes$dt$node.id), ]


  out[is.na(loose.cn.left), loose.cn.left := 0]
  out[is.na(loose.cn.right), loose.cn.right := 0]

  if (any(ix <- (out$loose.cn.left<0 | out$loose.cn.right<0)))
    stop(sprintf('some nodes (%s) have a higher incoming or outgoing edge copy number than the number of node copies', paste(ix %>% which, collapse = ',')))

  gg$nodes$mark(loose.cn.left = out$loose.cn.left,
                loose.cn.right = out$loose.cn.right)

  nodes = gg$nodes;
  nodes$loose.left = out$loose.cn.left > 0
  nodes$loose.right = out$loose.cn.right > 0

  return(gg)
}


#' @name peel
#' @title peel
#' @description
#'
#' Greedily "peels" walks off a genome graph with node and edge field "cn"
#' by finding successive flows that maximizes some edge metadata field (specified as field)
#' which if not specified will be a logical field type == 'ALT' specifying whether that
#' junction is an ALT edge
#'
#' @param gg gGraph with field cn
#' @param field edge metadata field to use to rank walk solutions (default edge field type == 'ALT')s
#' @param verbose flag = 1 regular verbosity, 2 = dump out Rcplex traces
#' @param embed.loops logical flag (FALSE) if TRUE will embed all the loops in the output into an arbitrary linear (path) walk
#' @author Marcin Imielinski
#' @return collection of gWalks annotated with cn on the original that when added will give the marginal copy profile on inputted nodes and edges
#' @export
peel = function(gg, field = NULL, embed.loops = FALSE, verbose = FALSE)
{
  gg = refresh(loosefix(gg))
  gg.og = refresh(gg)

  ## mini function to compute the cn of a walk in a graph by finding the min cn across its associated
  ## nodes, edges, and loose ends
  .mincn = function(walks)
  {
    source = walks$eval(node = snode.id[1])
    sink = walks$eval(node = snode.id[.N])

    source.loose = walks$eval(node = ifelse(strand[1] == '+', loose.cn.left[1], loose.cn.right[1]))
    sink.loose =  walks$eval(node = ifelse(strand[.N] == '+', loose.cn.right[.N], loose.cn.left[.N]))

    edge.min = Inf
    if (walks$edges %>% length)
      edge.min = walks$eval(edge = data.table(cn, id = abs(sedge.id))[, .(CN = cn[1]/.N), by = id][, min(floor(CN),  na.rm = TRUE)])

    node.min = walks$eval(node = data.table(cn, id = abs(snode.id))[, .(CN = cn[1]/.N), by = id][, min(floor(CN),  na.rm = TRUE)])

    pmin(
      ifelse(walks$circular, Inf, ## if circular no loose end capacity constraints
             ## fold back edge case where we mihgt overestimate our loose end capacity
      ifelse(source == -sink, floor(source.loose/2), 
             pmin(source.loose, sink.loose))), # otherwise we are limited by the loose end cn on each walk side
      edge.min, ## and the internal edge cn
      node.min ## and the internal node cn
    )
  }
   

  if (is.null(field))
  {
    gg$edges$mark(alt = gg$edges$dt$type == 'ALT')
    field = 'alt'
  }

  ## .check = function(gw, gg)
  ## {
  ##   gg2 = gW(grl = rep(gw$grl, gw$dt$cn), circular = rep(gw$circular, gw$dt$cn))$graph;
  ##   gg2$nodes$mark(cn = 1);
  ##   gg2$edges$mark(cn = 1);
  ##   gg2 = c(gg$copy, gg2)$disjoin()
  ##   gg2$nodes$mark(CN = gg2$nodes$dt$cn)
  ##   return(gg2)
  ## }

  out = gW(graph = gg)    

  ## loop next batch of max flow walks until no more walks
  while (length(walks <- gg$maxflow(walk = TRUE, path.only = FALSE, multi = TRUE, efield = field, cfield = 'cn', verbose = verbose)))
  {
    ## gg2 = .check(out, gg)
    ## cc =  gr2dt(gg.og$nodes$gr %*% gg2$nodes$gr)[CN != cn, .(node.id, cn, CN)]
    ## ee = merge(gg.og$junctions, gg2$junctions)$dt[cn != cn.1, .(edge.id, edge.id, cn, cn.1)]
    ## if (nrow(cc))
    ##   browser()

    ## first check for any walks that contain
    ## a node from another walk (edge cases involving reverse complements)
    walks = walks[rev(order(walks$dt$wid))]
    dups = dunlist(walks$snode.id)[, .(count = length(unique(listid)), listid = unique(listid)),
                              by = .(id = abs(V1))][count>1, ]

    if (walks$edges %>% length)
    {
      dups = rbind(dups, dunlist(walks$sedge.id)[, .(count = length(unique(listid)), listid = unique(listid)),
                              by = .(id = paste0('e', abs(V1)))][count>1, ])
    }
    
    if (nrow(dups))
    {
      keep = !(1:length(walks) %in% dups[duplicated(id), listid %>% unique])
      walks = walks[keep]
    }


    ## cn of walk is the min cn of edges, correcting for edges that get hit multiple times
    walks$set(cn = .mincn(walks))

    if (verbose)
      message('Peeling off ', length(walks), ' walks with max width ',
              max(walks$dt$wid), ' and max copy ', max(walks$dt$cn, na.rm = TRUE))

    ## add batch to output
    out = c(walks, out)
    
    ## peel off graph
    ## annotate walk with min cn
    while (length(walks <- walks[cn>0]))
    {
      ## figure out how much to decrement
      ## aggregate in case both strands of node is in a walk
#      gg.cache = gGnome::refresh(gg)

      ndec = dunlist(walks$snode.id)[, wcn := walks$dt$cn[listid]][, .(dec = sum(wcn)), keyby = .(nid = abs(V1))]
      ndec$cn = gg$nodes[ndec$nid]$dt$cn - ndec$dec
      gg$nodes[ndec$nid]$mark(cn = ndec$cn)

      if (walks$edges %>% length)
      {
        edec = dunlist(walks$sedge.id)[, wcn := walks$dt$cn[listid]][!is.na(V1), .(dec = sum(wcn)), keyby = .(eid = abs(V1))]
        edec$cn = gg$edges[edec$eid]$dt$cn - edec$dec
        gg$edges[edec$eid]$mark(cn = edec$cn)
      }
            
      ## to recompute walk cn have to take into account edges that get hit more than once
      gg = loosefix(gg)
      walks$set(cn = .mincn(walks))
#      bc = bcheck(gg)[balanced==FALSE, ]      

      ## if (gg$nodes$dt[, any(loose.cn.left<0 | loose.cn.right<0)])
      ##   browser()

      ## if (nrow(bc))
      ##   browser()
    }

    if (verbose)
    {
      ploidy = gg$nodes$dt[, sum(cn*width, na.rm = TRUE)/sum((1+0*cn)*width, na.rm = TRUE)]
        message('... remaining ploidy ', round(ploidy,2), ' across ', sum(gg$nodes$dt$cn>0, na.rm = TRUE), ' nodes and ', sum(gg$edges$dt$cn>0), ' edges with nonzero CN remaining in graph with total ', length(out), ' walks in collection')
      }
  }

  ## recast walks around original (input) graph
  walks = gW(snode.id = out$snode.id, circular = out$circular, meta = out$dt, graph = gg.og)
  walks = walks[rev(order(wid))]

  if (embed.loops)
    walks = embedloops(walks[circular == TRUE], walks[circular == FALSE], verbose = verbose>1)

  return(walks)
}


#' @name embedloops
#' @description embedloops
#'
#' Attempts to embed / tarnsplant a set of circular walks (loops)
#' into a set of recipients, both defined on the same gGraph, and 
#' both optionally having the metadata $cn.  Returns a gWalk with as
#' many of the loops embedded into the recipients as possible.
#'
#' A loop needs to intersect some node in order to be successfully embedded.
#' (note: that node may be in another loop, if that loop gets embedded)
#' (As a result though, we don't guarantee that every loop will be enbedded)
#'
#' Noter: loops with cn>1 will be embedded in tandem.  
#'
#' @param loops circular gWalk object
#' @param recipient set of gWalks defined on the same object
#' @return gWalk with as many loops embedded into recipients as possible
#' @export
#' @author Marcin Imielinski
embedloops = function(loops, recipients, verbose = FALSE)
{
  if (length(loops)==0)
    return(recipients)

  if (length(recipients)==0)
    return(loops)

  if (!all(loops$circular))
    stop('Loops must be circular walks')

  ## if cn not set then set to 1
  if (is.null(loops$dt$cn))
    loops$set(cn = 1)

  ## if cn not set then set to 1
  if (is.null(recipients$dt$cn))
    recipients$set(cn = 1)

  ## subroutine to embed
  .embed = function(donor, recipients)
  {
    entry = dunlist(recipients$snode.id)[abs(V1) %in% abs(donor$snode.id[[1]]), ][1,]
    
    if (is.na(entry$V1))
      return(list(embedded = FALSE, recipients = recipients))
    
    if (entry$V1 %in% -donor$snode.id[[1]])
      donor = donor[-1] ## reverse complement
    
    ## pivot donor around candidate node
    snid = donor$snode.id[[1]]
    ix = match(entry$V1, snid)
    snidp = snid[c(ix:length(snid))]
    snidp = c(snidp, setdiff(snid, snidp))
    
    this.recipient = recipients[entry$listid %>% as.integer]
    snidr = this.recipient$snode.id[[1]]

    ix = match(entry$V1, snidr)

    ## insert (several copies) of pivoted donor cycle just upstream of the entry point
    snid.new = c(snidr[seq_len(ix-1)], rep(snidp, donor$dt$cn), snidr[ix:length(snidr)])

    ## make new walk
    gw.new = gW(snode.id = list(snid.new), graph = this.recipient$graph, circular = this.recipient$circular, meta = this.recipient$meta)

    ## if recipient walk has more than one copy then we will need to preserve additional
    ## "original" copies of this.recipient walk
    if (this.recipient$dt$cn>1)
    {
      this.recipient$set(cn = this.recipient$dt$cn - 1)
      gw.new$set(cn = 1)
      gw.new = c(gw.new, this.recipient)
    }
   
    recipients = c(recipients[setdiff(1:length(recipients), entry$listid)], gw.new)

    recipients$set(numalt = recipients$eval(edge = sum(type == 'ALT')))
    recipients = recipients[rev(order(numalt))]
                   
    if (verbose)
    {
      message('Embedded loop ', paste(snidp, collapse = '->'), ' into recipient ', paste(snidr, collapse = '->'), ' giving ', paste(snid.new, collapse = '->'))
    }

    return(list(embedded = TRUE, recipients = recipients))
  }
  
  embedded = done = rep(FALSE, length(loops))

  while (!all(done))
  {
    old.recipients = copy(recipients)
    for (i in which(!done))
    {
      res = .embed(loops[i], recipients)
      embedded[i] = done[i] = res$embedded
      recipients = res$recipients
    }

    ## we give up if a round of embedding does not change recipients
    if (identical(old.recipients$snode.id, recipients$snode.id))
      done = rep(TRUE, length(loops))
  }

  ## return (modified) recipients and any leftover loops
  return(c(recipients, loops[which(!embedded)]))
}

#' @name binstats
#' @title binstats
#' @description
#'
#' Given GRanges of binned eg read depth data with field $cn, crosses 
#' nodes in graph and annotates graph nodes with
#' 'cn' and 'weight'.  Done upstream of balance. 
#'
#' If "by" field(s) specified and these fields exist in both
#' bins and graph nodes, then will use these in the overlaps query
#'
#' If field, purity, and ploidy provided then will
#' also transform read depth data in bin column "field"
#' using purity and ploidy to generate
#' @param gg gGraph
#' @param bins GRanges with field $cn or field field
#' @param by optional character vector specifying metadata field(s) shared by gg$nodes and bins to which bin / node overlaps will be limited to
#' @param field character vector field of bins to convert to (NULL)
#' @param purity purity parameter either specified together with field or embedded in gg$meta, must be specified if field is not NULL
#' @param ploidy ploidy parameter either specified together with field or embedded in gg$meta, must be specified if field is not NULL
#' @param min.bins minimum number of bins to use for intra segment variance computation (3)
#' @param loess logical flag whether to smooth / fit variance using loess (FALSE)
#' @param min.var minimal allowable per segment bin variance, which will ignore segments with very low variance due to all 0 or other reasons (0.1)
#' @param lp (logical) return weights consistent with LP optimization
#' 
#' @return gGraph whose nodes are annotated with $cn and $weight field
#' @export
#' @author Marcin Imielinski
binstats = function(gg, bins, by = NULL, field = NULL, purity = gg$meta$purity, ploidy = gg$meta$ploidy, loess = TRUE, min.bins = 3, verbose = TRUE, min.var = 0.1, lp = FALSE)
{
  gg = gg$copy

  if (!is.null(field) & !is.null(purity) & !is.null(ploidy) && is.numeric(purity) && is.numeric(ploidy))
  {
    if (verbose)
      message('Converting ', field, ' to cn using purity ', purity, ' and ploidy ', ploidy)

    bins$cn = rel2abs(bins, field = field, purity = purity, ploidy = ploidy)
  }

  if (is.null(bins$cn))
    stop('bins must have field cn or a field, purity, and ploidy must be specified where field is a column in bins')

  if (verbose)
    message('crossing nodes and bins via gr.findoverlaps')
  ov = gr.findoverlaps(gg$nodes$gr, bins, by = by, scol = names(values(bins)), return.type = 'data.table')
  if (verbose)
    message('aggregating bin stats per node')
  dt = ov[!is.na(cn), .(mean = mean(cn, na.rm = TRUE), var = var(cn, na.rm = TRUE), nbins = .N), keyby = query.id][.(1:length(gg$nodes)), ]
  dt[nbins<min.bins, var := NA]

  if (loess)
  {
    min.var = pmax(min.var, min(dt$var, na.rm = TRUE)) ## min allowable var
    loe = dt[!is.na(var) & !is.na(mean), loess(var ~ mean, weights = nbins)] ## loe = tmp[, loess(var ~ mean)] ##                  
    dt$var.obs = dt$var
    dt$var = pmax(min.var, predict(loe, dt$mean))
  }


  if (verbose)
    message('computing weights and returning')
  if (lp) {
      dt$weight = dt$nbins/(sqrt(dt$var) / sqrt(2))
  } else {
      dt$weight = dt$nbins/(2*dt$var)
  }


  if (any(is.infinite(dt$weight), na.rm = TRUE))
    warning('variance computation yielded infinite weight, consider setting nbins higher or using loess fit')

  gg$nodes$mark(cn = dt$mean, weight = dt$weight)
  return(gg)
}

#' @name phased.postprocess
#' @title phased.postprocess
#' @description
#'
#' Postprocess junction-balanced phased graph and creates unphased regions
#' This identifies regions without allelic CN imbalance
#'
#' @param gg junction-balanced phased gGraph. each node must have associated og node.id
#' @param phase.blocks (GRanges) granges of phase blocks from linked reads. default = NULL
#' @param mc.cores (int) number of cores. default = 8.
#' @param verbose (bool) verbose > 0 prints stuff. default 1.
#'
#' @export
phased.postprocess = function(gg, phase.blocks = NULL, mc.cores = 8, verbose = 1)
{
    ## check that gg nodes and edges have og node
    if (!("og.node.id" %in% colnames(gg$nodes$dt)) | !("og.edge.id" %in% colnames(gg$edges$dt))) {
        stop("gGraph must have og.node.id and og.edge.id node/edge metadata columns")
    }
    ## check that graph has been balanced (need cn.old and cn)
    if (!("cn" %in% colnames(gg$nodes$dt)) | !("cn.old" %in% colnames(gg$nodes$dt))) {
        stop("run balance to populate nodes with cn and cn.old")
    }

    ## make a copy of balanced graph to prevent mutation
    if (verbose) {
        message("Making a copy of input gGraph")
    }
    gg = gg$copy

    ## identify nodes without CN imbalance
    if (verbose) {
        message("Identifying nodes without CN imbalance")
    }
    og.node.balance = gg$nodes$dt[, .(og.node.id, allele, cn)] %>%
        dcast.data.table(og.node.id ~ allele, value.var = "cn")

    og.node.balance[, cn.imbalance := (major != minor)]
    og.node.balance[, cn.total := (major + minor)]

    og.node.balance[, phased := ifelse(cn.imbalance == TRUE | cn.total == 0, TRUE, FALSE)]

    if (verbose) {
        message("Number of potentially unphased nodes: ", sum(og.node.balance$phased == FALSE))
    }

    if (!is.null(phase.blocks)) {
        ## need to edit this to distinguish between phase blocks and CBS blocks
        
        if (verbose) {
            message("Identifying nodes lying within one phase block")
        }
        ## identify nodes lying entirely within one phase block
        n.pblocks = gg$nodes$gr %N% phase.blocks

        ## add number of phase blocks
        og.node.balance[, nblocks := n.pblocks[match(og.node.id, gg$nodes$dt$og.node.id)]]

        ## if within a single block then phased
        og.node.balance[nblocks == 1, phased := TRUE]

        if (verbose) {
            message("Number of potentially unphased nodes after considering phase blocks: ",
                    sum(og.node.balance$phased == FALSE))
        }
    }

    ## identify unphased og nodes
    unphased.og.nodes = og.node.balance[phased == FALSE, og.node.id]

    ## identify corresponding minor/major nodes
    unphased.minor.nodes = gg$nodes$dt[og.node.id %in% unphased.og.nodes & allele == "minor", node.id]
    unphased.major.nodes = gg$nodes$dt[og.node.id %in% unphased.og.nodes & allele == "major", node.id]

    ## create new data.table for nodes
    new.nodes.dt = gg$nodes$dt[!(node.id %in% unphased.minor.nodes),]

    ## mark major nodes as unphased
    new.nodes.dt[node.id %in% unphased.major.nodes, allele := "unphased"]

    ## reset CN to total CN
    new.nodes.dt[node.id %in% unphased.major.nodes,
                 cn := og.node.balance$cn.total[match(og.node.id, og.node.balance$og.node.id)]]

    ## fix the CN of all of these nodes
    new.nodes.dt[, fix := 1]

    ## create new data.table for edges
    new.edges.dt = gg$edges$dt[!(n1 %in% unphased.minor.nodes) | !(n2 %in% unphased.minor.nodes),]

    ## reformat nodes
    new.nodes.dt[allele == "unphased", col := alpha("gray", 0.5)]

    ## get nodes as GRanges
    new.nodes.gr = dt2gr(new.nodes.dt[, .(seqnames, start, end,
                                          og.node.id, marginal.cn, allele,
                                          var, nbins, weight, index, col,
                                          cn.old, cn, fix, ywid,
                                          old.node.id = node.id)]) %>% gr.sort


    ## reset edge endpoints
    dt = gg$nodes$dt[og.node.id %in% unphased.og.nodes, .(og.node.id, allele, node.id)] %>%
        dcast.data.table(og.node.id ~ allele, value.var = "node.id")
    new.edges.dt[(n1 %in% unphased.minor.nodes), n1 := dt$major[match(n1, dt$minor)]]
    new.edges.dt[(n2 %in% unphased.minor.nodes), n2 := dt$major[match(n2, dt$minor)]]

    ## reset all edge endpoints to new node.ids
    new.edges.dt[, n1 := match(n1, new.nodes.gr$old.node.id)]
    new.edges.dt[, n2 := match(n2, new.nodes.gr$old.node.id)]

    new.edges.dt = new.edges.dt[cn > 0,]

    ## remove edge CN and fix
    if ("fix" %in% colnames(new.edges.dt)) {
        new.edges.dt$fix = NULL
    }

    ## if ("cn" %in% colnames(new.edges.dt)) {
    ##     new.edges.dt$cn = NULL
    ## }
    

    if (verbose) {
        message("Creating new gGraph")
    }
    ## postprocessed.gg = balance(gG(nodes = new.nodes.gr, edges = new.edges.dt),
    ##                            M = 1e3, ism = FALSE, verbose = verbose, epgap = 1e-4,
    ##                            marginal = NULL)

    new.nodes.gr = inferLoose(new.nodes.gr, new.edges.dt)

    postprocessed.gg = gG(nodes = new.nodes.gr, edges = new.edges.dt)
    postprocessed.gg$set(y.field = "cn")
    return(postprocessed.gg)
}


#' @name phased.binstats
#' @title phased.binstats
#' @description
#'
#' Given GRanges containing major/minor allele counts and a balanced but unphased gGraph,
#' prepares phased gGraph input to balance.
#' 
#' @param gg gGraph
#' @param bins GRanges with:
#' @param purity (numeric)
#' @param ploidy (numeric)
#' @param count.field (str) field containing allele read counts (default count)
#' @param allele.field (str) field for containing allele label for read counts (default allele)
#' @param phase.blocks (GRanges) GRanges containing phase blocks (e.g. from HAPCUT). default NULL.
#' @param edge.phase.dt (data.table) with columns n1.major, n2.major, n1.minor, n2.minor and edge.id providing major/minor allele counts
#' @param vbase.count.thres (int) number of variant base counts required to phase edges (default 5)
#' @param vbase.prop.thres (float) proportion of allele excess required to phase edges (default 0.9)
#' @param min.bins (numeric) minimum number of bins for intra segment variance (default 3)
#' @param min.var (numeric) min allowable variance (default 0.1)
#' @param verbose (bool) default TRUE for debugging
#' @param mc.cores (int) number of cores
#' @return gGraph whose nodes are annotated with $cn.major, $cn.minor, $haplotype, and $weight fields
#' @export
phased.binstats = function(gg, bins = NULL, purity = NULL, ploidy = NULL,
                           count.field = "count", allele.field = "allele",
                           phase.blocks = NULL,
                           edge.phase.dt = NULL,
                           vbase.count.thres = 5, vbase.prop.thres = 0.9,
                           min.bins = 3, min.var = 1e-3,
                           verbose = TRUE, mc.cores = 8)
{
    if (verbose) {
        message("Checking inputs")
    }
    if (is.null(purity)) {
        warning("Purity not provided, setting to 1.0")
        purity = 1
    }
    if (is.null(ploidy)) {
        warning("Ploidy not provided, setting to 2.0")
        ploidy = 2
    }
    if (!inherits(bins, 'GRanges')) {
        stop("bins must be GRanges")
    }
    if (!(count.field %in% names(values(bins)))) {
        stop("count.field not found in bins metadata")
    }
    if (!(allele.field %in% names(values(bins)))) {
        stop("allele.field not found in bins metadata")
    }
    allele.values = unique(values(bins)[[allele.field]])
    if (!("major" %in% allele.values) | !("minor" %in% allele.values)) {
        stop("allele.field must contain labels 'major' and 'minor'")
    }
    if (!is.null(phase.blocks)) {
        if (!inherits(phase.blocks, 'GRanges')) {
            stop("phase.blocks must be GRanges")
        }
    }
    if (!is.null(edge.phase.dt)) {
        if (!is.data.table(edge.phase.dt)) {
            warning("edge.phase.dt must be data.table. ignoring this input.")
            edge.phase.dt = NULL
        }
        if (!all(c("edge.id", "n1.major", "n2.major", "n1.minor", "n2.minor") %in% colnames(edge.phase.dt))) {
            warning("edge.phase.dt does not have the required columns")
            edge.phase.dt = NULL
        }
    }
    
    ## helper function for mean allele-specific read count
    reads.to.allele.cn = function(bins, count.field, purity, ploidy) {
        y = values(bins)[[count.field]]
        y.bar = 2 * mean(y, na.rm = TRUE)

        ## purity and ploidy
        alpha = purity
        tau = ploidy

        ## linear equation
        denom = alpha * tau + 2 * (1 - alpha)
        beta = (y.bar * alpha) / denom
        gamma =(y.bar * (1 - alpha)) / denom

        cn = (y - gamma) / beta
        return(cn)
    }

    if (verbose) {
        message("Preparing phased gGraph nodes")
    }
    ## get original nodes from unphased gGraph
    og.nodes = gg$nodes$gr[,c()]

    ## keep original node ID annotation
    og.nodes$og.node.id = gg$nodes$dt$node.id

    if ("cn" %in% colnames(gg$nodes$dt)) {
        og.nodes$marginal.cn = gg$nodes$dt$cn
    }

    ## add phase block information to nodes
    if (!is.null(phase.blocks)) {
        if (verbose) {
            message("Adding phase block information")
        }
        pblocks = phase.blocks[, c()]
        pblocks$pblock = 1:length(pblocks) %>% as.character
        og.nodes = og.nodes %$% pblocks[, "pblock"]
    }

    phased.gg.nodes = c(og.nodes, og.nodes)
    phased.gg.nodes$allele = c(rep("major", length(og.nodes)),
                               rep("minor", length(og.nodes)))

    if (verbose) {
        message("Computing allele CN")
    }
    
    allele.cn.bins = bins[, c()]
    allele.cn.bins$cn = reads.to.allele.cn(bins, count.field, purity, ploidy)
    allele.cn.bins$allele = values(bins)[[allele.field]]

    ## overlap nodes with bins
    ov = gr.findoverlaps(phased.gg.nodes, allele.cn.bins,
                         by = c("allele"),
                         scol = c("cn"),
                         return.type = "data.table",
                         mc.cores = mc.cores)

    dt = ov[, .(mean = mean(cn, na.rm = T), var = var(cn, na.rm = T), nbins = .N), by = query.id]
    ov.match = match(1:length(phased.gg.nodes), dt$query.id)

    ## add information to nodes
    phased.gg.nodes$cn = dt$mean[ov.match]
    phased.gg.nodes$var = dt$var[ov.match]
    phased.gg.nodes$nbins = dt$nbins[ov.match]

    ## set weight to inverse variance
    phased.gg.nodes$weight = ifelse((phased.gg.nodes$var > 0) &
                                    (phased.gg.nodes$nbins > min.bins) &
                                    (phased.gg.nodes$var > min.var),
                                    (phased.gg.nodes$nbins / sqrt(phased.gg.nodes$var)),
                                    NA)

    if (verbose) {
        message("Preparing phased gGraph edges")
    }

    phased.gg.edges = rbind(
        gg$edges$dt[, .(n1, n2, n1.side, n2.side, type,
                        og.edge.id = edge.id,
                        n1.allele = "major",
                        n2.allele = "major")],
        gg$edges$dt[, .(n1 = n1 + length(og.nodes), n2 = n2 + length(og.nodes), type,
                        n1.side, n2.side,
                        og.edge.id = edge.id,
                        n1.allele = "minor",
                        n2.allele = "minor")],
        gg$edges$dt[, .(n1, n2 = n2 + length(og.nodes), type,
                        n1.side, n2.side,
                        og.edge.id = edge.id,
                        n1.allele = "major",
                        n2.allele = "minor")],
        gg$edges$dt[, .(n1 = n1 + length(og.nodes), n2, type,
                        n1.side, n2.side,
                        og.edge.id = edge.id,
                        n1.allele = "minor",
                        n2.allele = "major")]
    )
    
    ## add n1/n2 chromosome information
    phased.gg.edges[, ":="(n1.chr = seqnames(phased.gg.nodes)[n1] %>% as.character,
                           n2.chr = seqnames(phased.gg.nodes)[n2] %>% as.character)]

    ## add edge connection type (straight/cross)
    phased.gg.edges[n1.chr == n2.chr & n1.allele == n2.allele, connection := "straight"]
    phased.gg.edges[n1.chr == n2.chr & n1.allele != n2.allele, connection := "cross"]

    ## add phase block information to edges (for linked reads)
    if (!is.null(phase.blocks)) {
        phased.gg.edges[, ":="(n1.pblock = phased.gg.nodes$pblock[n1],
                           n2.pblock = phased.gg.nodes$pblock[n2])]

        ## fix cross REF edges to zero within phase blocks
        phased.gg.edges[(n1.pblock == n2.pblock) & type == "REF" & connection == "cross",
                    ":="(cn = 0, fix = 1)]
        if (verbose) {
            message("Number of REF cross edges within phased blocks: ",
                    nrow(phased.gg.edges[(n1.pblock == n2.pblock) &
                                         type == "REF" &
                                         connection == "cross"]))
        }
    }

    ## identify phased edges (for linked reads)
    if (!is.null(edge.phase.dt)) {

        ## compute totals
        ephase = edge.phase.dt[, .(edge.id, n1.major, n2.major, n1.minor, n2.minor,
                                   n1.total = n1.major + n1.minor,
                                   n2.total = n2.major + n2.minor)][
                                       (n1.total > vbase.count.thres) | (n2.total > vbase.count.thres)]

        ## count fraction of reads corresponding with each allele
        ephase[, n1.major.frac := n1.major / n1.total]
        ephase[, n2.major.frac := n2.major / n2.total]
        ephase[, n1.minor.frac := n1.minor / n1.total]
        ephase[, n2.minor.frac := n2.minor / n1.total]

        ## set phase if passing proportion threshold (vbase.prop.thres)
        ephase[n1.major.frac > vbase.prop.thres, n1.phase := "major"]
        ephase[n1.minor.frac > vbase.prop.thres, n1.phase := "minor"]
        ephase[n2.major.frac > vbase.prop.thres, n2.phase := "major"]
        ephase[n2.minor.frac > vbase.prop.thres, n2.phase := "minor"]

        ## add phase information to edges data frame
        phased.gg.edges[, n1.phase := ephase$n1.phase[match(og.edge.id, ephase$edge.id)]]
        phased.gg.edges[, n2.phase := ephase$n2.phase[match(og.edge.id, ephase$edge.id)]]

        ## fix things to zero
        phased.gg.edges[n1.phase == "major" & n1.allele == "minor", ":="(fix = 1, cn = 0)]
        phased.gg.edges[n2.phase == "major" & n2.allele == "minor", ":="(fix = 1, cn = 0)]
        phased.gg.edges[n1.phase == "minor" & n1.allele == "major", ":="(fix = 1, cn = 0)]
        phased.gg.edges[n2.phase == "minor" & n2.allele == "major", ":="(fix = 1, cn = 0)]

        if (verbose) {
            message("Number of ALT edges with n1 side fixed: ", sum(!is.na(phased.gg.edges$n1.phase)))
            message("Number of ALT edges with n2 side fixed: ", sum(!is.na(phased.gg.edges$n2.phase)))
        }
    }

    if (verbose) {
        message("Creating gGraph")
        message("Number of nodes: ", length(phased.gg.nodes))
        message("Number of REF edges: ", nrow(phased.gg.edges[type == "REF"]))
        message("Number of ALT edges: ", nrow(phased.gg.edges[type == "ALT"]))
    }

    phased.gg = gG(nodes = phased.gg.nodes, edges = phased.gg.edges)

    if (verbose) {
        message("Formatting gGraph")
    }

    ref.edge.col = alpha("blue", 0.3)
    alt.edge.col = alpha("red", 0.3)
    ref.edge.lwd = 0.5
    alt.edge.lwd = 1.0
    phased.gg$edges$mark(col = ifelse(phased.gg$edges$dt$type == "REF", ref.edge.col, alt.edge.col),
                         lwd = ifelse(phased.gg$edges$dt$type == "REF", ref.edge.lwd, alt.edge.lwd))

    major.node.col = alpha("red", 0.5)
    minor.node.col = alpha("blue", 0.5)
    phased.gg$nodes$mark(col = ifelse(phased.gg$nodes$dt$allele == "major", major.node.col, minor.node.col),
                         ywid = 0.8)

    return(phased.gg)
}

#' @name fitcn
#' @title fitcn
#' @author Julie Behr, Xiaotong Yao
#'
#' @param gw input gWalks
#' @param cn.field character column names of each graph's CN data
#' 
#' @export
fitcn = function (gw, cn.field = "cn", trim = TRUE, weight = NULL, obs.mat = NULL, verbose = TRUE, 
                  min.alt = TRUE, edgeonly = FALSE, evolve = FALSE, n.sol = 2, return.gw = TRUE,
                  sep = "_")
{
    ## gw = self$copy
    gw = gw$copy
    gg = gw$graph$copy
    stopifnot(all(cn.field %in% colnames(gg$nodes$dt)) & all(cn.field %in% colnames(gg$edges$dt)))
    ## if (is.null(gw$graph$nodes$dt$cn) | is.null(gw$graph$edges$dt$cn)) {
    ##     stop("cn field is missing from node and edge metadata")
    ## }
    gg$nodes$mark(cn = rowSums(as.matrix(gg$nodes$dt[, cn.field, with = FALSE])))
    rcn = sapply(cn.field, function(colnm) gg$nodes$eval(sum(cn)))
    lcn = sapply(cn.field, function(colnm) gg$nodes$eval(sum(cn), right = FALSE))
    rcn = gg$nodes$eval(sum(cn))
    lcn = gg$nodes$eval(sum(cn), right = FALSE)
    jbal.left = all(
    (lcn == rowSums(as.matrix(gg$nodes$dt[, cn.field, with = FALSE])) & !gg$nodes$dt$loose.left) | 
    (lcn <= rowSums(as.matrix(gg$nodes$dt[, cn.field, with = FALSE])) & gg$nodes$dt$loose.left),
    na.rm = TRUE)
    jbal.right = all(
    (rcn == rowSums(as.matrix(gg$nodes$dt[, cn.field, with = FALSE])) & !gg$nodes$dt$loose.right) |
    (rcn <= rowSums(as.matrix(gg$nodes$dt[, cn.field, with = FALSE])) & gg$nodes$dt$loose.right),
    na.rm = TRUE)
    if (!jbal.left | !jbal.right) 
        warning("graph does not appear to be junction balanced, please check inputs")
    ## helper function to get the problem right
    constrain.evolution = function(K, gw, A, b, sense){
        h = K[gw$edges[type == "ALT"]$dt[!duplicated(edge.id), 
                                         edge.id], ]
        A = rbind(A, cbind(sparseMatrix(1, 1, x = 0, dims = dim(h)), 
                           h))
        b = c(b, rep(1, nrow(h)))
        sense = c(sense, rep("L", nrow(h)))
        return(list(A = A, b = b, sense = sense))
    }
    constrain.observations = function(obs.mat, A, b, cvec, sense, 
                                      vtype) {
        if (!(ncol(obs.mat) * 2) == ncol(A)) 
            stop("input obs.mat contains the wrong number of columns; should match length of gw")
        p = nrow(obs.mat)
        w = ncol(obs.mat)
        Zero = sparseMatrix(1, 1, x = 0, dims = c(2 * w * p, 
                                                  2 * w * p))
        A0 = Zero[rep(1, nrow(A)), 1:(2 * p)]
        Ap = cbind(Zero[rep(1, p), 1:w], sign(obs.mat), diag(rep(-1, 
                                                                 p)), Zero[rep(1, p), 1:p])
        Mpub = cbind(Zero[rep(1, p), 1:(2 * w)], diag(rep(1, 
                                                          p)), diag(rep(-1e+07, p)))
        Mplb = cbind(Zero[rep(1, p), 1:(2 * w)], diag(rep(1, 
                                                          p)), diag(rep(-0.1, p)))
        Amp = rbind(cbind(A, A0), Ap, Mpub, Mplb)
        b = c(b, rep(0, 3 * p))
        cvec = c(cvec, rep(0, p), -1 * rowMax(obs.mat))
        sense = c(sense, rep("E", p), rep("L", p), rep("G", p))
        vtype = c(vtype, rep("I", p), rep("B", p))
        return(list(A = Amp, b = b, c = cvec, sense = sense, 
                    vtype = vtype))
    }
    generate.Ke = function(gw) {
        dt = gw$edgesdt[, c("walk.id", "sedge.id")][, `:=`(edge.id, abs(sedge.id))]
        dt$listid = factor(dt$walk.id, 1:length(gw))
        dt$edge.id = factor(dt$edge.id, gg$edgesdt$edge.id)
        cdt = dcast(dt[!is.na(sedge.id), ], listid ~ edge.id, 
                    fun.aggregate = length, value.var = "edge.id", drop = FALSE)
        mat = cdt[, -1]
        rownames(mat) = cdt$listid
        return(t(mat))
    }
    generate.Kn = function(gw) {
        dt = gw$nodesdt[, c("walk.id", "snode.id")][, `:=`(node.id, 
                                                           abs(snode.id))]
        dt$listid = factor(dt$walk.id, 1:length(gw))
        dt$node.id = factor(dt$node.id, gg$nodes$dt$node.id)
        cdt = dcast(dt[!is.na(snode.id), ], listid ~ node.id, 
                    fun.aggregate = length, value.var = "node.id", drop = FALSE)
        mat = cdt[, -1]
        rownames(mat) = cdt$listid
        return(t(mat))
    }
    generate.Amat = function(K, nblock = 1) {
        M = 1e+07
        K = as(K, "sparseMatrix")
        w = ncol(K)
        Zero = sparseMatrix(1, 1, x = 0, dims = c(2 * w * nblock, 2 * w * nblock))
        Amub = cbind(
            do.call(`cbind`, lapply(seq_len(nblock), function(i) diag(rep(1, w)))), diag(rep(-M, w))
        )
        Amlb = cbind(
            do.call(`cbind`, lapply(seq_len(nblock), function(i) diag(rep(1, w)))), diag(rep(-0.1, w))
        )
        A = rbind(
            ## cbind(K, Zero[rep(1, nrow(K)), (w + 1:w)]),
            cbind(Reduce(`diagc`, lapply(seq_len(nblock), function(i) K)), Zero[rep(1, nrow(K) * nblock), (w + 1:w)]),
            Amub,
            Amlb
        )
        return(A)
    }
    generate.bvec = function(e, K) {
        w = ncol(K)
        bvec = c(c(e), rep(0, 2 * w))
        return(bvec)
    }
    generate.cvec = function(K, weight, min.alt, gw, nblock = 1) {
        if (!is.null(weight) | min.alt) 
            weight = prep.weight(K, weight, min.alt, gw)
        w = ncol(K)
        if (is.null(weight)) 
            weight = rep(1, w)
        cvec = c(rep(0, w * nblock), weight)
        return(cvec)
    }
    prep.weight = function(K, weight, min.alt, gw) {
        if (!is.null(weight)) {
            if (length(weight) == 1 & is.character(weight) & 
                weight %in% colnames(gw$dt)) 
                weight = gw$dt[, weight, with = F]
            if (!(is.numeric(weight))) {
                stop("weight must either be numeric vector of same length as gw or the name of a single numeric annotation in gw")
            }
        }
        if (min.alt) {
            if (!is.null(weight)) {
                warning("modifying input weight to satisfy min.alt=TRUE")
            }
            else weight = rep(1, ncol(K))
            numalt = gw$eval(edge = sum(type == "ALT"))
            weight[is.na(numalt) | numalt == 0] = 0
        }
        return(weight)
    }
    generate.vtype = function(K, nblock = 1) {
        w = ncol(K)
        vtype = c(rep("I", w * nblock), rep("B", w))
        return(vtype)
    }
    generate.sense = function(K, nblock = 1) {
        w = ncol(K)
        r = nrow(K)
        sense = c(rep("E", r * nblock), rep("L", w), rep("G", w))
        return(sense)
    }
    diagc = function(mat1, mat2){
        out = matrix(0, nrow = nrow(mat1) + nrow(mat2), ncol = ncol(mat1) + ncol(mat2))
        el1 = which(as.matrix(mat1)!=0, arr.ind = TRUE)
        el2 = el2.ori = which(as.matrix(mat2)!=0, arr.ind = TRUE)
        el2[, "row"] = el2.ori[, "row"] + nrow(mat1)
        el2[, "col"] = el2.ori[, "col"] + ncol(mat1)
        out[rbind(el1, el2)] = c(mat1[el1], mat2[el2.ori])
        return(out)
    }
    ## start building the problem
    if (any(!cn.field %in% colnames(gw$graph$edges$dt))) {
        stop("cn field must be populated in the input graph node and edges metadata")
        ## already stopped, what is this for??
        ## if (trim) {
        ##     e = rep(1, length(unique(abs(unlist(gw$sedge.id)))))
        ##     e2 = rep(1, length(unique(abs(unlist(gw$snode.id)))))
        ## }
        ## else {
        ##     e = rep(1, length(gw$graph$edges))
        ##     e2 = rep(1, length(gw$graph$nodes))
        ## }
    }
    else {
        if (trim) {
            e = as.matrix(
                gw$graph$edges[sedge.id %in% abs(unlist(gw$sedge.id))]$dt[, cn.field, with = FALSE])
            e2 = as.matrix(
                gw$graph$nodes[snode.id %in% abs(unlist(gw$snode.id))]$dt[, cn.field, with = FALSE])
        }
        else {
            e = as.matrix(gw$graph$edges$dt[, cn.field, with = FALSE])
            e2 = as.matrix(gw$graph$nodes$dt[, cn.field, with = FALSE])
        }
    }
    if (edgeonly) {
        K = generate.Ke(gw)
    }
    else {
        K = rbind(generate.Ke(gw), generate.Kn(gw))
        e = rbind(e, e2)
    }
    if (nrow(K) != nrow(e)) {
        stop("Mismatch between size of A matrix and length of b vector. Some edges in gw$graph are not covered by gw. If this was intended, try trim=TRUE")
     }
    K.unit = K
    K = Reduce(`diagc`, lapply(seq_len(ncol(e)), function(i) K))
    ## a bit of cheating, doing the old procedure column by column as vectors
    ## tmp = lapply(seq_len(ncol(e)), function(i){
    ##     this.e = e[, i, drop = TRUE]
    ##     A = generate.Amat(unname(K))
    ##     b = generate.bvec(this.e, K)        
    ##     sense = generate.sense(K)
    ##     return(list(A = A, b = b, sense = sense, vtype = vtype))
    ## })
    A = generate.Amat(unname(K.unit), ncol(e))
    b = generate.bvec(e, K.unit)
    sense = generate.sense(K.unit, ncol(e))
    vtype = generate.vtype(K.unit, ncol(e))
    ## restructure the objective function
    c = generate.cvec(K.unit, weight, min.alt, gw, ncol(e))
    ## row bind them into final constraints
    ## browser()
    if (evolve) {
        ll = constrain.evolution(K, gw, A, b, sense)
        A = ll$A
        b = ll$b
        sense = ll$sense
    }
    if (!is.null(obs.mat)) {
        ll = constrain.observations(obs.mat, A, b, c, sense, 
                                    vtype)
        A = ll$A
        b = ll$b
        c = ll$c
        sense = ll$sense
        vtype = ll$vtype
    }
    ## browser()
    ## load customized lb and ub if present ## TODO make ub lb specific to samples
    if (any(grepl("lb", colnames(gw$dt)))){
        if (identical(grep("lb", colnames(gw$dt), value = TRUE), "lb")){
            lb = c(rep(gw$dt$lb, ncol(e)), rep(0, ncol(K.unit)))
        } else {
            lb.cols = gsub("cn", "lb", cn.field)
            if (!all(lb.cols %in% colnames(gw$dt))){
                lb = 0
            }
            lb = do.call("c", c(lapply(lb.cols, function(x) gw$dt[[x]]), rep(0, length(gw))))
        }
    } else {
        lb = 0
    }
    if (any(grepl("ub", colnames(gw$dt)))){
        if (identical(grep("ub", colnames(gw$dt), value = TRUE), "ub")){
            ub = c(rep(gw$dt$ub, ncol(e)), rep(Inf, ncol(K.unit)))
        } else {
            ub.cols = gsub("cn", "ub", cn.field)
            if (!all(ub.cols %in% colnames(gw$dt))){
                ub = Inf
            }
            ub = do.call("c", c(lapply(ub.cols, function(x) gw$dt[[x]]), rep(Inf, length(gw))))
        }
    } else {
        ub = Inf
    }
    ## if (len(lb) != len(ub)){
    ##     lb = rep(lb, length.out = pmax(len(lb), len(ub)))
    ##     ub = rep(ub, length.out = pmax(len(lb), len(ub)))
    ## }
    ## TODO: implement lb and ub of walk CNs
    sol = Rcplex::Rcplex(
        cvec = c,
        Amat = A,
        bvec = b,
        sense = sense, 
        Qmat = NULL,
        lb = lb,
        ub = ub,
        n = n.sol,
        objsense = "min", 
        vtype = vtype,
        control = list(
            trace = ifelse(verbose >= 1, 1, 0),
            tilim = 100,
            epgap = 1))
    
    if (!is.null(sol$xopt)) {
        sol = list(sol)
    }
    if (length(sol) == 0) {
        stop("No solutions found satisfying given constraints")
    }
    ## what is this rerun for???
    ## to exhaust solutions
    rerun = T
    while (rerun) {
        z = sign(vtype == "B")
        P = do.call(rbind, lapply(sol, function(x) x$xopt * z))
        p = rowSums(P) - 1
        Ahat = rbind(A, P)
        bhat = c(b, p)
        sensehat = c(sense, rep("L", length(p)))
        sol.new = Rcplex::Rcplex(cvec = c, Amat = Ahat, bvec = bhat, 
            sense = sensehat, Qmat = NULL, lb = lb, ub = ub, 
            n = n.sol, objsense = "min", vtype = vtype, control = list(trace = ifelse(verbose >= 
                1, 1, 0), tilim = 100, epgap = 1))
        if (length(sol.new) == 0) {
            rerun = F
        }
        else {
            sol = c(sol, sol.new)
            if (length(sol) >= n.sol) {
                sol = sol[1:n.sol]
                rerun = F
            }
        }
    }
    ## return(sol)
    ## for (cnm in cn.field){
    ##     gw$set(eval(paste0("cn.", i)) = xmat[, cnm])
    ## }
    ## if (length(sol) > 1) {

    if (return.gw){
        for (i in seq_along(sol)){
            this.sol = sol[[i]]
            this.x = this.sol$xopt
            this.cnf = rep(c(cn.field, "indicator"), each = length(gw))
            this.ls = split(this.x, this.cnf)
            names(this.ls) = paste(names(this.ls), i, sep = sep)
            do.call(gw$set, this.ls)
        }
        return(invisible(gw))
    } else {
        return(sol)
    }
}
