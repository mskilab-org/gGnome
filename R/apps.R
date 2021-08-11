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
#' @param force.major (logical) force major allele CN to be >= minor allele CN (default FALSE)
#' @param force.alt (logical) force incorporation of ALT edges, only applicable for phasing (default TRUE)
#' @param cnloh (logical) allow CN LOH? only relevant if phasing = TRUE. default FALSE.
#' @param lp (logical) solve as linear program using abs value (default TRUE)
#' @param M  (numeric) big M constraint for L0 norm loose end penalty (default 1e3)
#' @param verbose (integer)scalar specifying whether to do verbose output, value 2 will spit out MIP (1)
#' @param tilim (numeric) time limit on MIP in seconds (10)
#' @param epgap (numeric) relative optimality gap threshhold between 0 and 1 (default 1e-3)

#' @param trelim (numeric) max size of uncompressed tree in MB (default 32e3)
#' @param nodefileind (numeric) one of 0 (no node file) 1 (in memory compressed) 2 (on disk uncompressed) 3 (on disk compressed) default 1
#' @param debug (logical) returns list with names gg and sol. sol contains full RCPLEX solution. (default FALSE)
#' @param gurobi (logical) use gurobi if TRUE uses gurobi else CPLEX default FALSE
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
                   force.major = FALSE,
                   force.alt = TRUE,
                   cnloh = FALSE,
                   lp = TRUE,
                   verbose = 1,
                   tilim = 10,
                   trelim = 32e3,
                   nodefileind = 1,
                   epgap = 1e-3,
                   debug = FALSE)
{
    if (verbose) {
        message("creating copy of input gGraph")
    }

    gg = gg$copy

    if (verbose) {
        message("Checking inputs")
    }

    if (nodefileind) {
        if (!(nodefileind %in% c(0,1,2,3))) {
            warning("Invalid choice for nodefileind, resetting to default 1")
            nodefileind = 1
        }
    }
    nodefileind = as.integer(nodefileind)

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
    if (verbose) {
      message("adding l0 penalty indicator")
    }

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
        ## vars[type == "edge", emarginal.id := junction.map[abs(sedge.id), seen.by.emarginal]]
        vars[type == "edge", emarginal.id := junction.map[abs(sedge.id), subject.id]]
        ## add weight and target total CN
        emarginal = merge.data.table(unique(
            vars[type == "edge" & !is.na(emarginal.id),][, type := "emresidual"][, .(emarginal.id, sedge.id, lb = -M, ub = M, gid, type, vtype = "C", from, to)],
            by = "emarginal.id"),
            junction.map[, .(subject.id, weight, cn, fix)],
            by.x = "emarginal.id",
            by.y = "subject.id")
        vars = rbind(vars, emarginal, fill = TRUE)
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
                vars[type == "emresidual", .(emarginal.id, weight, vtype, type = "emdelta.plus")][, gid := emarginal.id],
                vars[type == "emresidual", .(emarginal.id, weight, vtype, type = "emdelta.minus")][, gid := emarginal.id]
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

        ## add telomeric annotation
        qtips = gr.end(si2gr(seqlengths(gg$nodes))) ## location of q arm tips
        term.in = c(which(start(gg$nodes$gr) == 1), ## beginning of chromosome
                    -which(gg$nodes$gr %^% qtips)) ## flip side of chromosome end
        term.out = -term.in ## out is reciprocal of in

        ## annotate loose indicators with this
        vars[!is.na(snode.id), telomeric := ifelse(snode.id %in% term.in | snode.id %in% term.out, TRUE, FALSE)]
        
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
            vars = merge.data.table(vars,
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
        #' zchoo Tuesday, Jun 15, 2021 11:53:15 AM
        #' this constraint appears to be valid even if running phasing.
        ## if (!phased) {
            ## extremity exclusivity (relevant for ALL graphs)
        loose.constraints = rbind(
            vars[type == "loose.in.indicator" & sign(snode.id) == 1 & telomeric == FALSE,
                 .(value = 1, id, cid = paste("extremity.exclusivity", ee.id))],
            vars[type == "loose.out.indicator" & sign(snode.id) == 1 & telomeric == FALSE,
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
            vars[type == "loose.in.indicator" & sign(snode.id) == 1 & ee.id %in% edge.ee.ids & telomeric == FALSE,
                 .(value = 1, id, cid = paste("extremity.exclusivity", ee.id))],
            vars[type == "loose.out.indicator" & sign(snode.id) == 1 & ee.id %in% edge.ee.ids & telomeric == FALSE,
                 .(value = 1, id, cid = paste("extremity.exclusivity", ee.id))]
        )

        loose.zeros.rhs = unique(loose.zeros[, .(cid, value = 0, sense = "E")], by = "cid")

        constraints = rbind(constraints, loose.zeros, fill = TRUE)
        b = rbind(b, loose.zeros.rhs, fill = TRUE)
        ## }

        if (phased) {
            ## homologous extremity exclusivity
            ## this is actually redundant with previous constraints
            loose.constraints = rbind(
                vars[type == "loose.in.indicator" & sign(snode.id)==1 & telomeric == FALSE,
                     .(value = 1, id, cid = paste("homol.extremity.exclusivity", hee.id))],
                vars[type == "loose.out.indicator" & sign(snode.id)==1 & telomeric == FALSE,
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
                vars[type == "loose.in.indicator" & sign(snode.id)==1 & telomeric == FALSE,
                     .(value = 1, sense = "L", cid = paste("homol.extremity.exclusivity", hee.id))],
                vars[type == "loose.out.indicator" & sign(snode.id)==1 & telomeric == FALSE,
                     .(value = 1, sense = "L", cid = paste("homol.extremity.exclusivity", hee.id))]
            ), by = "cid")
            
            b = rbind(b, rhs, fill = TRUE)

            if (verbose) {
                message("Number of homologous extremity exclusivity constraints: ",
                        nrow(rhs))
            }

            ## grab node ids associated with ALT edges on the left
            left.og.node.ids = c(gg$edges$dt[n1.side == "left" & type == "ALT", n1],
                                 gg$edges$dt[n2.side == "left" & type == "ALT", n2])
            right.og.node.ids = c(gg$edges$dt[n1.side == "right" & type == "ALT", n1],
                                   gg$edges$dt[n2.side == "right" & type == "ALT", n2])

            ## fix loose ends for these nodes to zero
            vars[type == "loose.in.indicator" & (snode.id %in% left.og.node.ids),
                 ":="(lb = 0, ub = 0)]
            vars[type == "loose.out.indicator" & (snode.id %in% right.og.node.ids),
                 ":="(lb = 0, ub = 0)]
            vars[type == "loose.in" & (snode.id %in% left.og.node.ids),
                 ":="(lb = 0, ub = 0)]
            vars[type == "loose.out" & (snode.id %in% right.og.node.ids),
                 ":="(lb = 0, ub = 0)]

            if (verbose) {
                message("Number of homologous loose ends: ",
                        length(left.og.node.ids) + length(right.og.node.ids))
            }
        
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
                vars[(!cnloh == TRUE) & type == "edge.indicator" & sedge.id > 0 & ref.or.alt == "ALT" & !is.na(s1),
                     .(value = 1, id, cid = paste("rhee", s1))],
                vars[(!cnloh == TRUE) & type == "edge.indicator" & sedge.id > 0 & ref.or.alt == "ALT" & !is.na(s2),
                     .(value = 1, id, cid = paste("rhee", s2))],
                vars[(!cnloh == TRUE) & type == "edge.indicator" & sedge.id > 0 & ref.or.alt == "ALT" & !is.na(s3),
                     .(value = 1, id, cid = paste("rhee", s3))],
                vars[(!cnloh == TRUE) & type == "edge.indicator" & sedge.id > 0 & ref.or.alt == "ALT" & !is.na(s4),
                     .(value = 1, id, cid = paste("rhee", s4))],
                vars[(!cnloh == TRUE) & type == "edge.indicator" & sedge.id > 0 & ref.or.alt == "ALT" & !is.na(c1),
                     .(value = 1, id, cid = paste("rhee", c1))],
                vars[(!cnloh == TRUE) & type == "edge.indicator" & sedge.id > 0 & ref.or.alt == "ALT" & !is.na(c2),
                     .(value = 1, id, cid = paste("rhee", c2))],
                vars[(!cnloh == TRUE) & type == "edge.indicator" & sedge.id > 0 & ref.or.alt == "ALT" & !is.na(c3),
                     .(value = 1, id, cid = paste("rhee", c3))],
                vars[(!cnloh == TRUE) & type == "edge.indicator" & sedge.id > 0 & ref.or.alt == "ALT" & !is.na(c4),
                     .(value = 1, id, cid = paste("rhee", c4))],

                ## loose indicators
                vars[type == "loose.in.indicator" & snode.id > 0 & !is.na(s) & telomeric == FALSE,
                     .(value = 1, id, cid = paste("rhee", s))],
                vars[type == "loose.in.indicator" & snode.id > 0 & !is.na(c) & telomeric == FALSE,
                     .(value = 1, id, cid = paste("rhee", c))],
                vars[type == "loose.out.indicator" & snode.id > 0 & !is.na(s) & telomeric == FALSE,
                     .(value = 1, id, cid = paste("rhee", s))],
                vars[type == "loose.out.indicator" & snode.id > 0 & !is.na(c) & telomeric == FALSE,
                     .(value = 1, id, cid = paste("rhee", c))]
                )

            ## filter constraints to only include things with >= 4 entries (e.g. must have an ALT edge)
            ## rhomol.constraints[, n.entries := .N, by = cid]
            ## remove this filter! due to some loose end violations!
            ## rhomol.constraints = rhomol.constraints[n.entries > 3, .(value, id, cid)]

            rhs = unique(rhomol.constraints[, .(value = 2, sense = "L", cid)], by = "cid")

            if (verbose) {
                message("Number of reciprocal homologous constraints: ", nrow(rhs))
            }

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

        ## force major allele to have higher CN than minor allele
        ## may not work for phased blocks
        if (force.major) {

            iconstraints = rbind(
                vars[type == "node" & allele == "major" & snode.id > 0,
                     .(value = 1, id, cid = paste("force.major", og.node.id))],
                vars[type == "node" & allele == "minor" & snode.id > 0,
                     .(value = -1, id, cid = paste("force.major", og.node.id))])

            rhs = unique(vars[type == "node" & snode.id > 0 & allele == "major",
                              .(value = 0, sense = "G", cid = paste("force.major", og.node.id))],
                         by = "cid")

            constraints = rbind(constraints, iconstraints, fill = TRUE)
            b = rbind(b, rhs, fill = TRUE)
            
        }
                

        ## force nonzero CN for ALT edges (because these have nonzero CN in original JaBbA output)
        ## can become infeasible if original graph is not compatible with ISM
        if (force.alt) {

            if (ism) {
                warning("Forcing ALT edges while running ISM can make some problems infeasible!")
            }
            
            iconstraints = unique(
                vars[type == "edge.indicator" & ref.or.alt == "ALT" & cnloh != TRUE,
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
                vars[type == "edge.indicator" & ref.or.alt == "ALT" & cnloh != TRUE,
                     .(value = 1, sense = "G", cid = paste("edge.indicator.sum.lb", og.edge.id))],
                by = "cid"
            )

            b = rbind(b, edge.indicator.b, fill = TRUE)
        }

        ## REF edge configuration constraint (added by default basically)
        ## only add this if there are no unphased nodes
        if (cnloh) {
            
            ## if allow CNLOH, the sum of edge indicators corresponding with og edge id is LEQ 2
            ## this is only allowed in constant CN regions and if breakpoint is not shared with any ALT edges

            ## penalize CNLOH edges

            if (!is.null(gg$edges$dt$cnloh)) {
                cnloh.edges = gg$edges$dt[cnloh == TRUE & type == "ALT", edge.id] %>% unique
                if (verbose) {
                    message("Number of marked CNLOH edges: ", length(cnloh.edges))
                }

                ## add CNLOH annotation to variables
                ## browser()
                vars[, cnloh := FALSE]
                vars[(type == "edge.indicator" | type == "edge" | type == "eresidual") &
                     ref.or.alt == "ALT" & (abs(sedge.id) %in% cnloh.edges),
                     ":="(cnloh = TRUE)]
                
            } else {
                warning("CNLOH not specified on edges. Disallowing!")
                cnloh.og.edges = c()
                vars[, cnloh := FALSE]
            }
        } else {

            cnloh.og.edges = c()
            vars[, cnloh := FALSE]
            
        }

        ## add CNLOH constraints for applicable edges
        
        ## iconstraints = unique(
        ##     vars[type == "edge.indicator" & ref.or.alt == "REF" & og.edge.id %in% cnloh.og.edges,
        ##          .(value = 1, id, edge.id = abs(sedge.id),
        ##            cid = paste("ref.configuration.constraint.cnloh", og.edge.id))],
        ##     by = "edge.id"
        ## )

        ## rhs = unique(
        ##     vars[type == "edge.indicator" & ref.or.alt == "REF" & og.edge.id %in% cnloh.og.edges,
        ##          .(value = 2, sense = "L",
        ##            cid = paste("ref.configuration.constraint.cnloh", og.edge.id))],
        ##     by = "cid"
        ## )

        ## constraints = rbind(constraints,
        ##                     iconstraints[, .(value, id, cid)],
        ##                     fill = TRUE)
        ## b = rbind(b, rhs, fill = TRUE)

        ## add ISM constraints for ALL REF edges (as CNLOH is now marked as ALT)
        
        iconstraints.from = unique(
            vars[type == "edge.indicator" & ref.or.alt == "REF", ##& !(og.edge.id %in% cnloh.og.edges),
                 .(value = 1, id,
                   edge.id = abs(sedge.id),
                   snode.id = from, ## this is actually a misleading name because from is the row in gg$dt
                   cid = paste("ref.configuration.constraint.from", from))],
            by = "edge.id"
        )

        iconstraints.to = unique(
            vars[type == "edge.indicator" & ref.or.alt == "REF", ##& !(og.edge.id %in% cnloh.og.edges),
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
                merge.data.table(vars[type == 'node', !"rid"], ov, by = 'snode.id')[, .(value = 1, id , cid = paste('mresidual', rid))],
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

  ## message("CVEC: ", length(cvec))

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

    ## browser()
    if (cnloh) {

        if ("cnloh" %in% colnames(vars)) {
            indices = which(vars$type == "edge.indicator" & !is.na(vars$cnloh) & vars$cnloh == TRUE)
            cvec[indices] = lambda

            message("Number of penalized CNLOH edges: ", length(indices))
        }
    }

    ## check constraints of CNLOH
    ## browser()
    ## vars[type == "edge.indicator" & cnloh == TRUE]
    ## vars[type == "edge.indicator" & cnloh == TRUE, .N, by = og.edge.id]
  

  lb = vars$lb
  ub = vars$ub

  control = list(trace = ifelse(verbose>=2, 1, 0), tilim = tilim, epgap = epgap, round = 1, trelim = trelim, nodefileind = nodefileind)
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
  ppc = function(x) (x %>% merge.data.table(vars, by = 'id') %>% merge.data.table(b, by = 'cid.char'))[, paste(paste(round(value.x, 1), '*', paste(type, gid, sep=  '_'), '(', signif(x, 2), ')', collapse = ' + '), ifelse(sense[1] == 'E', '=', ifelse(sense[1] == 'G', '>=', '<=')), round(value.y[1],2)), by = cid.char]
  
  ppv = function(x) {tmp = x %>% merge.data.table(constraints, by = 'id'); constraints[cid %in% tmp$cid, ] %>% ppc}

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
    ref.edge.col = alpha("blue", 0.2)
    alt.edge.col = alpha("red", 0.4)
    ref.edge.lwd = 0.5
    alt.edge.lwd = 1.0
    edge.col = ifelse(gg$edges$dt$type == "REF", ref.edge.col, alt.edge.col)
    edge.lwd = ifelse(gg$edges$dt$type == "REF", ref.edge.lwd, alt.edge.lwd)
    gg$edges$mark(col = edge.col, lwd = edge.lwd)

    ## mark zero cn edges
    zero.cn.col = alpha("gray", 0)
    zero.cn.lwd = 0.5
    zero.cn.edges = which(gg$edges$dt$cn == 0)
    gg$edges[zero.cn.edges]$mark(col = zero.cn.col, lwd = zero.cn.lwd)
  }
    if (debug) {
        return(list(gg = gg, sol = sol))
    }    
    return(gg)

}

## #' @name jbaLP
## #' @title jbaLP
## #' @description jbaLP
## #'
## #' Simple (probably temporary) wrapper around balance for JaBbA LP
## #' Takes karyograph as input and balances it.
## #'
## #' @param kag.file (character)
## #' @param kag (karyograph object)
## #' @param cn.field (character) column in karyograph with CN guess, default cnmle
## #' @param var.field (character) column in karyograph with node variance estimate, default loess.var
## #' @param bins.field (character) column in karyograph containing number of bins, default nbins
## #' @param min.var (numeric) min allowable variance default 1e-3
## #' @param min.bins (numeric) min allowable bins default 5
## #' @param lambda (numeric) slack penalty, default 100
## #' @param L0 (logical) default TRUE
## #' @param loose.collapse (logical) default FALSE
## #' @param M (numeric) max CN (default 1e3)
## #' @param verbose (numeric) 0 (nothing) 1 (everything  MIP) 2 (print MIP), default 2 print MIP
## #' @param tilim (numeric) default 1e3
## #' @param ism (logical) add infinite site assumption constraints? default TRUE
## #' @param epgap (numeric) default 1e-3
## #'
## #' @return
## #' karyograph with modified segstats/adj. Adds fields epgap, cl, ecn.in, ecn.out, eslack.in, eslack.out to $segstats and edge CNs to $adj
## #' 
## #' @author Marcin Imielinski, Zi-Ning Choo
## #' @export
## jbaLP = function(kag.file = NULL,
##                  kag = NULL,
##                  cn.field = "cnmle",
##                  var.field = "loess.var",
##                  bins.field = "nbins",
##                  min.var = 1e-3,
##                  min.bins = 1,
##                  lambda = 100,
##                  L0 = TRUE,
##                  loose.collapse = FALSE,
##                  M = 1e3,
##                  verbose = 2,
##                  tilim = 1e3,
##                  ism = TRUE,
##                  epgap = 1e-3)
## {
##     if (is.null(kag.file) & is.null(kag)) {
##         stop("one of kag or kag.file must be supplied")
##     }
##     if (!is.null(kag.file) & !is.null(kag)) {
##         warning("both kag.file and kag supplied. using kag.")
##     }
##     if (!is.null(kag)) {
##         if (verbose) {
##             message("using supplied karyograph")
##         }
##     } else {
##         if (file.exists(kag.file)) {
##             if (verbose) {
##                 message("reading karyograph from file")
##             }
##             kag = readRDS(kag.file)
##         } else {
##             stop("kag.file does not exist and kag not supplied")
##         }
##     }
##     kag.gg = gG(jabba = kag)

##     if (verbose) {
##         message("Marking nodes with cn contained in column: ", cn.field)
##     }
    
##     if (is.null(values(kag.gg$nodes$gr)[[cn.field]])) {
##         stop("karyograph must have field specified in cn.field")
##     }
##     kag.gg$nodes$mark(cn  = values(kag.gg$nodes$gr)[[cn.field]])

##     if (verbose) {
##         message("Computing node weights using variance contained in column: ", var.field)
##     }
    
##     if (is.null(values(kag.gg$nodes$gr)[[var.field]]) | is.null(values(kag.gg$nodes$gr)[[bins.field]])) {
##         warning("karyograph missing var.field. setting weights to node widths")
##         wts = width(kag.gg$nodes$gr)
##     } else {
        
##         ## process variances
##         vars = values(kag.gg$nodes$gr)[[var.field]]
##         vars = ifelse(vars < min.var, NA, vars) ## filter negative variances
##         sd = sqrt(vars) * kag$beta ## rel2abs the standard deviation

##         ## process bins
##         bins = values(kag.gg$nodes$gr)[[bins.field]]
##         bins = ifelse(bins < min.bins, NA, bins)

##         ## compute node weights
##         wts = bins / (sd / sqrt(2)) ## for consistency with Laplace distribution
##         wts = ifelse(is.infinite(wts) | is.na(wts) | wts < 0, NA, wts)
##     }
##     kag.gg$nodes$mark(weight = wts)
    
##     ## no edge CNs
##     kag.gg$edges$mark(cn = NULL)
##     kag.gg$nodes[cn > M]$mark(cn = NA, weight = NA)

##     if (verbose) {
##         message("Starting LP balance on gGraph with...")
##         message("Number of nodes: ", length(kag.gg$nodes))
##         message("Number of edges: ", length(kag.gg$edges))
##     }

##     res = balance(kag.gg,
##                   debug = TRUE,
##                   lambda = lambda,
##                   L0 = TRUE,
##                   verbose = verbose,
##                   tilim = tilim,
##                   epgap = epgap,
##                   lp = TRUE,
##                   ism = ism)
    
##     bal.gg = res$gg
##     sol = res$sol
    
##     ## just replace things in the outputs
##     ## this can create weird errors if the order of kag and bal.gg isn't the same
##     out = copy(kag)
##     new.segstats = bal.gg$gr
##     nnodes = length(new.segstats)

##     new.segstats$cl = 1 ## everything same cluster
##     new.segstats$epgap = sol$epgap ## add epgap from genome-side opt
    
##     ## weighted adjacency
##     adj = sparseMatrix(i = bal.gg$sedgesdt$from, j = bal.gg$sedgesdt$to,
##                        x = bal.gg$sedgesdt$cn, dims = c(nnodes, nnodes))
##     ## add the necessary columns
##     new.segstats$ecn.in = Matrix::colSums(adj)
##     new.segstats$ecn.out = Matrix::rowSums(adj)
##     target.less = (Matrix::rowSums(adj, na.rm = T) == 0)
##     source.less = (Matrix::colSums(adj, na.rm = T) == 0)
##     new.segstats$eslack.out[!target.less] = new.segstats$cn[!target.less] - Matrix::rowSums(adj)[!target.less]
##     new.segstats$eslack.in[!source.less] =  new.segstats$cn[!source.less] - Matrix::colSums(adj)[!source.less]
##     out$adj = adj

##     ## add metadata
##     out$segstats = new.segstats
##     out$status = sol$status
##     out$epgap = sol$epgap
##     return(out)
## }

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

#' @name jbaLP
#' @description jbaLP
#'
#' Simple (probably temporary) wrapper around balance for JaBbA LP
#' Reads karyograph.rds file and balances it using cnmle as CN estimate
#'
#' @param kag.file (character)
#' @param kag (karyograph object)
#' @param cn.field (character) column in karyograph with CN guess, default cnmle
#' @param var.field (character) column in karyograph with node weight guess, default
#' @param lambda (numeric) slack penalty, default 10
#' @param L0 (logical) default TRUE
#' @param loose.collapse (logical) default FALSE
#' @param M (numeric) max CN
#' @param verbose (numeric) 0 (nothing) 1 (everything but MIP) 2 (print the MIP), default 1
#' @param tilim (numeric) default 1e3
#' @param epgap (numeric) default 1e-3
#' @export
jbaLP = function(kag.file = NULL,
                 kag = NULL,
                 cn.field = "cnmle",
                 var.field = "raw.var",
                 min.var = 1e-3,
                 lambda = 10,
                 L0 = TRUE,
                 loose.collapse = FALSE,
                 M = 1e3,
                 verbose = 1,
                 tilim = 1e3,
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
    if (is.null(values(kag.gg$nodes$gr)[[cn.field]])) {
        stop("karyograph must have field specified in cn.field")
    }
    kag.gg$nodes$mark(cn  = values(kag.gg$nodes$gr)[[cn.field]])
    if (is.null(values(kag.gg$nodes$gr)[[var.field]])) {
        warning("karyograph missing var.field. setting weights to node widths")
        wts = width(kag.gg$nodes$gr)
    } else {
        vars = values(kag.gg$nodes$gr)[[var.field]]
        wts = ifelse(vars > min.var, 1/vars, NA)
    }
    kag.gg$nodes$mark(weight = wts)
    ## no edge CNs
    kag.gg$edges$mark(cn = NA)
    ## NA all the really big nodes, otherwise possibly feasibility issues
    kag.gg$nodes[cn > M]$mark(cn = NA, weight = NA)
    if (verbose) {
        message("Starting LP balance")
    }
    ## empirical lambda?
    res = balance(kag.gg, lambda = lambda, L0 = L0, loose.collapse = loose.collapse,
                     M = M, verbose = verbose, tilim = tilim, epgap = epgap, lp = TRUE,
                  ref.config = FALSE, phased = FALSE, marginal = NULL, debug = TRUE)
    bal.gg = res$gg
    sol = res$sol
    ## just replace things in the output
    out = copy(kag)
    new.segstats = bal.gg$gr
    nnodes = length(out$segstats)
    ## check if converged or just ran out of time
    if (sol$status == 1) {
        eg = epgap
    } else {
        eg = NA ## not sure how to extract optimality gap unfortunately
    }
    new.segstats$epgap = eg
    new.segstats$cl = 1 ## everything same cluster
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
    out$segstats = new.segstats
    return(out)
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
#' Attempts to embed / transplant  a set of circular walks (loops)
#' into a set of recipients, both defined on the same gGraph, and 
#' both optionally having the metadata $cn.  Returns a gWalk with as
#' many of the loops embedded into the recipients as possible.
#'
#' A loop needs to intersect some node in order to be successfully embedded.
#' (note: that node may be in another loop, if that loop gets embedded)
#' (As a result though, we don't guarantee that every loop will be enbedded)
#'
#' Note: loops will be embedded in order and loops with cn>1 will be embedded in tandem unless random = TRUE
#' in which case loops will be embedded randomly .. 
#' 
#'
#' @param loops gWalk object that may contain one or more circular walks and some linear ones as well
#' @param recipient set of gWalks defined on the same object (can also be circular) 
#' @return gWalk with as many loops embedded into recipients as possible
#' @export
#' @author Marcin Imielinski
embedloops = function(loops, recipients = loops[c()], random = FALSE, verbose = FALSE)
{
  if (length(loops)==0)
    return(recipients)

  ## move all non circular gWalks in loops into recipients 
  ix = !loops$dt$circular
  if (any(ix))
    {
      recipients = c(loops[ix], recipients)
      loops = loops[!ix]
    }
  
  if (all(!loops$dt$circular))
    return(recipients)

  ## if cn not set then set to 1
  if (is.null(loops$dt$cn))
    loops$set(cn = 1)

  ## if cn not set then set to 1
  if (is.null(recipients$dt$cn))
    recipients$set(cn = 1)

  if (random)
  {
    ix = rep(1:length(loops), loops$dt$cn) %>% sample
    loops = loops[ix]
    loops$set(cn = 1)
  }

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

#' @name find_na_ranges
#' @title find_na_ranges
#'
#' @description
#'
#' Identify NA ranges in a phased gGraph that actually cannot be phased
#'
#' @param gg (gGraph) junction-balanced allelic graph
#' @param min.bins (numeric) default 1
#' @param verbose (logical) default FALSE
#'
#' @return GRanges representing NA ranges that phased graph can be disjoined against in phased.postprocess
find_na_ranges = function(gg, min.bins = 1, verbose = FALSE) {

    if (!inherits(gg, "gGraph")) {
        stop("gg is not gGraph")
    }

    if (is.null(gg$nodes$dt$nbins)) {
        stop("gg nodes must have metadata $nbins")
    }
    
    gg.nodes.gr = gg$nodes$gr[, c("nbins")]
    gg.nodes.gr$na.node = is.na(gg.nodes.gr$nbins) | gg.nodes.gr$nbins < min.bins

    ## reduce NA nodes
    na.nodes.gr = gg.nodes.gr %Q% (na.node == TRUE)

    ## if there are magically not any of these, return
    if (length(na.nodes.gr) == 0) {
        return(GRanges())
    }
    
    na.nodes.gr = gr.reduce(na.nodes.gr, by = "na.node")

    ## grab left and right endpoints (use OG node id)
    na.nodes.gr$left.node.id = gr.findoverlaps(gr.start(na.nodes.gr), gg$nodes$gr[, "og.node.id"], first = TRUE, scol = "og.node.id")$og.node.id
    na.nodes.gr$right.node.id = gr.findoverlaps(gr.end(na.nodes.gr), gg$nodes$gr[, "og.node.id"], first = TRUE, scol = "og.node.id")$og.node.id

    ## check if these ranges are flanked by a rearrangement
    og.loose.dt = gg$nodes$dt[, .(loose.left = any(loose.left == TRUE, na.rm = TRUE), loose.right = any(loose.right == TRUE, na.rm = TRUE)), by = og.node.id]
    og.alt.dt = rbind(gg$edges$dt[type == "ALT", .(og.node.id = gg$nodes$dt$og.node.id[n1], side = n1.side)], gg$edges$dt[type == "ALT", .(og.node.id = gg$nodes$dt$og.node.id[n2], side = n2.side)])[, .(alt.left = any(side == "left", na.rm = TRUE), alt.right = any(side == "right", na.rm = TRUE)), by = og.node.id]

    ## annotate loose and ALTs
    na.nodes.gr$loose.left = na.nodes.gr$left.node.id %in% og.loose.dt[loose.left == TRUE, og.node.id]
    na.nodes.gr$loose.right = na.nodes.gr$right.node.id %in% og.loose.dt[loose.right == TRUE, og.node.id]
    na.nodes.gr$alt.left = na.nodes.gr$left.node.id %in% og.alt.dt[alt.left == TRUE, og.node.id]
    na.nodes.gr$alt.right = na.nodes.gr$right.node.id %in% og.alt.dt[alt.right == TRUE, og.node.id]

    ## next check if ranges are neighboring a rearrangement
    ## create a data table with left and right neighbor of every og node
    neighbors.dt = unique(gg$edges$dt[type == "REF", .(left.neighbor = gg$nodes$dt$og.node.id[n1], right.neighbor = gg$nodes$dt$og.node.id[n2])], by = "left.neighbor")

    ## annotate with og.node.id of left and right neighbor
    na.nodes.gr$left.neighbor = neighbors.dt$left.neighbor[match(na.nodes.gr$left.node.id, neighbors.dt$right.neighbor)]
    na.nodes.gr$right.neighbor = neighbors.dt$right.neighbor[match(na.nodes.gr$right.node.id, neighbors.dt$left.neighbor)]

    ## check if the neighbors of each node are flanked by a loose end
    na.nodes.gr$loose.left.neighbor = na.nodes.gr$left.neighbor %in% og.loose.dt[loose.right == TRUE, og.node.id]
    na.nodes.gr$loose.right.neighbor = na.nodes.gr$right.neighbor  %in% og.loose.dt[loose.left == TRUE, og.node.id]

    ## check if the neighbors of each node are flanked by an ALT edge
    na.nodes.gr$alt.left.neighbor = na.nodes.gr$left.neighbor %in% og.alt.dt[alt.right == TRUE, og.node.id]
    na.nodes.gr$alt.right.neighbor = na.nodes.gr$right.neighbor %in% og.alt.dt[alt.left == TRUE, og.node.id]

    ## get marginal copy number of left and right neighbor
    marginal.cn.dt = gg$nodes$dt[allele == "major", .(og.node.id, marginal.cn, major.cn = cn)]
    na.nodes.gr$left.neighbor.marginal.cn = marginal.cn.dt$marginal.cn[match(na.nodes.gr$left.neighbor, marginal.cn.dt$og.node.id)]
    na.nodes.gr$right.neighbor.marginal.cn = marginal.cn.dt$marginal.cn[match(na.nodes.gr$right.neighbor, marginal.cn.dt$og.node.id)]
    na.nodes.gr$left.neighbor.major.cn = marginal.cn.dt$major.cn[match(na.nodes.gr$left.neighbor, marginal.cn.dt$og.node.id)]
    na.nodes.gr$right.neighbor.major.cn = marginal.cn.dt$major.cn[match(na.nodes.gr$right.neighbor, marginal.cn.dt$og.node.id)]

    ## get major and marginal copy number of left and right endpoints
    na.nodes.gr$left.marginal.cn = marginal.cn.dt$marginal.cn[match(na.nodes.gr$left.node.id, marginal.cn.dt$og.node.id)]
    na.nodes.gr$right.marginal.cn = marginal.cn.dt$marginal.cn[match(na.nodes.gr$right.node.id, marginal.cn.dt$og.node.id)]
    na.nodes.gr$left.major.cn = marginal.cn.dt$major.cn[match(na.nodes.gr$left.node.id, marginal.cn.dt$og.node.id)]
    na.nodes.gr$right.major.cn = marginal.cn.dt$major.cn[match(na.nodes.gr$right.node.id, marginal.cn.dt$og.node.id)]

    ## easier manipulation
    na.nodes.dt = as.data.table(na.nodes.gr)

    ## filter by marginal CN (e.g. if left and right neighbors have to have the same marginal)
    #' zchoo Friday, Jul 30, 2021 11:28:55 AM
    ## removed this filter
    ## na.nodes.dt = na.nodes.dt[(left.neighbor.marginal.cn == right.neighbor.marginal.cn) | is.na(left.marginal.cn) | is.na(right.marginal.cn),]

    if (!nrow(na.nodes.dt)) {
        return(GRanges())
    }

    ## remove LOH ranges
    na.nodes.dt = na.nodes.dt[(left.neighbor.marginal.cn != left.neighbor.major.cn) |
                              (left.marginal.cn != left.major.cn) |
                              (right.neighbor.marginal.cn != right.neighbor.major.cn) |
                              (right.marginal.cn != right.major.cn),]
    
    if (!nrow(na.nodes.dt)) {
        return(GRanges())
    }

    ## check if left/right are telomeric
    na.nodes.dt[, left.telomeric := (start == 1)]
    na.nodes.dt[, sl := seqlengths(gg$nodes$gr)[match(seqnames, names(seqlengths(gg$nodes$gr)))]]
    na.nodes.dt[, right.telomeric := (end == sl)]

    ## resize
    na.nodes.dt[right.telomeric == FALSE & (alt.right.neighbor == TRUE | loose.right.neighbor == TRUE), end := end + 1]
    na.nodes.dt[left.telomeric == FALSE & (alt.left.neighbor == TRUE | loose.left.neighbor == TRUE), start := start - 1]

    ## filter by major CN (e.g. if left and right neighbors have LOH)
    return(dt2gr(na.nodes.dt[, .(seqnames, start, end)], seqlengths = seqlengths(gg$nodes$gr)))
}



#' @name phased.postprocess
#' @title phased.postprocess
#' @description
#'
#' Postprocess junction-balanced phased graph and creates unphased regions
#' This identifies regions without allelic CN imbalance
#'
#' @param gg junction-balanced phased gGraph. each node must have metadata og.nodes.id, allele, nbins, cn
#' @param min.bins (numeric) minimum number of bins to be marked as an NA node. default 1
#' @param phase.blocks (GRanges) granges of phase blocks from linked reads. default = NULL
#' @param mc.cores (int) number of cores. default = 8.
#' @param verbose (bool) verbose > 0 prints stuff. default 1.
#'
#' @return balanced gGraph with unphased nodes marked and compressed
#'
#' @export
phased.postprocess = function(gg, min.bins = 1, phase.blocks = NULL, mc.cores = 8, verbose = 1)
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

    ## identify unphased nodes to disjoin against
    ## browser()

    if (verbose) {
        message("Disjoining input graph against unphased node GRanges")
    }
    
    seed.dt = merge.data.table(gg$nodes$dt[allele == "major", .(seqnames, start, end, width,
                                                                major.cn = cn, og.node.id,
                                                                major.loose.left = loose.left,
                                                                major.loose.right = loose.right,
                                                                major.node.id = node.id)],
                               gg$nodes$dt[allele == "minor", .(og.node.id, minor.cn = cn,
                                                                minor.loose.left = loose.left,
                                                                minor.loose.right = loose.right,
                                                                minor.node.id = node.id)],
                               by = "og.node.id")

    ## check whether there is an allelic CN change on either the left or right
    ## sort nodes by start and end. nodes should be sorted, but just in case
    seed.dt = seed.dt %>% split(seed.dt$seqnames) %>% lapply(function(dt) {dt[order(start),]}) %>% rbindlist

    ## check for major or minor allele CN change on the left/right
    seed.dt[, major.prev := data.table::shift(major.cn, n = 1, type = "lag")]
    seed.dt[, major.next := data.table::shift(major.cn, n = 1, type = "lead")]
    seed.dt[, minor.prev := data.table::shift(minor.cn, n = 1, type = "lag")]
    seed.dt[, minor.next := data.table::shift(minor.cn, n = 1, type = "lead")]

    seed.dt[, major.left.cn := major.cn != major.prev]
    seed.dt[, major.right.cn := major.cn != major.next]
    seed.dt[, minor.left.cn := minor.cn != minor.prev]
    seed.dt[, minor.right.cn := minor.cn != minor.next]

    ## drop nodes with CN imbalance
    seed.dt = seed.dt[major.cn == minor.cn,]

    ## identify whether the major or minor node is joined to an ALT edge with non-zero CN
    nonzero.alt.left = c(gg$edges$dt[type == "ALT" & cn > 0 & n1.side == "left", n1],
                         gg$edges$dt[type == "ALT" & cn > 0 & n2.side == "left", n2])
    nonzero.alt.right = c(gg$edges$dt[type == "ALT" & cn > 0 & n1.side == "right", n1],
                          gg$edges$dt[type == "ALT" & cn > 0 & n2.side == "right", n2])

    seed.dt[, ":="(major.alt.left = major.node.id %in% nonzero.alt.left,
                   major.alt.right = major.node.id %in% nonzero.alt.right,
                   minor.alt.left = minor.node.id %in% nonzero.alt.left,
                   minor.alt.right = minor.node.id %in% nonzero.alt.right)]

    ## label telomeric
    sl = seqlengths(gg$nodes$gr)
    seed.dt[, left.telomeric := start == 1]
    seed.dt[, right.telomeric := end == sl[as.character(seqnames)]]

    ## check whether the left and right sides are EITHER loose
    seed.dt[, left.alt := (major.alt.left == TRUE | minor.alt.left == TRUE |
                           major.loose.left == TRUE | minor.loose.left == TRUE) &
                  (left.telomeric == FALSE)]
            ## (major.left.cn == TRUE | minor.left.cn == TRUE) & (left.telomeric == FALSE)]

    seed.dt[, right.alt := (major.alt.right == TRUE | minor.alt.right == TRUE |
                           major.loose.right == TRUE | minor.loose.right == TRUE) &
                  (right.telomeric == FALSE)]
            ## (major.right.cn == TRUE | minor.right.cn == TRUE) & (right.telomeric == FALSE)]

    ## shift the end points
    seed.dt[left.alt == TRUE & width > 1, start := start + 1]
    seed.dt[right.alt == TRUE & width > 1, end := end - 1]

    ## create GRanges
    seed.gr = dt2gr(seed.dt[, .(seqnames, start, end)], seqlengths = seqlengths(gg$nodes$gr), seqinfo = seqinfo(gg$nodes$gr))

    ## merge with NA ranges
    if (verbose) {
        message("Identifying NA ranges")
        }
    na.gr = find_na_ranges(gg, min.bins = min.bins)
    all.seed.gr = gr.reduce(c(seed.gr, na.gr))

    ## disjoin gGraph against this GRanges
    gg = gg$disjoin(all.seed.gr, collapse = FALSE)

    ## any new edges introduced have to be straight
    gg$edges[is.na(connection)]$mark(connection = "straight")

    ## fill in other metadata
    n1 = gg$edges$dt[, n1]
    n2 = gg$edges$dt[, n2]

    ## borrow CN from surrounding nodes
    n1.na = gg$edges$dt[is.na(cn), n1]
    gg$edges[is.na(cn)]$mark(cn = gg$nodes$dt$cn[match(n1.na, gg$nodes$dt$node.id)])

    ## fix loose end CNs
    gg = gGnome:::loosefix(gg)

    ## label n1/n2 allele and chromosome
    gg$edges$mark(n1.allele = gg$nodes$dt$allele[match(n1, gg$nodes$dt$node.id)])
    gg$edges$mark(n2.allele = gg$nodes$dt$allele[match(n2, gg$nodes$dt$node.id)])
    gg$edges$mark(n1.chr = gg$nodes$dt$seqnames[match(n1, gg$nodes$dt$node.id)])
    gg$edges$mark(n2.chr = gg$nodes$dt$seqnames[match(n2, gg$nodes$dt$node.id)])

    ## reset og.node.ids and og.edge.ids
    node.id.key = gg$nodes$dt[, .(seqnames, start, end, rg = paste0(seqnames, ":", start, "-", end), node.id)]
    node.id.key[, rg := as.integer(factor(rg))]
    gg$nodes$mark(og.node.id = node.id.key[, rg])

    ## reset og edge ids
    edge.id.key = gg$edges$dt[, .(n1, n1.side, n2, n2.side, type)]
    edge.id.key[, ":="(n1.og = node.id.key$rg[match(n1, node.id.key$node.id)],
                       n2.og = node.id.key$rg[match(n2, node.id.key$node.id)])]
    edge.id.key[, rg := paste(n1.og, n2.og, n1.side, n2.side, type)]
    edge.id.key[, rg := as.integer(factor(rg))]
    gg$edges$mark(og.edge.id = edge.id.key[, rg])


    ## identify nodes without CN imbalance
    if (verbose) {
        message("Identifying nodes without CN imbalance")
    }
    ## browser()
    og.node.balance = dcast.data.table(gg$nodes$dt[, .(og.node.id, allele, cn)], og.node.id ~ allele, value.var = "cn") %>% merge.data.table(gg$nodes$dt[, .(og.node.id, width)], by = "og.node.id", all.x = TRUE)

    og.node.balance[, cn.imbalance := (major != minor)]
    og.node.balance[, cn.total := (major + minor)]

    og.node.balance[, phased := ifelse(cn.imbalance == TRUE & width > 1, TRUE, FALSE)]

    ## mark these specifically as being allele-balanced
    og.node.balance[, ab := phased == FALSE]

    if (verbose) {
        message("Identifying nodes in NA stretches")
    }
    ## browser()
    na.node.ids = as.data.table(gg$nodes$gr[, c("og.node.id", "node.id")] %&% na.gr)[, og.node.id]

    ## mark na nodes as unphased
    og.node.balance[og.node.id %in% na.node.ids, phased := FALSE]

    ## annotate

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

    ## identify specifically allele-balanced nodes (vs. NA nodes)
    ab.og.nodes = og.node.balance[ab == TRUE, og.node.id]
    ab.major.nodes = gg$nodes$dt[og.node.id %in% ab.og.nodes & allele == "major", node.id]

    ## create new data.table for nodes
    new.nodes.dt = gg$nodes$dt[!(node.id %in% unphased.minor.nodes),]

    ## mark major nodes as unphased
    new.nodes.dt[node.id %in% unphased.major.nodes, allele := "unphased"]

    ## mark allele-balanced nodes specifically
    new.nodes.dt[node.id %in% ab.major.nodes, ab := TRUE]

    ## reset CN to total CN
    new.nodes.dt[node.id %in% unphased.major.nodes, cn := og.node.balance$cn.total[match(og.node.id, og.node.balance$og.node.id)]]

    ## fix the CN of all of these nodes
    ## new.nodes.dt[, fix := 1]

    ## reformat nodes
    new.nodes.dt[allele == "unphased", col := alpha("gray", 0.5)]

    ## get nodes as GRanges
    new.nodes.dt = new.nodes.dt %>% split(new.nodes.dt$seqnames) %>%
        lapply(function(dt) {dt[order(start),]}) %>% rbindlist
    new.nodes.gr = dt2gr(new.nodes.dt[, .(seqnames, start, end,
                                          og.node.id, marginal.cn, allele,
                                          var, nbins, weight, index, col,
                                          cn.old, cn, fix, ywid, ab,
                                          old.node.id = node.id)],
                         seqinfo = seqinfo(gg$nodes$gr),
                         seqlengths = seqlengths(gg$nodes$gr))

    ## create new data.table for edges
    new.edges.dt = gg$edges$dt
    
    ## reset edge endpoints
    ## browser()
    dt = gg$nodes$dt[, .(og.node.id, allele, node.id)] %>% dcast.data.table(og.node.id ~ allele, value.var = "node.id")
    new.edges.dt[(n1 %in% unphased.minor.nodes), n1 := dt$major[match(n1, dt$minor)]]
    new.edges.dt[(n2 %in% unphased.minor.nodes), n2 := dt$major[match(n2, dt$minor)]]

    ## reset all edge endpoints to new node.ids
    new.edges.dt[, n1 := match(n1, new.nodes.gr$old.node.id)]
    new.edges.dt[, n2 := match(n2, new.nodes.gr$old.node.id)]

    ## label REF edges as straight or cross based on og.node.id
    ## browser()
    ## new.edges.dt[type == "REF" & cn > 0, length(unique(connection)), by = og.edge.id] %>% summary
    new.edges.dt[type == "REF", orientation := .SD$connection[which(.SD$cn > 0)][1], by = og.edge.id]

    ## only keep REF edges in the correct orientation (regardless of CN)
    ## only keep ALT edges with CN > 0
    new.edges.dt = new.edges.dt[(type == "REF" & connection == orientation) |
                                (type == "ALT" & cn > 0),]

    ## deduplicate edges
    new.edges.dt = new.edges.dt[, .(connection = connection[1], type = type[1], cn = sum(cn, na.rm = TRUE), col = col[1]), by = .(n1, n1.side, n2, n2.side)]

    if (verbose) {
        message("Creating new gGraph")
    }
    new.nodes.gr = gGnome:::inferLoose(new.nodes.gr, new.edges.dt)

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
#' @param allele.field (str) field for containing major/minor allele label (default allele)
#' @param haplotype.field (str) field containing h1/h2 label (default haplotype)
#' @param phase.blocks (GRanges) GRanges containing phase blocks (e.g. from HAPCUT). default NULL.
#' @param breaks (GRanges) extra breakpoints to introduce. default NULL.
#' @param edge.phase.dt (data.table) with columns n1.major, n2.major, n1.minor, n2.minor and edge.id providing major/minor allele counts
#' @param vbase.count.thres (int) number of variant base counts required to phase edges (default 5)
#' @param vbase.prop.thres (float) proportion of allele excess required to phase edges (default 0.9)
#' @param min.bins (numeric) minimum number of bins for intra segment variance (default 3)
#' @param min.var (numeric) min allowable variance (default 0.1)
#' @param max.span (numeric) max span before penalizing CNLOH
#' @param verbose (bool) default TRUE for debugging
#' @param min.width (numeric) min allowable width for cnloh-adjacent node. default 1 Mbp
#' @param mc.cores (int) number of cores
#' @return gGraph whose nodes are annotated with $cn.major, $cn.minor, $haplotype, and $weight fields
#' @export
phased.binstats = function(gg, bins = NULL, purity = NULL, ploidy = NULL,
                           count.field = "count",
                           allele.field = "allele",
                           haplotype.field = "haplotype",
                           phase.blocks = NULL,
                           breaks = NULL,
                           edge.phase.dt = NULL,
                           vbase.count.thres = 5, vbase.prop.thres = 0.9,
                           min.bins = 3, min.var = 1e-3,
                           max.span = 1e6,
                           verbose = TRUE, min.width = 1e6, mc.cores = 8)
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
    if (!is.null(breaks)) {
        if (!inherits(breaks, 'GRanges')) {
            stop("breaks must be GRanges")
        }
    }
    if (!is.null(edge.phase.dt)) {
        if (!is.data.table(edge.phase.dt)) {
            warning("edge.phase.dt must be data.table. ignoring this input.")
            edge.phase.dt = NULL
        }
        if (!all(c("edge.id", "n1.h1", "n2.h1", "n1.h2", "n2.h2") %in% colnames(edge.phase.dt))) {
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
        message("Checking for extra breakpoints")
    }

    if (!is.null(breaks) && length(breaks) > 0) {
        ## get just the starting point and remove strand
        br = gr.stripstrand(gr.start(breaks))

        if (verbose) {
            message("Number of new breaks: ", length(br))
            message("Creating new gGraph with incorporated breaks...")
        }

        ## save old node IDs
        old.nodes = gg$nodes$gr[, "node.id"]
        old.nodes$unsplit.id = old.nodes$node.id

        ## create new gGraph with additional breakpoints
        new.nodes = gr.breaks(bps = br, query = gg$nodes$gr)
        gg = gG(breaks = new.nodes, junctions = gg$junctions[type == "ALT"])

        ## create node map to find which nodes were split up
        node.map = as.data.table(gg$nodes$gr[, "node.id"] %$% old.nodes[, "unsplit.id"])

        ## create a list of short nodes (CNLOH not allowed here!)
        short.nodes = as.data.table(gg$nodes$gr)[width < min.width, node.id]

        ## identify which REF edges are internal to an OG node
        edge.map = gg$edges$dt
        edge.map[, n1.unsplit := node.map$unsplit.id[match(n1, node.map$node.id)]]
        edge.map[, n2.unsplit := node.map$unsplit.id[match(n2, node.map$node.id)]]

        internal.edges = edge.map[type == "REF" &
                                  n1.unsplit == n2.unsplit &
                                  !(n1 %in% short.nodes) &
                                  !(n2 %in% short.nodes), edge.id]

        ## mark CNLOH in gg
        gg$edges$mark(cnloh = FALSE)
        gg$edges[internal.edges]$mark(cnloh = TRUE)

        if (verbose) {
            message("Number of internal edges marked in parent graph: ", length(internal.edges))
        }
    } else {
        gg$edges$mark(cnloh = FALSE)
    }

    if (verbose) {
        message("Marking pseudo-CNLOH")
    }

    ## pull intra-chromosomal ALT edges, because these are cnloh candidates
    pseudo.cnloh.junctions = gg$junctions[class %in% c("DEL-like", "INV-like", "DUP-like")]
    pseudo.cnloh.junctions.dt = pseudo.cnloh.junctions$dt

    if (nrow(pseudo.cnloh.junctions.dt)) {

        ## grab span - limit to below max.span
        pseudo.cnloh.junctions.dt[, span := pseudo.cnloh.junctions$span]

        ## compute node overlap with shadow
        node.overlap = pseudo.cnloh.junctions$shadow %N% gg$nodes$gr
        pseudo.cnloh.junctions.dt[, node.overlap.count := node.overlap]

        ## mark candidates
        ## pseudo.cnloh.edges = c(pseudo.cnloh.junctions.dt[span < max.span &
        ##                                                  class == "DEL-like" &
        ##                                                  node.overlap.count <= 3, edge.id],
        ##                        pseudo.cnloh.junctions.dt[span < max.span &
        ##                                                  class == "INV-like" &
        ##                                                  node.overlap.count <= 2, edge.id],
        ##                        pseudo.cnloh.junctions.dt[span < max.span &
        ##                                                  class == "DUP-like" &
        ##                                                  node.overlap.count <= 1, edge.id])

        #' zchoo Friday, Jul 30, 2021 11:31:11 AM
        ## changed this to include any junction with span under max.span
        pseudo.cnloh.edges = c(pseudo.cnloh.junctions.dt[span < max.span &
                                                         class == "DEL-like", edge.id],
                               pseudo.cnloh.junctions.dt[span < max.span &
                                                         class == "INV-like", edge.id],
                               pseudo.cnloh.junctions.dt[span < max.span &
                                                         class == "DUP-like", edge.id])
        
        ## pseudo.cnloh.edges = pseudo.cnloh[span < max.span & type == "ALT", edge.id]
        gg$edges[pseudo.cnloh.edges]$mark(cnloh = TRUE)

        if (verbose) {
            message("Number of pseudo-CNLOH edges marked:", length(pseudo.cnloh.edges))
        }

    } else {

        gg$edges$mark(cnloh = FALSE)

        if (verbose) {
            message("No pseudo-CNLOH edges detected.")
        }

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

        ## find index of the right and left phase blocks
        left.pblock = gr.match(gr.start(og.nodes), pblocks)
        right.pblock = gr.match(gr.end(og.nodes), pblocks)

        ## annotate GRanges with L and R phase block indices
        og.nodes$left.pblock = left.pblock
        og.nodes$right.pblock = right.pblock
    }

    if (haplotype.field %in% names(values(bins))) {
        if (verbose) {
            message("Adding haplotype information")
        }

        ## check that the entries in this field are valid
        if (!all(values(bins)[[haplotype.field]] %in% c("h1", "h2"))) {
            stop("values in haplotype.field must be either h1 or h2")
        }

        ## identify phase block of the left side of the node
        pblock.map = gr.match(pblocks, bins[values(bins)[[allele.field]] == "major"])
        pblock.major.haplotype = values(bins)[[haplotype.field]][pblock.map]
        node.left.major.haplotype = pblock.major.haplotype[og.nodes$left.pblock]
        node.right.major.haplotype = pblock.major.haplotype[og.nodes$right.pblock]

        major.nodes = copy(og.nodes)
        major.nodes$left.haplotype = node.left.major.haplotype
        major.nodes$right.haplotype = node.right.major.haplotype
        minor.nodes = copy(og.nodes)
        minor.nodes$left.haplotype = ifelse(node.left.major.haplotype == "h1",
                                            "h2",
                                            "h1")
        minor.nodes$right.haplotype = ifelse(node.right.major.haplotype == "h1",
                                             "h2",
                                             "h1")
        
        ## major.map = gr.match(og.nodes, bins[values(bins)[[allele.field]] == "major"])
        ## major.nodes = copy(og.nodes)
        ## major.nodes$haplotype = values(bins)[[haplotype.field]][major.map]
        ## minor.nodes = copy(og.nodes)
        ## minor.nodes$haplotype = ifelse(major.nodes$haplotype == "h1", "h2", "h1")

        major.nodes$allele = "major"
        minor.nodes$allele = "minor"
        phased.gg.nodes = c(major.nodes, minor.nodes)
    } else {
        phased.gg.nodes = c(og.nodes, og.nodes)
        phased.gg.nodes$allele = c(rep("major", length(og.nodes)),
                                   rep("minor", length(og.nodes)))
    }

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
                                    (phased.gg.nodes$nbins / (2 * sqrt(phased.gg.nodes$var))),
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

        ## add haplotype information of edge endpoints
        phased.gg.edges[, ":="(n1.haplotype = ifelse(n1.side == "left",
                                                     phased.gg.nodes$left.haplotype[n1],
                                                     phased.gg.nodes$right.haplotype[n1]),
                               n2.haplotype = ifelse(n2.side == "left",
                                                     phased.gg.nodes$left.haplotype[n2],
                                                     phased.gg.nodes$right.haplotype[n2]))]

        ## add phase block infromation of edge end points
        phased.gg.edges[, ":="(n1.pblock = ifelse(n1.side == "left",
                                                  phased.gg.nodes$left.pblock[n1],
                                                  phased.gg.nodes$right.pblock[n1]),
                               n2.pblock = ifelse(n2.side == "left",
                                                  phased.gg.nodes$left.pblock[n2],
                                                  phased.gg.nodes$right.pblock[n2]))]
        
        ## fix cross REF edges to zero within phase blocks
        phased.gg.edges[(n1.pblock == n2.pblock) & type == "REF" & (n1.haplotype != n2.haplotype),
                    ":="(ub = 0, lb = 0)]
        
        if (verbose) {
            message("Number of REF cross edges within phased blocks: ",
                    phased.gg.edges[(n1.pblock == n2.pblock) &
                                    type == "REF" &
                                    (n1.haplotype != n2.haplotype), .N])
        }
    }

    ## mark CNLOH edges
    if (!is.null(gg$edges$dt$cnloh)) {
        og.cnloh.edges = gg$edges$dt[cnloh == TRUE & type == "REF", edge.id]

        ## identify CNLOH cross edges and add extra edges that are ALT
        phased.cnloh.edges = phased.gg.edges[og.edge.id %in% og.cnloh.edges & (n1.allele != n2.allele),]
        phased.cnloh.edges[, ":="(cnloh = TRUE, type = "ALT", class = "CNLOH")]

        ## phased.gg.edges[og.edge.id %in% og.cnloh.edges, cnloh := TRUE]
        phased.gg.edges = rbind(phased.gg.edges, phased.cnloh.edges, fill = TRUE)

        if (verbose) {
            message("Number of CNLOH ALT edges added: ", phased.cnloh.edges[, .N])
        }

        ## mark pseudo-cnloh edges in child graph and make sure these are all zero
        og.pseudo.cnloh.edges = gg$edges$dt[cnloh == TRUE & type == "ALT", edge.id]
        phased.gg.edges[og.edge.id %in% og.pseudo.cnloh.edges & (n1.allele != n2.allele), cnloh := TRUE]

        if (verbose) {
            message("Number of pseudo-CNLOH ALT edges marked: ",
                    phased.gg.edges[og.edge.id %in% og.pseudo.cnloh.edges & (n1.allele != n2.allele), .N])
        }

    }

    ## identify phased edges (for linked reads)
    if (!is.null(edge.phase.dt)) {

        ## compute totals
        ephase = edge.phase.dt[, .(edge.id, n1.h1, n2.h1, n1.h2, n2.h2,
                                   n1.total = n1.h1 + n1.h2,
                                   n2.total = n2.h1 + n2.h2)][
                                       (n1.total > vbase.count.thres) | (n2.total > vbase.count.thres)]

        ## count fraction of reads corresponding with each allele
        ephase[, n1.h1.frac := n1.h1 / n1.total]
        ephase[, n2.h1.frac := n2.h1 / n2.total]
        ephase[, n1.h2.frac := n1.h2 / n1.total]
        ephase[, n2.h2.frac := n2.h2 / n1.total]

        ## set phase if passing proportion threshold (vbase.prop.thres)
        ephase[n1.h1.frac > vbase.prop.thres, n1.phase := "h1"]
        ephase[n1.h2.frac > vbase.prop.thres, n1.phase := "h2"]
        ephase[n2.h1.frac > vbase.prop.thres, n2.phase := "h1"]
        ephase[n2.h2.frac > vbase.prop.thres, n2.phase := "h2"]

        ## add phase information to edges data frame
        phased.gg.edges[, n1.phase := ephase$n1.phase[match(og.edge.id, ephase$edge.id)]]
        phased.gg.edges[, n2.phase := ephase$n2.phase[match(og.edge.id, ephase$edge.id)]]

        ## fix things to zero
        phased.gg.edges[n1.phase == "h1" & n1.haplotype == "h2", ":="(ub = 0, lb = 0)]
        phased.gg.edges[n2.phase == "h1" & n2.haplotype == "h2", ":="(ub = 0, lb = 0)]
        phased.gg.edges[n1.phase == "h2" & n1.haplotype == "h1", ":="(ub = 0, lb = 0)]
        phased.gg.edges[n2.phase == "h2" & n2.haplotype == "h1", ":="(ub = 0, lb = 0)]

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

    sol = Rcplex2(cvec = c,
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

        sol.new = Rcplex2(cvec = c, Amat = Ahat, bvec = bhat, 
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

#' @name parental
#' @title parental
#'
#' @description
#'
#' Converts an input unphased gGraph to a potential parental haplotype graph by randomly assigning ALT edges to a parental haplotype
#'
#' @param gg (gGraph) input gGraph. if desired can present haplotype field on edge metadata
#' @param fix (logical) fix marginal in balance? default TRUE
#' @param fix.emarginal (logical) fix edge marginal? default FALSE
#' @param force.alt (logical) force incorporation of all junctions? default TRUE
#' @param verbose (logical) default FALSE
#' @param lambda (numeric) default 10
#' @param eweight (numeric) edge weight default 1e3
#' @param epgap (numeric) default 1e-4
#' @param tilim (numeric) default 60
#' @param ... additional inputs to balance (e.g. epgap and whatnot)
#'
#' @return phased, balanced gGraph with og.node.id and allele annotation on nodes and og.edge.id annotation on edges
parental = function(gg,
                    haplotype.frac = 0.5,
                    fix = 1,
                    fix.emarginal = 0,
                    force.alt = TRUE,
                    verbose = FALSE,
                    lambda = 10,
                    eweight = 1e3,
                    epgap = 1e-4,
                    tilim = 60,
                    ...) {

    gg = gg$copy
    
    ## if (!("haplotype" %in% colnames(gg$edges$dt))) {

    ##     if (verbose) {
    ##         message("Assigning haplotypes with fraction ", haplotype.frac)
    ##     }

    ##     ## get number of ref and alt edges
    ##     n.alt = gg$edges$dt[type == "ALT", .N]
    ##     n.h1 = round(haplotype.frac * n.alt)
    ##     n.h2 = n.alt - n.h1
    ##     ht = sample(c(rep("h1", n.h1), rep("h2", n.h2)), size = n.alt, replace = FALSE)

    ##     ## mark haplotypes
    ##     gg$edges[type == "ALT"]$mark(haplotype = ht)
    ## } else {
    ##     if (verbose) {
    ##         message("using pre-assigned haplotypes")
    ##     }
    ## }

    n.og.nodes = nrow(gg$nodes$dt)
    new.nodes.dt = rbind(
        gg$nodes$dt[, .(og.node.id = node.id, haplotype = "h1", seqnames, start, end)],
        gg$nodes$dt[, .(og.node.id = node.id, haplotype = "h2", seqnames, start, end, cn = 1, weight = 1)],
        fill = TRUE
    )

    new.nodes.dt[, node.id := 1:.N]
    new.nodes.dt[, allele := "unphased"]

    new.edges.dt = rbind(
        gg$edges$dt[type == "REF", .(og.edge.id = edge.id, n1, n1.side, n2, n2.side, type)],
        gg$edges$dt[type == "REF", .(og.edge.id = edge.id,
                                     n1 = n1 + n.og.nodes, n1.side,
                                     n2 = n2 + n.og.nodes, n2.side, type)],
        gg$edges$dt[type == "ALT",
                    .(og.edge.id = edge.id,
                      n1, n1.side,
                      n2, n2.side, type)],
        gg$edges$dt[type == "ALT",
                    .(og.edge.id = edge.id,
                      n1 = n1 + n.og.nodes,
                      n1.side,
                      n2 = n2 + n.og.nodes,
                      n2.side,
                      type)],
        fill = TRUE
    )

    new.edges.dt[, connection := "straight"]

    haplotype.gg = gG(nodes = dt2gr(new.nodes.dt), edges = new.edges.dt)

    ## grab marginals...
    marginal.gr = gg$nodes$gr[, "cn"]
    marginal.gr$fix = fix

    ## grab edge marginals
    emarginals = this.complex$junctions[type == "ALT"]
    emarginals$set(fix = fix.emarginal)
    emarginals$set(weight = eweight)

    if (verbose) {
        message("Starting balance")
    }
    bal.gg = balance(haplotype.gg,
                     marginal = marginal.gr,
                     emarginal = emarginals,
                     phased = TRUE,
                     lp = TRUE,
                     tilim = tilim,
                     epgap = epgap,
                     verbose = verbose,
                     lambda = lambda,
                     ism = TRUE,
                     force.alt = force.alt)

    ## fix allele annotations
    if (verbose) {
        message("Formatting output graph and relabeling alleles")
    }
    bal.nodes.dt = bal.gg$nodes$dt
    bal.nodes.dt[, which.major := .SD$haplotype[which.max(.SD$cn)], by = og.node.id]
    bal.nodes.dt[, allele := ifelse(haplotype == which.major, "major", "minor")]
    bal.nodes.dt[, col := ifelse(allele == "major", alpha("red", 0.5), alpha("blue", 0.5))]

    bal.gg$nodes$mark(allele = bal.nodes.dt$allele, col = bal.nodes.dt$col)
    return(bal.gg)
}
    
        
#' @name simulate.hets
#' @title simulate.hets
#'
#' @description
#'
#' takes purity/ploidy 
#' 
#' @param gg (gGraph) phased balanced gGraph, such as from output of parental
#' @param bins (GRanges) locations of heterozygous sites with metadata columns allele and count
#' @param purity (numeric) default 1
#' @param ploidy (numeric) default 2
#' @param depth (numeric) (mean number of reads per site) default 50
#' @param theta (numeric) NB parameter, positive, infinite gives Poisson, default 1. recommend on the same order of magnitude as depth.
#'
#' @return GRanges with metadata fields count and allele representing simulated read counts given supplied graph and parameters
simulate.hets = function(gg,
                         bins,
                         purity = 1,
                         ploidy = 2,
                         depth = 50,
                         theta = 50) {

    if (!all(c("cn", "allele", "og.node.id") %in% colnames(gg$nodes$dt))) {
        stop("gg nodes missing metadata 'cn' and 'allele'")
    }
    if (!all(c("count", "allele") %in% names(values(bins)))) {
        stop("bins missing fields 'count' and 'allele'")
    }
    require(MASS)

    unique.bins = unique(bins[, c()])
    unique.bins$id = 1:length(unique.bins)

    ## create data.table with absolute dosage at all snp sites
    new.bins = c(unique.bins %$% gg$nodes[allele == "major"]$gr[, c("cn", "og.node.id", "allele")],
                 unique.bins %$% gg$nodes[allele == "minor"]$gr[, c("cn", "og.node.id", "allele")]) %>%
        as.data.table

    ## calculate slope and intercept from purity, ploidy, depth
    denom = 2 * (1 - purity) + purity * ploidy
    beta = depth * purity / denom
    gamma = depth * (1 - purity) / denom

    ## inverse rel2abs transformation
    new.bins[, mu := beta * cn + gamma]
    new.bins[, count := rnegbin(mu, theta = theta)]

    ## readjust so that major is always bigger than minor
    new.bins[, which.major := .SD$allele[which.max(.SD$count)], by = id]
    new.bins[, allele := ifelse(allele == which.major, "major", "minor")]
    
    return(dt2gr(new.bins[, .(seqnames, start, end, count, allele)]))
}
    

    
        
