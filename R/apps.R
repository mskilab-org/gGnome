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
#' @param lambda positive number specifying default loose end penalty (100), note if gg$node metadata contain $lambda field then this lambda will be multiplied by the node level lambda
#' @param marginal GRanges with field $cn and optional $weight field will be used to fit the summed values at each base of the genome to optimally fit the marginal value, optional field $fix will actually constrain the marginal to be the provided value
#' @param tight indices or epxression on node metadata specifying at which nodes to disallow loose ensd
#' @param nfix indices or expression on node metadata specifying which node cn to fix
#' @param efix indices or expression on edge metadata specifying which edge cn to fix
#' @param nrelax indices or expression on node metadata specifying which nodes cn to relax
#' @param erelax  indices or expression on edge metadata specifying which edges cn to relax
#' @param L0  flag whether to apply loose end penalty as L1 (TRUE)
#' @param loose.collapse (parameter only relevant if L0 = TRUE) will count all unique (by coordinate) instances of loose ends in the graph as the loose end penalty, rather than each instance alone ... useful for fitting a metagenome graph   (FALSE)
#' @param phased (bool) indicates whether to run phased/unphased. default = FALSE
#' @param M  big M constraint for L0 norm loose end penalty, should be >1000
#' @param verbose integer scalar specifying whether to do verbose output, value 2 will spit out MIP (1)
#' @param tilim time limit on MIP in seconds (10)
#' @param epgap relative optimality gap threshhold between 0 and 1 (0.01)
#' @return balanced gGraph maximally resembling input gg in CN while minimizing loose end penalty lambda.
#' @author Marcin Imielinski
#' @export 
balance = function(gg,
                   lambda = 0.1,
                   marginal = NULL,
                   tight = NULL,
                   nfix = NULL, efix = NULL, nrelax = NULL, erelax = NULL,
                   L0 = TRUE,
                   loose.collapse = FALSE,
                   M = 1e2,
                   phased = FALSE,
                   verbose = 1,
                   tilim = 10,
                   epgap = 0.01)
{
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

  gg = gg$copy
  
  ## default local lambda lambda is node width
  if (!('lambda' %in% names(gg$nodes$dt)))
    gg$nodes$mark(lambda = 1)
#    gg$nodes$mark(lambda = width(gg$nodes$gr))

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

  if (phased) {
    if (verbose) {
      message("adding indicator variables for edge CN")
    }

    ## add og.edge.id information for each edge
    sedge.to.og.dt = gg$edges$dt[,
                                 .(sedge.id, og.edge.id, ## map sedge.id to og.edge.id
                                   ref.or.alt = type, ## get REF or ALT annotations (important for constraints)
                                   connection) ## get straight/cross annotations
                                 ]
    setkey(sedge.to.og.dt, "sedge.id")

    ## add binary indicator variables for each edge
    edge.indicator.vars = vars[type == "edge",][, type := "edge.indicator"][, vtype := "B"][, gid := sedge.id]
    setkey(edge.indicator.vars, "sedge.id")

    ## use sedge.id as a key to join with edge metadata
    edge.indicator.vars = edge.indicator.vars[sedge.to.og.dt]

    ## add one indicator sum variable per og edge ID to vars table
    vars = rbind(vars, edge.indicator.vars, fill = TRUE)

    if (verbose) {
      message("adding major/minor allele CN and og.node.id to vars")
    }

    ## idea here is to add variables that force major allele CN to be at least as large as minor CN
    ## create data table where keys are node ids, and og.node.id/allele can be easily found
    snode.to.og.dt = gg$nodes$dt[, .(snode.id, og.node.id, allele)]
    setkey(snode.to.og.dt, "snode.id")

    ## merge og.node.id and allele information into vars
    vars = merge(vars, snode.to.og.dt, by="snode.id", all.x = TRUE)
  }

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
  
  setkeyv(vars, c('type', 'gid'))

  ## add marginal copy number residual if specified
  vars$mfix = NA
  if (!is.null(marginal))
  {
    if (verbose) {
      message("adding marginal CN variables")
    }

    if (!inherits(marginal, 'GRanges') || is.null(marginal$cn))
    {
      stop('marginal must be a GRanges with field $cn')
    }
    
    if (is.null(marginal$weight))
      marginal$weight = 1

    if (is.null(marginal$fix))
      marginal$fix = FALSE
    
    ## first disjoin marginal against the nodes
    ## ie wee ned to create a separate residual variable for every unique
    ## disjoint overlap of marginal with the nodes
    dmarginal = gg$nodes$gr %>% gr.stripstrand %*% grbind(marginal %>% gr.stripstrand) %>% disjoin %$% marginal[, c('cn', 'weight', 'fix')] %Q% (!is.na(cn)) %Q% (!is.na(weight)) %Q% (!is.infinite(weight))

    vars = rbind(vars,
                 gr2dt(dmarginal)[, .(cn, weight, mfix = fix>0, rid = 1:.N, type = 'mresidual', vtype = 'C')],
                 fill = TRUE
                 )
    message("Done adding marginal vars")
  }

  vars[, id := 1:.N] ## set id in the optimization
  vars[is.na(lb), lb := -Inf]
  vars[is.na(ub), ub := Inf]
  vars[, relax := FALSE][, fix := FALSE]
  vars[type == 'mresidual' & mfix == TRUE, ":="(lb = 0, ub = 0)]
  vars[type %in% c('node', 'edge'), lb := pmax(lb, 0, na.rm = TRUE)]
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
  vars[type %in% c('loose.in', 'loose.in.indicator') & snode.id %in% term.in, terminal := TRUE]
  vars[type %in% c('loose.out', 'loose.out.indicator') & snode.id %in% term.out, terminal := TRUE]

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

  if (phased)
  {
    #'#########################
    ## add constraints forcing major CN to be larger than minor CN
    #'#########################
    if (!("allele" %in% colnames(vars)) | !("og.node.id" %in% colnames(vars))) {
      stop("allele field needs to be added to vars")
    }

    ## major allele coefficient is 1, minor allele coefficient is -1
    ## make sure that there's only one per node.id (abs of snode.id)
    allele.constraints = rbind(
      unique(
        vars[type == "node" & allele == "major",
             .(value = 1, id, node.id = abs(snode.id),
               cid = paste("allele.constraint", og.node.id))],
        by = "node.id"),
      unique(
        vars[type == "node" & allele == "minor",
             .(value = -1, id, node.id = abs(snode.id),
               cid = paste("allele.constraint", og.node.id))],
        by = "node.id"),
      fill = TRUE)

    ## add these constraints
    constraints = rbind(constraints,
                        allele.constraints[, .(value, id, cid)],
                        fill = TRUE)

    ## RHS: force (major CN - minor CN) to be >= 0
    allele.rhs = unique(
      vars[type == "node",
           .(value = 0, sense = "G", cid = paste("allele.constraint", og.node.id))],
      by = "cid")

    b = rbind(b, allele.rhs, fill = TRUE)

    #'#########################
    ## add constraints that force indicators to be 1 if edge CN > 0
    ## 
    #'#########################

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
    ## we may want to change this later to JaBbA CN instead of (potentially large) M for stability?
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
    iconstraints = vars[type == "edge", .(value = 1, id, sedge.id, cid = paste("edge.indicator.lb", sedge.id))]

    ## add matching indicator variables for LB
    iconstraints = rbind(
      iconstraints,
      vars[type == "edge.indicator", ][
        sedge.id %in% iconstraints$sedge.id, .(value = -0.1, id, cid = iconstraints$cid, sedge.id)],
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

    ###################
    ## add the edge indicator sum constraints
    ###################

    ## ALT edges: only one of four edges can have nonzero CN
    ## set upper bound (no need to set lower bound because these are binary variables and cannot be negative)
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

    ## REF edge configuration constraint
    ## iconstraints.from = unique(
    ##   vars[type == "edge.indicator" & ref.or.alt == "REF",
    ##        .(value = 1, id,
    ##          edge.id = abs(sedge.id),
    ##          cid = paste("ref.configuration.constraint.from", from))],
    ##   by = "edge.id"
    ## )

    ## iconstraints.to = unique(
    ##   vars[type == "edge.indicator" & ref.or.alt == "REF",
    ##        .(value = 1, id,
    ##          edge.id = abs(sedge.id),
    ##          cid = paste("ref.configuration.constraint.to", to))],
    ##   by = "edge.id"
    ## )

    ## iconstraints = rbind(iconstraints.from, iconstraints.to)
    ## constraints = rbind(
    ##   constraints,
    ##   iconstraints[, .(value, id, cid)],
    ##   fill = TRUE)

    ## ## sum to at most 1
    ## edge.indicator.b = unique(
    ##   iconstraints[, .(value = 1, sense = "L", cid)],
    ##   by = "cid"
    ## )

    ## ## add to b
    ## b = rbind(b, edge.indicator.b, fill = TRUE)


    ## ## REF edges: up to two of four edges can have nonzero CN (easiest to implement...)
    iconstraints = unique(
      vars[type == "edge.indicator" & ref.or.alt == "REF",
           .(value = 1, id,
             edge.id = abs(sedge.id),
             cid = paste("edge.indicator.sum.ub", og.edge.id))],
      by = "edge.id"
    )

    constraints = rbind(
      constraints,
      iconstraints[, .(value, id, cid)],
      fill = TRUE)

    edge.indicator.b = unique(
      vars[type == "edge.indicator" & ref.or.alt == "REF",
           .(value = 2, sense = "L", cid = paste("edge.indicator.sum.ub", og.edge.id))],
      by = "cid"
    )

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

  ########
  ## MAKE MATRICES
  ########

  ## now Rcplex time
  ## remove any rows with b = NA

  b = b[!is.na(value), ]
  constraints = constraints[cid %in% b$cid, ]

  ## convert constraints to integers
  ucid = unique(b$cid)
  b[, cid.char := cid]
  b[, cid := cid %>% factor(ucid) %>% as.integer]
  constraints[, cid.char := cid]
  constraints[, cid := cid %>% factor(ucid) %>% as.integer]
  setkey(b, cid)

  ## create constraint matrix, Qmat, and cobj, lb, ub from vars and constraints  lambda = 10
  Amat = sparseMatrix(constraints$cid, constraints$id, x = constraints$value, dims = c(length(ucid), nrow(vars)))
  vars[is.na(weight), weight := 0]

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
          cvec = lambda*(vars[, lambda*(type %in% c('loose.in.indicator.sum.indicator', 'loose.out.indicator.sum.indicator', 'loose.in.indicator', 'loose.out.indicator') & !terminal)] %>% as.numeric)
        }
      else
        {
          cvec = lambda*(vars[, lambda*(type %in% c('loose.in.indicator', 'loose.out.indicator') & !terminal)] %>% as.numeric)
        }
    }
  else
    cvec = lambda*(vars[, lambda*(type %in% c('loose.in', 'loose.out') & !terminal)] %>% as.numeric)

  ## implement reward if provided
  if (length(ix <- which(vars$reward!=0)))
  {
    if (verbose)
      message('Applying reward')
    cvec[ix] = -vars$reward[ix]
  }

  lb = vars$lb
  ub = vars$ub
  bvec = b[.(1:nrow(Amat)), value]
  sense = b[.(1:nrow(Amat)), sense]

  control = list(trace = ifelse(verbose>=2, 1, 0), tilim = tilim, epgap = epgap, round = 1)
  sol = Rcplex::Rcplex(cvec = cvec, Amat = Amat, bvec = bvec, Qmat = Qmat, lb = lb, ub = ub, sense = sense, vtype = vars$vtype, objsense = 'min', control = control)
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

  return(gg)
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
#' @return gGraph whose nodes are annotated with $cn and $weight field
#' @export
#' @author Marcin Imielinski
binstats = function(gg, bins, by = NULL, field = NULL, purity = gg$meta$purity, ploidy = gg$meta$ploidy, loess = TRUE, min.bins = 3, verbose = TRUE, min.var = 0.1)
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
  ov = gr.indoverlaps(gg$nodes$gr, bins, by = by, scol = names(values(bins)), return.type = 'data.table')
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
  dt$weight = dt$nbins/(2*dt$var)


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
#' Identifies regions without allelic CN imbalance
#'
#' @param gg junction-balanced phased gGraph. each node must have associated og node.id
#' @param mc.cores (int) number of cores. default = 8.
#' @param verbose (bool) verbose = TRUE prints stuff. default TRUE.
phased.postprocess = function(gg, mc.cores = 8, verbose = TRUE)
{
  ## check that gg nodes and edges have og node id
  if (!("og.node.id" %in% colnames(gg$nodes$dt)) | !("og.edge.id" %in% colnames(gg$edges$dt))) {
    stop("gGraph must have og.node.id and og.edge.id node/edge metadata columns")
  }

  ## check that graph has been balanced (need cn.old and cn)
  if (!("cn" %in% colnames(gg$nodes$dt)) | !("cn.old" %in% colnames(gg$nodes$dt))) {
    stop("run balance to populate nodes with cn and cn.old")
  }


  if (verbose) {
    message("Identifying node pairs without CN imbalance")
  }

  ## identify nodes without allelic CN imbalance
  og.node.balance = gg$nodes$dt[, .(og.node.id, allele, cn)] %>%
    dcast(og.node.id ~ allele, value.var = "cn") ## data.table keyed by og.node.id with major/minor CN

  ## identify og.node.id without cn.imbalance (add a column to indicate if there is CN imbalance)
  og.node.balance[, ":="(cn.imbalance = (major != minor),
                         total.cn = major + minor)]
  unphased.og.nodes = og.node.balance[cn.imbalance == FALSE, og.node.id]

  ## identify minor nodes that will be deleted (map og.node.id to major and minor nodes)
  deleted.minor.nodes = gg$nodes$dt[og.node.id %in% unphased.og.nodes, .(og.node.id, allele, node.id)] %>%
    dcast(og.node.id ~ allele, value.var = "node.id") ## data.table keyed by og.node.id in unphased regions

  ## modify major allele of unphased nodes so that CN is now total CN
  new.major.dt = gg$nodes$dt[allele == "major" & og.node.id %in% unphased.og.nodes,]
  new.major.dt$cn = NULL
  new.major.dt = new.major.dt %>%
    merge(og.node.balance[, .(og.node.id, cn = total.cn)], by = "og.node.id")
  ## mark these nodes
  new.major.dt[, ":="(col = "purple",
                      allele = "unphased")]

  ## prepare GRanges for new nodes
  new.nodes.dt = rbind(new.major.dt,
                       gg$nodes$dt[!(og.node.id %in% unphased.og.nodes),])
  

  if (verbose) {
    message("Identifying NA valued nodes")
  }


  ## tbh this could easily be a preprocessing step instead of post-processing. consider adding to binstats.
  ## need to check cn.old for NA value and map og.node.id to major/minor node.ids
  ## na.node.dt = gg$nodes$dt[is.na(cn.old), .(og.node.id, allele, node.id)] %>%
  ##   dcast(og.node.id ~ allele, value.var = "node.id")

  ## ## just mark for now! don't remove/merge these.
  ## new.nodes.dt[node.id %in% 

  ## ## mark major nodes
  ## ## new.nodes.dt[node.id %in% na.node.dt$major, ":="(allele = "na.node", col = "black")]

  ## ## remove minor nodes
  ## new.nodes.dt = new.nodes.dt[!(node.id %in% na.node.dt$minor),]

  if (verbose) {
    message("Processing edges and reindexing")
  }

  ## lol just delete edges for now and rebalance
  new.edges.dt = gg$edges$dt[(n1 %in% new.nodes.dt$node.id) & (n2 %in% new.nodes.dt$node.id),]

  ## in balance, we should allow the number of alt edges associated with these nodes to be up to two
  ## otherwise there may be feasibility issues with the marginal constraint :(
  ## but they would need to be both straight or both cross
  
  ## reindex nodes and edges
  new.nodes.gr = dt2gr(new.nodes.dt)
  reindex = 1:length(new.nodes.gr)
  names(reindex) = as.character(new.nodes.gr$node.id)
  new.edges.dt[, ":="(n1 = reindex[as.character(n1)],
                      n2 = reindex[as.character(n2)])]

  return(gG(nodes = new.nodes.gr, edges = new.edges.dt[, .(n1, n2, n1.side, n2.side,
                                                           cn = cn.old, og.edge.id, connection)]))
}


#' @name phased.binstats
#' @title phased.binstats
#' @description
#'
#' Given GRanges containing major/minor allele counts and a balanced but unphased gGraph,
#' prepares phased gGraph input to balance.
#' 
#' @param gg gGraph
#' @param bins GRanges with: (($allele | type type) & ($cn | field field))
#' @param fix.dels (bool) remove obviously deleted alleles. default TRUE
#' @param fix.hets (bool) remove phasing if allele CNs are equal? default TRUE
#' @param field (str) field for allele read counts
#' @param type (str) field for specifying whether the listed CN corresponds with major or minor allele
#' @param purity purity parameter either specified together with field or embedded in gg$meta, must be specified if field is not NULL
#' @param ploidy ploidy parameter either specified together with field or embedded in gg$meta, must be specified if field is not NULL
#' @param min.bins minimum number of bins to use for intra segment variance computation (3)
#' @param verbose (bool) default TRUE for debugging
#' @param min.var minimal allowable per segment bin variance, which will ignore segments with very low variance due to all 0 or other reasons (0.1)
#' @param mc.cores (int) number of cores
#' @return gGraph whose nodes are annotated with $cn, $allele, and $weight field
phased.binstats = function(gg,
                           bins,
                           fix.del = FALSE,
                           fix.het = FALSE,
                           field = "cn",
                           type = "allele",
                           purity = gg$meta$purity,
                           ploidy = gg$meta$ploidy,
                           min.bins = 3,
                           verbose = TRUE,
                           min.var = 0.1,
                           mc.cores = 8)
{
  reads.to.allelic.cn = function(bins, field, purity = 1.0, ploidy = 2.0) {
    #' params
    #' bins (GRanges): should contain metadata column specified by field containing allele-specific reads
    #' field (string): metadata column name. values in this column are allele-specific reads
    #' purity (float)
    #' ploidy (float)
    #' returns
    #' numeric vector corresponding to the allele-specific CN on each position in bins
    #' check that every entry is doubled in bins?

    ## check that expected metadata columns are present
    if (!(field %in% names(values(bins)))) {
      stop("field specified by field is not a metadata column in bins")
    }

    ## mean allele-specific read count across all heterozygous SNPs
    y = values(bins)[[field]]
    y.bar = mean(y, na.rm = TRUE) * 2 ## careful! assumes a perfect GRanges with exactly duplicated rows. might want to check for this in the future!

    ## purity and ploidy for notational consistency
    alpha = purity
    tau = ploidy

    ## 2x inter-peak space
    denom = alpha * tau + 2.0 * (1 - alpha) ## non-allelic CN given purity
    beta = (y.bar * alpha) / denom ## 1/2 gap between non-alleleic peaks
    gamma = y.bar * (1 - alpha) / denom ## first jump from zero for non-allelic CN. differs from paper, missing factor of two
    if (verbose) {
      message("Y-bar (mean allelic reads): ", y.bar)
      message("Beta (slope): ", beta)
      message("Gamma (intercept): ", gamma)
    }

    ## allele-specific CN
    cn = (y - gamma) / beta
    return(cn)
  }

  #' prepare skeleton for phased gGraph (to be populated with CN estimates)
  if (verbose == TRUE) {
    message("Preparing phased gGraph...")
  }

  #' create GRanges corresponding to nodes of major and minor allele graphs
  n.nodes = length(gg$nodes) ## get number of nodes in the original unphased graph
  major.nodes.gr = gg$nodes$gr[,c("node.id")]
  minor.nodes.gr = gg$nodes$gr[,c("node.id")]
  #' assign unique node.id and store original node id as og.node.id
  names(values(major.nodes.gr)) = c("og.node.id") ## store original node ID
  names(values(minor.nodes.gr)) = c("og.node.id")
  major.nodes.gr$node.id = major.nodes.gr$og.node.id ## keep node.id of major allele
  minor.nodes.gr$node.id = minor.nodes.gr$og.node.id + n.nodes ## shift node.id of minor allele
  #' label whether node belongs to major or minor allele
  major.nodes.gr$allele = "major"
  minor.nodes.gr$allele = "minor"

  #' create data.tables corresponding to  edges that go straight across
  major.edges.dt = gg$edges$dt[, .(og.edge.id = edge.id,
                                   n1.side, n2.side,
                                   n1, n2,
                                   connection = "straight")] ## indicate connection type
  minor.edges.dt = gg$edges$dt[, .(og.edge.id = edge.id,
                                   n1.side, n2.side,
                                   n1 = n1 + n.nodes, ## convert n1, n2 node.id to minor allele counterparts
                                   n2 = n2 + n.nodes,
                                   connection = "straight")] ## indicate connection type
  #' get data.table for edges that cross from major to minor allele
  cross.edges.dt = mclapply(1:nrow(major.edges.dt),
                            function(ix) {
                              row = major.edges.dt[ix,]
                              n1.side = row$n1.side
                              n2.side = row$n2.side
                              new.n1 = row$n1
                              new.n2 = row$n2 + n.nodes ## covert n2 to minor allele node index
                              new.row = data.table(
                                og.edge.id = row$og.edge.id,
                                n1.side = n1.side,
                                n2.side = n2.side,
                                n1 = c(row$n1, row$n1 + n.nodes),
                                n2 = c(row$n2 + n.nodes, row$n2),
                                connection = "cross" ## indicate connection type
                              )
                              return(new.row)
                            },
                            mc.cores = mc.cores) %>% rbindlist()

  #' create new gGraph
  phased.nodes = c(major.nodes.gr, minor.nodes.gr)
  phased.edges = list(major.edges.dt, minor.edges.dt, cross.edges.dt) %>% rbindlist()
  phased.gg = gG(nodes = phased.nodes, edges = phased.edges)

  #' update edge colors for plotting
  phased.gg$edges[connection == "cross" & type == "REF"]$mark(col = "light blue")
  phased.gg$edges[connection == "cross" & type == "ALT"]$mark(col = "pink")
  phased.gg$edges[connection == "straight" & type == "REF"]$mark(col = "blue")
  phased.gg$edges[connection == "straight" & type == "ALT"]$mark(col = "red")

  #' update node colors for plotting
  phased.gg$nodes[allele == "major"]$mark(col = "red")
  phased.gg$nodes[allele == "minor"]$mark(col = "blue")

  if (!is.null(field) & !is.null(purity) & !is.null(ploidy) && is.numeric(purity) && is.numeric(ploidy))
  {
   if (verbose)
      message('Converting ', field, ' to cn using purity ', purity, ' and ploidy ', ploidy)
    bins$cn = reads.to.allelic.cn(bins, field, purity = purity, ploidy = ploidy)
  }

  #' check that bins has required fields for CN
  if (is.null(bins$cn)) {
    stop("bins must have metadata column $cn if field is not specified")
  }

    #' ensure that allele column is present or type column is present
  if (is.null(bins$allele)) {
    stop("bins must have metadata column $allele")
  }

  #' make sure that major and minor are the only entries in allele
  allele.entries = unique(bins$allele) %>% .[order(.)]
  if (all(allele.entries != c("major", "minor"))) {
    stop("metadata column $allele must contain only values 'major' and 'minor'")
  }

   #' overlap major/minor alleles with bins separately
  if (verbose) {
    message("crossing nodes and bins via gr.findoverlaps")
  }

  ov = gr.findoverlaps(phased.nodes, ## concatenated GRanges for phased gGraph
                       bins, 
                       by = c("allele"), ## only find overlaps if alleles field is matching
                       qcol = c("node.id"), ## columns to keep from the nodes GRanges
                       scol = c("allele", "cn"), ## columns to keep from bins
                       return.type = "data.table")

  #' compute bin stats per node
  if (verbose) {
    message("aggregating bin stats per node")
  }
  #' get node information
  dt = ov[, .(mean = mean(cn, na.rm = TRUE),
              var = var(cn, na.rm = TRUE),
              nbins = .N),
          by = node.id]

  message(colnames(dt))

  allele.info = data.table(
    node.id = phased.nodes$node.id,
    allele = phased.nodes$allele
  )

  message(colnames(allele.info))

  dt = merge(dt, allele.info, by = "node.id", all.y = TRUE) %>% .[order(node.id),]

  message(colnames(dt))

  #' identifying deleted segments on minor allele
  if (fix.del) {
    if (verbose) {
      message("identifying deleted alleles")
    }
    ## mark nodes on minor allele with zero CN and low variance
    mean.thres = 0.1
    var.thres = 0.1
    dt[allele == "minor" & mean < mean.thres & var < var.thres,
       ":="(phasing = "del")]
    ## mark nodes connected by straight ref edges to deleted nodes
    ## idea here is to remove any short nodes/NA mean nodes for cleaner deletions
    del.nodes = dt[phasing == "del", node.id]
    adj.nodes = union(phased.gg$edges$dt[n1 %in% del.nodes & connection == "straight" & type == "REF", n2],
                      phased.gg$edges$dt[n2 %in% del.nodes & connection == "straight" & type == "REF", n1])
    dt[node.id %in% adj.nodes & (is.na(mean) | mean < mean.thres),
       ":="(phasing = "del")]
  }

  #' identifying segments where minor and major allele CN are equivalent
  ## maximally stupid/easy thing where just means are compared
  if (fix.het) {
    if (verbose) {
      message("identifying het regions")
    }
    ## identify het og.nodes
    tmp = phased.gg$nodes$dt[, .(og.node.id, node.id, allele)] %>%
      dcast(og.node.id ~ allele, value.var = "node.id") %>%
      as.data.table()
    tmp[match(dt[allele == "minor", node.id], minor), ":="(minor.cn = dt[allele == "minor", mean])]
    tmp[match(dt[allele == "major", node.id], major), ":="(major.cn = dt[allele == "major", mean])]
    tmp[, ":="(diff = major.cn - minor.cn)]
    ## identify het nodes (using thres)
    ## match back diff to dt
    dt[match(tmp$minor, node.id), ":="(diff = tmp$diff)]
    dt[match(tmp$major, node.id), ":="(diff = tmp$diff)]
    diff.thres = 0.5
    dt[diff < 0.5, ":="(phasing = "het")]
    ## mark adjacent nodes
    het.nodes = dt[phasing == "het", node.id]
    adj.nodes = union(phased.gg$edges$dt[n1 %in% het.nodes & connection == "straight" & type == "REF", n2],
                      phased.gg$edges$dt[n2 %in% het.nodes & connection == "straight" & type == "REF", n1])
    dt[node.id %in% adj.nodes & is.na(diff),
       ":="(phasing = "het")]
    if (verbose) {
      message("found ", length(dt[phasing=="het", node.id]), " het regions")
    }
  }

  #' set variance to NA if number of bins is less than specificied minimum
  dt[nbins < min.bins, var := NA]

  #' compute weights (nbins / variance)
  ## for now adding a jitter
  dt[, ":="(weight = nbins / (2 * var + 1e-2))]

  #' add cn (dt$mean) and weight to phased gGraph
  phased.gg$nodes$mark(cn = dt$mean, weight = dt$weight, fixed = dt$phasing)

  if (any(is.infinite(dt$weight), na.rm = TRUE)) {
    warning('variance computation yielded infinite weight, consider setting min.bins higher or using loess fit')
  }

  #' remove deleted nodes
  if (fix.del) {
    phased.gg$nodes[fixed == "del"]$mark(col = "black") ## mark for now with different color
  }

  if (fix.het) {
    phased.gg$nodes[fixed == "het"]$mark(col = "purple") ## mark for now with different color
  }
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
            lb = do.call("c", c(lapply(lb.cols, function(x) gw$dt[[x]]), rep(0, len(gw))))
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
            ub = do.call("c", c(lapply(ub.cols, function(x) gw$dt[[x]]), rep(Inf, len(gw))))
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
            this.cnf = rep(c(cn.field, "indicator"), each = len(gw))
            this.ls = split(this.x, this.cnf)
            names(this.ls) = paste(names(this.ls), i, sep = sep)
            do.call(gw$set, this.ls)
        }
        return(invisible(gw))
    } else {
        return(sol)
    }
}
