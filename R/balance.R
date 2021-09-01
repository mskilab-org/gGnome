#' @name balance2
#' @title balance2
#'
#' @description
#'
#' beta reimplementation of balance
#'
#' @param gg (gGraph) with node fields
#' - $cn (numeric) CN guess
#' - $weight (numeric)
#' @param lambda
#' @param marginal
#' @param phased (logical) default FALSE
#' @param ism (logical) default TRUE
#' @param M (numeric) default 1000
#' @param epgap (numeric) default 0.1
#' @param trelim (numeric) max size of uncompressed tree in GB
#' @param nodefileind (numeric) one of 0 (no node file) 1 (in memory compressed) 2 (on disk uncompressed) 3 (on disk compressed) default 1
#' @return junction-balanced gGraph
balance2 = function(gg = NULL,
                    lambda = 100,
                    marginal = NULL,
                    phased = FALSE,
                    ism = TRUE,
                    M = 1000,
                    epgap = 0.1,
                    trelim = 16,
                    nodefileind = 1) {
    
}

#' @name check_balance_inputs
#' @title check_balance_inputs
#'
#' @description
#'
#' Validate input gGraph for balance
#'
#' @param gg (gGraph)
check_balance_inputs = function(gg) {
}

#' @name create_node_variables
#' @title create_node_variables
#'
#' @description
#' create node CN, node residual, and helper variables for nodes
#' 
#' @param gg
#' @return data.table
create_node_variables = function(gg) {
}

#' @name create_edge_variables
#' @title create_edge_variables
#'
#' @description
#' edge CN, edge residual, and helper variables for edges
#'
#' @param gg
#' @return data.table
create_edge_variables = function(gg) {
}
