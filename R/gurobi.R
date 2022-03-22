#' @name run_gurobi
#' @title run_gurobi
#'
#' @description
#' Rcplex2-like wrapper for calling Gurobi for MIP optimization
#' 
#' @name cvec (numeric) linear part of objective function
#' @name Amat (Matrix) sparse matrix of constraints
#' @name bvec (numeric) RHS of constraints
#' @name Qmat (Matrix) sparse matrix, quadratic
#' @name lb (numeric) lower bounds
#' @name ub (numeric) upper bounds
#' @name sense (character) directionality of constraints (e.g. E, G, L)
#' @name objsense (character) optimization type (one of min, max)
#' @name control (list) (MIP control parameters)
#' @name threads (numeric) number of threads for parallelization
run_gurobi = function(cvec = NULL,
                      Amat = NULL,
                      bvec = NULL,
                      Qmat = NULL,
                      lb = NULL,
                      ub = NULL,
                      sense = NULL,
                      vtype = NULL,
                      objsense = 'min',
                      control = list(epgap = 1e-2, tilim = 360, trace = 2),
                      threads = 32)
{

    ## check that gurobi is installed
    if (!requireNamespace("gurobi", quietly = TRUE)) {
        stop("The gurobi package must be installed to use this functionality")
    }
    ## build model
    model = list(
        obj = cvec,
        A = Amat,
        rhs = bvec,
        Q = Qmat,
        lb = lb,
        ub = ub,
        vtype = vtype,
        sense = c("E"="=", "G"=">", "L"="<")[sense], ## inequalities are leq, geq (e.g. not strict)
        modelsense = objsense)

    ## params
    params = list()
    if (!is.null(control$epgap)) {
        params$MIPGap = control$epgap
    }
    if (!is.null(control$tilim)) {
        params$TimeLimit = control$tilim
    }
    if (!is.null(control$trace)) {
        params$LogToConsole = ifelse(control$trace > 0, 1, 0)
    }

    ## TODO: set up env list for running on compute cluster

    ## run gurobi
    sol = gurobi::gurobi(model = model, params = params)

    ## make solution consistent with Rcplex output
    ## but return all the things
    sol$epgap = sol$mipgap

    return(sol)
}
