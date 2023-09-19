#' @name testCPLEX
#' @description
#' Tests if CPLEX functionality has been added to gGnome
#'
#' 
#' @return returns a message whether the check was successful, else does a fresh installation of gGnome 
#' @author Tanubrata Dey
#'  Allows for subsetting of the gEdge Object using bracket notation


testOptimizationFunction <- function() {
                                        # Define your input data and parameters for testing
                                        # Replace these with your actual input data and parameters
  cvec <- c(1, 2, 3)
  Amat <- matrix(c(1, 2, 3, 4, 5, 6), ncol = 3)
  bvec <- c(10, 20)
  
                                        # Define a default matrix for Qmat (replace this with your default matrix)
  defaultQmat <- matrix(0, ncol = 3, nrow = 3)
  
  tryCatch({
                                        # Call the optimization function with the default matrix
    gGnome::Rcplex2(cvec, Amat, bvec, Qmat = defaultQmat)
    message("Optimization successful.\n")
                                        # You can add additional checks on the result here if needed
  }, error = function(e) {
                                        # If there is an error, install the gGnome package
    message("CPLEX is not wired with gGnome, Make sure CPLEXDIR exists in your .bashrc to add CPLEX functionality in gGnome.")
    cplex.dir = Sys.getenv("CPLEX_DIR")
    if (is.null(cplex.dir)){
      stop("CPLEX_DIR environment variable not found!")
    } else if (!file.exists(paste0(cplex.dir, "/cplex"))) {
      stop("${CPLEX_DIR}/cplex not found")
    } else if (!file.exists(paste0(cplex.dir, "/cplex/include")) ||
               !file.exists(paste0(cplex.dir, "/cplex/lib"))){
      stop("${CPLEX_DIR}/cplex/[(include)|(lib)] do not both exist")
    } else {
      message("Force re-installing gGnome to add CPLEX...")
    # devtools::install_github("mskilab-org/gGnome", force = TRUE)
      install.packages("~/git/gGnome_0.1.tar.gz", force = TRUE)
    }
  })
}







