\docType{class}
\name{gGraph}
\alias{gGraph}
\title{gGraph-class}
\description{
  The main work horse of this package. Rearrangement graph G=(V, E), where V is a set of
  strand-specific \code{GRanges} that both strand of any range must be present, and E is
  a set of directed edges connecting adjacent nodes stored in the form of \code{data.table}
  with two required columns \code{from} and \code{to} that matches the
  node's index inV.

  Every gGraph must be defined on a reference genome, and that is stored in the \code{seqinfo}
  of V. Optional metadata is allowed and appended as extra columns in V or E, some of which
  are required by the descendant classes like \code{bGraph}.

  Nodes and edges are necessary and sufficient to define a \code{gGraph}
  instance, while optional metadata fields like copy numbers, edge
  attributes can be extended. In the following examples \code{gg} is a gGraph object.
}

\section{Details}{
  \describe{
    \item{\strong{Constructors}}{
      \describe{
	\item{\strong{Usage}}{
	  \preformatted{
  gGraph$new(tile = NULL, junctions = NULL, cn = FALSE,
	     jabba = NULL,
	     weaver = NULL,
	     prego = NULL,
	     segs = NULL, es = NULL, ploidy = NULL, purity = NULL)
	  }
	}
	
	\item{\strong{Arguments}}{
	  \describe{
	    \item{tile}{the \code{GRanges} genome segmentation}
	    
	    \item{junctions}{the \code{GRangesList} of SV junctions}
	    
	    \item{cn}{\code{logical}, if \code{TRUE} will consider the copy
	      number of the input \code{tile} and \code{junction}. Otherwise,
	      create a \code{gGraph} without copy number specified.
	    }
	  }
	}
	
	\item{\strong{Value}}{
	  A \code{gGraph} object
	}
	
	\item{\strong{Details}}{
	  Direct use of constructors is only intended for development
	  purposes. We suggest using \code{\link{gread}}.

	  If either/both \code{tile} or \code{junctions} is provided,
	  the constrctor will use \code{\link{gg$karyograph}} method to
	  construct the output. \code{cn} will also be considered here.

	  If \code{jabba} provided as a \code{list} output object from
	  \href{http://github.com/mskilab/JaBbA}{JaBbA} or a path to
	  serialized such object,the constructor will create an output
	  based on the JaBbA output.
	}
      }
    }


    \item{\strong{Public Fields}}{
      In the following examples, \code{gg} is a \code{gGraph} instance.
      
      \describe{
	\item{gg$segstats}{The \code{GRanges} containing nodes}
	\item{gg$edges}{The \code{data.table} containing edges}
      }
    }


    \item{\strong{Public methods}}{
      \describe{
	\item{gg$addJuncs}{
	  Adding new junctions into the graph.
	  \describe{
	    \item{junc}{a \code{GRangesList} of input junctions}
	  }
	}
      }
    }  
  }
}

\examples{

}


\seealso{

}
