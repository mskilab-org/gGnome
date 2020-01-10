## appease R CMD check vs data.table 
sid=side1=side2=side_1=side_2=silent=snid=splice.variant=splicevar=str1=str2=strand_1=strand_2=subject.id=suffix=tag=threep.cc=threep.coord=threep.exon=threep.frame=threep.pc=threep.sc=threep.sc.frame=to=transcript.id.x=transcript.id.y=transcript_associated=transcript_id=twidth=tx.cc=tx.ec=tx_strand=tx_strand.x=tx_strand.y=txid=type=uid=uids=val=walk.id=walk.iid=walk.iid.x=walk.iid.y=wkid=NULL

#' @name convex.basis
#' @description
#'
#' Outputs a matrix K of the convex basis of matrix A
#'
#' i.e. each column x = K[,i] is a minimal solution (with respect to sparsity) to
#' Ax = 0, x>=0
#'
#' exclude.basis =  0, 1 matrix of dimension k x ncol(A) specifying k sparsity patterns that we would
#' like to exclude from the convex.basis.  This can speed up computation since any non-negative
#' combinations of vectors that satisfy an exclusion property will also be excludable, and thus
#' we can remove such vectors as soon as we detect them..
#'
#' exclude.range = 9, 1 matrix of dimension k x nrow(A) specifying k sparsity patterns that we would like
#' exclude, but these are specified in terms of the range of abs(A) .. i.e. we want to exclude all
#' basis vectors v such that nz(exclude.ranges[i, ]) C  nz(abs(A)*v)) for some pattern i.  Again
#' any non-neg linear comb of any intermediate-basis vector that satisfies this property will satisfy it,
#' as a result we can exclude these vectors when we see them.
#'
#' @param A nxm signed incidence matrix of directed graph of n nodes and m edges
#' @param interval chunks with which to do all vs all basis comparison
#' @param chunksize  chunksize with which to proceed through higher level iteration
#' @param verbose logical flag whether to output
#' @return m x k matrix of elementary paths and cycles comprising the non-negative convex basis of the' non-negative column span of A
#' @author Marcin Imielinski
convex.basis = function(A,
                        interval = 80,
                        chunksize = 100,
                        exclude.basis = NULL,
                        exclude.range = NULL,
                        maxchunks = Inf,
                        verbose = F){
    ZERO = 1e-8;
    remaining = 1:nrow(A);
    iter = 0;
    i = 0;
                                        #    order = c()
    numelmos = c()
    K_i = I = as(diag(rep(1, ncol(A))), 'sparseMatrix');
                                        #    A_i = as(A %*% K_i, 'sparseMatrix');
    K_i = I = diag(rep(1, ncol(A)))
    A_i = A %*% K_i

    if (!is.null(exclude.basis)){
        exclude.basis = sign(exclude.basis)
        exclude.basis = exclude.basis[rowSums(exclude.basis)>0, ]
        if (nrow(exclude.basis) == 0){
            exclude.basis = NULL
        }
    }

    if (!is.null(exclude.range)){
        exclude.range = sign(exclude.range)
        exclude.range = exclude.range[rowSums(exclude.range)>0, ]
        if (nrow(exclude.range) == 0){
            exclude.range = NULL
        }
    }

                                        # vector to help rescale matrix (avoid numerical issues)
    mp  = apply(abs(A), 1, min); # minimum value of each column
    mp[mp[ZERO]] = 1; # columns with zero minimum get scale "1"

    st = Sys.time()
                                        # iterate through rows of A, "canceling" them out
    while (length(remaining)>0){   
        ## TODO figure out why we have to check this so many times
        if (nrow(K_i)==0 | ncol(K_i)==0){
            return(matrix())
        }

        iter = iter+1;
        K_last = K_i;

        if (verbose){
            print(Sys.time() - st)
        }

        if (verbose){
            cat('Iter ', iter, '(of',  nrow(A_i),  ') Num basis vectors: ', nrow(K_i), " Num active components: ", sum(Matrix::rowSums(K_i!=0)), "\n")
        }

        i = remaining[which.min(Matrix::rowSums(A_i[remaining,, drop = FALSE]>=ZERO)*Matrix::rowSums(A_i[remaining,, drop = FALSE]<=(-ZERO)))]  # chose "cheapest" rows

        remaining = setdiff(remaining, i);
                                        #        order = c(order, i);

        zero_elements = which(abs(A_i[i, ]) <= ZERO);
        K_i1 = K_last[zero_elements, , drop = FALSE];  ## K_i1 = rows of K_last that are already orthogonal to row i of A
        K_i2 = NULL; ## K_i1 = will store positive combs of K_last rows that are orthogonal to row i of A (will compute these below)

        pos_elements = which(A_i[i, ]>ZERO)
        neg_elements = which(A_i[i, ]<(-ZERO))

        if (verbose){
            cat('Iter ', iter, " Row ", i, ":", length(zero_elements), " zero elements ", length(pos_elements), " pos elements ", length(neg_elements), " neg elements \n")
        }

        if (length(pos_elements)>0 & length(neg_elements)>0)
            for (m in seq(1, length(pos_elements), interval))
                for (l in seq(1, length(neg_elements), interval)){
                    ind_pos = c(m:min(c(m+interval, length(pos_elements))))
                    ind_neg = c(l:min(c(l+interval, length(neg_elements))))

                    indpairs = cbind(rep(pos_elements[ind_pos], length(ind_neg)),
                                     rep(neg_elements[ind_neg], each = length(ind_pos))); # cartesian product of ind_pos and ind_neg
                    pix = rep(1:nrow(indpairs), 2)
                    ix = c(indpairs[,1], indpairs[,2])
                                        #                coeff = c(-A_i[i, indpairs[,2]], A_i[i, indpairs[,1]])  ## dealing with Matrix ghost
                    coeff = c(-A_i[i, ][indpairs[,2]], A_i[i, ][indpairs[,1]])  ##
                    combs = Matrix::sparseMatrix(pix, ix, x = coeff, dims = c(nrow(indpairs), nrow(K_last)))
                    combs[cbind(pix, ix)] = coeff;

                    H = combs %*% K_last;

                                        # remove duplicated rows in H (with respect to sparsity)
                    H = H[!duplicated(as.matrix(H)>ZERO), ];

                                        # remove rows in H that have subsets in H (with respect to sparsity) ..
                    if ((as.numeric(nrow(H))*as.numeric(nrow(H)))>maxchunks){
                        print('Exceeding maximum number of chunks in convex.basis computation')
                        stop('Exceeding maximum number of chunks in convex.basis computation')
                    }
                    keep = which(Matrix::colSums(sparse_subset(abs(H)>ZERO, abs(H)>ZERO, chunksize = chunksize, quiet = !verbose))<=1) # <=1 since every H is its own subset
                    H = H[keep, , drop = FALSE]

                                        # remove rows in H that have subsets in K_i2
                    if (!is.null(K_i2))
                        if (nrow(K_i2)>0){
                            if ((as.numeric(nrow(K_i2))*as.numeric(nrow(H)))>maxchunks){
                                print('Exceeding maximum number of chunks in convex.basis computation')
                                stop('Exceeding maximum number of chunks in convex.basis computation')
                            }
                            keep = which(Matrix::colSums(sparse_subset(abs(K_i2)>ZERO, abs(H)>ZERO, chunksize = chunksize, quiet = !verbose))==0)
                            H = H[keep, , drop = FALSE]
                        }

                                        # remove rows in H that have subsets in K_i1
                    if (!is.null(K_i1))
                        if (nrow(K_i1)>0){
                            if ((as.numeric(nrow(K_i1))*as.numeric(nrow(H)))>maxchunks)
                            {
                                print('Exceeding maximum number of chunks in convex.basis computation')
                                stop('Exceeding maximum number of chunks in convex.basis computation')
                            }
                            keep = which(Matrix::colSums(sparse_subset(abs(K_i1)>ZERO, abs(H)>ZERO, chunksize = chunksize, quiet = !verbose))==0)
                            H = H[keep, , drop = FALSE]
                        }

                                        # maintain numerical stability
                    if ((iter %% 10)==0){
                        H = diag(1/apply(abs(H), 1, max)) %*% H

                                        #                K_i2 = rBind(K_i2, H)
                    }
                    K_i2 = rbind(K_i2, as.matrix(H))
                }

                                        #        K_i = rBind(K_i1, K_i2)
        K_i = rbind(K_i1, K_i2) ## new basis set

        if (nrow(K_i)==0){
            return(matrix())
        }
        
        ## only keep vectors that fail to intersect all vectors "exclude" in matrix
        if (!is.null(exclude.basis)) {
            if ((as.numeric(nrow(exclude.basis))*as.numeric(nrow(K_i)))>maxchunks){
                print('Exceeding maximum number of chunks in convex.basis computation')
                stop('Exceeding maximum number of chunks in convex.basis computation')
            }
            keep = Matrix::colSums(sparse_subset(exclude.basis>0, K_i>ZERO))==0
            if (verbose){
                cat('Applying basis exclusion and removing', sum(keep==0), 'basis vectors\n')
            }
            K_i = K_i[keep, , drop = F]
        }
        
        ## only keep vectors that fail to intersect all vectors "exclude" in matrix
        if (!is.null(exclude.range)){
            A_i_abs = abs(A) %*% t(K_i)
            if ((as.numeric(nrow(exclude.range))*as.numeric*ncol(A_i_abs))>maxchunks){
                print('Exceeding maximum number of chunks in convex.basis computation')
                stop('Exceeding maximum number of chunks in convex.basis computation')
            }
            keep = Matrix::colSums(sparse_subset(exclude.range>0, t(A_i_abs), quiet = !verbose))==0
            if (verbose){
                cat('Applying range exclusion and removing', sum(keep==0), 'basis vectors\n')
            }
            K_i = K_i[keep, , drop = F]
        }

        A_i = A %*% t(K_i)
    }

    return(t(K_i))
}

#' @name all.paths
#' @title all.paths
#' @description
#'
#'
#' Low level function to enumerate all elementary paths and cycles through graph
#'
#' takes directed graph represented by n x n binary adjacency matrix  A and outputs all cycles and paths between source.vertices, sink.vertices
#'
#'
#' @param A nxn adjacency matrix
#' @param all logical flag, if all = T, will include all sources (parentless vertices) and sinks (childless vertices) in path computati
#' @param ALL logical flag, if ALL = T, will also include vertices without outgoing and incoming edges in paths
#' @param sources graph indices to treat as sources (by default is empty)
#' @param sinks graph indices to treat as sinks (by default is empty)
#' @param verbose logical flag
#' @return list of integer vectors corresponding to indices in A (i.e. vertices)
#' $paths = paths indices
#' $cycles = cycle indices
#' @keywords internal
#' @author Marcin Imielinski
#' @noRd
all.paths = function(A, all = F, ALL = F, sources = c(), sinks = c(), source.vertices = sources, sink.vertices = sinks,
                     exclude = NULL, ## specifies illegal subpaths, all such paths / cycles and
                     ## their supersets will be excluded, specified as k x nrow(A) matrix of vertex sets
                     verbose = FALSE,...)
{
    blank.vertices = Matrix::which(Matrix::rowSums(A)==0 & Matrix::colSums(A)==0)

    if (ALL)
        all = T

    if (all)
    {
        source.vertices = Matrix::which(Matrix::rowSums(A)>0 & Matrix::colSums(A)==0)
        sink.vertices = Matrix::which(Matrix::colSums(A)>0 & Matrix::rowSums(A)==0)
    }

    out = list(cycles = NULL, paths = NULL)

    node.ix = which(Matrix::rowSums(A!=0)>0 | Matrix::colSums(A!=0)>0)
    if (length(node.ix)==0)
        return(out)

    A = A[node.ix, node.ix]

    if (!is.null(exclude))
        exclude = sign(abs(exclude[, node.ix]))

    ij = Matrix::which(A!=0, arr.ind = T)
    B = Matrix::sparseMatrix(c(ij[,1], ij[,2]), rep(1:nrow(ij), 2), x = rep(c(-1, 1), each = nrow(ij)), dims = c(nrow(A), nrow(ij)))
    I = diag(rep(1, nrow(A)))

    source.vertices = setdiff(match(source.vertices, node.ix), NA)
    sink.vertices = setdiff(match(sink.vertices, node.ix), NA)

    B2 = cbind(B, I[, source.vertices, drop = FALSE], -I[, sink.vertices, drop = FALSE])

    if (verbose)
        cat(sprintf('Computing paths for %s vertices and %s edges\n', nrow(B2), ncol(B2)))

    K = convex.basis(B2, verbose = verbose, exclude.range = exclude, ...)

    if (all(is.na(K)))
        return(out)

    K = K[, Matrix::colSums(K[1:ncol(B), ,drop = FALSE])!=0, drop = FALSE] ## remove any pure source to sink paths

    is.cyc = Matrix::colSums(B %*% K[1:ncol(B), ,drop = FALSE]!=0)==0


    out$cycles = lapply(which(is.cyc),
                        function(i)
                        {
                            k = which(K[1:ncol(B), i]!=0)
                            v.all = unique(as.vector(ij[k, , drop = FALSE]))
                            sG = graph.edgelist(ij[k, , drop = FALSE])
                            tmp.v = v.all[c(1,length(v.all))]
                            p.fwd = get.shortest.paths(sG, tmp.v[1], tmp.v[2])
                            p.bwd = get.shortest.paths(sG, tmp.v[2], tmp.v[1])
                            return(node.ix[unique(unlist(c(p.fwd, p.bwd)))])
                        })

    out$paths = lapply(which(!is.cyc),
                       function(i)
                       {
                           k = K[1:ncol(B), i]
                           eix = which(k!=0)
                           v.all = unique(as.vector(ij[eix, , drop = FALSE]))
                           sG = graph.edgelist(ij[eix, , drop = FALSE])
                           io = B %*% k
                           v.in = Matrix::which(io<0)[1]
                           v.out = Matrix::which(io>0)[1]
                           return(node.ix[unlist(get.shortest.paths(sG, v.in, v.out))])
                       })

    if (length(out$cycles)>0)
    {
        tmp.cix = cbind(unlist(lapply(1:length(out$cycles), function(x) rep(x, length(out$cycles[[x]])))), unlist(out$cycles))
        out$cycles = out$cycles[!duplicated(as.matrix(Matrix::sparseMatrix(tmp.cix[,1], tmp.cix[,2], x = 1)))]
    }

    if (length(out$paths)>0)
    {
        tmp.pix = cbind(unlist(lapply(1:length(out$paths), function(x) rep(x, length(out$paths[[x]])))), unlist(out$paths))
        out$paths = out$paths[!duplicated(as.matrix(Matrix::sparseMatrix(tmp.pix[,1], tmp.pix[,2], x = 1)))]
    }

    if (ALL & length(blank.vertices)>0)
        out$paths = c(out$paths, lapply(blank.vertices, identity))

    return(out)
}


#' @name skrub
#' @title skrub
#' @description
#'
#' Converts data.table columns to standard types or replaces with NA and throws a warning
#'
#' Also converts factor columns to character columns in data.table, making
#' everyone's life easier. 
#'
#' @author Marcin Imielinski
#' @param dt data.table or data.frame
#' @keywords internal
#' @noRd
skrub = function(dt)
{
    cl = lapply(names(dt), function(x) class(dt[[x]]))
    names(cl) = names(dt)
    for (nm in names(cl)[cl=='factor'])
        dt[[nm]] = as.character(dt[[nm]])

    for (nm in names(cl)[cl=='integer'])
        dt[[nm]] = as.numeric(dt[[nm]])

    ## clean up any additional weird types before aggregating
    ALLOWED.CLASSES = c('integer', 'numeric', 'logical', 'character', 'list')
    if (any(tofix <- !sapply(dt, class) %in% ALLOWED.CLASSES))
    {
        warning(sprintf('found non-standard data types among one or more gEdge metadata columns (%s): converting to character before aggregating.  Consider manually converting these columns to one of the standard types: %s',
                        paste(names(dt)[tofix], collapse = ', '),
                        paste(ALLOWED.CLASSES, collapse = ', ')))
        
        for (fix in which(tofix))
        {
            replace = tryCatch(as.character(dt[[fix]]), error = function(e) NULL)
            
            if (is.null(replace))
            {
                warning(sprintf('Conversion of column character failed for column %s, replacing values with NA', names(dt)[fix]))
                replace = NA
            }
            dt[[fix]] = replace
        }
    }  
    return(dt)
}


#' @name sparse_subset
#' @title sparse_subset
#' @description
#'
#' given k1 x n matrix A and k2 x n matrix B
#' returns k1 x k2 matrix C whose entries ij = 1 if the set of nonzero components of row i of A is
#' a (+/- strict) subset of the nonzero components of row j of B
#'
sparse_subset = function(A, B, strict = FALSE, chunksize = 100, quiet = FALSE)
{
    nz = Matrix::colSums(as.matrix(A)!=0, 1)>0

    if (is.null(dim(A)) | is.null(dim(B)))
        return(NULL)

    C = Matrix::sparseMatrix(i = c(), j = c(), dims = c(nrow(A), nrow(B)))

    for (i in seq(1, nrow(A), chunksize))
    {
        ixA = i:min(nrow(A), i+chunksize-1)
        for (j in seq(1, nrow(B), chunksize))
        {
            ixB = j:min(nrow(B), j+chunksize-1)

            if (length(ixA)>0 & length(ixB)>0 & !quiet)
                cat(sprintf('\t interval A %s to %s (%d) \t interval B %d to %d (%d)\n', ixA[1], ixA[length(ixA)], nrow(A), ixB[1], ixB[length(ixB)], nrow(B)))
            if (strict)
                C[ixA, ixB] = (sign((A[ixA, , drop = FALSE]!=0)) %*% sign(t(B[ixB, , drop = FALSE]!=0))) * (sign((A[ixA, , drop = FALSE]==0)) %*% sign(t(B[ixB, , drop = FALSE]!=0))>0)
            else
                C[ixA, ixB] = (sign(A[ixA, nz, drop = FALSE]!=0) %*% sign(t(B[ixB, nz, drop = FALSE]==0)))==0
        }
    }

    return(C)
}


#' @name label.runs
#' @title label.runs
#' @description
#'
#' For logical input labels all instances of "TRUE" with a unique label and everything else as false
#'
#' For non-logical (e.g. character) input labels, labels each contiguous runs of the same value with a unique label
#' (note: even subsequent runs of an earlier used value in the vector will be given a new unique label)
#' 
#' 
#' @author Marcin Imielinski
#' @export
label.runs = function(x)
{
    if (!is.logical(x))
    {
        cumsum(abs(diff(as.numeric(c(0, as.integer(factor(x))))))>0)
    }
    else ## note will label all runs of FALSE with NA
    {
        as.integer(ifelse(x, cumsum(diff(as.numeric(c(FALSE, x)))>0), NA))
    }
}

#' @name dunlist
#' @title dunlist
#'
#' @description
#' unlists a list of vectors, matrices, data.tables into a data.table indexed by the list id
#' $listid
#'
#' does fill = TRUE in case the data.tables inside the list do not have compatible column names 
#' 
#' @param x list of vectors, matrices, or data frames
#' @return data.frame of concatenated input data with additional fields $ix and $iix specifying the list item and within-list index from which the given row originated from
#' @author Marcin Imielinski
#' @export
#' @keywords internal
#' @noRd 
#############################################################
dunlist = function(x)
{

    if (length(x)==0)
        return(data.table())

    if (is.null(names(x)))
        names(x) = seq_along(x)
    tmp = lapply(x, as.data.table)
    
    out = cbind(data.table(listid = rep(names(x), elementNROWS(x)), rbindlist(tmp, fill = TRUE)))
    setkey(out, listid)
    return(out)
}



#' @name read_vcf
#' @title parses VCF into GRanges or data.table
#'
#' @description
#'
#' Wrapper around Bioconductor VariantAnnotation. Reads VCF into GRanges or data.table format
#'
#' @param fn argument to parse via bcftools
#' @param gr GRanges input GRanges (default = NULL)
#' @param hg string Human reference genome (default = 'hg19')
#' @param geno boolean Flag whether to pull the genotype information information in the GENO vcf fields (default = NULL)  
#' @param swap.header string Pathn to another VCF file (in case of VCF with malformed header)(default = NULL)   
#' @param verbose boolean Flag (default = FALSE)
#' @param add.path boolean Flag to add the path of the current VCF file to the output (default = FALSE)
#' @param tmp.dir string Path to directory for temporary files (default = '~/temp/.tmpvcf')
#' @param ... extra parameters
#' @author Marcin Imielinski
#' @keywords internal
#' @noRd
read_vcf = function(fn, gr = NULL, hg = 'hg19', geno = NULL, swap.header = NULL, verbose = FALSE, add.path = FALSE, tmp.dir = '~/temp/.tmpvcf', ...)
{

    in.fn = fn

    if (verbose){
        cat('Loading', fn, '\n')
    }

    if (!is.null(gr)){

        tmp.slice.fn = paste(tmp.dir, '/vcf_tmp', gsub('0\\.', '', as.character(runif(1))), '.vcf', sep = '')
        cmd = sprintf('bcftools view %s %s > %s', fn,  paste(gr.string(gr.stripstrand(gr)), collapse = ' '), tmp.slice.fn)

        if (verbose){
            cat('Running', cmd, '\n')
        }
        system(cmd)
        fn = tmp.slice.fn
    }

    if (!is.null(swap.header)){

        if (!file.exists(swap.header)){
            stop(sprintf('Error: Swap header file %s does not exist\n', swap.header))
        }

        system(paste('mkdir -p', tmp.dir))
        tmp.name = paste(tmp.dir, '/vcf_tmp', gsub('0\\.', '', as.character(runif(1))), '.vcf', sep = '')
        if (grepl('gz$', fn)){
            system(sprintf("zcat %s | grep '^[^#]' > %s.body", fn, tmp.name))
        }
        else{
            system(sprintf("grep '^[^#]' %s > %s.body", fn, tmp.name))
        }

        if (grepl('gz$', swap.header)){
            system(sprintf("zcat %s | grep '^[#]' > %s.header", swap.header, tmp.name))
        }
        else{
            system(sprintf("grep '^[#]' %s > %s.header", swap.header, tmp.name))
        }

        system(sprintf("cat %s.header %s.body > %s", tmp.name, tmp.name, tmp.name))
        vcf = readVcf(tmp.name, hg, ...)
        system(sprintf("rm %s %s.body %s.header", tmp.name, tmp.name, tmp.name))

    }
    else{
        vcf = readVcf(fn, hg, ...)
    }

    out = granges(vcf)

    if (!is.null(values(out))){
        values(out) = cbind(values(out), info(vcf))
    }
    else{
        values(out) = info(vcf)
    }

    if (!is.null(geno)){

        if (!is.logical(geno)){
            geno = TRUE
        }

        if (geno){
            for (g in  names(geno(vcf))){
                geno = names(geno(vcf))
                warning(sprintf('Warning: Loading all geno fields:\n\t%s', paste(geno, collapse = ',')))
            }
        }

        gt = NULL

        if (length(g) > 0){

            for (g in geno){
                m = as.data.frame(geno(vcf)[[g]])
                names(m) = paste(g, names(m), sep = '_')
                if (is.null(gt)){
                    gt = m
                }
                else{
                    gt = cbind(gt, m)
                }
            }
            
            values(out) = cbind(values(out), as(gt, 'DataFrame'))
        }
    }

    if (!is.null(gr)){
        system(paste('rm', tmp.slice.fn))
    }

    if (add.path){
        values(out)$path = in.fn
    }

    return(out)
}




#' @name file.url.exists
#' @title Check if a file or url exists
#' @param f File or url
#' @return TRUE or FALSE
#' @importFrom RCurl url.exists
#' @noRd
file.url.exists <- function(f) {
    return(file.exists(f) || RCurl::url.exists(f))
}

#' @name read.rds.url
#' @title Checks if path is URL or file, then reads RDS
#' @param f File or url
#' @return data
#' @noRd
read.rds.url <- function(f) {
    if (grepl("^http",f))
        return(readRDS(gzcon(url(f))))
    return(readRDS(f))
}




#' @name ra.overlaps
#' @title ra.overlaps
#' @description
#'
#' Determines overlaps between two piles of rearrangement junctions ra1 and ra2 (each GRangesLists of signed locus pairs)
#' against each other, returning a sparseMatrix that is T at entry ij if junction i overlaps junction j.
#'
#' if argument pad = 0 (default) then only perfect overlap will validate, otherwise if pad>0 is given, then
#' padded overlap is allowed
#'
#' strand matters, though we test overlap of both ra1[i] vs ra2[j] and gr.flipstrand(ra2[j])
#'
#' @param ra1 \code{GRangesList} with rearrangement set 1
#' @param ra2 \code{GRangesList} with rearrangement set 2
#' @param pad Amount to pad the overlaps by. Larger is more permissive. Default is exact (0)
#' @param arr.ind Default TRUE
#' @param ignore.strand Ignore rearrangement orientation when doing overlaps. Default FALSE
#' @param ... params to be sent to \code{\link{gr.findoverlaps}}
#' @name ra.overlaps
#' @keywords internal
#' @noRd
ra.overlaps = function(ra1, ra2, pad = 0, arr.ind = TRUE, ignore.strand=FALSE, ...)
{
    bp1 = grl.unlist(ra1) + pad
    bp2 = grl.unlist(ra2) + pad
    ix = gr.findoverlaps(bp1, bp2, ignore.strand = ignore.strand, ...)

    .make_matches = function(ix, bp1, bp2)
    {
        if (length(ix) == 0){
            return(NULL)
        }
        tmp.match = cbind(bp1$grl.ix[ix$query.id], bp1$grl.iix[ix$query.id], bp2$grl.ix[ix$subject.id], bp2$grl.iix[ix$subject.id])
        tmp.match.l = lapply(split(1:nrow(tmp.match), paste(tmp.match[,1], tmp.match[,3])), function(x) tmp.match[x, , drop = F])

        ## match only occurs if each range in a ra1 junction matches a different range in the ra2 junction
        matched.l = sapply(tmp.match.l, function(x) all(c('11','22') %in% paste(x[,2], x[,4], sep = '')) | all(c('12','21') %in% paste(x[,2], x[,4], sep = '')))

        return(do.call('rbind', lapply(tmp.match.l[matched.l], function(x) cbind(x[,1], x[,3])[!duplicated(paste(x[,1], x[,3])), , drop = F])))
    }

    tmp = .make_matches(ix, bp1, bp2)

    if (is.null(tmp)){
        if (arr.ind){

            return(as.matrix(data.table(ra1.ix = as.numeric(NA), ra2.ix = as.numeric(NA))))
        }
        else{
            return(Matrix::sparseMatrix(length(ra1), length(ra2), x = 0))
        }
    }

    rownames(tmp) = NULL

    colnames(tmp) = c('ra1.ix', 'ra2.ix')

    if (arr.ind) {
        ro = tmp[order(tmp[,1], tmp[,2]), , drop = FALSE]
        if (class(ro)=='integer'){
            ro <- matrix(ro, ncol=2, nrow=1, dimnames=list(c(), c('ra1.ix', 'ra2.ix')))
        }
        return(ro)
    } else {
        ro = Matrix::sparseMatrix(tmp[,1], tmp[,2], x = 1, dims = c(length(ra1), length(ra2)))
        return(ro)
    }
}

#' @name gt.gencode
#' @description
#'
#' internal function to format transcript annotations for fusion output
#'
#' @param gencode GRanges output of rtracklayer::import of GENCODE gtf
#' @param bg.col character representing color to put in colormap
#' @param cds.col color for CDS
#' @param utr.col color for UTR
#' @param st.col scalar character representing color of CDS start
#' @param en.col scalar character representing color of CDS end
#' @keywords internal
#' @noRd
#' @return gTrack object of gencode output
gt.gencode = function(gencode, bg.col = alpha('blue', 0.1), cds.col = alpha('blue', 0.6), utr.col = alpha('purple', 0.4), st.col = 'green',
                      en.col = 'red')  
{
    if (length(gencode)==0)
        return(gTrack())

    tx = gencode[gencode$type =='transcript']
    genes = gencode[gencode$type =='gene']
    exons = gencode[gencode$type == 'exon']
    utr = gencode[gencode$type == 'UTR']
    ## ut = unlist(utr$tag)
    ## utix = rep(1:length(utr), sapply(utr$tag, length))
    ## utr5 = utr[unique(utix[grep('5_UTR',ut)])]
    ## utr3 = utr[unique(utix[grep('3_UTR',ut)])]
    ## utr5$type = 'UTR5'
    ## utr3$type = 'UTR3'
    startcodon = gencode[gencode$type == 'start_codon']
    stopcodon = gencode[gencode$type == 'stop_codon']
    OUT.COLS = c('gene_name', 'transcript_name', 'transcript_id', 'type', 'exon_number', 'type')
    tmp = c(genes, tx, exons, utr, startcodon, stopcodon)[, OUT.COLS]
    
    ## compute tx ord of intervals
    ord.ix = order(tmp$transcript_id, match(tmp$type, c('gene', 'transcript', 'exon', 'UTR', 'start_codon','stop_codon')))
    tmp.rle = rle(tmp$transcript_id[ord.ix])
    tmp$tx.ord[ord.ix] = unlist(lapply(tmp.rle$lengths, function(x) 1:x))
    tmp = tmp[rev(order(match(tmp$type, c('gene', 'transcript', 'exon', 'UTR', 'start_codon','stop_codon'))))] 
    tmp.g = tmp[tmp$type != 'transcript']
    cmap = list(type = c(gene = bg.col, transcript = bg.col, exon = cds.col, start_codon = st.col, stop_codon = en.col, UTR = utr.col))
    tmp.g = gr.disjoin(gr.stripstrand(tmp.g))
    return(gTrack(tmp.g[, c('type', 'gene_name')], colormaps = cmap))
}


#' @name alpha
#' @title alpha
#' @description
#' Give transparency value to colors
#'
#' Takes provided colors and gives them the specified alpha (ie transparency) value
#'
#' @author Marcin Imielinski
#' @param col RGB color
#' @keywords internal
#' @noRd
alpha = function(col, alpha)
{
    col.rgb = grDevices::col2rgb(col)
    out = grDevices::rgb(red = col.rgb['red', ]/255, green = col.rgb['green', ]/255, blue = col.rgb['blue', ]/255, alpha = alpha)
    names(out) = names(col)
    return(out)
}


#' @name dunlist
#'  @title dunlist
#'
#' @description
#' unlists a list of vectors, matrices, data.tables into a data.table indexed by the list id
#' $listid
#'
#' does fill = TRUE in case the data.tables inside the list do not have compatible column names
#'
#' @param x list of vectors, matrices, or data frames
#' @return data.frame of concatenated input data with additional fields $ix and $iix specifying the list item and within-list index from which the given row originated from
#' @author Marcin Imielinski
#' @keywords internal
#' @noRd
dunlist = function(x)
{
    listid = rep(1:length(x), elementNROWS(x))

    if (!is.null(names(x))) ## slows things down
        listid = names(x)[listid]
    
    xu = unlist(x, use.names = FALSE)  

    if (is.null(xu))
    {
      return(as.data.table(list(listid = c(), V1 = c())))
    }
    
    if (!(inherits(xu, 'data.frame')) | inherits(xu, 'data.table'))
        xu = data.table(V1 = xu)
    
    
    out = cbind(data.table(listid = listid), xu)
    setkey(out, listid)
    return(out)  
}


#' @name pdist
#' @title pdist
#'
#' @description
#'
#' Given two GRanges gr1 and gr2 each of the same length, returns the reference
#' distance between them, subject to ignore.strand = TRUE
#' 
#' @param gr1 GRanges
#' @param gr2 GRAnges
#' @return vector of 
#' @author Marcin Imielinski
#' @keywords internal
#' @noRd
pdist = function(gr1, gr2, ignore.strand = TRUE)
{
    if (length(gr1) != length(gr2))
        stop('arguments have to be the same length')

    d = ifelse(as.logical(seqnames(gr1) != seqnames(gr2)), Inf,
               pmin(abs(start(gr1)-start(gr2)),
                    abs(start(gr1)-end(gr2)),
                    abs(end(gr1)-start(gr2)),
                    abs(end(gr1)-end(gr2))))

    if (!ignore.strand && any(ix <- strand(gr1) != strand(gr2)))
        d[as.logical(ix)] = Inf

    return(d)
}


#' @name ra.duplicated
#' @title ra.duplicated
#' @description
#'
#' Show if junctions are Deduplicated
#'
#' Determines overlaps between two or more piles of rearrangement junctions (as named or numbered arguments) +/- padding
#' and will merge those that overlap into single junctions in the output, and then keep track for each output junction which
#' of the input junctions it was "seen in" using logical flag  meta data fields prefixed by "seen.by." and then the argument name
#' (or "seen.by.ra" and the argument number)
#'
#' @author Xiaotong Yao
#' @param grl GRangesList representing rearrangements to be merged
#' @param pad non-negative integer specifying padding
#' @param ignore.strand whether to ignore strand (implies all strand information will be ignored, use at your own risk)
#' @return \code{GRangesList} of merged junctions with meta data fields specifying which of the inputs each outputted junction was "seen.by"
#' @name ra.duplicated
ra.duplicated = function(grl, pad=500, ignore.strand=FALSE){

    if (!is(grl, "GRangesList")){
        stop("Error: Input must be GRangesList!")
    }

    ##if (any(elementNROWS(grl)!=2)){
    ##    stop("Error: Each element must be length 2!")
    ##}

    if (length(grl)==0){
        return(logical(0))
    }

    if (length(grl)==1){
        return(FALSE)
    }

    if (length(grl)>1){

        ix.pair = as.data.table(ra.overlaps(grl, grl, pad=pad, ignore.strand = ignore.strand))[ra1.ix!=ra2.ix]

        if (nrow(ix.pair)==0){
            return(rep(FALSE, length(grl)))
        }
        else {
            ##           dup.ix = unique(rowMax(as.matrix(ix.pair)))
            dup.ix = unique(apply(as.matrix(ix.pair), 1, max))
            return(seq_along(grl) %in% dup.ix)
        }
    }
}



#' Merges rearrangements represented by \code{GRangesList} objects
#'
#' Determines overlaps between two or more piles of rearrangement junctions (as named or numbered arguments) +/- padding
#' and will merge those that overlap into single junctions in the output, and then keep track for each output junction which
#' of the input junctions it was "seen in" using logical flag  meta data fields prefixed by "seen.by." and then the argument name
#' (or "seen.by.ra" and the argument number)
#'
#' @param ... GRangesList representing rearrangements to be merged
#' @param pad non-negative integer specifying padding
#' @param ind  logical flag (default FALSE) specifying whether the "seen.by" fields should contain indices of inputs (rather than logical flags) and NA if the given junction is missing
#' @param ignore.strand whether to ignore strand (implies all strand information will be ignored, use at your own risk)
#' @return \code{GRangesList} of merged junctions with meta data fields specifying which of the inputs each outputted junction was "seen.by"
#' @name ra.merge
#' @export
#' @examples
#'
#' # generate some junctions
#' gr1 <- GRanges(1, IRanges(1:10, width = 1), strand = rep(c('+', '-'), 5))
#' gr2 <- GRanges(1, IRanges(4 + 1:10, width = 1), strand = rep(c('+', '-'), 5))
#' ra1 = split(gr1, rep(1:5, each = 2))
#' ra2 = split(gr2, rep(1:5, each = 2))
#'
#' ram = ra.merge(ra1, ra2)
#' values(ram) # shows the metadata with TRUE / FALSE flags
#'
#' ram2 = ra.merge(ra1, ra2, pad = 5) # more inexact matching results in more merging
#' values(ram2)
#'
#' ram3 = ra.merge(ra1, ra2, ind = TRUE) #indices instead of flags
#' values(ram3)
ra.merge = function(..., pad = 0, ind = FALSE, ignore.strand = FALSE){
    ra = list(...)
    ra = ra[which(!sapply(ra, is.null))]
    nm = names(ra)
    if (is.null(nm)){
        nm = paste('ra', 1:length(ra), sep = '')
    }
    nm = paste('seen.by', nm, sep = '.')
    if (length(nm)==0){
        return(NULL)
    }
    out = ra[[1]]
    values(out) = cbind(as.data.frame(matrix(FALSE, nrow = length(out), ncol = length(nm), dimnames = list(NULL, nm))), values(out))

    if (!ind){
        values(out)[, nm[1]] = TRUE
    } else{
        values(out)[, nm[1]] = 1:length(out)
    }

    if (length(ra)>1){
        for (i in 2:length(ra)){
            this.ra = ra[[i]]
            if (length(this.ra)>0){
                values(this.ra) = cbind(as.data.frame(matrix(FALSE, nrow = length(this.ra), ncol = length(nm), dimnames = list(NULL, nm))), values(this.ra))
                ovix = ra.overlaps(out, this.ra, pad = pad, ignore.strand = ignore.strand)

                if (!ind){
                    values(this.ra)[[nm[i]]] = TRUE
                } else{
                    values(this.ra)[[nm[i]]] = 1:length(this.ra)
                }

                if (!ind){
                    if (!all(is.na(ovix))){
                        values(out)[, nm[i]][ovix[,1]] = TRUE
                    }
                } else{
                    values(out)[, nm[i]] = NA
                    if (!all(is.na(ovix))){
                        values(out)[, nm[i]][ovix[,1]] = ovix[,1]
                    }
                }
                ## which are new ranges not already present in out, we will add these
                if (!all(is.na(ovix))){
                    nix = setdiff(1:length(this.ra), ovix[,2])
                } else{
                    nix = 1:length(this.ra)
                }

                if (length(nix)>0){
                    val1 = values(out)
                    val2 = values(this.ra)
                    if (ind){
                        val2[, nm[1:(i-1)]] = NA
                    }
                    else{
                        val2[, nm[1:(i-1)]] = FALSE
                    }
                    values(out) = NULL
                    values(this.ra) = NULL
                    out = grl.bind(out, this.ra[nix])
                    values(out) = rrbind(val1, val2[nix, ])
                }
            }
        }
    }
    return(out)
}

#' @name tstamp
#' @title tstamp
#' @description
#' Timestamp used to check for staleness of gGraph and other objects 
#' @keywords internal
#' @noRd 
tstamp = function()
{
    return(paste(as.character(Sys.time()), runif(1)))
}

#' @name dodo.call
#' @title dodo.call
#' @description
#' do.call implemented using eval parse for those pesky (e.g. S4) case when do.call does not work
#' @keywords internal
#' @noRd 
dodo.call = function(FUN, args)
{
    if (!is.character(FUN))
        FUN = substitute(FUN)
    cmd = paste(FUN, '(', paste('args[[', 1:length(args), ']]', collapse = ','), ')', sep = '')
    return(eval(parse(text = cmd)))
}

#' @name dedup
#' @title dedup
#'
#' @description
#' relabels duplicates in a character vector with .1, .2, .3
#' (where "." can be replaced by any user specified suffix)
#'
#' @param x input vector to dedup
#' @param suffix suffix separator to use before adding integer for dups in x
#' @return length(x) vector of input + suffix separator + integer for dups and no suffix for "originals"
#' @author Marcin Imielinski
#' @noRd
dedup = function(x, suffix = '.')
{
    dup = duplicated(x);
    udup = setdiff(unique(x[dup]), NA)
    udup.ix = lapply(udup, function(y) which(x==y))
    udup.suffices = lapply(udup.ix, function(y) c('', paste(suffix, 2:length(y), sep = '')))
    out = x;
    out[unlist(udup.ix)] = paste(out[unlist(udup.ix)], unlist(udup.suffices), sep = '');
    return(out)
}


#' @name spmelt
#' @title spMelt
#' @description
#' Melts sparse matrix into data.table
#' 
#' @param A 
#' @return data.table of all non
#' @author Marcin Imielinski
spmelt = function(A, baseval = 0) {
    if (!is.null(baseval))
    {
        ij = Matrix::which(A!=baseval, arr.ind = TRUE)
    }
    else ## take everything
    {
        ij = as.matrix(expand.grid(1:nrow(A), 1:ncol(A)))
    }
    dt = data.table(i = ij[,1], j = ij[,2], val = A[ij])
}






#' @name gstat
#'
#' @export
gstat = function(gg,
                 INF = max(seqlengths(gg)) + 1,
                 thresh = 1e6){
    if (!inherits(gg, "gGraph")){
        stop("Input is not a gGraph object.")
    }
    
    if (length(gg$nodes)==0 |
        length(gg$edges)==0){
        return(NULL)
    }
    if (!is.element("cn", colnames(gg$nodes$dt)) |
        !is.element("cn", colnames(gg$edges$dt)) |
        !any(gg$edges$dt[, type=="ALT"])){
        return(NULL)
    }
    if (is.na(INF)){
        INF = 1e9
    }
    ## needed "fb" in the edge table
    if (!is.element("fb", colnames(gg$edges$dt))){
        gg$annotate("fb",
                    data=gg$junctions$span<=1e4 & gg$junctions$sign==1,
                    id=gg$junctions$dt$edge.id,
                    class="edge")
    }
    ## needed "term" field in the node table
    if (!is.element("term", colnames(gg$nodes$dt)) |
        !is.element("fused", colnames(gg$nodes$dt))){
        n2e = rbind(
            gg$edges$dt[, .(nid = n1,
                            type,
                            side = n1.side)],
            gg$edges$dt[, .(nid = n2,
                            type,
                            side = n2.side)])
        n2e[, term := length(unique(side))<2, by=nid]
        gg$annotate("term",
                    n2e[, term],
                    n2e[, nid],
                    "node")
        fused = gg$nodes$dt$node.id %in% n2e[type=="ALT", unique(nid)]
        gg$nodes$mark(fused = fused)
    }
    edt = gg$edges$dt
    ndt = gg$nodes$dt
    ## four traditional classes, plus fold back
    n.del = edt[, sum(class=="DEL-like", na.rm = TRUE)]
    n.dup = edt[, sum(class=="DUP-like", na.rm = TRUE)]
    n.inv = edt[, sum(class=="INV-like", na.rm = TRUE)]
    n.tra = edt[, sum(class=="TRA-like", na.rm = TRUE)]
    n.fb = edt[, sum(fb==TRUE, na.rm = TRUE)]
    ## highest CN of a tra, dup, fb
    max.cn.del = edt[class=="DEL-like", pmax(max(cn, na.rm=T), 0)]
    max.cn.dup = edt[class=="DUP-like", pmax(max(cn, na.rm=T), 0)]
    max.cn.inv = edt[class=="INV-like", pmax(max(cn, na.rm=T), 0)]
    max.cn.tra = edt[class=="TRA-like", pmax(max(cn, na.rm=T), 0)]
    max.cn.fb = edt[fb==TRUE, pmax(max(cn, na.rm=T), 0)]
    ## edge wise features
    this.alt.sg = gg[, type=="ALT"]
    diam = this.alt.sg$diameter
    alt.on.diam = length(diam$edges)
    ## walk the whole ALT graph
    alt.wks = this.alt.sg$walks()
    circ.ix = alt.wks$dt[circular==TRUE, walk.id]
    n.circ = length(circ.ix)
    if (n.circ == 0){
        max.len.circ = 0
        max.cn.circ = 0
    } else {
        max.len.circ = alt.wks$dt[circular==TRUE, max(length)]
        max.cn.circ = max(sapply(circ.ix,
                                 function(y) {
                                     alt.wks[y]$edges$dt[, min(cn)]
                                 }),
                          na.rm=T)
    }
    ## number of junctions
    n.junc = length(gg$edges)
    jcn.tab = table(edt[type=="ALT", cn])
    ## junction copy numbers
    n.jcn1 = ifelse(is.na(jcn.tab['1']), 0, jcn.tab['1'])
    n.jcn1.prop = n.jcn1/n.junc
    n.jcn2 = ifelse(is.na(jcn.tab['2']), 0, jcn.tab['2'])
    n.jcn2.prop = n.jcn2/n.junc
    n.jcn3p = n.junc - n.jcn1 - n.jcn2
    max.jcn = edt[type=="ALT", max(cn, na.rm=T)]
    ## diam cn1
    if (n.jcn1>0){
        this.alt.1.sg = this.alt.sg[, cn==1]
        alt.1.diam = this.alt.1.sg$diameter
        alt.1.on.diam = length(alt.1.diam$edges)
    } else {
        alt.1.on.diam = 0
    }
    ## diam cn2
    if (n.jcn2>0){
        this.alt.2.sg = this.alt.sg[, cn==2]
        alt.2.diam = this.alt.2.sg$diameter
        alt.2.on.diam = length(alt.2.diam$edges)
    } else {
        alt.2.on.diam = 0
    }
    ## node wise features
    n.fused = ndt[, sum(fused==TRUE, na.rm=T)]
    n.unfus = ndt[, sum(fused==FALSE, na.rm=T)]
    ## fused sized and unfused sizes, without terminal node
    fused.size.med = median(ndt[term==FALSE & fused==TRUE, width])
    fused.size.mean = mean(ndt[term==FALSE & fused==TRUE, width])
    unfus.size.med = median(ndt[term==FALSE & fused==FALSE, width])
    if (is.na(unfus.size.med)){unfus.size.med = 0}
    unfus.size.mean = mean(ndt[term==FALSE & fused==FALSE, width])
    if (is.na(unfus.size.mean)){unfus.size.mean = 0}
    ## foot print
    n.chr = length(unique(as.character(seqnames(gg$gr))))
    ## excluding terminal nodes
    footprint = gg[is.na(term) | term==FALSE]$footprint
    if (length(footprint)>1){
        footprint.ds = gr.dist(footprint, footprint)
        footprint.ds[is.na(footprint.ds)] = INF
        footprint.ds[is.infinite(footprint.ds)] = INF
        footprint.cl = hclust(as.dist(footprint.ds), "single")
        footprint.gp = cutree(footprint.cl, h = thresh)
        n.fp.gp = length(unique(footprint.gp))
    } else {
        n.fp.gp = 1
    }
    ## node CN
    cn.min = ndt[, min(cn, na.rm=T)]
    cn.max = ndt[, max(cn, na.rm=T)]
    cn.tab = ndt[, table(cn)]
    cn.mode = as.numeric(names(which.max(cn.tab)))
    cn.mode.prop = ndt[, sum(cn==cn.mode)/.N]
    cn.states = ndt[, length(unique(cn))]
    ## fused CN
    cn.fused.mode = as.numeric(
        names(which.max(ndt[fused==TRUE, table(cn)])))
    cn.fused.mode.prop = ndt[, sum(fused==TRUE & cn==cn.fused.mode)/sum(fused)]
    ## unfused CN
    cn.unfus.mode = as.numeric(
        names(which.max(ndt[fused==FALSE, table(cn)]))
    )
    if (length(cn.unfus.mode)==0){
        cn.unfus.mode = 0
        cn.unfus.mode.prop = 0
    } else {
        cn.unfus.mode.prop = ndt[, sum(fused==FALSE & cn==cn.unfus.mode)/sum(!fused)]
    }
    out =
        data.table(alt.all = n.junc,
                   n.del,
                   n.dup,
                   n.inv,
                   n.tra,
                   max.cn.del,
                   max.cn.dup,
                   max.cn.inv,
                   max.cn.tra,
                   max.cn.fb,
                   n.jcn1,
                   n.jcn2,
                   n.jcn3p,
                   n.jcn1.prop,
                   n.jcn2.prop,
                   alt.on.diam,
                   alt.on.diam.prop = alt.on.diam/n.junc,
                   alt.1.on.diam,
                   alt.1.on.diam.prop = ifelse(n.jcn1==0, 0, alt.1.on.diam/n.jcn1),
                   alt.2.on.diam,
                   alt.2.on.diam.prop = ifelse(n.jcn2==0, 0, alt.2.on.diam/n.jcn2),
                   n.fused, fused.size.med, fused.size.mean,
                   n.unfus, unfus.size.med, unfus.size.mean,
                   n.chr,
                   n.fp = length(footprint),
                   n.fp.gp,
                   cn.mode, cn.mode.prop,
                   cn.fused.mode, cn.fused.mode.prop,
                   cn.unfus.mode, cn.unfus.mode.prop,
                   cn.states,
                   cn.max, cn.min, max.jcn,
                   max.len.circ, max.cn.circ)
    return(out)
}

#' @name cov2csv
#' @description
#' prepare csv file for gGnome.js
#' Col1 x
#' Col2 y
#' Col3 cumulative x
cov2csv = function(x,
        field = "ratio",
        id = "data",
        outdir = "./coverage/",
        meta.js = NULL)
{            
    ## first x must be either a GRanges,
    ## a RDS file containing a GRanges,
    ## or a TXT file that can be read as a GRanges
    if (inherits(x, "character")){
        if (!file.exists(x)){
            stop("Input file not found.")
        }
        fn = x
        if (grepl("rds$", fn)){
            x = readRDS(fn)
        } else if (grepl("[(bed)|(bw)|(wig)]$", fn)){
            x = rtracklayer::import(fn)
        } else if (grepl("[(txt)|(tsv)]$", fn)) {
            x = dt2gr(fread(fn))
        } else {
            stop("Input file not in valid format: txt, tsv, bed, bw, wig, rds")
        }
    }

    ## respect the seqlengths in meta.js
    if (is.character(meta.js) && file.exists(meta.js)){
        meta = rbindlist(jsonlite::read_json(meta.js)$metadata)
        sl = meta[, setNames(endPoint, chromosome)]
    }

    if (inherits(x, "GRanges")){
        if (!exists("sl")){
            sl = seqlengths(x)
        }
        if (all(is.na(sl))){
            stop("No seqlengths in the input.")
        }
    } else {
        stop("input is not a GRanges")
    }
    ## build the cumulative coordinates
    dt = data.table(seqlevels = names(sl),
                    seqlengths = as.double(sl),
                    cstart = c(1, cumsum(shift(as.double(sl))[-1])))
    ## build the data.table
    dat = merge(gr2dt(x), dt,
                by.x = "seqnames",
                by.y = "seqlevels",
                all.x = TRUE)
    if (!is.element(field, colnames(dat))){
        stop("'field' is not in the input data")
    }

    ## write to output file
    dat[, new.start := start + cstart - 1]

    if (!dir.exists(outdir)){
        dir.create(outdir)
    }

    dir.create(paste0(outdir, "/", id))

    dat[!is.na(get(field)) & !is.infinite(get(field)), {
        this.fn = paste0(outdir, "/",
                         id, "/",
                         id, ".",
                         seqnames, ".csv")
        write.table(.SD[, .(x = start, y = get(field), chromosome = seqnames, place = -1e-4 + new.start)],
                    file = this.fn,
                    quote = FALSE,
                    col.names = TRUE,
                    row.names = FALSE,
                    sep = ",")
    }, by = seqnames]
    return(paste0(outdir, "/", id))
}

#' @name j.dist
#' @description
#' @param j1
#' @param j2
#' @export
j.dist = function(j1, j2 = NULL){
    if (is.null(j2)){
        j2 = j1
    }
    if (!inherits(j1, "Junction")){
        j1 = Junction$new(grl = j1)
    }
    if (!inherits(j2, "Junction")){
        j2 = Junction$new(grl = j2)
    }
    ij = data.table(expand.grid(
        list(i = seq_along(j1),
             j = seq_along(j2))))[i<j]
    ij[, ":="(i1)]
}

#' @name jab2json
#' @description a wrapper function to dump JaBbA results run with Flow to gGnome.js viz
#' @export
jab2json = function(fn = "./jabba.simple.rds",
                    gGnome.js.dir = "~/git/gGnome.js"){

}
