#' @name gen_PGV_instance
#' @rdname gen_js_instance
#' @description internal
#'
#' Generate a PGV instance
#'
#' Takes a table with paths to gGraphs and coverage files (optional) and generates an instance of a gGnome.js directory that is ready to visualize using gGnome.js
#' 
#' @param data either a path to a TSV/CSV or a data.table
#' @param outdir the path where to save the files. This path should not exist, unless you want to add more files to an existing directory in which case you must use --append
#' 
#' @export
gen_PGV_instance = function(data,
                           name.col = 'pair',
                           outdir = './gGnome.js',
                           cov.col = 'cov',
                           gg.col = 'complex',
                           append = FALSE,
                           cov.field = 'ratio',
                           cov.field.col = NA,
                           cov.bin.width = 1e4,
                           color.field = NULL,
                           dataset_name = NA,
                           ref.name = NA,
                           overwrite = FALSE,
                           annotation = c('simple', 'bfb', 'chromoplexy',
                                       'chromothripsis', 'del', 'dm', 'dup',
                                       'pyrgo', 'qrdel', 'qrdup', 'qrp', 'rigma',
                                       'tic', 'tyfonas'),
                           mc.cores = 1
                     ){
    return(gen_js_instance(data = data,
                           name.col = name.col,
                           outdir = outdir,
                           cov.col = cov.col,
                           gg.col = gg.col,
                           append = append,
                           js.type = 'PGV',
                           cov.field = cov.field,
                           cov.field.col = cov.field.col,
                           cov.bin.width = cov.bin.width,
                           color.field = color.field,
                           dataset_name = dataset_name,
                           ref.name = ref.name,
                           overwrite = overwrite,
                           annotation = annotation,
                           mc.cores = mc.cores))
}

#' @name gen_gGnomejs_instance
#' @rdname gen_js_instance
#' @description internal
#'
#' Generate a gGnome.js instance
#'
#' Takes a table with paths to gGraphs and coverage files (optional) and generates an instance of a gGnome.js directory that is ready to visualize using gGnome.js
#' 
#' @param data either a path to a TSV/CSV or a data.table
#' @param outdir the path where to save the files. This path should not exist, unless you want to add more files to an existing directory in which case you must use --append
#' 
#' @export
gen_gGnomejs_instance = function(data,
                           name.col = 'pair',
                           outdir = './gGnome.js',
                           cov.col = 'cov',
                           gg.col = 'complex',
                           append = FALSE,
                           cov.field = 'ratio',
                           cov.field.col = NA,
                           cov.bin.width = 1e4,
                           ref.name = NA, #TODO: we need 
                           overwrite = FALSE,
                           annotation = c('simple', 'bfb', 'chromoplexy',
                                       'chromothripsis', 'del', 'dm', 'dup',
                                       'pyrgo', 'qrdel', 'qrdup', 'qrp', 'rigma',
                                       'tic', 'tyfonas'),
                           mc.cores = 1
                     ){
    return(gen_js_instance(data = data,
                           name.col = name.col,
                           outdir = outdir,
                           cov.col = cov.col,
                           gg.col = gg.col,
                           append = append,
                           js.type = 'gGnome.js',
                           cov.field = cov.field,
                           cov.field.col = cov.field.col,
                           cov.bin.width = cov.bin.width,
                           ref.name = ref.name,
                           overwrite = overwrite,
                           annotation = annotation,
                           mc.cores = mc.cores))
}


#' @name gen_js_instance
#' @description internal
#'
#' Generate a gGnome.js instance
#'
#' Takes a table with paths to gGraphs and coverage files (optional) and generates an instance of a gGnome.js directory that is ready to visualize using gGnome.js
#' 
#' @param data either a path to a TSV/CSV or a data.table
#' @param outdir the path where to save the files. This path should not exist, unless you want to add more files to an existing directory in which case you must use --append
#' 
#' @export
gen_js_instance = function(data,
                           name.col = 'pair',
                           outdir = './gGnome.js',
                           cov.col = 'cov',
                           gg.col = 'complex',
                           append = FALSE,
                           js.type = 'gGnome.js',
                           cov.field = 'ratio',
                           cov.field.col = NA,
                           cov.bin.width = 1e4,
                           color.field = NULL,
                           dataset_name = NA,
                           ref.name = NA,
                           overwrite = FALSE,
                           annotation = c('simple', 'bfb', 'chromoplexy',
                                       'chromothripsis', 'del', 'dm', 'dup',
                                       'pyrgo', 'qrdel', 'qrdup', 'qrp', 'rigma',
                                       'tic', 'tyfonas'),
                           mc.cores = 1
                     ){
    # check the path and make a clone of the github repo if needed
    outdir = js_path(outdir, js.type = js.type, append = append)

    # get the path to the metadata file
    meta.js = get_path_to_meta_js(outdir, js.type = js.type)

    # read and check the input data
    data = read.js.input.data(data, name.col = name.col)

    # generate coverage files
    message('Generating coverage files')
    coverage_files = gen_js_coverage_files(data, outdir, name.col = name.col, overwrite = overwrite, cov.col = cov.col,
                          js.type = js.type, cov.field = cov.field,
                          cov.field.col = cov.field.col,
                          bin.width = cov.bin.width, dataset_name = dataset_name,
                          ref.name = ref.name, color.field = color.field,
                          meta.js = meta.js, mc.cores = mc.cores)

    data$coverage = coverage_files

    message('Generating json files')
    gg.js.files = gen_gg_json_files(data, outdir, meta.js = meta.js, name.col = name.col, gg.col = gg.col,
                                    js.type = js.type, dataset_name = dataset_name, ref.name = ref.name,
                                    overwrite = overwrite, annotation = annotation)

    data$gg.js = gg.js.files

    # generate the datafiles
    dfile = gen_js_datafiles(data, outdir, js.type, name.col = name.col, ref.name = ref.name, dataset_name = dataset_name)
}

gen_js_datafiles = function(data, outdir, js.type, name.col = NA, meta_col = NA, ref.name = NA, dataset_name = NA){
    dfile = get_js_datafiles_path(outdir, js.type)

    if (is.na(meta_col)){
        # if there is no metadata column then add empty metadata values in a new "meta_col" column
        meta_col = 'description'
        data$description = ''
    }

    message(paste0('Writing description file to: ', dfile))
    if (js.type == 'gGnome.js'){
        if (file.exists(dfile)){
            datafiles = fread(dfile)
            # if some of the samples that we are adding are already in the datafiles then we want to override these
            jsons = data[, paste0(get(name.col), '.json')]
            datafiles_trimmed = datafiles[!(datafile %in% jsons)]
            datafiles = rbind(datafiles_trimmed, data[, .(datafile = paste0(get(name.col), '.json'), description)])
        } else {
            datafiles = data[, .(datafile = paste0(get(name.col), '.json'), description)]
        }
        fwrite(datafiles, dfile)
    }
    if (js.type == 'PGV'){
        if (is.na(dataset_name)){
            stop('dataset_name must be provided for PGV.')
        }
        if (is.na(ref.name)){
            stop('ref.name must be provided for PGV.')
        }

        if (!('visible' %in% names(data))){
            data$visible = TRUE
        }

        plots = lapply(1:data[,.N], function(idx){
                     gg.js = data[idx, gg.js]
                     cov.fn = data[idx, coverage]
                     gg.track = NULL
                     cov.track = NULL
                     if (!is.na(gg.js)){
                         if (file.exists(gg.js)){
                             gg.track = list('type' = 'genome',
                                             'source' = paste0(data[idx, get(name.col)], '.json'),
                                             'title' = data[idx, get(name.col)],
                                             'visible' = ifelse(data[idx, visible] == TRUE, TRUE, FALSE))
                         }
                     }
                     if (!is.na(cov.fn)){
                         if (file.exists(cov.fn)){
                             cov.track = list('type' = 'scatterplot',
                                             'source' = paste0(data[idx, get(name.col)], '-coverage.arrow'),
                                             'title' = paste0(data[idx, get(name.col)], ' Coverage Distribution'),
                                             'visible' = FALSE) # coverage tracks will always be set to not visible on load
                         }
                     }
                     tracks = list(gg.track, cov.track)
                     return(tracks)
        })

        plots = do.call(c, plots)

        item = list(filename = paste0(dataset_name, '.json')) # TODO: once this trello task is resolved then we need to update this line: https://trello.com/c/VZALB1we
        item$description = list(paste0('dataset=', dataset_name)) # TODO: we need to figure out the purpose of the description in PGV and update this accordingly
        item$reference = ref.name
        item$plots = plots

        if (file.exists(dfile)){
            # there is already a file and we want to extend/update it
            library(jsonlite)
            # if the data data.table has a "visible" column then we will use it to determine which tracks will be visible on load
            datafiles = c(jsonlite::read_json(dfile), list(item))
        } else {
            datafiles = list(item)
        }
        jsonlite::write_json(datafiles, dfile,
                             pretty=TRUE, auto_unbox=TRUE, digits=4)
    }
    return(dfile)
}

get_js_datafiles_path = function(outdir, js.type){
    is.acceptable.js.type(js.type)
    if (js.type == 'gGnome.js'){
        dfile = paste0(outdir, '/datafiles.csv')
    }
    if (js.type == 'PGV'){
        dfile = paste0(outdir, '/public/datafiles.json')
    }
    return(dfile)
}

gen_gg_json_files = function(data, outdir, meta.js, name.col = 'pair', gg.col = 'complex', js.type = 'gGnome.js',
                             dataset_name = NA, ref.name = NULL, overwrite = FALSE, annotation = NULL){
    json_dir = get_gg_json_dir_path(outdir, js.type, dataset_name)
    json_files = lapply(1:data[, .N], function(idx){
        gg.js = get_gg_json_path(data[idx, get(name.col)], json_dir)
        if (!file.exists(gg.js) | overwrite){
            # TODO: at some point we need to do a sanity check to see that a valid rds of gGraph was provided
            gg = readRDS(data[idx, get(gg.col)])
            sl = parse.js.seqlenghts(meta.js, js.type = js.type, ref.name = ref.name)
            gg.js = refresh(gg[seqnames %in% names(sl)])$json(filename = gg.js,
                        verbose = TRUE,
                        annotation = annotation)
        } else {
            message(gg.js, ' found. Will not overwrite it.')
        }
        return(normalizePath(gg.js))
    })
    return(unlist(json_files))
}

is_git_lfs_available = function(raise = TRUE){
    if (length(readLines(pipe('command -v git-lfs'))) == 0){
        if (raise){
            stop('git-lfs is not installed, please install git-lfs (https://git-lfs.github.com/)')
        }
        return(FALSE)
    }
    return(TRUE)
}

is.dir.a.PGV.instance = function(outdir){
    if (!file.exists(paste0(outdir, '/public/settings.json'))){
        stop(outdir, ' does not seem to be a proper clone of the PGV github repository.')
    }
}

is.dir.a.gGnome.js.instance = function(outdir){
    if (!file.exists(paste0(outdir, '/public/metadata.json'))){
        stop(outdir, ' does not seem to be a proper clone of the gGnome.js github repository.')
    }
}

is.dir.a.js.instance = function(outdir, js.type){
    if (js.type == 'gGnome.js'){
        is.dir.a.gGnome.js.instance(outdir)
    }
    if (js.type == 'PGV'){
        is.dir.a.PGV.instance(outdir)
    }
}

js_path = function(outdir, append = FALSE, js.type = 'gGnome.js'){

    is.acceptable.js.type(js.type)
    outdir = suppressWarnings(normalizePath(outdir))

    if (file.exists(outdir) & !dir.exists(outdir)){
        stop('The output directory must be a valid path for a diretory, but you provided a path of a file that already exists:', outdir)
    }

    if (dir.exists(outdir) & !append){
        # if the folder exists and there is no append flag then throw error
        stop('The output directory already exists. If you wish to generate a new gGnome.js isntance, please provide a path for a new directory. If you wish to add more file to an existing instance of gGnome.js then use --append.')
    }

    if (!append){
        # clone the repository from github
        message('Cloning the ', js.type, ' repository from github.')
        if (js.type == 'gGnome.js'){
            system(paste0('git clone https://github.com/mskilab/gGnome.js.git ', outdir))
        } else {
            system(paste0('git clone https://github.com/mskilab/PGV.git ', outdir))
        }
    }
    # normalize path one more time to make sure we return the absolute path
    outdir = suppressWarnings(normalizePath(outdir))
    is.dir.a.js.instance(outdir, js.type)
    return(outdir)
}

acceptable.js.types = c('PGV', 'gGnome.js')
is.acceptable.js.type = function(js.type){
    if (!(js.type %in% acceptable.js.types)){
        stop('Invalid js.type. The only js types familiar to us are: ', acceptable.js.types)
    }
}

gen_js_coverage_files = function(data, outdir, name.col = 'pair', overwrite = FALSE, cov.col = 'cov',
                                 js.type = 'gGnome.js', cov.field = 'ratio', cov.field.col = NA,
                                 bin.width = 1e4, dataset_name = NA, ref.name = 'hg19',
                                 color.field = NULL, meta.js = NULL, mc.cores = 1){
    if (!is.na(cov.field.col)){
        if (!(cov.field.col %in% names(data))){
            stop(paste0('You provided the following invalid column name for cov.field.col: ', cov.field.col))
        }
        message('cov.field.col provided: "', cov.field.col, '". Will be reading coverage field name from this column.')
    }

    if (!(cov.col %in% names(data))){
        stop('Invalid cov.col. There is no column "', cov.col, '" in your data.')
    }

    cov_dir = get_js_cov_dir_path(outdir, js.type, dataset_name)
    cov_files = mclapply(1:data[, .N], function(idx){
        skip_cov = FALSE
        covfn = get_js_cov_path(data[idx, get(name.col)], cov_dir, js.type)
        if (!is.na(cov.field.col)){
            cov.field = data[idx, get(cov.field.col)] 
            message('cov.field: ', cov.field)
        }
        if (!file.exists(covfn) | overwrite){
            if (is.na(cov.field)){
                warning(paste0('No coverage field was provided for ', data[idx, get(name.col)], ' so no coverage will be generated.'))
                skip_cov = TRUE
            } else {
                if (is.na(cov.col)){
                    warning(paste0('No coverage data was provided for ', data[idx, get(name.col)], ' so no coverage will be generated.'))
                    skip_cov = TRUE
                } else {
                    cov_input_file = data[idx, get(cov.col)]
                    if (is.na(cov_input_file)){
                        warning(paste0('No coverage file was provided for ', data[idx, get(name.col)], ' so no coverage will be generated.'))
                        skip_cov = TRUE
                    } else {
                        if (!file.exists(cov_input_file)){
                            warning(paste0('No coverage file was provided for ', data[idx, get(name.col)], ' so no coverage will be generated.'))
                            skip_cov = TRUE
            }}}}
            if (skip_cov){
                return(NA)
            } else {
                if (js.type == 'gGnome.js'){
                    cov2csv(cov_input_file, field = cov.field,
                            output_file = covfn)
                } else {
                    cov2arrow(cov_input_file, field = cov.field,
                              output_file = covfn, ref.name = ref.name,
                              color.field = color.field, overwrite = overwrite,
                              meta.js = meta.js, bin.width = bin.width)
                }
            }
        } else {
            message(covfn, ' found. Will not overwrite it.')
        }
        return(normalizePath(covfn))
    }, mc.cores = mc.cores)
    return(unlist(cov_files))
}

read.js.input.data = function(data, name.col = 'pair'){
    if (!inherits(data, 'data.table')){
        if (!is.character(data)){
            stop('Invalid input data of class: "', class(data), '". Expected data.table or path to CSV/TSV file.')
        }
        if (!file.exists(data)){
            stop('Invalid input data. The input data must be a data.table or a path to a CSV/TSV file')
        }
        data = fread(data)
    }
    if (!(name.col %in% names(data))){
        stop('Invalid name.col provided: "', name.col, '".')
    }
    if (any(duplicated(data[, get(name.col)]))){
        stop('The name.col must hold non-redundant values, but the name.col you provided has duplicates. Here is an example for a value with duplicates: ', data[duplicated(get(name.col)), get(name.col)])
    }
    return(data)
}

get_gg_json_path = function(nm, gg_json_dir){
    gg.js = paste0(gg_json_dir, "/", nm, ".json")
    return(gg.js)
}

get_gg_json_dir_path = function(outdir, js.type, dataset_name = NA){
    if (js.type == 'gGnome.js'){
        gg_json_dir = paste0(outdir, '/json')
    } else {
        if (js.type == 'PGV'){
            gg_json_dir = get_pgv_data_dir(outdir, dataset_name = dataset_name)
        }
    }
    return(gg_json_dir)
}

get_js_cov_path = function(nm, cov_dir, js.type){
    if (js.type == 'gGnome.js'){
        covfn = paste0(cov_dir, "/", nm, ".csv")
        return(covfn)
    } else {
        if (js.type == 'PGV'){
            covfn = paste0(cov_dir, '/', nm, '-coverage.arrow')
        }
    }
    return(covfn)
}

get_js_cov_dir_path = function(outdir, js.type, dataset_name = NA){
    if (js.type == 'gGnome.js'){
        cov_dir = paste0(outdir, '/scatterPlot')
    } else {
        if (js.type == 'PGV'){
            cov_dir = get_pgv_data_dir(outdir, dataset_name)
        }
    }
    return(cov_dir)
}

get_pgv_data_dir = function(outdir, dataset_name = NA){
    if (is.na(dataset_name)){
        stop('dataset_name must be provided for PGV.')
    }
    data_dir = paste0(outdir, '/public/data/', dataset_name, '/')
    # make sure the directory exists
    if (!dir.exists(data_dir)){
        message('Creating a directory for the PGV data files here: ', data_dir)
        dir.create(data_dir, recursive = TRUE)
    }
    return(data_dir)
}

#' @export
#' @name cov2cov.js
#' @description
#' Takes a GRanges with coverage data and converts it to a data.table with the info needed for gGnome.js and PGV
#'
#' if bin.width is specified then coverage data will also be rebinned
#' if convert.to.cn == TRUE then rel2abs will be applied
#' @param cov coverage GRanges or path to file with coverage data
#' @param meta.js
cov2cov.js = function(cov, meta.js = NULL, js.type = 'gGnome.js', field = 'ratio',
                      bin.width = NA, ref.name = NULL, color.field = NULL){
    ## respect the seqlengths in meta.js
    if (is.character(meta.js) && file.exists(meta.js)){
        sl = parse.js.seqlenghts(meta.js, js.type = js.type, ref.name = ref.name)
    }

    x = readCov(cov)

    if (!exists("sl")){
        sl = seqlengths(x)
    }
    if (all(is.na(sl))){
        stop("No seqlengths in the input.")
    }

    if (!is.element(field, names(mcols(x)))){
        stop("The provided field '", field, "' is not in the input data")
    }

    fields = field

    if (!is.null(color.field)){
        if (!is.element(color.field, names(mcols(x)))){
            stop("The color.field '", color.field, "' is not in the input data")
        }
        fields = c(field, color.field)
    }

    if (!is.na(bin.width)){
        message('Rebinning coverage with bin.width=', bin.width)
        if (!is.numeric(bin.width)){
            stop('bin.width must be numeric')
        }
        x.rebin = rebin(x, bin.width, field, FUN = median)
        if (!is.null(color.field)){
            # we will take the median value of numeric values and the a single value (by majority votes) for any other type
            # this is intended so that we can keep the colors if they existed
            my_cool_fn = function(value, width, na.rm){
                ifelse(is.numeric(value), median(value),
                                   names(sort(table(value), decreasing = T)[1]))
            }
            x = gr.val(x.rebin, x,
                       val = fields, FUN = my_cool_fn)
        } else {
            x = x.rebin
        }
        message('Done rebinning')
    }

    # TODO: add option to do rel2abs

    ## build the cumulative coordinates
    dt = data.table(seqlevels = names(sl),
                    seqlengths = as.double(sl),
                    cstart = c(1, 1 + cumsum(as.double(sl))[-length(sl)]))

    overlap.seqnames = intersect(seqlevels(x), names(sl))
    if (length(overlap.seqnames) == 0){
        stop('The names of sequences in the input coverage and in the reference don\'t match. This is an example seqname from the ref: "',
             names(sl)[1],
        '". And here is an example from the coverage file: "', seqnames(x)[1], '".')
    }
    # make sure that seqnames overlap between coverage and reference
    invalid.seqnames = setdiff(seqnames(x), names(sl))
    if (length(invalid.seqnames) > 0){
        warning(sprintf('The coverage input includes sequence names that are not in the specified reference, and hence these sequence names will be excluded. These are the excluded sequences: %s', paste(as.character(invalid.seqnames), collapse = ', ')))
    }

    ## build the data.table
    dat = as.data.table(merge(gr2dt(x), dt,
                by.x = "seqnames",
                by.y = "seqlevels",
                all.x = TRUE))
    dat = dat[seqnames %in% overlap.seqnames] # only keep seqnames that are in the reference
    dat[, new.start := start + cstart - 1]
    # convert Inf to NA
    dat[get(field) == Inf, (field) := NA]

    return(dat)
}

#' @name cov2csv
#' @description
#' prepare csv file for gGnome.js
#' Col1 x
#' Col2 y
#' Col3 cumulative x
cov2csv = function(cov,
        field = "ratio",
        output_file = "coverage.csv",
        ...)
{

    cov = readCov(cov)

    if (!is.element(field, names(mcols((cov))))){
        stop('"', field, '" is not in the input data')
    }

    dat = cov2cov.js(cov, ...)

    outdir = dirname(output_file)
    ## make sure the path to the file exists
    if (!dir.exists(outdir)){
        dir.create(outdir)
    }

    cov.out = dat[, .(x = start, y = get(field), chromosome = as.character(seqnames))][!is.na(y)]
    message('Writing coverage data to file: ', output_file)
    fwrite(cov.out, output_file,  sep = ",")
    return(output_file)
}


#' @name cov2arrow
#' @description
#'
#' Prepares an scatter plot arrow file with coverage info for PGV (https://github.com/mskilab/pgv)
#'
#' @param cov input coverage data (GRanges)
#' @param field which field of the input data to use for the Y axis
#' @param output_file output file path.
#' @param ref.name the name of the reference to use. If not provided, then the default reference that is defined in the meta.js file will be loaded.
#' @param color.field a field in the input GRanges object to use to determine the color of each point
#' @param overwrite (logical) by default, if the output path already exists, it will not be overwritten.
#' @param meta.js path to JSON file with metadata for PGV (should be located in "public/settings.json" inside the repository)
#' @param bin.width (integer) bin width for rebinning the coverage (default: 1e4)
#' @author Alon Shaiber
#' @export
cov2arrow = function(cov,
        field = "ratio",
        output_file = 'coverage.arrow',
        ref.name = 'hg19',
        color.field = NULL,
        overwrite = FALSE,
        meta.js = NULL,
        ...){

    outdir = dirname(output_file)
    dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

    if (!file.exists(output_file) | overwrite){
        if (!requireNamespace("arrow", quietly = TRUE)) {
            stop('You must have the package "arrow" installed in order for this function to work. Please install it.')
        }



        message('Converting coverage format')
        dat = cov2cov.js(cov, meta.js = meta.js, js.type = 'PGV', field = field,
                         ref.name = ref.name, color.field = color.field, ...)
        message('Done converting coverage format')


        if (!is.null(color.field)){
            dat[, color := color2numeric(get(color.field))]
        } else {
            if (!is.null(meta.js)){
                ref_meta = get_ref_metadata_from_PGV_json(meta.js, ref.name)
                setkey(ref_meta, 'chromosome')
                dat$color = color2numeric(ref_meta[dat$seqnames]$color)
            } else {
                # no color field and no meta.js so set all colors to black
                dat$color = 0
            }
        }

        outdt = dat[, .(x = new.start, y = get(field), color)]

        # if there are any NAs for colors then set those to black
        outdt[is.na(color), color := 0]

        # remove NAs
        outdt = outdt[!is.na(y)]

        message('Writing arrow file (using write_feather)')
        arrow_table = arrow::Table$create(outdt, schema = arrow::schema(x = arrow::float32(), y = arrow::float32(), color = arrow::float32()))
        arrow::write_feather(arrow_table, output_file)
    } else {
        message('arrow file, "', output_file, '" already exists.')
    }
    return(output_file)
}

#' @name color2hex
#' @description
#'
#' Takes a vector of colors and returns the hex color code for the colors
#'
#' Any color that could be parsed by col2rgb is acceptable. Missing or invalid values are assigned a default color. The color names could be a mix of hex color codes and names (e.g. "black")
#'
#' @param x vector of colors names
#' @param default_color the color to default to for NAs and invalid values
#' @return vector of hex color codes
#' @author Alon Shaiber
color2hex = function(x, default_color = '#000000'){
    cols = lapply(x, function(y){
       tryCatch(col2rgb(y),
                error = function(e) 'default')
    })
    out = sapply(cols, function(y){
       tryCatch(rgb(y[1], y[2], y[3], maxColorValue = 255),
                error = function(e) 'default')
    })

    default_pos = sum(out == 'default')
    if (default_pos > 0){
        warning(sprintf('There were %s entries with missing or invalid colors and these were set to the default color: %s', default_pos, default_color))
        out[which(out == 'default')] = default_color
    }
    return(out)
}

#' @name colorhex2numeric
#' @description
#'
#' Takes a vector of colors hex code and returns a vector with integers corresponding to each hex number
#'
#' @param x vector of colors hex codes
#' @return numeric vector
#' @author Alon Shaiber
colorhex2numeric = function(x){
     return(strtoi(gsub('\\#', '0x', x)))
}

#' @name color2numeric
#' @description
#'
#' Takes a vector of colors and returns a numeric vector as expected by PGV
#'
#' Any color that could be parsed by col2rgb is acceptable. Missing or invalid values are assigned a default color. The color names could be a mix of hex color codes and names (e.g. "black")
#'
#' @param x vector of colors names
#' @return numeric vector
#' @author Alon Shaiber
color2numeric = function(x, default_color = '#000000'){
     return(colorhex2numeric(color2hex(x, default_color = default_color)))
}


#####################################################
#' @name gtf2json
#' @description Turning a GTF format gene annotation into JSON
#'
#' @param gtf path to GTF input file.
#' @param gtf.rds path to rds file which includes a data.table holding the GTF information.
#' @param gtf.gr.rds path to rds file which includes a GRanges holding the GTF information.
#' @param metadata.filename metadata JSON output file name (./metadata.json).
#' @param genes.filename genes JSON output file name (./genes.json).
#' @param genes
#' @param gene_weights table with weights for genes. The first column is the gene name and the second column is the numeric weight. Either a data.frame, data.table, or a path to a file need to be provided. Genes with weights above 10 would be prioritized to show when zoomed out in gGnome.js.
#' @param grep
#' @param grepe
#' @param chrom.sizes if not provided then the default hg19 chromosome lengths will be used.
#' @param include.chr chromosomes to include in the output. If not provided then all chromosomes in the reference are included.
#' @param gene.collapse
#' @param verbose
#' @author Xiaotong Yao, Alon Shaiber
#' @return file_list list containing the paths of the metadata and genes JSON-formatted output files.
#' @export
####################################################
gtf2json = function(gtf=NULL,
                    gtf.rds=NULL,
                    gtf.gr.rds=NULL,
                    metadata.filename="./metadata.json",
                    genes.filename="./genes.json",
                    genes=NULL,
                    gene_weights=NULL,
                    grep=NULL,
                    grepe=NULL,
                    chrom.sizes=NULL,
                    include.chr=NULL,
                    gene.collapse=TRUE,
                    verbose = TRUE){
    require(data.table)
    require(gUtils)
    require(rtracklayer)

    if (!is.null(gtf.gr.rds)){
        message("Using GRanges from rds file.")
        infile = gtf.gr.rds
        gr = readRDS(gtf.gr.rds)
        dt = gr2dt(gr)
    } else if (!is.null(gtf.rds)){
        message("Using GTF data.table from rds file.")
        infile = gtf.rds
        dt = as.data.table(readRDS(gtf.rds))
    } else if (!is.null(gtf)){
        message("Using raw GTF file.")
        infile = gtf

        gr = rtracklayer::import.gff(gtf)
        dt = gr2dt(gr)

    } else {
        warning("No input gene annotation. Use the built-in GENCODE v19 in gUtils package")
        require(skidb)
        gr = read_gencode()
        infile = "default"
        dt = gr2dt(gr)
    }

    if (verbose){
        message("Finished reading raw data, start processing.")
    }

    ## get seqlengths
    if (is.null(chrom.sizes)){
        message("No ref genome seqlengths given, use default.")
        ## chrom.sizes = system.file("extdata", "hg19.regularChr.chrom.sizes", package="gGnome")
        ## system.file("extdata", "hg19.regularChr.chrom.sizes", package="gGnome")
        Sys.setenv(DEFAULT_BSGENOME=system.file("extdata", "hg19.regularChr.chrom.sizes", package="gUtils"))
    } else {
                Sys.setenv(DEFAULT_BSGENOME=chrom.sizes)
    }

    sl = hg_seqlengths(include.junk=TRUE)

    if (!is.null(include.chr)){
        sl = sl[include.chr]
    }
    chrs = data.table(seqnames = names(sl), seqlengths=sl)

    ## meta data field
    require(RColorBrewer)
    qw = function(x) paste0('"', x, '"') ## quote

    meta.json =paste0(paste0("{", qw("metadata"),': [\n'),
                     chrs[, paste("\t\t{",
                                  qw("chromosome"),": ", qw(seqnames),
                                  ", ", qw("startPoint"),": ", 1,
                                  ", ", qw("endPoint"), ": ", seqlengths,
                                  ", ", qw("color"),
                                  ": ", qw(substr(tolower(brewer.master( max(.I), 'BrBG' )), 1, 7)), " }",
                                  collapse=",\n",
                                  sep="")],
                     '\n  ],\n',
                     paste(
                     paste0(qw("sequences"), ": {", qw("T"), ": ", qw("#E6E431"), ", ", qw("A"), ": ", qw("#5157FB"), ", ", qw("G"), ": ", qw("#1DBE21"), ", ", qw("C"), ": ",qw("#DE0A17"), ", ", qw("backbone"), ": ", qw("#AD26FA"), "}"),
                     paste0(qw("coveragePointsThreshold"), ":  30000"),
                     paste0(qw("scatterPlot"), ": {", qw("title"), ": ", qw("Coverage"), "}"),
                     paste0(qw("barPlot"), ":  {", qw("title"), ":  ", qw("RPKM"), "}"),
                     paste0(qw("intervalsPanelHeightRatio"), ": 0.6"),
                     sep = ",\n")
                    )


    if (verbose){
        message("Metadata fields done.")
    }

    ## reduce columns: seqnames, start, end, strand, type, gene_id, gene_name, gene_type, transcript_id
    ## reduce rows: gene_status, "KNOWN"; gene_type, not "pseudo", not "processed transcript"
    dtr = dt
    if ('gene_status' %in% names(dt)){
        dtr = dt[gene_status=="KNOWN"]
    }
    dtr = dtr[!grepl("pseudo", gene_type) &
             gene_type != "processed_transcript",
             .(chromosome=seqnames, startPoint=start, endPoint=end, strand,
               title = gene_name, gene_name, type, gene_id, gene_type,
               transcript_id, transcript_name)]

    if(!is.null(include.chr)){
           dtr = dtr[chromosome %in% include.chr]
    }
    if (!is.null(genes)){
        dtr = dtr[title %in% genes]
    } else {
            if (!is.null(grep) | !is.null(grepe)) {
                if (!is.null(grep)){
                dtr = dtr[grepl(grep, title)]
            }
            if (!is.null(grepe)){
                dtr = dtr[!grepl(grepe, title)]
            }
        }
    }

    if (nrow(dtr)==0){
        stop("Error: No more data to present.")
    }

    if (gene.collapse){
        ## collapse by gene
        dtr[, hasCds := is.element("CDS", type), by=gene_id]
        dtr = rbind(dtr[hasCds==TRUE][type %in% c("CDS","UTR","gene")],
                    dtr[hasCds==FALSE][type %in% c("exon", "gene")])
        ## dedup
        dtr = dtr[!duplicated(paste(chromosome, startPoint, endPoint, gene_id))]
        dtr[, title := gene_name]
        dtr = dtr[type != "transcript"]

        ## group id
        dtr[, gid := as.numeric(as.factor(gene_id))]
        if (verbose){
            message("Intervals collapsed to gene level.")
        }
    } else {
        ## collapse by transcript
        dtr[, hasCds := is.element("CDS", type), by=transcript_id]
        dtr = rbind(dtr[hasCds==TRUE][type %in% c("CDS","UTR","transcript")],
                    dtr[hasCds==FALSE][type %in% c("exon","transcript")])
        ## dedup
        dtr = dtr[!duplicated(paste(chromosome, startPoint, endPoint, transcript_id))]
        dtr[, title := transcript_name]
        dtr = dtr[type != "gene"]

        ## group id
        dtr[, gid := as.numeric(as.factor(transcript_id))]
        if (verbose){
            message("Intervals collapsed to transcript level.")
        }
    }

    dtr[, iid := 1:nrow(dtr)]

    #' incorporate gene_weights 
    if (!is.null(gene_weights)){
        # check that this is a data.table
        if (inherits(gene_weights, 'data.frame')){
            gene_weights = gene_weights %>% as.data.table
        } else {
            if (!file.exists(gene_weights)){
                stop('Gene weights must be provided either as a dataframe or as a path to a file with a tabular text format (with no header).')
            }
            gene_weights = fread(gene_weights, header = FALSE)
        }

        if (dim(gene_weights)[2] != 2){
            stop('gene_weights must be a table with just two columns.')
        }
        setnames(gene_weights, names(gene_weights), c('gene_name', 'weight'))

        # make sure that weights are numeric
        gene_weights[, weight := as.numeric(weight)]

        if (gene_weights[is.na(weight), .N] > 0){
            print('Some weights provided in gene_weights are either not-valid or missing and would be set to the default value (1).')
        }

        # check names of genes
        genes_in_gene_weights_but_not_in_dtr = setdiff(gene_weights$gene_name, dtr$gene_name)
        if (length(genes_in_gene_weights_but_not_in_dtr) > 0){
            print(sprintf('Warning: the following gene names appear in the provided gene_weights, but do not match any of the genes in the reference genome (and hence will be ignored): %s', genes_in_gene_weights_but_not_in_dtr))
        }
        dtr = merge(dtr, gene_weights, by = 'gene_name', all.x = TRUE)

        #' set all missing weights to 1
        dtr[is.na(weight), weight := 1]
    }

    ## processing genes
    genes.json = dtr[, paste0(
        c(paste0('{', qw("genes"),": ["),
          paste(
              "\t{",
              qw("iid"), ": ", iid,
              ", ", qw("chromosome"), ": ", qw(chromosome),
              ", ", qw("startPoint"), ": ", startPoint,
              ", ", qw("endPoint"), ": ", endPoint,
              ", ", qw("y"), ": ", 0,
              ", ", qw("title"), ": ", qw(title),
              ", ", qw("group_id"), ": ", qw(gid),
              ", ", qw("type"), ": ", qw(type),
              ", ", qw("strand"), ": ", qw(strand),
              ", ", qw("weight"), ": ", weight,
              "}",
              sep = "",
              collapse = ',\n'),
          "]"),
        collapse = '\n')
        ]


    ## assembling the JSON
    out_meta = paste(c(meta.json, "}"),
                     sep = "")

    writeLines(out_meta, metadata.filename)
    message(sprintf('Wrote JSON file of %s to %s', infile, metadata.filename))


    out_genes = paste(c(genes.json, "}"),
                     sep = "")
    writeLines(out_genes, genes.filename)
    message(sprintf('Wrote JSON file of %s to %s', infile, genes.filename))

    return(list(metadata.filename = metadata.filename, genes.filename = genes.filename))
}


#' @name jab2json
#' @description a wrapper function to dump JaBbA results run with Flow to gGnome.js viz
#' @export
jab2json = function(fn = "./jabba.simple.rds",
                    gGnome.js.dir = "~/git/gGnome.js"){

}


#' @name parse.js.seqlenghts
#' @description
#' Takes a settings JSON file from either gGnome.js or PGV and parses it into a data.table
#' @param meta.js input settings JSON file
#' @param js.type either 'gGnome.js' or 'PGV' to determine the format of the JSON file
#' @param ref.name the name of the reference to load (only relevant for PGV). If not provided, then the default reference will be loaded.
#' @author Alon Shaiber
parse.js.seqlenghts = function(meta.js, js.type = 'gGnome.js', ref.name = NULL){
    if (!(js.type %in% c('gGnome.js', 'PGV'))){
        stop('js.type must be either gGnome.js or PGV')
    }
     if (js.type == 'gGnome.js'){
        message('Getting seqlenghts for gGnome.js')
        settings = jsonlite::read_json(meta.js)
        if (!('metadata' %in% names(settings))){
            stop('Input JSON file is not a valid gGnome.js settings JSON. Please check the required format.')
        }
        meta = rbindlist(settings$metadata)
        sl = meta[, setNames(endPoint, chromosome)]
     } else {
        message('Getting seqlenghts for PGV')
        ref_meta = get_ref_metadata_from_PGV_json(meta.js, ref.name)
        sl = ref_meta[, setNames(endPoint, chromosome)]
    }
    return(sl)
}

#' @export
#' @name get_ref_metadata_from_PGV_json
#' @description
#' get a data.table with the metadata for the reference (columns: chromosome, startPoint, endPoint, color)
#' @param ref.name
#' @param meta.js
get_ref_metadata_from_PGV_json = function(meta.js, ref.name = NULL){
    if (!is.character(ref.name)){
        stop('Invalid ref.name: ', ref.name)
    }
    meta = jsonlite::read_json(meta.js)
    if (!('coordinates' %in% names(meta))){
        stop('Input meta file is not a proper PGV settings.json format.')
    }
    coord = meta$coordinates
    if (is.null(ref.name)){
        # if no ref was provided then take the default
        if (!('default' %in% names(coord))){
            stop('Invalid meta file. The meta file, ', meta.js, ', is missing a default coordinates value.')
        }
        ref.name = coord$default
        message('No reference name provided so using default: "', ref.name, '".')
    }
    if (!('sets' %in% names(coord))){
        stop('Input meta file is not a proper PGV settings.json format.')
    }
    sets = coord$sets
    if (!(ref.name %in% names(sets))){
        stop('Invalid ref.name: ', ref.name, '. The ref.name does not appear to be described in ', meta.js)
    }
    seq.info = rbindlist(sets[[ref.name]])
    return(seq.info)
}

get_path_to_meta_js = function(outdir, js.type = js.type){
    is.acceptable.js.type(js.type)
    if (!dir.exists(outdir)){
        stop('No such directory: ', outdir)
    }
    if (js.type == 'gGnome.js'){
        meta.js = suppressWarnings(normalizePath(paste0(outdir, '/public/metadata.json')))
    }
    if (js.type == 'PGV'){
        meta.js = suppressWarnings(normalizePath(paste0(outdir, '/public/settings.json')))
    }
    if (!file.exists(meta.js)){
        stop('We could not find a metadata file where we expected it to be ("',
             meta.js,
             '"). Something must have gone wrong, please check that this is a valid clone of the ',
             js.type, ' github repository.')
    }
    return(meta.js)
}

