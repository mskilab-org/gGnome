#' @name pgv
#' @rdname gen_js_instance
#' @description
#'
#' Generate a PGV instance
#'
#' Takes a table with paths to gGraphs and coverage files (optional) and generates an instance of a gGnome.js directory that is ready to visualize using gGnome.js
#' 
#' @param data either a path to a TSV/CSV or a data.table
#' @param dataset_name the name of the dataset. This should be the name of the project that all the samples belong to. You must provide a name since PGV stores all the data under a folder matching your dataset name. This allows a single PGV instance to include multiple datasets which could be browsed by going to the "Data Selection" page in the browser
#' @param name.col column name in the input data table containing the sample names (default: "sample")
#' @param outdir the path where to save the files. This path should not exist, unless you want to add more files to an existing directory in which case you must use append = TRUE
#' @param cov.col column name in the input data table containing the paths to coverage files
#' @param gg.col column name in the input data table containing the paths to RDS files containing the gGnome objects
#' @param append use this flag if you already have an existing PGV folder that you want to add more files to
#' @param cov.field the name of the field in the coverage GRanges that should be used (default: "ratio")
#' @param cov.field.col column name in the input data table containing the name of the field in the coverage GRanges that should be used. If this is supplied then it overrides the value in "cov.field". Use this if some of your coverage files differ in the field used.
#' @param cov.bin.width bin width to use when rebinning the coverage data (default: 1e4). If you don't want rebinning to be performed then set to NA.
#' @param cov.color.field field in the coverage GRanges to use in order to set the color of coverage data points. If nothing is supplied then default colors are used for each seqname (namely chromosome) by reading the colors that are defined in the settings.json file for the specific reference that is being used for this dataset.
#' @param ref.name the genome reference name used for this dataset. This reference name must be defined in the settings.json file. By default PGV accepts one of the following: hg19, hg38, covid19. If you are using a different reference then you must first add it to the settings.json file.
#' @param overwrite by default only files that are missing will be created. If set to TRUE then existing coverage arrow files and gGraph JSON files will be overwritten
#' @param annotation which node/edge annotation fields to add to the gGraph JSON file. By default we assume that gGnome::events has been executed and we add the following SV annotations: 'simple', 'bfb', 'chromoplexy', 'chromothripsis', 'del', 'dm', 'dup', 'pyrgo', 'qrdel', 'qrdup', 'qrp', 'rigma', 'tic', 'tyfonas'
#' @param tree.path path to newick file containing a tree to incorporate with the dataset. IF provided then the tree is added to datafiles.json and will be visualized by PGV. If the names of leaves of the tree match the names defined in the name.col then PGV will automatically assocaited these leaves with the samples and hence upon clicking a leaf of the tree the browser will scroll down to the corresponding genome graph track
#' @param mc.cores how many cores to use
#' 
#' @export
pgv = function(data,
                    dataset_name = NA,
                    name.col = 'sample',
                    outdir = './gGnome.js',
                    cov.col = 'coverage',
                    gg.col = 'graph',
                    append = FALSE,
                    cov.field = 'ratio',
                    cov.field.col = NA,
                    cov.bin.width = 1e4,
                    cov.color.field = NULL,
                    ref.name = NA,
                    overwrite = FALSE,
                    annotation = c('simple', 'bfb', 'chromoplexy',
                                'chromothripsis', 'del', 'dm', 'dup',
                                'pyrgo', 'qrdel', 'qrdup', 'qrp', 'rigma',
                                'tic', 'tyfonas'),
                    tree.path = NA,
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
                           cov.color.field = cov.color.field,
                           dataset_name = dataset_name,
                           ref.name = ref.name,
                           overwrite = overwrite,
                           annotation = annotation,
                           tree.path = tree.path,
                           mc.cores = mc.cores))
}

#' @name gGnome.js
#' @rdname gen_js_instance
#' @description
#'
#' Generate a gGnome.js instance
#'
#' Takes a table with paths to gGraphs and coverage files (optional) and generates an instance of a gGnome.js directory that is ready to visualize using gGnome.js
#' 
#' @param data either a path to a TSV/CSV or a data.table
#' @param name.col column name in the input data table containing the sample names (default: "sample")
#' @param outdir the path where to save the files. This path should not exist, unless you want to add more files to an existing directory in which case you must use append = TRUE
#' @param cov.col column name in the input data table containing the paths to coverage files
#' @param gg.col column name in the input data table containing the paths to RDS files containing the gGnome objects
#' @param append use this flag if you already have an existing PGV folder that you want to add more files to
#' @param cov.field the name of the field in the coverage GRanges that should be used (default: "ratio")
#' @param cov.field.col column name in the input data table containing the name of the field in the coverage GRanges that should be used. If this is supplied then it overrides the value in "cov.field". Use this if some of your coverage files differ in the field used.
#' @param cov.bin.width bin width to use when rebinning the coverage data (default: 1e4). If you don't want rebinning to be performed then set to NA.
#' @param ref.name the genome reference name used for this dataset. This reference name must be defined in the settings.json file. By default PGV accepts one of the following: hg19, hg38, covid19. If you are using a different reference then you must first add it to the settings.json file.
#' @param overwrite by default only files that are missing will be created. If set to TRUE then existing coverage arrow files and gGraph JSON files will be overwritten
#' @param annotation which node/edge annotation fields to add to the gGraph JSON file. By default we assume that gGnome::events has been executed and we add the following SV annotations: 'simple', 'bfb', 'chromoplexy', 'chromothripsis', 'del', 'dm', 'dup', 'pyrgo', 'qrdel', 'qrdup', 'qrp', 'rigma', 'tic', 'tyfonas'
#' @param mc.cores how many cores to use
#' 
#' @export
gGnome.js = function(data,
                           name.col = 'sample',
                           outdir = './gGnome.js',
                           cov.col = 'coverage',
                           gg.col = 'graph',
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
#' @description
#'
#' Generate a gGnome.js instance
#'
#' Takes a table with paths to gGraphs and coverage files (optional) and generates an instance of a gGnome.js directory that is ready to visualize using gGnome.js
#' 
#' @param data either a path to a TSV/CSV or a data.table
#' @param name.col column name in the input data table containing the sample names (default: "sample")
#' @param outdir the path where to save the files. This path should not exist, unless you want to add more files to an existing directory in which case you must use append = TRUE
#' @param cov.col column name in the input data table containing the paths to coverage files
#' @param gg.col column name in the input data table containing the paths to RDS files containing the gGnome objects
#' @param append use this flag if you already have an existing folder that you want to add more files to
#' @param js.type either "PGV" or "gGnome.js"
#' @param cov.field the name of the field in the coverage GRanges that should be used (default: "ratio")
#' @param cov.field.col column name in the input data table containing the name of the field in the coverage GRanges that should be used. If this is supplied then it overrides the value in "cov.field". Use this if some of your coverage files differ in the field used.
#' @param cov.bin.width bin width to use when rebinning the coverage data (default: 1e4). If you don't want rebinning to be performed then set to NA.
#' @param cov.color.field field in the coverage GRanges to use in order to set the color of coverage data points. If nothing is supplied then default colors are used for each seqname (namely chromosome) by reading the colors that are defined in the settings.json file for the specific reference that is being used for this dataset.
#' @param dataset_name the name of the dataset. Only relevant for PGV. This should be the name of the project that all the samples belong to. You must provide a name since PGV stores all the data under a folder matching your dataset name. This allows a single PGV instance to include multiple datasets which could be browsed by going to the "Data Selection" page in the browser
#' @param ref.name the genome reference name used for this dataset. For specific behaviour refer to the PGV/gGnome.js wrappers
#' @param overwrite by default only files that are missing will be created. If set to TRUE then existing coverage arrow files and gGraph JSON files will be overwritten
#' @param annotation which node/edge annotation fields to add to the gGraph JSON file. By default we assume that gGnome::events has been executed and we add the following SV annotations: 'simple', 'bfb', 'chromoplexy', 'chromothripsis', 'del', 'dm', 'dup', 'pyrgo', 'qrdel', 'qrdup', 'qrp', 'rigma', 'tic', 'tyfonas'
#' @param tree.path path to newick file containing a tree to incorporate with the dataset. Only relevant for PGV. IF provided then the tree is added to datafiles.json and will be visualized by PGV. If the names of leaves of the tree match the names defined in the name.col then PGV will automatically assocaited these leaves with the samples and hence upon clicking a leaf of the tree the browser will scroll down to the corresponding genome graph track
#' @param mc.cores how many cores to use
#' 
#' @export
gen_js_instance = function(data,
                           name.col = 'sample',
                           outdir = './gGnome.js',
                           cov.col = 'coverage',
                           gg.col = 'graph',
                           append = FALSE,
                           js.type = 'gGnome.js',
                           cov.field = 'ratio',
                           cov.field.col = NA,
                           cov.bin.width = 1e4,
                           cov.color.field = NULL,
                           dataset_name = NA,
                           ref.name = NA,
                           overwrite = FALSE,
                           annotation = c('simple', 'bfb', 'chromoplexy',
                                       'chromothripsis', 'del', 'dm', 'dup',
                                       'pyrgo', 'qrdel', 'qrdup', 'qrp', 'rigma',
                                       'tic', 'tyfonas'),
                           tree.path = NA,
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
                          ref.name = ref.name, cov.color.field = cov.color.field,
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

#' @name gen_js_datafiles
#' @description internal
#'
#' Generate the datafiles object for a PGV or gGnome.js instance
#'
#' @param data either a path to a TSV/CSV or a data.table
#' @param outdir the path to the PGV/gGnome.js repository clone
#' @param js.type either "PGV" or "gGnome.js"
#' @param name.col column name in the input data table containing the sample names (default: "sample")
#' @param meta_col column in the input data table containing the description of each sample. A single string is expected in which each description term is separated by a semicolon and space ("; "). For example: "ATCC; 2014; Luciferase; PTEN-; ESR1-""
#' @param ref.name the genome reference name used for this dataset. Only relevant for PGV
#' @param dataset_name the name of the dataset. Only relevant for PGV. This should be the name of the project that all the samples belong to. You must provide a name since PGV stores all the data under a folder matching your dataset name. This allows a single PGV instance to include multiple datasets which could be browsed by going to the "Data Selection" page in the browser
#' 
#' @export
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
                     nm = data[idx, get(name.col)]
                     if (!is.na(gg.js)){
                         if (file.exists(gg.js)){
                             gg.track = list('sample' = nm,
                                             'type' = 'genome',
                                             'source' = paste0(nm, '.json'),
                                             'title' = nm,
                                             'visible' = ifelse(data[idx, visible] == TRUE, TRUE, FALSE))
                         }
                     }
                     if (!is.na(cov.fn)){
                         if (file.exists(cov.fn)){
                             cov.track = list('sample' = nm,
                                              'type' = 'scatterplot',
                                              'source' = paste0(nm, '-coverage.arrow'),
                                              'title' = paste0(nm, ' Coverage Distribution'),
                                              'visible' = FALSE) # coverage tracks will always be set to not visible on load
                         }
                     }
                     tracks = list(gg.track, cov.track)
                     return(tracks)
        })

        plots = do.call(c, plots)

        item = list()
        item$description = list(paste0('dataset=', dataset_name)) # TODO: we need to figure out the purpose of the description in PGV and update this accordingly
        item$reference = ref.name
        item$plots = plots

        if (file.exists(dfile)){
            # there is already a file and we want to extend/update it
            library(jsonlite)
            datafiles = jsonlite::read_json(dfile)
            if (dataset_name %in% names(datafiles)){
                warning('Notice that an entry for "', dataset_name, '" previously existed  in your datafiles.json and will now be override.')
            }
        } else {
            datafiles = list()
        }
        datafiles[[dataset_name]] = item
        jsonlite::write_json(datafiles, dfile,
                             pretty=TRUE, auto_unbox=TRUE, digits=4)
    }
    return(dfile)
}

#' @name get_js_datafiles_path
#' @description
#'
#' Get the path to the datafiles (CSV for gGnome.js, JSON for PGV) inside the clone of the repository
#'
#' @param outdir the path to the PGV/gGnome.js repository clone
#' @param js.type either "PGV" or "gGnome.js"
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

#' @name gen_gg_json_files
#' @description internal
#'
#' Generate a gGnome.js instance
#'
#' Takes a table with paths to gGraphs and coverage files (optional) and generates an instance of a gGnome.js directory that is ready to visualize using gGnome.js
#' 
#' @param data either a path to a TSV/CSV or a data.table
#' @param outdir the path to the PGV/gGnome.js repository clone
#' @param meta.js path to JSON file with metadata (for PGV should be located in "public/settings.json" inside the repository and for gGnome.js should be in public/genes/metadata.json)
#' @param name.col column name in the input data table containing the sample names (default: "sample")
#' @param gg.col column name in the input data table containing the paths to RDS files containing the gGnome objects
#' @param js.type either "PGV" or "gGnome.js"
#' @param dataset_name the name of the dataset. Only relevant for PGV. This should be the name of the project that all the samples belong to. You must provide a name since PGV stores all the data under a folder matching your dataset name. This allows a single PGV instance to include multiple datasets which could be browsed by going to the "Data Selection" page in the browser
#' @param ref.name the genome reference name used for this dataset. For specific behaviour refer to the PGV/gGnome.js wrappers
#' @param overwrite by default only files that are missing will be created. If set to TRUE then existing coverage arrow files and gGraph JSON files will be overwritten
#' @param annotation which node/edge annotation fields to add to the gGraph JSON file. By default we assume that gGnome::events has been executed and we add the following SV annotations: 'simple', 'bfb', 'chromoplexy', 'chromothripsis', 'del', 'dm', 'dup', 'pyrgo', 'qrdel', 'qrdup', 'qrp', 'rigma', 'tic', 'tyfonas'
gen_gg_json_files = function(data, outdir, meta.js, name.col = 'sample', gg.col = 'graph', js.type = 'gGnome.js',
                             dataset_name = NA, ref.name = NULL, overwrite = FALSE, annotation = NULL){
    json_dir = get_gg_json_dir_path(outdir, js.type, dataset_name)
    json_files = lapply(1:data[, .N], function(idx){
        gg.js = get_gg_json_path(data[idx, get(name.col)], json_dir)
        if (!file.exists(gg.js) | overwrite){
            # TODO: at some point we need to do a sanity check to see that a valid rds of gGraph was provided
            gg = readRDS(data[idx, get(gg.col)])
            sl = parse.js.seqlenghts(meta.js, js.type = js.type, ref.name = ref.name)
            # check for overlap in sequence names
            gg.reduced = gg[seqnames %in% names(sl)]
            if (length(gg.reduced) == 0){
                stop(sprintf('There is no overlap between the sequence names in the reference used by gGnome.js and the sequences in your gGraph. Here is an example sequence from your gGraph: "%s". And here is an example sequence from the reference used by gGnome.js: "%s"', seqlevels(gg$nodes$gr)[1], names(sl)[1]))
            }
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

#' @name is_git_lfs_available
#' @description internal
#'
#' Check if git-lfs command is available
#'
#' @param raise by default if git-lfs is not available then an error will be raised. Set raise = FALSE if you don't want an error to occur but just want to know if git-lfs is available
is_git_lfs_available = function(raise = TRUE){
    conn = pipe('command -v git-lfs')
    available = length(readLines(conn)) > 0
    close(conn)
    if (!available){
        if (raise){
            stop('git-lfs is not installed, please install git-lfs (https://git-lfs.github.com/)')
        }
        return(FALSE)
    }
    return(TRUE)
}

#' @name is.dir.a.PGV.instance
#' @description internal
#'
#' Check if a path matches something that looks like a clone of the PGV github repository
#'
#' This is done by checking if the folder contains the subdirectory "public" and within it "settings.json".
#' If the file is not found then an error is raised
#' 
#' @param outdir path to directory
is.dir.a.PGV.instance = function(outdir){
    if (!file.exists(paste0(outdir, '/public/settings.json'))){
        stop(outdir, ' does not seem to be a proper clone of the PGV github repository.')
    }
}

#' @name is.dir.a.gGnome.js.instance
#' @description internal
#'
#' Check if a path matches something that looks like a clone of the gGnome.js github repository
#'
#' This is done by checking if the folder contains the subdirectory "public" and within it "metadata.json".
#' If the file is not found then an error is raised
#' 
#' @param outdir path to directory
is.dir.a.gGnome.js.instance = function(outdir){
    if (!file.exists(paste0(outdir, '/public/metadata.json'))){
        stop(outdir, ' does not seem to be a proper clone of the gGnome.js github repository.')
    }
}

#' @name is.dir.a.js.instance
#' @description internal
#'
#' Check if a path matches something that looks like a clone of the PGV or gGnome.js github repositories
#'
#' This is done by checking if the folder contains the subdirectory "public" and within it "settings.json".
#' If the file is not found then an error is raised
#' 
#' @param outdir path to directory
#' @param js.type either "PGV" or "gGnome.js"
is.dir.a.js.instance = function(outdir, js.type){
    if (js.type == 'gGnome.js'){
        is.dir.a.gGnome.js.instance(outdir)
    }
    if (js.type == 'PGV'){
        is.dir.a.PGV.instance(outdir)
    }
}

#' @name js_path
#' @description internal
#'
#' Takes a path and checks if it is a valid path to a gGnome.js/PGV directory. 
#'
#' If the directory does not exist then a clone from github is generated.
#' 
#' @param outdir path to directory
#' @param append if set to TRUE then the directory is expected to already exist
#' @param js.type either "PGV" or "gGnome.js"
js_path = function(outdir, append = FALSE, js.type = 'gGnome.js'){

    is.acceptable.js.type(js.type)
    outdir = suppressWarnings(normalizePath(outdir))

    if (file.exists(outdir) & !dir.exists(outdir)){
        stop('The output directory must be a valid path for a diretory, but you provided a path of a file that already exists:', outdir)
    }

    if (dir.exists(outdir) & !append){
        # if the folder exists and there is no append flag then throw error
        stop('The output directory already exists. If you wish to generate a new ', js.type, ' isntance, please provide a path for a new directory. If you wish to add more file to an existing instance of gGnome.js then use "append = TRUE".')
    }

    if (!append){
        # clone the repository from github
        message('Cloning the ', js.type, ' repository from github.')
        is_git_lfs_available()
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

#' @name is.acceptable.js.type
#' @description internal
#'
#' Checks that the provided js.type is valid
#' 
#' @param js.type either "PGV" or "gGnome.js"
is.acceptable.js.type = function(js.type){
    if (!(js.type %in% acceptable.js.types)){
        stop('Invalid js.type. The only js types familiar to us are: ', acceptable.js.types)
    }
}

#' @name gen_js_coverage_files
#' @description internal
#'
#' Generate  the CSV (for gGnome.js) or arrow (for PGV) coverage files
#'
#' accepts any data type that is acceptable for #' 
#' @param data either a path to a TSV/CSV or a data.table
#' @param outdir the path to the PGV/gGnome.js repository clone
#' @param name.col column name in the input data table containing the sample names (default: "sample")
#' @param overwrite by default only files that are missing will be created. If set to TRUE then existing coverage arrow files and gGraph JSON files will be overwritten
#' @param cov.col column name in the input data table containing the paths to coverage files
#' @param js.type either "PGV" or "gGnome.js"
#' @param cov.field the name of the field in the coverage GRanges that should be used (default: "ratio")
#' @param cov.field.col column name in the input data table containing the name of the field in the coverage GRanges that should be used. If this is supplied then it overrides the value in "cov.field". Use this if some of your coverage files differ in the field used.
#' @param cov.bin.width bin width to use when rebinning the coverage data (default: 1e4). If you don't want rebinning to be performed then set to NA.
#' @param dataset_name the name of the dataset. Only relevant for PGV. This should be the name of the project that all the samples belong to. You must provide a name since PGV stores all the data under a folder matching your dataset name. This allows a single PGV instance to include multiple datasets which could be browsed by going to the "Data Selection" page in the browser
#' @param ref.name the genome reference name used for this dataset. For specific behaviour refer to the PGV/gGnome.js wrappers
#' @param cov.color.field field in the coverage GRanges to use in order to set the color of coverage data points. If nothing is supplied then default colors are used for each seqname (namely chromosome) by reading the colors that are defined in the settings.json file for the specific reference that is being used for this dataset.
#' @param meta.js path to JSON file with metadata (for PGV should be located in "public/settings.json" inside the repository and for gGnome.js should be in public/genes/metadata.json)
#' @param mc.cores how many cores to use
#' 
#' @export
gen_js_coverage_files = function(data, outdir, name.col = 'sample', overwrite = FALSE, cov.col = 'cov',
                                 js.type = 'gGnome.js', cov.field = 'ratio', cov.field.col = NA,
                                 bin.width = 1e4, dataset_name = NA, ref.name = 'hg19',
                                 cov.color.field = NULL, meta.js = NULL, mc.cores = 1){
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
                              cov.color.field = cov.color.field, overwrite = overwrite,
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

#' @name read.js.input.data
#' @description internal
#'
#' Accepts a TSV/CSV or a data.table and checks that it is properly formatted
#'
#' @param data either a path to a TSV/CSV or a data.table
#' @param name.col column name in the input data table containing the sample names (default: "sample")
read.js.input.data = function(data, name.col = 'sample'){
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

#' @name get_gg_json_path
#' @description internal
#'
#' get the path to the gGraph JSON file inside a js directory
#'
#' @param nm name of the sample
#' @param gg_json_dir the path to the directory holding the gGraphs JSON files
get_gg_json_path = function(nm, gg_json_dir){
    gg.js = paste0(gg_json_dir, "/", nm, ".json")
    return(gg.js)
}

#' @name get_gg_json_dir_path
#' @description internal
#'
#' get the path to the gGraph JSON directory inside a js repository clone
#'
#' @param outdir the path to the PGV/gGnome.js repository clone
#' @param js.type either "PGV" or "gGnome.js"
#' @param dataset_name the name of the dataset. Only relevant for PGV. This should be the name of the project that all the samples belong to. You must provide a name since PGV stores all the data under a folder matching your dataset name. This allows a single PGV instance to include multiple datasets which could be browsed by going to the "Data Selection" page in the browser
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

#' @name get_js_cov_path
#' @description internal
#'
#' get the path to the gGraph JSON file inside a js directory
#'
#' @param nm name of the sample
#' @param cov_dir the path to the directory holding the coverage files
#' @param js.type either "PGV" or "gGnome.js"
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

#' @name get_js_cov_dir_path
#' @description internal
#'
#' get the path to the directory holding the coverage files inside a PGV/gGnome.js repository clone
#'
#' @param outdir the path to the PGV/gGnome.js repository clone
#' @param js.type either "PGV" or "gGnome.js"
#' @param dataset_name the name of the dataset. Only relevant for PGV. This should be the name of the project that all the samples belong to. You must provide a name since PGV stores all the data under a folder matching your dataset name. This allows a single PGV instance to include multiple datasets which could be browsed by going to the "Data Selection" page in the browser
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

#' @name get_pgv_data_dir
#' @description internal
#'
#' get the path to the dataset's data dir inside the PGV directory
#'
#' @param outdir the path to the PGV/gGnome.js repository clone
#' @param dataset_name the name of the dataset. Only relevant for PGV. This should be the name of the project that all the samples belong to. You must provide a name since PGV stores all the data under a folder matching your dataset name. This allows a single PGV instance to include multiple datasets which could be browsed by going to the "Data Selection" page in the browser
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

#' @name cov2cov.js
#' @description
#'
#' Takes a GRanges with coverage data and converts it to a data.table with the info needed for gGnome.js and PGV
#'
#' if bin.width is specified then coverage data will also be rebinned
#' if convert.to.cn == TRUE then rel2abs will be applied
#'
#' @param cov coverage GRanges or path to file with coverage data (see gGnome::readCov for details
#' @param meta.js path to JSON file with metadata (for PGV should be located in "public/settings.json" inside the repository and for gGnome.js should be in public/genes/metadata.json)
#' @param js.type either "PGV" or "gGnome.js"
#' @param field the name of the field in the coverage GRanges that should be used (default: "ratio")
#' @param bin.width bin width to use when rebinning the coverage data (default: 1e4). If you don't want rebinning to be performed then set to NA.
#' @param ref.name the genome reference name used for this dataset. For specific behaviour refer to the PGV/gGnome.js wrappers
#' @param cov.color.field field in the coverage GRanges to use in order to set the color of coverage data points. If nothing is supplied then default colors are used for each seqname (namely chromosome) by reading the colors that are defined in the settings.json file for the specific reference that is being used for this dataset.
#' @return data.table containing the coverage data formatted according to the format expected by PGV and gGnome.js
#' @export
cov2cov.js = function(cov, meta.js = NULL, js.type = 'gGnome.js', field = 'ratio',
                      bin.width = NA, ref.name = NULL, cov.color.field = NULL){
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
        stop("The provided field '", field, "' is not in the input coverage data")
    }

    fields = field

    if (!is.null(cov.color.field)){
        if (!is.element(cov.color.field, names(mcols(x)))){
            stop("The cov.color.field '", cov.color.field, "' is not in the input data")
        }
        fields = c(field, cov.color.field)
    }

    if (!is.na(bin.width)){
        message('Rebinning coverage with bin.width=', bin.width)
        if (!is.numeric(bin.width)){
            stop('bin.width must be numeric')
        }
        x.rebin = rebin(x, bin.width, field, FUN = median)
        if (!is.null(cov.color.field)){
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

#' @name cov2cov.js
#' @description
#'
#' Takes a GRanges with coverage data and converts it to a data.table with the info needed for gGnome.js and PGV
#'
#' if bin.width is specified then coverage data will also be rebinned
#' if convert.to.cn == TRUE then rel2abs will be applied
#'
#' @param cov coverage GRanges or path to file with coverage data
#' @param meta.js path to JSON file with metadata (for PGV should be located in "public/settings.json" inside the repository and for gGnome.js should be in public/genes/metadata.json)
#' @param js.type either "PGV" or "gGnome.js"
#' @param field the name of the field in the coverage GRanges that should be used (default: "ratio")
#' @param bin.width bin width to use when rebinning the coverage data (default: 1e4). If you don't want rebinning to be performed then set to NA.
#' @param ref.name the genome reference name used for this dataset. For specific behaviour refer to the PGV/gGnome.js wrappers
#' @param cov.color.field field in the coverage GRanges to use in order to set the color of coverage data points. If nothing is supplied then default colors are used for each seqname (namely chromosome) by reading the colors that are defined in the settings.json file for the specific reference that is being used for this dataset.
#' @return data.table containing the coverage data formatted according to the format expected by PGV and gGnome.js
#' @export
cov2csv = function(cov,
        field = "ratio",
        output_file = "coverage.csv",
        ...)
{

    dat = cov2cov.js(cov, field = field, ...)

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
#' @param cov.color.field a field in the input GRanges object to use to determine the color of each point
#' @param overwrite (logical) by default, if the output path already exists, it will not be overwritten.
#' @param meta.js path to JSON file with metadata for PGV (should be located in "public/settings.json" inside the repository)
#' @param bin.width (integer) bin width for rebinning the coverage (default: 1e4)
#' @author Alon Shaiber
#' @export
cov2arrow = function(cov,
        field = "ratio",
        output_file = 'coverage.arrow',
        ref.name = 'hg19',
        cov.color.field = NULL,
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
                         ref.name = ref.name, cov.color.field = cov.color.field, ...)
        message('Done converting coverage format')


        if (!is.null(cov.color.field)){
            dat[, color := color2numeric(get(cov.color.field))]
        } else {
            if (!is.null(meta.js)){
                ref_meta = get_ref_metadata_from_PGV_json(meta.js, ref.name)
                setkey(ref_meta, 'chromosome')
                dat$color = color2numeric(ref_meta[dat$seqnames]$color)
            } else {
                # no cov.color.field and no meta.js so set all colors to black
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
#' @param gr GRanges object with the GTF information
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
                    gr=NULL,
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
    } else if (!is.null(gr)){
        message("Using input GRanges.")
        infile = 'Input GRanges'
        dt = gr2dt(gr)
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
    message(sprintf('Wrote JSON metadata file of %s to %s', infile, metadata.filename))


    out_genes = paste(c(genes.json, "}"),
                     sep = "")
    writeLines(out_genes, genes.filename)
    message(sprintf('Wrote JSON genes file of %s to %s', infile, genes.filename))

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
#' @param meta.js path to JSON file with metadata (for PGV should be located in "public/settings.json" inside the repository and for gGnome.js should be in public/genes/metadata.json)
#' @param js.type either 'gGnome.js' or 'PGV' to determine the format of the JSON file
#' @param ref.name the name of the reference to load (only relevant for PGV). If not provided, then the default reference (which is set in the settings.json file) will be loaded.
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

#' @name get_ref_metadata_from_PGV_json
#' @description internal
#' get a data.table with the metadata for the reference (columns: chromosome, startPoint, endPoint, color)
#' @param meta.js path to JSON file with metadata (for PGV should be located in "public/settings.json" inside the repository and for gGnome.js should be in public/genes/metadata.json)
#' @param ref.name the name of the reference to load (only relevant for PGV). If not provided, then the default reference (which is set in the settings.json file) will be loaded.
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

#' @name get_path_to_meta_js
#' @description internal
#' Get the path to the meta.js file
#' @param outdir the path where to the PGV/gGnome.js folder.
#' @param js.type either "PGV" or "gGnome.js"
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

