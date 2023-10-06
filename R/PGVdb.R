# one create PGV instance 
# inside we have the add function to initialize the PGV datafile.json and get pointers to where all the files are
# also set up datafiles.json and get pointers to where the data folder and datafiles.json are being stored
# will need to set up datafiles.json and data folder along with backups for each
# probably would get bloat backing up data folder so that is a bit suspect.
set_up_dev_build_PGV = {}


# so the file locs across prod vs dev is just where the two files are:
# for prod: is flat with pgv
# for dev: stored with public/
# not sure what to do with this yet
# .pgv_type = function(pgv_dir = NULL){
#   if (file.exists(paste0(pgv_dir,"/datafiles.json")) &&
#                   dir.exists(paste0(pgv_dir, "/data"))){
#     message("production build found")
#     datafiles.json_path = normalizePath(paste0(pgv_dir,"/datafiles.json"))
#     datafolder_path = normalizePath(paste0(pgv_dir, "/data"))
#   } else if (file.exists(paste0(pgv_dir, "/public/datafiles.json")) &&
#              dir.exists(paste0(pgv_dir, "/public/data"))){
#     message("development build found")
#     datafiles.json_path = normalizePath(paste0(pgv_dir,"/public/datafiles.json"))
#     datafolder_path = normalizePath(paste0(pgv_dir, "/public/data"))
#   } else{
#     stop("pgv dir not set up as expected. For the development build we expect datafiles.json 
#          and data to be stored in the public folder in pgv. For Production, datafiles.json 
#          and data folder should be stored flat in the pgv folder.")
#   }
#   return(c(datafiles.json_path, datafolder_path))
# }

# function to check db validity
.valid_json_db = function(json_db, 
                          full_check = FALSE){
  # get patient IDs
  if (!file.exists(json_db$settings.js)){
    stop("stop settings.json path must be provided.")
  }
  # check reference
  if (!length(unique(json_db$references$patient.id)) == 
      length(json_db$references$patient.id)){
    stop("duplicate patient ids. found in references table")
  }
  # check description
  if (any(!(unique(json_db$descriptions$patient.id) %in% json_db$references$patient.id))){
    stop("patient ids were added to description tags that are not yet added. 
         Please run add_patients_PGV to add patients before adding their 
         tags.")
  }
  if (any(dim(unique(json_db$descriptions)) != dim(json_db$descriptions))){
    message( "non unique descriptions for patients in the json_db$descriptions found.
          Removing them before proceeding")
    json_db$descriptions = unique(json_db$descriptions) 
  }
  #check that there are graphs for each patient added
  if (any(!(unique(json_db$descriptions$patient.id) %in% unique(json_db$plots$patient.id)))){
    message("patient ids with no graphs were found. Dropping them.")
    ids_drop = json_db$references$patient.id[which(!(unique(json_db$descriptions$patient.id) %in% 
                                                       unique(json_db$plots$patient.id)))]
    if (any(json_db$descriptions$patient.id %in% ids_drop)){
      json_db$descriptions = json_db$descriptions[-which(patient.id %in% ids_drop)]
    }
    json_db$references = json_db$references[-which(patient.id %in% ids_drop)]
  }
  # check for graphs that do not have patient.id 
  
  # check file loc
  if (!file.exists(json_db$datafiles.json)){
    stop("file does not exist for json_db$datafiles.json")
  }
  # check write out location
  if (!dir.exists(json_db$data_folder)){
    stop("data folder does not exist")
  }
  # check that all patient graphs have source files exists.
  if (full_check){
    for (i in 1:nrow(json_db$plots)){
      if (!is.na(json_db$plots$source[i])){
        if (!file.exists(normalizePath(paste0(json_db$data_folder, "/",
                                              json_db$plots$patient.id[i], "/",
                                              json_db$plots$source[i])))){
          message("file not found for ", paste0(json_db$data_folder, "/",
                                               json_db$plots$patient.id[i], "/",
                                               json_db$plots$source[i]))
          warning("consider dropping row ", i, " in json_db$plots. ")
        }
      }
    }
  message("all files found and exists")
  }
  return(json_db)
}  

# function to check table add validity
.valid_add_table = function(table_add, descriptors){
  # expect columns: patient.id is essential.
  if (!("patient.id" %in% colnames(table_add))){
    stop("patient.id column must be in table_add and not NA")
  }
  if (!("name.col" %in% colnames(table_add))){
    message("name.col not provided, passing patient.id as name.col")
    table_add$name.col = table_add$patient.id
  }
  if (!any(descriptors %in% colnames(table_add))){
    # create a tag of sample name 
    table_add$tags = paste0("patient.id=", table_add$patient.id)
  }
  if (!("ref" %in% colnames(table_add))){
    stop("expects reference column to be added. Currently supporting hg19/hg38")
  } else {
    if (any(is.na(table_add$ref))){
      stop("No NAs allowed for ref. Ref must always be provided for all samples")
    }
  }
  return(table_add)
}
# add new participants/patient IDs to the database
# this will be used to populate the db initially as well
# table with expected columns
# descriptor cols if tags was not used
# push_PGV pushes JSON_db to datafiles.json
add_patients_PGV = function(json_db,
                            table_add = NULL, 
                            descriptors = 'tags',
                            # overwrite = F,
                            push_PGV = T,
                            cores = 10){
  # check table add for having right args in
  if (is.null(table_add)){
    # returns a table of any expected args
    message("Returning an add_table example format:
            patient.id: character vector of unique patient id.
            name.col: character vector of a specific sample for a patient. 
            If NA, no name for graphs added are given.
            gg.col: ggraph of sample
            gw.col: gwalk of sample
            cov.col: coverage of sample
            annotation: if passed events as graph what events to look for, 
            Default= c('simple', 'bfb', 'chromoplexy',
                              'chromothripsis', 'del', 'dm', 'dup',
                              'pyrgo', 'qrdel', 'qrdup', 'qrp', 'rigma',
                              'tic', 'tyfonas')
            tree: tree file for phylogeny of that patient sample
            ref: reference build of sample: hg19/hg38
            cov.field: 'ratio' or 'background' (dryclean) if NA passed, uses 1st column
            overwrite: FALSE, rewrites graphs if these graphs already exist
            tags: description tags of patient id to filter by")
    table_add = data.table(patient.id = NA,
                           name.col = NA,
                           gg.col = NA,
                           gw.col = NA,
                           cov.col = NA,
                           annotation = list(c('simple', 'bfb', 'chromoplexy',
                                          'chromothripsis', 'del', 'dm', 'dup',
                                          'pyrgo', 'qrdel', 'qrdup', 'qrp', 'rigma',
                                          'tic', 'tyfonas')),
                           tree = NA,
                           ref = 'hg19',
                           cov.field = 'ratio',
                           tags = NA)
    return(table_add)
  } else {
    table_add = .valid_add_table(table_add, descriptors)
    # add tags to our db
    json_db = add_tags_PGV(json_db, table_add, descriptors)
    # check if the graphs being added and sample name does not exist
    dirpaths = file.path(paste0(json_db$data_folder, 
                                    "/", table_add$patient.id))
    table_add$dirpaths = dirpaths
    if (any(dir.exists(dirpaths))){
      # some dir paths are found -> output to user to let them know they exist
      paths_exist = table_add$patient.id[which(dir.exists(dirpaths))]
      # and drop these ids.
      message(paths_exist, " path exists for this patient.id. To add the sample back, 
      manually add to the PGVdb object with associated tags, ref and sample names.")
      table_add = table_add[-which(dir.exists(dirpaths)),]
    } 
    
    #update json_db$references with given patient.ids.
    json_db_add = unique(table_add[,c("patient.id", "ref"), with=FALSE])
    colnames(json_db_add) = c("patient.id", "reference")
    json_db$references = rbind(json_db$references, json_db_add)
    
    # then for each patient we run this on an mclapply
    l_graphs = mclapply(1:nrow(table_add), function(i){
    # we create sample dir if it does not exist
    if (!dir.exists(file.path(paste0(json_db$data_folder, "/", 
                                         table_add$patient.id[[i]])))){
      dir.create(file.path(paste0(json_db$data_folder, "/", 
                                      table_add$patient.id[[i]])))
    }
    # note to MAX: we need to account for patient.id cols with multiple pairs
    # within the add_table and create unique labels probably
    new_graphs = add_graphs_PGV(json_db, 
                                table_row = table_add[i,])
  }, mc.cores = cores)
  l_graphs = l_graphs %>% 
      data.table::rbindlist(., fill = T)
    # from here we add this data.table to our original json_db$plots
  json_db$plots = rbind(json_db$plots, l_graphs, fill =T)
  if (push_PGV){
    push_PGV_db(json_db)
  }  
  return(json_db)
  }
}

# table add is a table of tags to add for patients
# for table it is patient.id, tags to add, tag, to add 2, etc.
# can be done multiple ways -> figure this out later
add_tags_PGV = function(json_db,
                        table_add, 
                        descriptors = "tags"){
  # read in json_db 
  # read in and error out if descriptors columns are not found
  if (any(descriptors %in% colnames(table_add))){
    if (any(!(descriptors %in% colnames(table_add)))){
      descriptors = descriptors[descriptors %in% colnames(table_add)] 
      message("dropping unused column: ", 
              descriptors[!(descriptors %in% colnames(table_add))])
    }
  } else {
    descriptors = "tags"
    if (!(descriptors %in% colnames(table_add))){
      stop("descriptor columns not found in table_add")
    }
  }
  # first grab deduplicated rows:
  if (any(duplicated(table_add$patient.id))){
    table_dedup = table_add[-which(duplicated(table_add$patient.id)),]
  } else {
    table_dedup = table_add
  }
  # if descriptors is a list of cols -> we convert to tags 
  if (length(descriptors) > 1){
    list_desc = lapply(1:nrow(table_dedup), function(x){
        out = table_dedup[x,(descriptors), with=F] %>% paste0(colnames(.), "=", .)
        data.table(patient.id = rep(table_dedup$patient.id[x], 
                                    length(out)),
                   tags = out)
      }) %>% data.table::rbindlist(.,fill = T)
  } else {
      # we expect it to be a tags list with either , or ; separated values
      list_desc = lapply(1:nrow(table_dedup), function(x){ 
        out = strsplit(table_dedup[,get(descriptors)][[x]], split = ";|,")
        data.frame(patient.id = rep(table_dedup$patient.id[x], 
                                    length(out[[1]])), 
                   tags = out[[1]])
      }) %>% data.table::rbindlist(.,fill = T)
    }
  # append the tags to json_db$description
  json_db$descriptions = rbind(json_db$descriptions, list_desc)
  # grab uniques and message out if we had any non unique rows added
  return(json_db)
}


# checks for type and then runs the gg_gen_... for that function
# table add format: patient.id, file path, 
# type i.e. ggraph, cov, gwalk, phylogeny, barplot, scatterplot, bigwig
add_graphs_PGV = function(json_db, table_row){
  # add graphs -> 1 row check for non NA entries in the 4 columns
  # dropping empty columns
  table_row[ , which(sapply(table_row, 
                            function(x) all(is.na(x)))) := NULL]
  # ggraph for this file
  if ("gg.col" %in% colnames(table_row)){
    gen_gg_json_PGV(table_row, json_db)
    # create a data.table
    gg.row = data.table(sample = table_row$name.col,
                        type = "genome",
                        source = paste0(table_row$name.col,".json"),
                        title = table_row$name.col,
                        visible = TRUE,
                        figure = NA,
                        server = NA,
                        uuid = NA,
                        patient.id = table_row$patient.id,
                        plot_id = NA)
  }
  # cov
  if ("cov.col" %in% colnames(table_row)){
    gen_js_cov_PGV(table_row, json_db)
    cov.row = data.table(sample = table_row$name.col,
                         type = "scatterplot",
                         source = paste0(table_row$name.col,
                                         "-coverage.arrow"),
                         title = paste0(table_row$name.col,
                                        ' Coverage Distribution'),
                         visible = TRUE,
                         figure = NA,
                         server = NA,
                         uuid = NA,
                         patient.id = table_row$patient.id,
                         plot_id = NA)
  }  
  # gwalk
  if ("gw.col" %in% colnames(table_row)){
    if (file.exists(table_row$gw.col)){
      gen_gw_json_PGV(table_row, json_db)
      gw.row = data.table(sample = table_row$name.col,
                        type = "walk",
                        source = paste0(table_row$name.col,".walks.json"),
                        title = paste0(table_row$name.col, " Walks"),
                        visible = TRUE,
                        figure = NA,
                        server = NA,
                        uuid = NA,
                        patient.id = table_row$patient.id,
                        plot_id = NA)
    } else {
      warning("gw file does not exist for ", table_row$name.col)
    }
  }
  # tree
  if ("tree" %in% colnames(table_row)){
    if (file.exists(table_row$tree)){
      file.copy(table_row$tree, paste0(json_db$data_folder,"/",
                                       table_row$patient.id,"/",
                                       table_row$patient.id, ".newick"))
    # copy from source dir to the end dir for this file
      tree.row = data.table(sample = NA,
                          type = "phylogeny",
                          source = paste0(table_row$patient.id, ".newick"),
                          title = paste0("Phylogenetic Information for ", 
                                         table_row$patient.id),
                          visible = TRUE,
                          figure = NA,
                          server = NA,
                          uuid = NA,
                          patient.id = table_row$patient.id,
                          plot_id = NA)
    } else {
      warning("tree file for ", table_row$patient.id, " did not exist.")
    }
  }
  # from our 1 row generate the data.table of row size up to 
  # 4 (graph, cov, gwalk, tree) for this row
  # check if variables exist and then rbindlist existing variables
  if (length(which(c(exists("gg.row"), exists("cov.row"),
        exists("gw.row"), exists("tree.row")))) > 0){
    dat <- data.table(sample=character(), type=character(), 
                      source=character(), title=character(),
                      visible=logical(), figure=character(),
                      server=character(), uuid=character(),
                      patient.id=character(), plot_id=character())
    dat <- rbind(dat, if(exists("gg.row")) gg.row)
    dat <- rbind(dat, if(exists("cov.row")) cov.row)
    dat <- rbind(dat, if(exists("gw.row")) gw.row)
    dat <- rbind(dat, if(exists("tree.row")) tree.row)
    return(dat)
  } else {
    warning("Nothing existed for this sample: ", table_row$name.col)
  }
}

# remove_patient removal of patient from the database. Hard remove or mask?
drop_patients_PGV = function(
  json_db,
  patient_ids,
  delete = FALSE){
  drop_ids = unique(patient_ids)
  if (any(drop_ids %in% json_db$references$patient.id)){
    if (any(!(drop_ids %in% json_db$references$patient.id))){
      warning(paste0(drop_ids[which(!(drop_ids %in% json_db$references$patient.id))], 
                     " did not match with any ids. "))
      drop_ids = drop_ids[-which(!(drop_ids %in% json_db$references$patient.id))]
    }
  }  else {
    stop("patient ids provided do not match any ids in json_db")
  }
  if (delete){
    message("deleting directories")
    # get all file paths to remove
    remove_dt = json_db$references[patient.id %in% drop_ids]
    patient_name = remove_dt$patient.id
    for (i in 1:length(patient_name)){
        message("deleting directory ", normalizePath(paste0(json_db$data_folder, "/", 
                                                       patient_name[i])))
        unlink(normalizePath(paste0(json_db$data_folder, "/",
                                         patient_name[i])), recursive = T)
    }
  }
  # clean up tables and drop vals
  json_db$descriptions = json_db$descriptions[!(patient.id %in% drop_ids)]
  json_db$references = json_db$references[!(patient.id  %in% drop_ids)]
  json_db$plots = json_db$plots[!(patient.id  %in% drop_ids)]
  return(json_db)
}

# drop_ids: a list of unique plot_ids to drop: 
# remove: fully remove from database? deletes object and value from data folder
drop_graphs_PGV = function(json_db, drop_ids, delete = FALSE){
  # check drop_ids to see which graphs to drop
  # get unique of these ids in case duplicates were provided
  drop_ids = unique(drop_ids)
  if(any(drop_ids %in% json_db$plots$plot_id)){
    if (any(!(drop_ids %in% json_db$plots$plot_id))){
      warning(paste0(drop_ids[which(!(drop_ids %in% json_db$plots$plot_id))], 
                     " did not match with any ids. "))
      drop_ids = drop_ids[-which(!(drop_ids %in% json_db$plots$plot_id))]
    }
  } else {
    stop("ids provided do not match any ids in json_db")
  }
  # validity check to see if these samples/graphs exist to drop
  if (delete){
    message("deleting files")
    # get all file paths to remove
    remove_dt = json_db$plots[plot_id %in% drop_ids]
    remove_files = remove_dt$source
    patient_name = remove_dt$patient.id
    for (i in 1:length(remove_files)){
      if (!is.na(remove_files[i])){
        message("deleting file ", normalizePath(paste0(json_db$data_folder, "/", 
                                                       patient_name[i], "/",
                                                       remove_files[i])))
        file.remove(normalizePath(paste0(json_db$data_folder, "/",
                                         patient_name[i], "/",
                                         remove_files[i])))
      }
    }
  }
  json_db$plots = json_db$plots[!(plot_id %in% drop_ids)]
  return(json_db)
}
# pass patient ids for subset, or if NA drop this tag across all
# drop_graphs removes graphs of subset
# pass patient ids & pair for graphs to be dropped

# return_PGV_db returns the datafiles.json in our object format.
# let user edit
return_PGV_db = function(datafiles.json,
                         data_folder,
                         PGV_public_dir){
  df_json = jsonlite::fromJSON(datafiles.json)
  df_description = lapply(names(df_json), function(x){
    data.table::data.table(patient.id = x, 
               tags = df_json[[x]]$description)
  }) %>% data.table::rbindlist(.,fill = T)
  df_reference = lapply(names(df_json), function(x){
    data.table::data.table(patient.id = x,  
                           reference = df_json[[x]]$reference)
  }) %>% data.table::rbindlist(.,fill = T)
  df_plots = lapply(names(df_json), function(x){
    plots = as.data.table(df_json[[x]]$plots)
    plots$patient.id = x
    # give plots unique ids for people to drop 
    plots$plot_id = paste0(x, "_", (1:nrow(df_json[[x]]$plots)))
    return(plots)
  }) %>% data.table::rbindlist(., fill = T)
  # shoddy fix to issue where we change NULL values to NA.. might have to keep
  # some of the weirdness to allow for the {} returns
  is.na(df_plots) <- df_plots == "NULL"
  
  # checks to make sure we have everything in the database server and that files
  # are as expected ->
  # dir.exists(data_folder){}
  
  # check for if the files listed in datafiles.json are there. 
  # If not spit warning 
  # for each file if expected file is missing
  if (file.exists(paste0(PGV_public_dir, "/settings.json"))){
    settings.js = file.path(paste0(PGV_public_dir, "/settings.json"))
  } else {
    warning("settings.json file not found in public dir.")
    settings.js = "settings json file was not found. Replace this with path to settings.json"
  }
  # add file.source to df_plots if we are generating new plots
  # df_plots$file.source = NA
  json_db = list(descriptions = df_description, references = df_reference, 
                  plots = df_plots, datafiles.json = datafiles.json,
                 data_folder = data_folder, 
                 settings.js = settings.js)
  return(json_db)
}


# push_PGV_db pushes the changes to datafiles.json & data folder
push_PGV_db = function(json_db, backup = TRUE){
  json_db = .valid_json_db(json_db)
  # clean up and remove added columns
  
  
  # plots
  json_db$plots[,plot_id:=NULL]
  # json_db$plots[,file.source:=NULL]
  #
  pt_ids = unique(json_db$references$patient.id)
  json_format = lapply(pt_ids, function(x){
    description = json_db$descriptions[patient.id == x]$tags
    plots_Df = json_db$plots[patient.id == x, !"patient.id"]
    return(list(description = description, 
           reference = json_db$references[patient.id == x]$reference, 
           plots = plots_Df))
  })
  names(json_format) = pt_ids
  
  # here we save a backup version with date & time for now
  file.rename(json_db$datafiles.json, 
              paste0(json_db$datafiles.json,
                     format(Sys.time(), "%Y%m%d_%H%M%S")))
  # re json
  
  jsonlite::write_json(json_format, json_db$datafiles.json,
                       pretty=TRUE, auto_unbox=TRUE, digits=4)
}


# this function will be default TRUE for append, update, remove_patient, 
# drop_tags, and drop_graphs
# push changes will check changes with respect to previous_PGV file 
# so read in, comp changes, generate files as needed and then push
# will need some stringent checking to make sure things are good to go



# revert_PGV_db reverts the version back based on link to file, otherwise
# default to most recent version
revert_PGV_db = function(current_json_path, json_db,
                         old_json_path=NULL){
  #if(!is.null(old_json_path)){
  file.rename(old_json_path, 
              current_json_path)
  #} else {
  # find most recent version saved in pgvdir
  # } 
}



### Gen graphs functions

#' @name gen_gg_json_PGV
#' @description internal
#'
#' Generate the json files that will represent your gGraphs
#'
#' @param data either a path to a TSV/CSV or a data.table
#' @param outdir the path to the PGV/gGnome.js repository clone
#' 
#' @details returns out a json to the save via R6 capabilities
gen_gg_json_PGV = function(table_row, json_db){
  gg.js = file.path(table_row$dirpaths, 
                    paste0(table_row$name.col,".json"))
  if (file.exists(gg.js)){
    warning("file ", gg.js, "already exists. Delete if you want to update.")
  } else {
    print(paste0("reading in ", table_row$gg.col))
    # TODO: at some point we need to do a sanity check to see that a valid rds of gGraph was provided
    if (grepl(table_row$gg.col, pattern = ".rds")){
      gg = readRDS(table_row$gg.col)
    } else{
      message("expected .rds ending for gGraph. Still attempting to read: ", table_row$gg.col)
      gg = readRDS(table_row$gg.col)
    }
    if (any(class(gg) == "gGraph")){  
        sl = parse.js.seqlengths(json_db$settings.js, 
                               js.type = "PGV", 
                               ref = table_row$ref)
        # check for overlap in sequence names
        gg.reduced = gg[seqnames %in% names(sl)]
        if (length(gg.reduced) == 0){
          stop(sprintf('There is no overlap between the sequence names in the reference 
                     used by PGV and the sequences in your gGraph. Here is an 
                     example sequence from your gGraph: "%s". And here is an 
                     example sequence from the reference used by gGnome.js: "%s"', 
                     seqlevels(gg$nodes$gr)[1], names(sl)[1]))
        }
        # sedge.id or other field
        if (exists("annotation")){
        # probably check for other cid.field names?
          #field = 'sedge.id'
          refresh(gg[seqnames %in% names(sl)])$json(filename = gg.js,
                                                  verbose = TRUE,
                                                  annotation = table_row$annotation[[1]])#,
                                                  #cid.field = field)
        } else {
          refresh(gg[seqnames %in% names(sl)])$json(filename = gg.js,
                                                  verbose = TRUE)
        }
    } else {
      warning(table_row$gg.col, " rds read was not a gGraph")
    }
  }
}


#' @name gen_js_cov_PGV
#' @description internal
#'
#' Generate arrow coverage files
#'
#' @param table_row single row in data.table for sample and patient.id
#' @param json_db deserialized json database for PGV.
#' 
#' @details creates a PGV ready arrow with respect to inputs in cov.col and stores it
#' in the json_db$data_folder / patient.id / name.col -coverage.arrow
#' @export
gen_js_cov_PGV = function(table_row, json_db){
  cov_dir = table_row$dirpaths
  covfn = file.path(paste0(table_row$dirpaths, "/",
                            table_row$name.col, "-coverage.arrow"))
  skip_cov = FALSE
  if (!file.exists(covfn)){
      if (is.na(table_row$cov.field)){
        warning(paste0('No coverage field was provided for ', 
                       table_row$cov.col, 
                       ' so no coverage will be generated.'))
        skip_cov = TRUE
      } else {
        if (is.na(table_row$cov.col)){
          warning(paste0('No coverage data was provided for ', 
                         table_row$name.col, 
                         ' so no coverage will be generated.'))
          skip_cov = TRUE
        } else {
          cov_input_file = table_row$cov.col
          if (is.na(cov_input_file)){
            warning(paste0('No coverage file was provided for ', table_row$name.col, 
                           ' so no coverage will be generated.'))
            skip_cov = TRUE
          } else {
            if (!file.exists(cov_input_file)){
              warning(paste0('Input coverage file does not exist for name: ', 
                             table_row$name.col, 
                             ' so no coverage will be generated.'))
              skip_cov = TRUE
            }}}
        }
      if (skip_cov){
        return(NA)
      } else {
          # load gGraph
          cov2arrowPGV(cov_input_file, 
                    field = table_row$cov.field,
                    meta.js = json_db$settings.js, #gg = gg, 
                    ref = table_row$ref,
                    output_file = covfn)
        }
    } else {
      message(covfn, ' found. Will not overwrite it.')
    }
}

#' @name cov2arrowPGV
#' @description
#'
#' Prepares an scatter plot arrow file with coverage info for PGV (https://github.com/mskilab/pgv)
#'
#' @param cov input coverage data (GRanges)
#' @param field which field of the input data to use for the Y axis
#' @param output_file output file path.
#' @param ref the name of the reference to use. If not provided, then the default reference that is defined in the meta.js file will be loaded.
#' @param cov.color.field a field in the input GRanges object to use to determine the color of each point
#' @param overwrite (logical) by default, if the output path already exists, it will not be overwritten.
#' @param meta.js path to JSON file with metadata for PGV (should be located in "public/settings.json" inside the repository)
#' @param bin.width (integer) bin width for rebinning the coverage (default: 1e4)
#' @author Alon Shaiber, Max Chao
#' @export
cov2arrowPGV = function(cov,
                     field = "ratio",
                     output_file = 'coverage.arrow',
                     ref = 'hg19',
                     meta.js = NULL,
                     ...){
  if (!file.exists(output_file)){
    if (!requireNamespace("arrow", quietly = TRUE)) {
      stop('You must have the package "arrow" installed in order for converting a
           coverage file to arrow file to work. Please install it.')
    }
    message('Converting coverage format')
    dat = cov2cov.js(cov, meta.js = meta.js, 
                     js.type = 'PGV', field = field,
                     ref = ref, ...)
    message('Done converting coverage format')
    if (!is.null(meta.js)){
        ref_meta = get_ref_metadata_from_PGV_json(meta.js, ref)
        setkey(ref_meta, 'chromosome')
        # create a map
        # 3.981s
        map_cols = data.table(color = unique(ref_meta$color),
                   numcolor = color2numeric(unique(ref_meta$color)))
        dat$color = merge(ref_meta[dat$seqnames], map_cols, 
                           by = "color", sort = F)$numcolor

    } else {
        # no cov.color.field and no meta.js so set all colors to black
        dat$color = 0
    }
    outdt = dat[, .(x = new.start, y = get(field), color)]
    # if there are any NAs for colors then set those to black
    outdt[is.na(color), color := 0]
    # remove NAs
    outdt = outdt[!is.na(y)]
    
    # sort according to x values (that is what PGV expects)
    outdt = outdt[order(x)]
    
    message('Writing arrow file (using write_feather)')
    arrow_table = arrow::Table$create(outdt, 
                                      schema = arrow::schema(x = arrow::float32(), 
                                                             y = arrow::float32(), 
                                                             color = arrow::float32()))
    arrow::write_feather(arrow_table, output_file)
  } else {
    message('arrow file, "', output_file, '" already exists.')
  }
  return(output_file)
}

#' @name gen_gw_json_PGV
#' @description internal
#'
#' Generate json files that will represent your gWalk objects
#'
#' @param table_row single row in data.table for sample and patient.id
#' @param json_db deserialized json database for PGV.
gen_gw_json_PGV= function(table_row, json_db){
  json_dir = table_row$dirpaths
  gw.js = file.path(json_dir, 
                    paste0(table_row$name.col, 
                           ".walks.json"))
  if (!file.exists(gw.js)){
    print(paste0("reading in ", table_row$gw.col))
      # TODO: at some point we need to do a sanity check to see that a valid rds of gWalk was provided
    gw = readRDS(table_row$gw.col) %>% 
      refresh
    if (gw$length == 0) {
      warning(sprintf("Zero walks in gWalk .rds file provided for sample %s! 
                      No walks json will be produced!", table_row$name.col))
        return(NA)
    }
    gw$json(filename = gw.js, verbose = TRUE,
            annotation = table_row$annotation[[1]],
            include.graph = FALSE)
    # sn.ref = parse.js.seqlengths(json_db$settings.js, 
    #                              js.type = "PGV", 
    #                              ref = table_row$ref) %>% names
    # sn.walks = seqlevels(gw)
    # sn.walks.only = sn.walks[!sn.walks %in% sn.ref]
    # gw.reduced = gw %&% sn.ref
    # if (length(sn.walks.only) > 0) { 
    #   gw.reduced = gw.reduced[gw.reduced %^% sn.walks.only == FALSE]
    # }
    # if (gw.reduced$length == 0){
    #   warning(sprintf('Provided gWalk .rds for sample %s contained walks, but they 
    #                   all involved sequences not contained in the chosen reference 
    #                   genome, so no walks json will be produced! Here is an example 
    #                   sequence name from your gWalks: "%s". And here is an example 
    #                   sequence from the reference used by PGV: "%s"', 
    #                   table_row$name.col, 
    #                   sn.walks.only[1],
    #                   sn.ref[1]))
    #     return(NA)
    #   }
    #   if (length(sn.walks.only) > 0) {
    #     gw.excluded = gw %&% sn.walks.only
    #   } else {
    #     gw.excluded = gW()
    #   }
    #   if (gw.excluded$length > 0) {
    #     warning(sprintf('%i walks excluded because they (fully or partially) 
    #                     fell outside of reference ranges.', 
    #                     gw.excluded$length))
    #   }
      # also.print.graph.to.json = ifelse(js.type == "PGV", FALSE, TRUE)
      # gw.js = gw.reduced$json(filename = gw.js, verbose = TRUE, 
      #                         annotation = table_row$annotation[[1]], 
      #                         include.graph = FALSE)
    
    } else {
      message(gw.js, ' found. Will not overwrite it.')
    }
}
