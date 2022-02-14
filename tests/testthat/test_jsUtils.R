library(testthat)
library(gUtils)

setDTthreads(1)
gg.jabba = gG(jabba = system.file('extdata/hcc1954', 'jabba.rds', package="gGnome"))

test_that('gtf2json', {
  # download light weight gtf
  system('wget http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/gencode.v19.long_noncoding_RNAs.gtf.gz')
  gtf = './gencode.v19.long_noncoding_RNAs.gtf.gz'
  # download chrom sizes from gGnome
  chrom.sizes = system.file("extdata", "human_g1k_v37.regular.chrom.sizes", package="gGnome")
  gr = rtracklayer::import.gff(gtf)
  jsons = gtf2json(gr = gr)
  expect_true(file.exists(jsons$metadata.filename) & file.exists(jsons$genes.filename))

  jsons = gtf2json(gtf = gtf, chrom.sizes = chrom.sizes)
  expect_true(file.exists(jsons$metadata.filename) & file.exists(jsons$genes.filename))

  expect_error(gtf2json(gr = gr, chrom.sizes = 4))
  expect_error(gtf2json(gr = gr, chrom.sizes = 'not a file'))
  expect_error(gtf2json(gr = gr, chrom.sizes = gtf)) # providing bad file for chrom.sizes

  saveRDS(gr, 'gtf.gr.rds')
  jsons = gtf2json(gtf.gr.rds = 'gtf.gr.rds', chrom.sizes = chrom.sizes)
  expect_true(file.exists(jsons$metadata.filename) & file.exists(jsons$genes.filename))

  saveRDS(gr2dt(gr), 'gtf.rds')
  jsons = gtf2json(gtf.rds = 'gtf.rds', chrom.sizes = chrom.sizes)
  expect_true(file.exists(jsons$metadata.filename) & file.exists(jsons$genes.filename))

  expect_error(gtf2json(verbose = TRUE, chrom.sizes = chrome.sizes))

  jsons = gtf2json(gr = gr, chrom.sizes = chrom.sizes, include.chr = c('chr1','chr2'), genes = 'MIR137HG')
  expect_true(file.exists(jsons$metadata.filename) & file.exists(jsons$genes.filename))

  # if we filter too much we get an error
  expect_error(gtf2json(gr = gr, chrom.sizes = chrom.sizes, include.chr = c('chr2'), genes = 'MIR137HG'))

  jsons = gtf2json(gr = gr, chrom.sizes = chrom.sizes, grep = 'MIR', grepe = 'MYR')
  expect_true(file.exists(jsons$metadata.filename) & file.exists(jsons$genes.filename))

  expect_error(gtf2json(gr = gr, chrom.sizes = chrom.sizes, grep = 'no_such_gene'))

  jsons =gtf2json(verbose = TRUE, chrom.sizes = chrom.sizes, gr = gr, gene.collapse = FALSE)
  expect_true(file.exists(jsons$metadata.filename) & file.exists(jsons$genes.filename))

  gene_weights = data.table(V1 = 'MIR137HG', V2 = 1e3)
  jsons = gtf2json(verbose = TRUE, chrom.sizes = chrom.sizes, gr = gr, grep = 'MIR', gene_weights = gene_weights)
  expect_true(file.exists(jsons$metadata.filename) & file.exists(jsons$genes.filename))

  expect_error(gtf2json(verbose = TRUE, chrom.sizes = chrom.sizes, gr = gr,
                        gene_weights = 'not.a.data.table'))
  expect_error(gtf2json(verbose = TRUE, chrom.sizes = chrom.sizes, gr = gr,
                        gene_weights = data.table(V1 = 'MIR1302-11', V2 = 1e3,
                                                  too_many_columns = 'error')))
  expect_error(gtf2json(verbose = TRUE, chrom.sizes = chrom.sizes, gr = gr,
                        gene_weights = data.table(V1 = 'MIR1302-11', V2 = 1e3,
                                                  too_many_columns = 'error')))

})

tmpdir = tempdir()
message('Running jsUtils tests in the temp dir: ', tmpdir)
gGnome.js.path = paste0(tmpdir, '/gGnome.js')
PGV.path = paste0(tmpdir, '/PGV')

message("Reading test coverage")
cov.fn = system.file("extdata/", "coverage.5k.txt", package = "gGnome")

message('Creating mock js.input.data')
# make a data.table as expected by gen_js_instance

gg.rds = paste0(tmpdir, '/gg.jabba.rds')
saveRDS(gg.jabba, gg.rds)

ncn.gr = copy(gg.jabba$gr)[,c()]
ncn.gr$ncn = 2

# make a fake karyograph
kag = list(segstats = ncn.gr)
kag.fn = paste0(tmpdir, '/fake-kag.rds')
saveRDS(kag, kag.fn)

bad.kag = list(not_segstats = '')
bad.kag.fn = paste0(tmpdir, '/bad-kag.rds')
saveRDS(bad.kag, bad.kag.fn)

js_data = data.table(sample = 'mypair', coverage = cov.fn, graph = gg.rds, kag = kag.fn, bad.kag = bad.kag.fn)
test_that('gen_js_instance', {

    system(paste0('rm -rf ', gGnome.js.path))
    gGnome.js(js_data,
                    outdir = paste0(tmpdir, '/gGnome.js'),
                    ncn.gr = ncn.gr,
                    annotation = NULL)

    system(paste0('rm -rf ', paste0(tmpdir, '/gGnome.js2')))
    gGnome.js(js_data,
                    outdir = paste0(tmpdir, '/gGnome.js2'),
                    reference = system.file('extdata/jsUtils', 'mock_ref_dir', package="gGnome"),
                    kag.col = 'bad.kag', # provide a bad karyograph file that does not contain segstats. Things should still work, but a warning is expected
                    annotation = NULL)

    # bad ref name
    system(paste0('rm -rf ', paste0(tmpdir, '/gGnome.js3')))
    expect_error(gGnome.js(js_data,
                    outdir = paste0(tmpdir, '/gGnome.js3'),
                    reference = 'no_such_ref',
                    annotation = NULL))

    # bad ref dir
    bad_ref_dir = paste0(tmpdir, '/bad_ref_dir')
    dir.create(bad_ref_dir, showWarnings = FALSE)
    system(paste0('rm -rf ', paste0(tmpdir, '/gGnome.js3')))
    expect_error(gGnome.js(js_data,
                    outdir = paste0(tmpdir, '/gGnome.js3'),
                    reference = bad_ref_dir,
                    annotation = NULL))

    expect_error(gGnome.js(data.table(sample = 'mypair2', coverage = cov.fn, graph = gg.rds),
                    outdir = paste0(tmpdir, '/gGnome.js'),
                    cov.field ='no.such.field',
                    append = TRUE,
                    annotation = NULL))

    system(paste0('rm -rf ', PGV.path))
    pgv(js_data,
                outdir = paste0(tmpdir, '/PGV'),
                ref = 'hg19',
                annotation = NULL,
                ncn.gr = ncn.gr,
                dataset_name = 'test',
                )

    # test what happens if an existing directory is provided but append is set to FALSE
    expect_error(pgv(data.table(sample = c('mypair', 'mypair2'),
                                  coverage = c(cov.fn, cov.fn),
                                  graph = c(gg.rds, gg.rds)),
                    outdir = paste0(tmpdir, '/PGV'),
                    ref = 'hg19',
                    cov.field = NA,
                    append = FALSE,
                    annotation = NULL,
                    dataset_name = 'test',
                    ))

    # test adding more data and also test what happens when the cov.field is NA (we expect a warning about skipping the coverage generation)
    expect_warning(pgv(data.table(sample = c('mypair', 'mypair2'),
                                  coverage = c(cov.fn, cov.fn),
                                  graph = c(gg.rds, gg.rds)),
                    outdir = paste0(tmpdir, '/PGV'),
                    ref = 'hg19',
                    cov.field = NA,
                    append = TRUE,
                    annotation = NULL,
                    dataset_name = 'test',
                    ))

    # re-run without adding a new samples, but with adding tree (so should just update the datafiles)
    # also adding connections.associations
    pgv(data.table(sample = 'mypair2', coverage = cov.fn, graph = gg.rds),
                    outdir = paste0(tmpdir, '/PGV'),
                    ref = 'hg19',
                    cov.field = NA,
                    append = TRUE,
                    annotation = NULL,
                    tree = system.file('extdata', 'phylogeny-pgv.newick', package="gGnome"),
                    connections.associations = TRUE,
                    dataset_name = 'test'
                    )

    # tree is not newick
    expect_warning(pgv(data.table(sample = 'mypair2', coverage = cov.fn, graph = gg.rds),
                    outdir = paste0(tmpdir, '/PGV'),
                    ref = 'hg19',
                    cov.field = NA,
                    append = TRUE,
                    annotation = NULL,
                    tree = '/dev/null',
                    connections.associations = TRUE,
                    dataset_name = 'test'
                    ))

    # test bad tree file. does not exist. should get a warning
    expect_warning(pgv(data.table(sample = 'mypair2', coverage = cov.fn, graph = gg.rds),
                    outdir = paste0(tmpdir, '/PGV'),
                    ref = 'hg19',
                    cov.field = NA,
                    append = TRUE,
                    annotation = NULL,
                    tree = 'nofilehere',
                    dataset_name = 'test'
                    ))

    expect_error(is.dir.a.PGV.instance('/dev/null'))
    expect_error(is.dir.a.gGnome.js.instance('/dev/null'))

    # provide an invalid ref
    expect_error(pgv(data.table(sample = 'mypair2', coverage = cov.fn, graph = gg.rds),
                    outdir = paste0(tmpdir, '/PGV_errored'),
                    ref = 'test_ref_error',
                    cov.field = NA,
                    annotation = NULL,
                    tree = 'nofilehere',
                    dataset_name = 'test'
                    ))
                    

})


test_that('cov2csv', {

    tmp.fn = tempfile()
    cov.csv = cov2csv(cov = cov.fn, field = 'ratio', output_file = tmp.fn, bin.width = NA)

    expect_true(file.exists(cov.csv))
    expect_error(cov2csv("non.existent"))
    expect_error(cov2csv(cov.fn, 'missing.field'))
})

test_that('cov2arrow', {

    tmp.fn = tempfile()
    cov.arrow = cov2arrow(cov = cov.fn, field = 'ratio', output_file = tmp.fn)

    expect_true(all(file.exists(cov.arrow)))
    expect_error(cov2arrow("non.existent"))
    expect_error(cov2arrow(cov.arrow)) # providing an invalid (existing) file as input
})

message('Download js metadata files')
gGnome.js.meta = tempfile()
# notice that if the repository structure will change then we will need to update this 
download.file('https://raw.githubusercontent.com/mskilab/gGnome.js/master/public/metadata.json', gGnome.js.meta)

PGV.meta = tempfile()
# notice that if the repository structure will change then we will need to update this 
download.file('https://raw.githubusercontent.com/mskilab/pgv/main/public/settings.json', PGV.meta)

test_that('parse.js.seqlenghts', {

    expect_is(parse.js.seqlenghts(meta.js = gGnome.js.meta, js.type = 'gGnome.js'), 'integer')
    expect_error(parse.js.seqlenghts(meta.js = gGnome.js.meta, js.type = 'PGV'))
    
    expect_is(parse.js.seqlenghts(meta.js = PGV.meta, js.type = 'PGV', ref = 'hg19'), 'integer')
    expect_error(parse.js.seqlenghts(meta.js = PGV.meta, js.type = 'gGnome.js'))
    expect_error(parse.js.seqlenghts(meta.js = PGV.meta, js.type = 'PGV'))
    expect_error(parse.js.seqlenghts(meta.js = PGV.meta, js.type = '/dev/null'))
})

test_that('cov2cov.js', {
    cov = readCov(cov.fn)
    cov_chr = gr.chr(cov)
    cov_combined = suppressWarnings(c(cov, cov_chr))
    expect_is(cov2cov.js(cov, meta.js = gGnome.js.meta), 'data.table')
    expect_error(cov2cov.js(cov_chr, meta.js = gGnome.js.meta))
    expect_error(cov2cov.js(cov, meta.js = PGV.meta, js.type = 'PGV'))
    expect_is(cov2cov.js(cov, meta.js = PGV.meta, js.type = 'PGV', ref = 'hg19'), 'data.table')
    expect_is(cov2cov.js(cov_chr, meta.js = PGV.meta, js.type = 'PGV', ref = 'hg38_chr'), 'data.table')
    expect_warning(cov2cov.js(cov_combined, meta.js = gGnome.js.meta))
})

test_that('js_path', {
    # now this directory already exists so we expect an error
    expect_error(js_path(PGV.path, js.type = 'PGV'))

    # providing a file name
        expect_error(js_path('/dev/null', js.type = 'PGV'))
})

test_that('get_path_to_meta_js', {

    expect_error(get_path_to_meta_js('non.existing.directory', js.type = 'gGnome.js'))

    meta.js.gGnome.js = get_path_to_meta_js(gGnome.js.path, js.type = 'gGnome.js')
    expect_true(file.exists(meta.js.gGnome.js))

    meta.js.PGV = get_path_to_meta_js(PGV.path, js.type = 'PGV')
    expect_true(file.exists(meta.js.PGV))

    expect_error(get_path_to_meta_js(PGV.path, js.type = 'noSuchJs'))

    # giving the wrong path so expecting not to find the metadata file
    expect_error(js_path(gGnome.js.path, js.type = 'PGV'))
    expect_error(js_path(PGV.path, js.type = 'gGnome.js'))
})

test_that('read.js.input.data', {
    expect_error(read.js.input.data(list()))
    expect_error(read.js.input.data('no.such.file'))

    js_test_data = read.js.input.data(js_data)
    expect_is(js_test_data, 'data.table')
})

test_that('get_js_cov_dir_path', {
    expect_true(dir.exists(get_js_cov_dir_path(gGnome.js.path, js.type = 'gGnome.js')))
    expect_error(get_js_cov_dir_path(PGV.path, js.type = 'PGV'))
    expect_true(dir.exists(get_js_cov_dir_path(PGV.path, js.type = 'PGV', dataset_name = 'test')))
})

test_that('gen_js_coverage_files', {

    expect_error(gen_js_coverage_files(js_data, outdir = gGnome.js.path, cov.col = 'no.such.column'))
    expect_error(gen_js_coverage_files(js_data, outdir = gGnome.js.path, cov.field.col = 'no.such.column'))

    fn =  gen_js_coverage_files(js_data, outdir = gGnome.js.path,
                          js.type = 'gGnome.js') 
    print(fn)
    expect_true(all(file.exists(unlist(fn))))

    fn =  gen_js_coverage_files(js_data, outdir = PGV.path,
                          js.type = 'PGV', dataset_name = 'test') 
    print(fn)
    expect_true(all(file.exists(unlist(fn))))

})

test_that('gen_gg_json_files', {

    fn =  gen_gg_json_files(js_data, outdir = gGnome.js.path,
                          gg.col = 'graph', js.type = 'gGnome.js') 
    print(fn)
    expect_true(all(file.exists(unlist(fn))))

    fn =  gen_gg_json_files(js_data, outdir = PGV.path,
                          gg.col = 'graph', js.type = 'PGV', dataset_name = 'test') 
    print(fn)
    expect_true(all(file.exists(unlist(fn))))

    # raise error when no dataset_name provided
    expect_error(fn =  gen_gg_json_files(js_data, outdir = PGV.path,
                          gg.col = 'graph', js.type = 'PGV')) 

})

test_that('get_cids', {
   cids = get_cids(gg.jabba, cid.field = 'sedge.id') 
   expect_warning(get_cids(gg.jabba, cid.field = 'invalid_field'))
   expect_warning(get_cids(gg.jabba, cid.field = 'n1.side')) # non-numeric
   expect_warning(get_cids(gg.jabba, cid.field = 'SUBN')) # contains NAs
})
