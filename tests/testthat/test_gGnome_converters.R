library(testthat)
library(gUtils)

setDTthreads(1)

test_that(desc = "test read.juncs for svaba vcf",
          code = {
              fn = system.file("extdata", "svaba.subset.vcf", package = "gGnome")
              svaba.all = gGnome:::read.juncs(fn)
              svaba.std = gGnome:::read.juncs(fn, standard.only = TRUE)
              svaba.chr = gGnome:::read.juncs(fn, standard.only = TRUE, chr.convert = TRUE)
              ## make sure correct number of junctions is read in
              expect_true(length(svaba.all) == 20)
              expect_true(length(svaba.std) == 20)
              expect_true(length(svaba.chr) == 20)
              ## make sure seqlengths only includes standard chr
              expect_true(all(names(seqlengths(svaba.std)) %in% paste0("chr", c(as.character(1:22), "X", "Y"))))
              ## make sure chr prefex is removed
              expect_true(!any(grepl("chr", names(seqlengths(svaba.chr)))))
          })

test_that(desc = "test ONT VCF",
          code = {
              fn = system.file("extdata", "ont.subset.vcf", package = "gGnome")
              ont.all = gGnome:::read.juncs(fn)
              ont.std = gGnome:::read.juncs(fn, standard.only = TRUE)
              expect_true(length(ont.all) == 25)
              expect_true(length(ont.std) < 25)
              expect_true(all(names(seqlengths(ont.std)) %in% c(as.character(1:22), "X", "Y")))
              bp1 = grl.pivot(ont.std)[[1]]
              bp2 = grl.pivot(ont.std)[[2]]
              expect_true(all(as.character(seqnames(bp1)) %in% c(as.character(1:22), "X", "Y")))
              expect_true(all(as.character(seqnames(bp2)) %in% c(as.character(1:22), "X", "Y")))
          })

## Keep synthetic VCFs local and remove them even when an expectation fails.
with.junction.vcf = function(rows, code, end.header = NULL,
                             mate.header = TRUE, parid.header = FALSE,
                             include.samples = TRUE) {
    fn = tempfile(fileext = ".vcf")
    on.exit(unlink(fn), add = TRUE)
    header = c('##fileformat=VCFv4.2',
               '##contig=<ID=chr1,length=10000>',
               '##contig=<ID=chr2,length=20000>',
               '##contig=<ID=chrUn,length=30000>',
               '##FILTER=<ID=lowSupport,Description="Synthetic low support">',
               '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="SV type">',
               '##INFO=<ID=LABEL,Number=1,Type=String,Description="Original row label">')
    if (!is.null(end.header)) {
        header = c(header, end.header)
    }
    if (mate.header) {
        header = c(header, '##INFO=<ID=MATEID,Number=.,Type=String,Description="Mate ID">')
    }
    if (parid.header) {
        header = c(header, '##INFO=<ID=PARID,Number=1,Type=String,Description="Partner ID">')
    }
    columns = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"
    if (include.samples) {
        header = c(header,
                   '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
                   '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">')
        columns = paste(columns, "FORMAT", "NORMAL", "TUMOR", sep = "\t")
    }
    writeLines(c(header, columns, rows), fn)
    code(fn)
}

junction.vcf.rows = function(id, pos, alt, svtype, chr = "chr1",
                             mate = NULL, end = NULL, parid = NULL,
                             filter = "PASS", normal = "0/0:10", tumor = "0/1:20") {
    info = paste0("SVTYPE=", svtype, ";LABEL=", id)
    if (!is.null(mate)) {
        info = paste0(info, ifelse(is.na(mate), "", paste0(";MATEID=", mate)))
    }
    if (!is.null(end)) {
        info = paste0(info, ifelse(is.na(end), "", paste0(";END=", end)))
    }
    if (!is.null(parid)) {
        info = paste0(info, ifelse(is.na(parid), "", paste0(";PARID=", parid)))
    }
    paste(chr, pos, id, "N", alt, "60", filter, info,
          "GT:DP", normal, tumor, sep = "\t")
}

expect.junction.coordinates = function(juncs, chr1, pos1, strand1, chr2, pos2, strand2) {
    expect_s4_class(juncs, "GRangesList")
    expect_length(juncs, length(pos1))
    expect_equal(as.integer(elementNROWS(juncs)), rep(2L, length(pos1)))
    bp = grl.pivot(juncs)
    expect_equal(as.character(seqnames(bp[[1]])), chr1)
    expect_equal(start(bp[[1]]), pos1)
    expect_equal(end(bp[[1]]), pos1)
    expect_equal(as.character(strand(bp[[1]])), strand1)
    expect_equal(as.character(seqnames(bp[[2]])), chr2)
    expect_equal(start(bp[[2]]), pos2)
    expect_equal(end(bp[[2]]), pos2)
    expect_equal(as.character(strand(bp[[2]])), strand2)
}

expect.junction.metadata = function(juncs, id, svtype, filter = "PASS",
                                   normal = "0/0", tumor = "0/1",
                                   normal.dp = 10L, tumor.dp = 20L) {
    meta = mcols(juncs)
    expect_equal(as.character(meta$LABEL), id)
    expect_equal(as.character(meta$SVTYPE), svtype)
    expect_equal(as.character(meta$FILTER), rep_len(filter, length(id)))
    expect_equal(as.character(meta$GT__S1), rep_len(normal, length(id)))
    expect_equal(as.character(meta$GT__S2), rep_len(tumor, length(id)))
    expect_equal(as.integer(meta$DP__S1), rep_len(normal.dp, length(id)))
    expect_equal(as.integer(meta$DP__S2), rep_len(tumor.dp, length(id)))
}

test_that("read.juncs routes typed paired ALT without END and retains original metadata", {
    ids = c("del-a", "del-b", "dup-a", "dup-b", "inv-a", "inv-b",
            "bnd-a", "bnd-b", "tra-a", "tra-b")
    types = rep(c("DEL", "DUP", "INV", "BND", "TRA"), each = 2)
    alts = c("N[chr1:200[", "]chr1:100]N", "]chr1:400]N", "N[chr1:300[",
             "N]chr2:600]", "N]chr1:500]", "[chr2:800[N", "[chr1:700[N",
             "N[chr2:1000[", "]chr1:900]N")
    rows = junction.vcf.rows(ids, seq(100, 1000, by = 100), alts, types,
                             chr = c(rep("chr1", 5), "chr2", "chr1", "chr2", "chr1", "chr2"),
                             mate = ids[c(2, 1, 4, 3, 6, 5, 8, 7, 10, 9)],
                             filter = rep(c("PASS", "lowSupport", "PASS", "lowSupport", "PASS"), each = 2),
                             normal = rep(c("0/0:11", "0/1:99"), 5),
                             tumor = rep(c("0/1:21", "1/1:88"), 5))
    with.junction.vcf(rows, function(fn) {
        juncs = read.juncs(fn)
        expect.junction.coordinates(juncs, rep("chr1", 5), c(100, 300, 500, 700, 900),
                                    c("-", "+", "-", "+", "-"),
                                    c("chr1", "chr1", "chr2", "chr2", "chr2"),
                                    c(200, 400, 600, 800, 1000), c("+", "-", "-", "+", "+"))
        expect.junction.metadata(juncs, ids[seq(1, 9, by = 2)],
                                 c("DEL", "DUP", "INV", "BND", "TRA"),
                                 c("PASS", "lowSupport", "PASS", "lowSupport", "PASS"),
                                 normal.dp = 11L, tumor.dp = 21L)
        expect_equal(as.character(unlist(mcols(juncs)$ALT)), alts[seq(1, 9, by = 2)])
        expect_equal(as.character(unlist(mcols(juncs)$MATEID)), ids[seq(2, 10, by = 2)])
    })
})

test_that("read.juncs mixes symbolic END events with paired DEL and INV", {
    ids = c("symbolic-del", "symbolic-dup", "symbolic-inv",
            "paired-del-a", "paired-del-b", "paired-inv-a", "paired-inv-b")
    rows = junction.vcf.rows(ids, c(1000, 1200, 1400, 2000, 2100, 2200, 2300),
                             c("<DEL>", "<DUP>", "<INV>", "N[chr1:2100[", "]chr1:2000]N",
                               "N]chr2:2300]", "N]chr1:2200]"),
                             c("DEL", "DUP", "INV", "DEL", "DEL", "INV", "INV"),
                             chr = c(rep("chr1", 6), "chr2"),
                             mate = c(NA, NA, NA, ids[5], ids[4], ids[7], ids[6]),
                             end = c(1100, 1300, 1500, 9000, NA, 9900, 9800))
    with.junction.vcf(rows, function(fn) {
        juncs = read.juncs(fn)
        expect.junction.coordinates(juncs, rep("chr1", 6),
                                    c(999, 1300, 1399, 1400, 2000, 2200),
                                    c("-", "-", "-", "+", "-", "-"),
                                    c(rep("chr1", 5), "chr2"),
                                    c(1101, 1200, 1500, 1501, 2100, 2300),
                                    c("+", "+", "-", "+", "+", "-"))
        expect.junction.metadata(juncs, ids[c(1, 2, 3, 3, 4, 6)],
                                 c("DEL", "DUP", "INV", "INV", "DEL", "INV"))
        expect_equal(as.integer(mcols(juncs)$END), c(1100L, 1300L, 1500L, 1500L, 9000L, 9900L))
    }, end.header = '##INFO=<ID=END,Number=1,Type=Integer,Description="Symbolic end">')
})

test_that("read.juncs skips only symbolic events when END is absent", {
    ids = c("no-end-del", "no-end-dup", "no-end-inv", "paired-a", "paired-b")
    rows = junction.vcf.rows(ids, c(100, 200, 300, 400, 500),
                             c("<DEL>", "<DUP>", "<INV>", "N[chr1:500[", "]chr1:400]N"),
                             c("DEL", "DUP", "INV", "DEL", "DEL"),
                             mate = c(NA, NA, NA, ids[5], ids[4]))
    with.junction.vcf(rows, function(fn) {
        juncs = suppressWarnings(read.juncs(fn))
        expect.junction.coordinates(juncs, "chr1", 400, "-", "chr1", 500, "+")
        expect.junction.metadata(juncs, "paired-a", "DEL")
    })
})

test_that("read.juncs validates missing nonnumeric and list END per symbolic row", {
    ids = c("valid-del", "missing-del", "text-dup", "multiple-inv", "zero-del",
            "negative-dup", "infinite-inv", "nan-del", "fractional-dup", "dot-inv",
            "valid-dup", "paired-a", "paired-b")
    types = c("DEL", "DEL", "DUP", "INV", "DEL", "DUP", "INV", "DEL", "DUP", "INV", "DUP", "DEL", "DEL")
    rows = junction.vcf.rows(ids, c(1000, 1200, 1300, 1400, 1500, 1600, 1700, 1800, 1900, 2000, 2100, 3000, 3100),
                             c(paste0("<", types[1:11], ">"), "N[chr1:3100[", "]chr1:3000]N"), types,
                             mate = c(rep(NA, 11), ids[13], ids[12]),
                             end = c("1100", NA, "not-an-end", "1450,1490", "0", "-1", "Inf", "NaN", "1950.5", ".", "2200", "ignored-text", NA))
    with.junction.vcf(rows, function(fn) {
        juncs = suppressWarnings(read.juncs(fn))
        expect.junction.coordinates(juncs, rep("chr1", 3), c(999, 2200, 3000), rep("-", 3),
                                    rep("chr1", 3), c(1101, 2100, 3100), rep("+", 3))
        expect.junction.metadata(juncs, ids[c(1, 11, 12)], c("DEL", "DUP", "DEL"))
        expect_equal(as.character(unlist(mcols(juncs)$END)), c("1100", "2200", "ignored-text"))
    }, end.header = '##INFO=<ID=END,Number=.,Type=String,Description="Possibly malformed symbolic end">')
})

test_that("read.juncs excludes INS SGL and malformed ALT without losing valid BND", {
    ids = c("insertion", "bracketed-ins-a", "bracketed-ins-b", "sgl-left", "sgl-right",
            "mixed-brackets", "zero-distal", "text-distal", "two-proximal",
            "multiple-distal", "symbolic-bnd", "sequence-bnd")
    rows = junction.vcf.rows(ids, seq(100, 1200, by = 100),
                             c("<INS>", "N[chr1:300[", "]chr1:200]N", "N.", ".N",
                               "N[chr2:200]", "N[chr2:0[", "N[chr2:bad[", "N[chr2:200[N",
                               "N[chr2:200[N[chr1:300[", "<BND>", "ACGT"),
                             c(rep("INS", 3), rep("SGL", 2), rep("BND", 7)))
    with.junction.vcf(rows, function(fn) {
        juncs = suppressWarnings(read.juncs(fn))
        expect_s4_class(juncs, "GRangesList")
        expect_length(juncs, 0)
    }, mate.header = FALSE)
    valid = junction.vcf.rows("valid-bnd", 1500, "N[chr2:1600[", "BND", filter = "lowSupport")
    with.junction.vcf(c(rows, valid), function(fn) {
        juncs = suppressWarnings(read.juncs(fn))
        expect.junction.coordinates(juncs, "chr1", 1500, "-", "chr2", 1600, "+")
        expect.junction.metadata(juncs, "valid-bnd", "BND", "lowSupport")
    }, mate.header = FALSE)
})

test_that("read.juncs accepts header-only VCF with and without sample columns", {
    for (samples in c(TRUE, FALSE)) {
        with.junction.vcf(character(), function(fn) {
            for (include.geno in c(TRUE, FALSE)) {
                juncs = read.juncs(fn, include.geno = include.geno)
                expect_s4_class(juncs, "GRangesList")
                expect_length(juncs, 0)
            }
        }, include.samples = samples)
    }
})

test_that("read.juncs uses original shuffled IDs and rejects inconsistent mates", {
    ids = c("right.Z9", "ignored.one", "alpha-Q", "lost-leaf", "chain.a", "chain.b", "chain.c",
            "wrong.a", "wrong.b", "orient.a", "orient.b", "left.X2", "omega-17", "to-skipped", "filtered-mate")
    rows = junction.vcf.rows(ids,
                             c(200, 50, 300, 500, 700, 800, 900, 1000, 1101, 1200, 1300, 100, 400, 1400, 1500),
                             c("]chr1:100]N", "<INS>", "N]chr2:400]", "N[chr2:600[",
                               "N[chr2:800[", "]chr1:700]N", "N[chr2:800[", "N[chr2:1100[",
                               "]chr1:1000]N", "N[chr2:1300[", "N]chr1:1200]", "N[chr2:200[",
                               "N]chr1:300]", "N[chr2:1500[", "]chr1:1400]N"),
                             c("DEL", "INS", "INV", rep("BND", 8), "DEL", "INV", "BND", "INS"),
                             chr = c("chr2", "chr1", "chr1", "chr1", "chr1", "chr2", "chr1",
                                     "chr1", "chr2", "chr1", "chr2", "chr1", "chr2", "chr1", "chr2"),
                             mate = c("left.X2", NA, "omega-17", "absent-id", "chain.b", "chain.c", "chain.b",
                                      "wrong.b", "wrong.a", "orient.b", "orient.a", "right.Z9", "alpha-Q", "filtered-mate", "to-skipped"),
                             filter = c("lowSupport", rep("PASS", 14)),
                             normal = paste0("0/0:", 100 + seq_along(ids)),
                             tumor = paste0("0/1:", 200 + seq_along(ids)))
    with.junction.vcf(rows, function(fn) {
        juncs = suppressWarnings(read.juncs(fn))
        expect.junction.coordinates(juncs, c("chr2", "chr1"), c(200, 300), c("+", "-"),
                                    c("chr1", "chr2"), c(100, 400), c("-", "-"))
        expect.junction.metadata(juncs, ids[c(1, 3)], c("DEL", "INV"), c("lowSupport", "PASS"),
                                 normal.dp = c(101L, 103L), tumor.dp = c(201L, 203L))
        expect_equal(as.character(unlist(mcols(juncs)$MATEID)), c("left.X2", "omega-17"))
    })
})

test_that("read.juncs uses PARID when MATEID is empty", {
    rows = junction.vcf.rows(c("arbitrary-one", "unrelated-two"), c(100, 200),
                             c("N]chr2:200]", "N]chr1:100]"), "BND", chr = c("chr1", "chr2"),
                             mate = c(".", "."), parid = c("unrelated-two", "arbitrary-one"))
    with.junction.vcf(rows, function(fn) {
        juncs = read.juncs(fn)
        expect.junction.coordinates(juncs, "chr1", 100, "-", "chr2", 200, "-")
        expect.junction.metadata(juncs, "arbitrary-one", "BND")
        expect_equal(as.character(mcols(juncs)$PARID), "unrelated-two")
    }, parid.header = TRUE)
})

test_that("read.juncs parses ordinary BND ALT without a MATEID header", {
    ids = c("random-A", "unrelated-Z", "singleton-9", "not-a-pair")
    rows = junction.vcf.rows(ids, c(100, 300, 500, 700),
                             c("N[chr2:200[", "]chr2:400]N", "N]chr2:600]", "[chr2:800[N"), "BND")
    with.junction.vcf(rows, function(fn) {
        juncs = read.juncs(fn)
        expect.junction.coordinates(juncs, rep("chr1", 4), c(100, 300, 500, 700), c("-", "+", "-", "+"),
                                    rep("chr2", 4), c(200, 400, 600, 800), c("+", "-", "-", "+"))
        expect.junction.metadata(juncs, ids, rep("BND", 4))
    }, mate.header = FALSE)
    with.junction.vcf(rows, function(fn) {
        juncs = suppressWarnings(read.juncs(fn))
        expect_s4_class(juncs, "GRangesList")
        expect_length(juncs, 0)
    })
})

test_that("read.juncs validates HMF suffix pairing when MATEID is empty", {
    ids = c("event.o", "bad.o", "event.h", "bad.h")
    rows = junction.vcf.rows(ids, c(100, 300, 200, 401),
                             c("N[chr2:200[", "N[chr2:400[", "]chr1:100]N", "]chr1:300]N"), "BND",
                             chr = c("chr1", "chr1", "chr2", "chr2"), mate = rep(".", 4))
    with.junction.vcf(rows, function(fn) {
        juncs = suppressWarnings(read.juncs(fn))
        expect.junction.coordinates(juncs, "chr1", 100, "-", "chr2", 200, "+")
        expect.junction.metadata(juncs, "event.o", "BND")
    })
})

test_that("read.juncs preserves chromosome strand sequence length and metadata options", {
    ids = c("standard-a", "standard-b", "nonstandard-first-a", "nonstandard-first-b",
            "nonstandard-second-a", "nonstandard-second-b")
    rows = junction.vcf.rows(ids, seq(100, 600, by = 100),
                             c("N[chr2:200[", "]chr1:100]N", "N]chr1:400]", "N]chrUn:300]",
                               "[chrUn:600[N", "[chr1:500[N"),
                             rep(c("DEL", "INV", "BND"), each = 2),
                             chr = c("chr1", "chr2", "chrUn", "chr1", "chr1", "chrUn"),
                             mate = ids[c(2, 1, 4, 3, 6, 5)], filter = "lowSupport")
    with.junction.vcf(rows, function(fn) {
        original = read.juncs(fn)
        expect.junction.coordinates(original, c("chr1", "chrUn", "chr1"), c(100, 300, 500), c("-", "-", "+"),
                                    c("chr2", "chr1", "chrUn"), c(200, 400, 600), c("+", "-", "+"))
        expect.junction.metadata(original, ids[c(1, 3, 5)], c("DEL", "INV", "BND"), "lowSupport")
        expect_equal(seqlengths(original), c(chr1 = 10000L, chr2 = 20000L, chrUn = 30000L))

        converted = read.juncs(fn, chr.convert = TRUE)
        expect.junction.coordinates(converted, c("1", "Un", "1"), c(100, 300, 500), c("-", "-", "+"),
                                    c("2", "1", "Un"), c(200, 400, 600), c("+", "-", "+"))
        expect_equal(seqlengths(converted), c("1" = 10000L, "2" = 20000L, "Un" = 30000L))
        expect_equal(mcols(converted), mcols(original))

        flipped = read.juncs(fn, flipstrand = TRUE)
        expect.junction.coordinates(flipped, c("chr1", "chrUn", "chr1"), c(100, 300, 500), c("+", "+", "-"),
                                    c("chr2", "chr1", "chrUn"), c(200, 400, 600), c("-", "+", "-"))
        expect_equal(mcols(flipped), mcols(original))

        standard = read.juncs(fn, standard.only = TRUE)
        expect.junction.coordinates(standard, "chr1", 100, "-", "chr2", 200, "+")
        expect_equal(seqlengths(standard), c(chr1 = 10000L, chr2 = 20000L))
        expect.junction.metadata(standard, "standard-a", "DEL", "lowSupport")

        lengths = c(chr1 = 12000L, chr2 = 22000L, chrUn = 32000L)
        supplied = read.juncs(fn, seqlengths = lengths)
        expect.junction.coordinates(supplied, c("chr1", "chrUn", "chr1"), c(100, 300, 500), c("-", "-", "+"),
                                    c("chr2", "chr1", "chrUn"), c(200, 400, 600), c("+", "-", "+"))
        expect_equal(seqlengths(supplied), lengths)
        expect_equal(mcols(supplied), mcols(original))

        no.geno = read.juncs(fn, include.geno = FALSE)
        expect.junction.coordinates(no.geno, c("chr1", "chrUn", "chr1"), c(100, 300, 500), c("-", "-", "+"),
                                    c("chr2", "chr1", "chrUn"), c(200, 400, 600), c("+", "-", "+"))
        geno.columns = c("GT__S1", "GT__S2", "DP__S1", "DP__S2")
        expect_false(any(geno.columns %in% names(mcols(no.geno))))
        info.columns = setdiff(names(mcols(original)), geno.columns)
        expect_equal(mcols(no.geno), mcols(original)[, info.columns, drop = FALSE])

        no.features = read.juncs(fn, keep.features = FALSE)
        expect.junction.coordinates(no.features, c("chr1", "chrUn", "chr1"), c(100, 300, 500), c("-", "-", "+"),
                                    c("chr2", "chr1", "chrUn"), c(200, 400, 600), c("+", "-", "+"))
        expect_equal(ncol(mcols(no.features)), 0L)

        combined = read.juncs(fn, chr.convert = TRUE, standard.only = TRUE, flipstrand = TRUE,
                             include.geno = FALSE, keep.features = FALSE)
        expect.junction.coordinates(combined, "1", 100, "+", "2", 200, "-")
        expect_equal(seqlengths(combined), c("1" = 10000L, "2" = 20000L))
        expect_equal(ncol(mcols(combined)), 0L)
    })
})

test_that("read.juncs keeps valid scalar END rows around missing and nonnumeric END", {
    ids = c("before-gap", "missing-end", "nonnumeric-end", "after-gap", "paired-first", "paired-last")
    rows = junction.vcf.rows(ids, c(100, 300, 500, 700, 900, 1000),
                             c("<DEL>", "<DUP>", "<INV>", "<DUP>", "N[chr2:1000[", "]chr1:900]N"),
                             c("DEL", "DUP", "INV", "DUP", "DEL", "DEL"),
                             chr = c(rep("chr1", 5), "chr2"),
                             mate = c(rep(NA, 4), ids[6], ids[5]),
                             end = c("200", NA, "unavailable", "800", NA, NA))
    with.junction.vcf(rows, function(fn) {
        juncs = suppressWarnings(read.juncs(fn))
        expect.junction.coordinates(juncs, rep("chr1", 3), c(99, 800, 900), rep("-", 3),
                                    c("chr1", "chr1", "chr2"), c(201, 700, 1000), rep("+", 3))
        expect.junction.metadata(juncs, ids[c(1, 4, 5)], c("DEL", "DUP", "DEL"))
    }, end.header = '##INFO=<ID=END,Number=1,Type=String,Description="Possibly missing symbolic end">')
})

test_that("read.juncs returns empty when every explicit mate is invalid", {
    rows = junction.vcf.rows(c("orphan", "wrong-left", "wrong-right"), c(100, 300, 401),
                             c("N[chr2:200[", "N[chr2:400[", "]chr1:300]N"), "BND",
                             chr = c("chr1", "chr1", "chr2"),
                             mate = c("absent", "wrong-right", "wrong-left"))
    with.junction.vcf(rows, function(fn) {
        juncs = suppressWarnings(read.juncs(fn))
        expect_s4_class(juncs, "GRangesList")
        expect_length(juncs, 0)
    })
})

test_that(desc = "test reading Junction object",
          code = {
              fn = system.file("extdata", "svaba.subset.junction.rds", package = "gGnome")
              svaba.all = gGnome:::read.juncs(fn)
              expect_true(length(svaba.all) == 20)
          })

test_that(desc = "test reading bedpe",
          code = {
              fn = system.file("extdata", "ont.subset.bedpe", package = "gGnome")
              ont.bpe = gGnome:::read.juncs(fn)
              ont.std = gGnome:::read.juncs(fn, standard.only = TRUE)
              expect_true(length(ont.bpe) == 25)
              expect_true(length(ont.std) < 25)
          })


test_that('jab2gg', {
    # test allelic annotation of jabba
    jabba = readRDS(system.file('extdata/hcc1954', 'jabba.rds', package="gGnome"))
    ccols = c("cn", "edges.in", "edges.out", "tile.id")
    # random values for cn.low and cn.high
    asegstats.low = gr2dt(jabba$segstats[, ccols])
    asegstats.low[, parent := .I]
    asegstats.high = asegstats.low
    asegstats.low[!is.na(cn), cn := floor(runif(1, 0, pmax(0, cn/2 - 1)))]
    asegstats.high$cn = jabba$segstats$cn - asegstats.low$cn
    asegstats.low[, type := 'low']
    asegstats.high[, type := 'high']
    # notice that this fake asegstats does not have all the fields that a real asegstats will have (for example "phased" is missing)
    jabba$asegstats = c(dt2gr(asegstats.low), dt2gr(asegstats.high))
    gg = gG(jabba = jabba)
    expect_true(gg$nodes$dt[!is.na(cn), all(cn == cn.low + cn.high)])
})

test_that('jab2gg with purity', {
    # test allelic annotation of jabba
    jabba = gG(jabba = system.file("extdata", "hcc1954.rigma.sg.with.pp.rds", package = "gGnome"))
    expect_true(!is.null(jabba$meta$purity))
    expect_true(!is.null(jabba$meta$ploidy))
    expect_true(!is.na(jabba$meta$purity))
    expect_true(!is.na(jabba$meta$ploidy))
})

test_that('cougar1gg', {
    # download cougar files from cougar repo

    system('wget https://github.com/compbio-UofT/CouGaR-viz/raw/master/TCGA.samples.all_Contigs.tar.gz')
    system('tar -zxf TCGA.samples.all_Contigs.tar.gz TCGA-LN-A4A3/')

    chrom.sizes = system.file("extdata", "human_g1k_v37.regular.chrom.sizes", package="gGnome")
    Sys.setenv(DEFAULT_GENOME = chrom.sizes)
    cougar_list = cougar2gg('TCGA-LN-A4A3')
    cougar_gg = gG(breaks = cougar_list$breaks, junctions = cougar_list$juncs)
    expect_true(inherits(cougar_gg, 'gGraph'))

    expect_error(cougar2gg('not a dir'))
    expect_error(cougar2gg('./')) # not a cougar output directory
})

test_that('rck2gg', {
    # download cougar files from cougar repo

    rck.tar = system.file('extdata', "rck.hcc1954.tar.gz", package = "gGnome")
    system(paste0('tar -xzvf ', rck.tar))
    rck.dir = 'rck.hcc1954'

    haploid.res = rck2gg(rck.dir, haploid = TRUE)
    haploid.gg = gG(nodes = haploid.res$nodes, edges = haploid.res$edges)
    expect_true(length(haploid.res$edges) == 10)
    expect_true(inherits(haploid.gg, 'gGraph'))

    haploid.res = rck2gg(rck.dir, haploid = FALSE)
    haploid.gg.not.haploid = gG(nodes = haploid.res$nodes, edges = haploid.res$edges)
    expect_true(length(haploid.res$edges) == 14)
    expect_true(inherits(haploid.gg.not.haploid, 'gGraph'))

    expect_error(rck2gg('not a dir'))
    expect_error(rck2gg('./')) # not a cougar output directory

    # make bad rck files
    tmp = paste0(tempdir(), '/bad.rck')
    dir.create(tmp)
    fwrite(data.table(V1 = 'bad', V2 = 'not better'), paste0(tmp, '/rck.scnt.tsv'))
    fwrite(data.table(V1 = 'bad', V2 = 'not better'), paste0(tmp, '/rck.acnt.tsv'))
    expect_error(rck2gg(tmp))

    # good scnt file, bad acnt file
    good_scnt = paste0(rck.dir, '/rck.scnt.tsv')
    system(paste0('cp ', good_scnt, ' ', tmp))
    expect_error(rck2gg(tmp))
})


test_that('haplograph', {
    setDTthreads(1)
    bionano.fn = system.file('extdata', 'H838_rare_variant_pipeline_EXP_REFINEFINAL1_bionano.xmap.gz', package = 'gGnome')
    # use utils.R read_xmap to read the xmap to GRanges and generate a gWalk
    bgw = bionano.fn %>% read_xmap(merge = FALSE) %>% gW(grl = .)
    # generate haplograph from the walks
    bgg = gG(walks = bgw)
    expect_true(inherits(bgg, 'gGraph'))
})

test_that('jab2gg', {
    setDTthreads(1)
    # almost everything is tested elsewhere. just testing some edge cases here
    expect_error(jab2gg(list()))
    expect_error(jab2gg('no.such.file.rds'))
    expect_error(jab2gg(123))
    empty_gg = jab2gg(gG()) # this should work
})

test_that('read.juncs', {
    setDTthreads(1)
    # almost everything is tested elsewhere. just testing some edge cases here
    expect_error(read.juncs(NA))
    expect_error(jab2gg('no.such.file.rds'))
    expect_error(jab2gg(123))
    empty_gg = jab2gg(gG()) # this should work
})

test_that('karyotype', {
    setDTthreads(1)
    # almost everything is tested elsewhere. just testing some edge cases here
    expect_true(inherits(karyotype(), 'gGraph'))
    expect_error(karyotype(karyo = 'not null')) # since this is not implemented yet we expect error
})

test_that('alignments2gg', {
    setDTthreads(1)
    # download lightweight (~40Mb) hg19 BAM file
    system('wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeUwRepliSeq/wgEncodeUwRepliSeqK562G1AlnRep1.bam')
    system('wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeUwRepliSeq/wgEncodeUwRepliSeqK562G1AlnRep1.bam.bai')
    reads = bamUtils::read.bam('wgEncodeUwRepliSeqK562G1AlnRep1.bam', intervals = 'chr1')
    expect_true(inherits(gG(alignments = reads), 'gGraph'))
    # cleanup
    system('rm wgEncodeUwRepliSeqK562G1AlnRep1.bam*')
})
