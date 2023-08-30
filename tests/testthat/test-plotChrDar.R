set.seed(230822)
data("chr1_genes")
foi <- sample(chr1_genes, 1)
features <- sample(chr1_genes, 20)
fl <- system.file("extdata", "chr1.vcf.bgz", package="darr")
genotypes <- readGenotypes(fl)
groups <- list(
    group1 = c("S2", "S7", "S9", "S10", "S19", "S20"),
    group2 = c("S3", "S6", "S11", "S12", "S15", "S16", "S18")
)
counts <- countAlleles(genotypes, groups)
props <- countsToProps(counts)
contrasts <- matrix(
    data = c(1, -1),
    dimnames = list(
        Levels = c("group1", "group2"),
        Contrasts = c("group1v2")
    )
)
dar <- dar(props, contrasts)$group1v2

test_that("plotChrDar_checks errors when passed incorrect arguments", {
    expect_error(plotChrDar(dar = NULL))
    expect_error(
        .plotChrDar_checks(dar, foi = NULL),
        "`foi` must be a GRanges object"
    )
    expect_error(
        .plotChrDar_checks(dar, features = NULL),
        "`features` must be a GRanges object"
    )
    expect_error(
        .plotChrDar_checks(dar, chr = "100"),
        "Chromosome '100' not found in ranges of `dar`"
    )
    expect_error(
        .plotChrDar_checks(dar, chr = 1),
        "`chr` must be a character"
    )
    expect_error(
        .plotChrDar_checks(dar, foi = foi, foi_anno = NULL),
        "Annotation column for `foi` must be a character"
    )
    expect_error(
        .plotChrDar_checks(dar, foi = foi, foi_anno = "doesntexist"),
        "Column 'doesntexist' not found in mcols of `foi`"
    )
    expect_error(
        .plotChrDar_checks(dar, features = features, features_anno = NULL),
        "Annotation column for `features` must be a character"
    )
    expect_error(
        .plotChrDar_checks(
            dar, features = features, features_anno = "doesntexist"
        ),
        "Column 'doesntexist' not found in mcols of `features`"
    )
    expect_error(
        .plotChrDar_checks(dar, foi_highlight = NULL),
        "`foi_highlight` must be logical"
    )
    expect_error(
        .plotChrDar_checks(dar, features_highlight = NULL),
        "`features_highlight` must be logical"
    )
    dar_2seqs <- dar
    seqlevels(dar_2seqs) <- c("1", "2")
    dar_2seqs <- c(
        dar_2seqs,
        GRanges(
            seqnames = "2", ranges = IRanges(start = 1, end = 1),
            seqinfo = seqinfo(dar_2seqs)
        )
    )
    expect_error(
        .plotChrDar_checks(dar_2seqs),
        "All ranges of `dar` must exist on the same chromosome"
    )
    foi_chr2 <- foi
    seqlevels(foi_chr2) <- "2"
    expect_error(
        .plotChrDar_checks(dar, foi = foi_chr2),
        "All ranges of supplied GRanges objects must exist on the same chromosome"
    )
})

test_that("dar_track is created with correct attributes", {
    ## Data
    dar_track <- .darTrack(
        dar = dar, darVal = "region", chr = "1",
        foi_highlight = TRUE, features_highlight = TRUE
    )
    expect_true(is(dar_track@trackList[[1]], "DataTrack"))
    expect_equal(dar_track@trackList[[1]]@name, "DAR")
    expect_equal(dar_track@trackList[[1]]@dp@pars$window, -1)
    expect_equal(dar_track@trackList[[1]]@dp@pars$windowSize, 1)
    expect_equal(dar_track@trackList[[1]]@dp@pars$cex, 0.4)
    expect_equal(dar_track@trackList[[1]]@dp@pars$col, "grey20")
    expect_equal(dar_track@trackList[[1]]@dp@pars$col.axis, "black")
    expect_equal(dar_track@trackList[[1]]@dp@pars$yTicksAt, seq(0, 1, 0.1))
    expect_equal(dar_track@trackList[[1]]@dp@pars$ylim, c(0, 1))
    expect_equal(dar_track@trackList[[1]]@dp@pars$rotation.title, 90)
    ## Highlight
    dar_track <- .darTrack(
        dar = dar, darVal = "region", chr = "1",
        foi = foi, foi_highlight = TRUE,
        features = features, features_highlight = TRUE
    )
    expect_true(is(dar_track, "HighlightTrack"))
    expect_equal(
        dar_track@dp@pars$col,
        c(rep("red", length(foi)), rep("#ffe0e0", length(features)))
    )
    expect_equal(
        dar_track@dp@pars$fill,
        c(rep("red", length(foi)), rep("#ffe0e0", length(features)))
    )
})

test_that("foi_track is created with correct attributes", {
    expect_null(.foiTrack())
    foi_track <- .foiTrack(chr = "1", foi = foi, foi_anno = "gene_name")
    expect_true(is(foi_track, "AnnotationTrack"))
    expect_equal(foi_track@name, "foi")
    expect_equal(foi_track@dp@pars$shape, "box")
    expect_equal(foi_track@dp@pars$col, "white")
    expect_equal(foi_track@dp@pars$fill, "white")
    expect_equal(foi_track@dp@pars$groupAnnotation, "group")
    expect_equal(foi_track@dp@pars$fontcolor.group, "red")
    expect_equal(foi_track@dp@pars$cex.group, 0.6)
    expect_equal(foi_track@dp@pars$just.group, "below")
})

test_that("axis_track is created with correct attributes", {
    ## Axis
    axis_track <- .axisTrack(chr = "1")
    expect_true(is(axis_track, "GenomeAxisTrack"))
    expect_equal(axis_track@name, "Chr1")
    expect_equal(axis_track@dp@pars$rotation.title, 0)
    expect_true(axis_track@dp@pars$add35)
    expect_true(axis_track@dp@pars$add53)
    ## Highlight
    axis_track <- .axisTrack(chr = "1", foi = foi)
    expect_true(is(axis_track, "HighlightTrack"))
    expect_true(is(axis_track@trackList[[1]], "GenomeAxisTrack"))
    expect_equal(axis_track@trackList[[1]]@name, "Chr1")
    expect_equal(axis_track@trackList[[1]]@dp@pars$rotation.title, 0)
    expect_true(axis_track@trackList[[1]]@dp@pars$add35)
    expect_true(axis_track@trackList[[1]]@dp@pars$add53)
})

test_that("features_track is created with correct attributes", {
    expect_null(.featuresTrack())
    features_track <- .featuresTrack(
        chr = "1", features = features, features_anno = "gene_name"
    )
    expect_true(is(features_track, "AnnotationTrack"))
    expect_equal(features_track@name, "Features")
    expect_equal(features_track@dp@pars$shape, "box")
    expect_equal(features_track@dp@pars$col, "darkgray")
    expect_equal(features_track@dp@pars$fill, "darkgray")
    expect_equal(features_track@dp@pars$groupAnnotation, "group")
    expect_equal(features_track@dp@pars$fontcolor.group, "black")
    expect_equal(features_track@dp@pars$cex.group, 0.6)
    expect_equal(features_track@dp@pars$rotation.title, 0)
})

test_that("Final plot is created as expected", {
    plot <- plotChrDar(
        dar = dar, darVal = "region", chr = "1",
        foi = foi, foi_anno = "gene_name", foi_highlight = TRUE,
        features = features, features_anno = "gene_name",
        features_highlight = TRUE,
        title = "Test plot"
    )
    expect_true(is(plot$foi, "AnnotationTrack"))
    expect_true(is(plot$Chr1, "GenomeAxisTrack"))
    expect_true(is(plot$Features, "AnnotationTrack"))
    expect_true(is(plot$DAR, "DataTrack"))
})
