fl <- system.file("extdata", "chr1.vcf.bgz", package="darr")
genotypes <- readGenotypes(fl)
groups <- list(
    group1 = paste0("sample", 1:6),
    group2 = paste0("sample", 7:13)
)
counts <- countAlleles(genotypes, groups)
counts_filt <- filterLoci(counts)
props <- countsToProps(counts_filt)
contrasts <- matrix(
    data = c(1, -1),
    dimnames = list(
        Levels = c("group1", "group2"),
        Contrasts = c("group1v2")
    )
)
dar <- dar(props, contrasts, win_size = 5)

test_that("flipRanges returns regions", {
    darRegions <- flipRanges(dar)
    widths <- width(darRegions$group1v2)
    expect_true(max(widths) > 1)
})

test_that("flipRanges extends ranges", {
    darRegions <- flipRanges(dar, extend_edges = TRUE)
    expect_equal(min(start(darRegions$group1v2)), 1)
    expect_equal(max(end(darRegions$group1v2)), 59578282)
})

test_that("flipRanges errors when no win_size in metadata", {
    dar <- endoapply(dar, function(x){
        metadata(x) <- metadata(x)[!(names(metadata(x)) %in% "win_size")]
        x
    })
    expect_error(flipRanges(dar), "No win_size detected\\.")
})

test_that("flipRanges errors when no seqlengths and extend_edges = TRUE", {
    dar <- endoapply(dar, function(x){
        seqinfo(x) <- Seqinfo(seqnames = "1")
        x
    })
    expect_error(
        flipRanges(dar, extend_edges = TRUE),
        "Cannot extend edges\\."
    )
})

test_that("flipRanges can revert back to output of dar()", {
    darRegions <- flipRanges(dar)
    expect_identical(flipRanges(darRegions), dar)
    darRegions <- flipRanges(dar, extend_edges = TRUE)
    expect_identical(flipRanges(darRegions), dar)
})
