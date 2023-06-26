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
dar <- dar(props, contrasts, winSize = 5)

test_that("convertRanges returns windows", {
    darWindows <- convertRanges(dar)
    widths <- width(darWindows$group1v2)
    expect_true(max(widths) > 1)
})

test_that("convertRanges extends ranges", {
    darWindows <- convertRanges(dar, extendEdges = TRUE)
    expect_equal(min(start(darWindows$group1v2)), 1)
    expect_equal(max(end(darWindows$group1v2)), 59578282)
})

test_that("convertRanges errors when no winSize in metadata", {
    dar <- endoapply(dar, function(x){
        metadata(x) <- metadata(x)[!(names(metadata(x)) %in% "winSize")]
        x
    })
    expect_error(convertRanges(dar), "No winSize detected\\.")
})

test_that("convertRanges errors when no seqlengths and extendEdges = TRUE", {
    dar <- endoapply(dar, function(x){
        seqinfo(x) <- Seqinfo(seqnames = "1")
        x
    })
    expect_error(
        convertRanges(dar, extendEdges = TRUE),
        "Cannot extend edges\\."
    )
})

test_that("convertRanges can revert back to output of dar()", {
    darWindows <- convertRanges(dar)
    expect_identical(convertRanges(darWindows), dar)
    darWindows <- convertRanges(dar, extendEdges = TRUE)
    expect_identical(convertRanges(darWindows), dar)
})
