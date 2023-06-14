fl <- system.file("extdata", "", package="darr")
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
dar <- dar(props, contrasts)
darSmooth <- smoothDar(dar, winSize = 5)

test_that("getWinRanges returns windows", {
    darWindows <- getWinRanges(darSmooth)
    widths <- width(darWindows$group1v2)
    expect_true(max(widths) > 1)
})

test_that("getWinRanges extends ranges", {
    darWindows <- getWinRanges(darSmooth, extendEdges = TRUE)
    expect_equal(min(start(darWindows$group1v2)), 1)
    expect_equal(max(end(darWindows$group1v2)), 59578282)
})

test_that("getWinRanges errors when no winSize in metadata", {
    expect_error(getWinRanges(dar), "No winSize detected\\.")
    darSmooth <- endoapply(darSmooth, function(x){
        metadata(x) <- list()
        x
    })
    expect_error(getWinRanges(darSmooth), "No winSize detected\\.")
})

test_that("getWinRanges errors when no seqlengths and extendEdges = TRUE", {
    darSmooth <- endoapply(darSmooth, function(x){
        seqinfo(x) <- Seqinfo(seqnames = "1")
        x
    })
    expect_error(
        getWinRanges(darSmooth, extendEdges = TRUE),
        "Cannot extend edges\\."
    )
})
