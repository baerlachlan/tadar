fl <- system.file("extdata", "chr1.vcf.gz", package="darr")
genotypes <- readGenotypes(fl)
groups <- list(
    group1 = c("S2", "S7", "S9", "S10", "S19", "S20"),
    group2 = c("S3", "S6", "S11", "S12", "S15", "S16", "S18")
)
counts <- countAlleles(genotypes, groups)

test_that("filterLoci returns expected output with default filter", {
    countsFilt <- filterLoci(counts)
    expect_equal(NROW(countsFilt$group1), 10701)
    expect_equal(NROW(countsFilt$group2), 12199)
})

test_that("filterLoci errors when metadata columns are named incorrectly", {
    testCounts <- lapply(counts, function(x) {
        mc <- mcols(x)
        names(mc) <- c("n_called", "n_missing", "n_1", "n_2", "n_3", "n_0")
        mcols(x) <- mc
        x
    })
    testCounts <- GRangesList(testCounts)
    expect_error(
        filterLoci(testCounts),
        'Names of metadata columns must equal c(.+)'
    )
})

test_that("filterLoci errors when incorrect filter expression is supplied", {
    ## Doesn't return logical vector
    expect_error(
        filterLoci(counts, n_called),
        "`filter` expression must return a logical vector of correct length"
    )
    ## Doesn't return vector of correct length
    expect_error(
        filterLoci(counts, (n_called > n_missing)[1]),
        "`filter` expression must return a logical vector of correct length"
    )
})
