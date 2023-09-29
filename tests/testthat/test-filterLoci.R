fl <- system.file("extdata", "chr1.vcf.bgz", package="tadar")
genotypes <- readGenotypes(fl)
groups <- list(
    group1 = paste0("sample", 1:6),
    group2 = paste0("sample", 7:13)
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
