fl <- system.file("extdata", "chr1.vcf.bgz", package="tadar")
genotypes <- readGenotypes(fl)
groups <- list(
    group1 = paste0("sample", 1:6),
    group2 = paste0("sample", 7:13)
)
counts <- countAlleles(genotypes, groups)

test_that("countAlleles returns GRangesList of same length as groups", {
    expect_equal(length(groups), length(counts))
})

test_that("countAlleles returns list elements with correct lengths", {
    expect_equal(length(counts[[1]]), length(counts[[2]]))
    expect_equal(length(counts[[1]]), length(genotypes))
})

test_that("countAlleles returns GRangesList with same names as groups", {
    expect_equal(names(counts), names(groups))
})

test_that("countAlleles returns mcol names in correct order", {
    expect_equal(
        unique(c(names(mcols(counts[[1]])), names(mcols(counts[[1]])))),
        c("n_called", "n_missing", "n_0", "n_1", "n_2", "n_3")
    )
})
