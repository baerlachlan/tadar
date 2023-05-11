fl <- system.file("extdata", "chr1.vcf.gz", package="darr")
geno <- readGenotypes(fl)
sampleGroups <- list(
    group1 = c("S2", "S7", "S9", "S10", "S19", "S20"),
    group2 = c("S3", "S6", "S11", "S12", "S15", "S16", "S18")
)
alleles <- countAlleles(geno, sampleGroups)

test_that("countAlleles returns GRangesList of same length as groups", {
    expect_equal(length(sampleGroups), length(alleles))
})

test_that("countAlleles returns list elements with correct lengths", {
    expect_equal(length(alleles[[1]]), length(alleles[[2]]))
    expect_equal(length(alleles[[1]]), length(geno))
})

test_that("countAlleles returns GRangesList with same names as groups", {
    expect_equal(names(alleles), names(sampleGroups))
})

test_that("countAlleles returns correct mcol names", {
    expect_equal(
        unique(c(names(mcols(alleles[[1]])), names(mcols(alleles[[1]])))),
        c("n_called", "n_nocall", "n_0", "n_1", "n_2", "n_3")
    )
})
