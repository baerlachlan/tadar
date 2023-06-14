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
dar <- dar(props, contrasts)

test_that("smoothDar adds an additional mcol", {
    darSmooth <- smoothDar(dar, winSize = 5)
    lenSmooth <- length(mcols(darSmooth[[1]]))
    lenDar <- length(mcols(dar[[1]]))
    expect_equal(lenSmooth, lenDar + 1)
})

test_that("smoothDar errors with incorrect winSize", {
    expect_error(
        smoothDar(dar, winSize = 0),
        "`winSize` must be an odd integer greater than 0"
    )
    expect_error(
        smoothDar(dar, winSize = 10),
        "`winSize` must be an odd integer greater than 0"
    )
    expect_error(
        smoothDar(dar, winSize = 10371),
        "`winSize` greater than number of ranges"
    )
})

test_that("smoothDar adds winSize to metadata", {
    winSize <- 5
    darSmooth <- smoothDar(dar, winSize = winSize)
    winMeta <- metadata(darSmooth[[1]])$winSize
    expect_equal(winSize, winMeta)
})
