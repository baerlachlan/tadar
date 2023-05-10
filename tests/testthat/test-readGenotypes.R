fl <- system.file("extdata", "chr1.vcf.gz", package="darr")

test_that("readGenotypes returns the expected names", {
    genotypes <- readGenotypes(fl)
    expect_equal(names(genotypes), NULL)
    expect_equal(
        names(mcols(genotypes)),
        c("S2", "S7", "S9", "S10", "S19", "S20",
          "S3", "S6", "S11", "S12", "S15", "S16", "S18")
    )
})

test_that("readGenotypes errors when unphase arg is not logical", {
    expect_error(readGenotypes(fl, unphase = c()), 'is.logical(.+) is not TRUE')
})
