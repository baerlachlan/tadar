fl <- system.file("extdata", "chr1.vcf.bgz", package="darr")

test_that("readGenotypes returns the expected names", {
    mcolNames <- c("S2", "S7", "S9", "S10", "S19", "S20",
      "S3", "S6", "S11", "S12", "S15", "S16", "S18")
    geno_fl <- readGenotypes(fl)
    expect_equal(names(geno_fl), NULL)
    expect_equal(names(mcols(geno_fl)), mcolNames)
    tbx <- TabixFile(fl)
    geno_tbx <- readGenotypes(tbx)
    expect_equal(names(geno_tbx), NULL)
    expect_equal(names(mcols(geno_tbx)), mcolNames)
})

test_that("readGenotypes errors when expected", {
    expect_error(readGenotypes(fl, unphase = c()), 'is.logical(.+) is not TRUE')
    expect_error(
        readGenotypes(fl, param = ScanVcfParam(geno = NA)),
        'any\\(svpChecks\\) is not TRUE')
})
