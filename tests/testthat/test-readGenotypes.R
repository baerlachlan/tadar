fl <- system.file("extdata", "chr22.vcf.gz", package="VariantAnnotation")

test_that("readGenotypes returns the expected names", {
    genos <- readGenotypes(fl)
    expect_equal(names(genos), NULL)
    expect_equal(
        names(mcols(genos)),
        c("HG00096", "HG00097", "HG00099", "HG00100", "HG00101")
    )
})

test_that("readGenotypes errors when unphase arg is not logical", {
    expect_error(readGenotypes(fl, unphase = c()), 'is.logical(.+) is not TRUE')
})
