data("chr1_genes")
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
darSmooth <- smoothDar(dar, winSize = 5)
darWindows <- getWinRanges(darSmooth)
darWindows_ex <- getWinRanges(darSmooth, extendEdges = TRUE)

test_that("assignFeatureDar returns the expected output", {
    geneDar_raw <- assignFeatureDar(chr1_genes, dar, darVal = "raw")
    expect_equal(length(geneDar_raw$group1v2), 820)
    geneDar_smooth <- assignFeatureDar(chr1_genes, darWindows, darVal = "smooth")
    expect_equal(length(geneDar_smooth$group1v2), 1455)
    geneDar_ex <- assignFeatureDar(chr1_genes, darWindows_ex, darVal = "smooth")
    expect_equal(length(geneDar_ex$group1v2), 1456)
})

test_that("assignFeatureDar errors whene expected", {
    dar <- endoapply(dar, function(x) x[,c()])
    expect_error(
        assignFeatureDar(chr1_genes, dar, darVal = "raw"),
        "No DAR values detected"
    )
    darWindows <- endoapply(darWindows, function(x) x[,c()])
    expect_error(
        assignFeatureDar(chr1_genes, darWindows, darVal = "smooth"),
        "No smoothed DAR values detected"
    )
})

test_that("assignFeatureDar gives warnings when expected", {
    expect_warning(
        assignFeatureDar(chr1_genes, darSmooth, darVal = "smooth"),
        "Range\\(s\\) detected with width == 1"
    )
    expect_warning(
        assignFeatureDar(chr1_genes, darWindows, darVal = "raw"),
        "Range\\(s\\) detected with width > 1"
    )
})
