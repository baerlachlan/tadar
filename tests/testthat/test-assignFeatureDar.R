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
darRegions <- flipRanges(dar)
darRegions_ex <- flipRanges(dar, extendEdges = TRUE)

test_that("assignFeatureDar returns the expected output", {
    geneDar_origin <- assignFeatureDar(chr1_genes, dar, darVal = "origin")
    expect_equal(length(geneDar_origin$group1v2), 820)
    geneDar_region <- assignFeatureDar(chr1_genes, darRegions, darVal = "region")
    expect_equal(length(geneDar_region$group1v2), 1455)
    geneDar_ex <- assignFeatureDar(chr1_genes, darRegions_ex, darVal = "region")
    expect_equal(length(geneDar_ex$group1v2), 1456)
})

test_that("assignFeatureDar errors when expected", {
    dar <- endoapply(dar, function(x) x[,c()])
    expect_error(
        assignFeatureDar(chr1_genes, dar, darVal = "origin"),
        "No dar_origin values detected"
    )
    darRegions <- endoapply(darRegions, function(x) x[,c()])
    expect_error(
        assignFeatureDar(chr1_genes, darRegions, darVal = "region"),
        "No dar_region values detected"
    )
})

test_that("assignFeatureDar gives warnings when expected", {
    expect_warning(
        assignFeatureDar(chr1_genes, dar, darVal = "region"),
        "Range\\(s\\) detected with width == 1"
    )
    expect_warning(
        assignFeatureDar(chr1_genes, darRegions, darVal = "origin"),
        "Range\\(s\\) detected with width > 1"
    )
})
