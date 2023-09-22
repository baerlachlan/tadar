data("chr1_genes")
fl <- system.file("extdata", "chr1.vcf.bgz", package="darr")
genotypes <- readGenotypes(fl)
groups <- list(
    group1 = paste0("sample", 1:6),
    group2 = paste0("sample", 7:13)
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
darRegions_ex <- flipRanges(dar, extend_edges = TRUE)

test_that("assignFeatureDar returns the expected output", {
    geneDar_origin <- assignFeatureDar(dar, chr1_genes, dar_val = "origin")
    expect_equal(length(geneDar_origin$group1v2), 820)
    geneDar_region <- assignFeatureDar(darRegions, chr1_genes, dar_val = "region")
    expect_equal(length(geneDar_region$group1v2), 1455)
    geneDar_ex <- assignFeatureDar(darRegions_ex, chr1_genes, dar_val = "region")
    expect_equal(length(geneDar_ex$group1v2), 1456)
})

test_that("assignFeatureDar errors when expected", {
    dar <- endoapply(dar, function(x) x[,c()])
    expect_error(
        assignFeatureDar(dar, chr1_genes, dar_val = "origin"),
        "No dar_origin values detected"
    )
    darRegions <- endoapply(darRegions, function(x) x[,c()])
    expect_error(
        assignFeatureDar(darRegions, chr1_genes, dar_val = "region"),
        "No dar_region values detected"
    )
})

test_that("assignFeatureDar gives warnings when expected", {
    expect_warning(
        assignFeatureDar(dar, chr1_genes, dar_val = "region"),
        "Range\\(s\\) detected with width == 1"
    )
    expect_warning(
        assignFeatureDar(darRegions, chr1_genes, dar_val = "origin"),
        "Range\\(s\\) detected with width > 1"
    )
})
