data("chr1_genes")
fl <- system.file("extdata", "chr1.vcf.bgz", package="tadar")
genotypes <- readGenotypes(fl)
groups <- list(
    group1 = paste0("sample", 1:6),
    group2 = paste0("sample", 7:13)
)
counts <- countAlleles(genotypes, groups)
counts_filt <- filterLoci(counts)
props <- countsToProps(counts_filt)
contrasts <- matrix(
    data = c(1, -1),
    dimnames = list(
        Levels = c("group1", "group2"),
        Contrasts = c("group1v2")
    )
)
dar <- dar(props, contrasts, region_loci = 5)
dar_regions <- flipRanges(dar)
dar_regions_ex <- flipRanges(dar, extend_edges = TRUE)

test_that("assignFeatureDar returns the expected output", {
    gene_dar_origin <- assignFeatureDar(dar, chr1_genes, dar_val = "origin")
    expect_equal(length(gene_dar_origin$group1v2), 1456)
    expect_equal(sum(!is.na(gene_dar_origin$group1v2$dar)), 820)
    gene_dar_region <- assignFeatureDar(
        dar_regions, chr1_genes, dar_val = "region"
    )
    expect_equal(length(gene_dar_region$group1v2), 1456)
    expect_equal(sum(!is.na(gene_dar_region$group1v2$dar)), 1455)
    gene_dar_ex <- assignFeatureDar(
        dar_regions_ex, chr1_genes, dar_val = "region"
    )
    expect_equal(length(gene_dar_ex$group1v2), 1456)
    expect_equal(sum(!is.na(gene_dar_ex$group1v2$dar)), 1456)
})

test_that("assignFeatureDar errors when expected", {
    dar <- endoapply(dar, function(x) x[,c()])
    expect_error(
        assignFeatureDar(dar, chr1_genes, dar_val = "origin"),
        "No dar_origin values detected"
    )
    dar_regions <- endoapply(dar_regions, function(x) x[,c()])
    expect_error(
        assignFeatureDar(dar_regions, chr1_genes, dar_val = "region"),
        "No dar_region values detected"
    )
})

test_that("assignFeatureDar gives warnings when expected", {
    expect_warning(
        assignFeatureDar(dar, chr1_genes, dar_val = "region"),
        "Range\\(s\\) detected with width == 1"
    )
    expect_warning(
        assignFeatureDar(dar_regions, chr1_genes, dar_val = "origin"),
        "Range\\(s\\) detected with width > 1"
    )
})
