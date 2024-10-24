data("chr1_genes")
data("chr1_tt")
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
dar_regions <- flipRanges(dar, extend_edges = TRUE)
gene_dar <- assignFeatureDar(dar_regions, chr1_genes, dar_val = "region")
chr1_tt <- merge(chr1_tt, mcols(gene_dar$group1v2), sort = FALSE)
modP(chr1_tt$PValue, chr1_tt$dar)

test_that("modP errors when lengths don't match", {
    expect_error(
        modP(chr1_tt$PValue[-1], chr1_tt$dar),
        "pvals and dar objects must be of same length"
    )
})

test_that("modP returns numeric of same length", {
    expect_equal(
        length(modP(chr1_tt$PValue, chr1_tt$dar)),
        length(chr1_tt$PValue)
    )
})
