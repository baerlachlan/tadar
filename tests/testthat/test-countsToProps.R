fl <- system.file("extdata", "chr1.vcf.bgz", package="darr")
genotypes <- readGenotypes(fl)
groups <- list(
    group1 = paste0("sample", 1:6),
    group2 = paste0("sample", 7:13)
)
counts <- countAlleles(genotypes, groups)
counts_filt <- filterLoci(counts)

test_that("countsToProps returns proportions", {
    props <- countsToProps(counts_filt)
    rowsums <- lapply(props, function(x){
        x <- mcols(x)
        x <- as.matrix(x)
        sums <- rowSums(x)
        unique(sums)
    })
    expect_equal(unique(unlist(rowsums)), 1)
})

test_that("countsToProps errors when counts are not filtered", {
    expect_error(
        countsToProps(counts),
        "Detected range\\(s\\) with no counts.+"
    )
})

test_that("countsToProps errors when metadata columns are named incorrectly", {
    counts_filt <- lapply(counts_filt, function(x) {
        mc <- mcols(x)
        names(mc) <- c("n_called", "n_missing", "n_1", "n_2", "n_3", "n_0")
        mcols(x) <- mc
        x
    })
    counts_filt <- GRangesList(counts_filt)
    expect_error(
        countsToProps(counts_filt),
        'Names of metadata columns must equal c(.+)'
    )
})
