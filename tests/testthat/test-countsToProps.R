fl <- system.file("extdata", "", package="darr")
genotypes <- readGenotypes(fl)
groups <- list(
    group1 = c("S2", "S7", "S9", "S10", "S19", "S20"),
    group2 = c("S3", "S6", "S11", "S12", "S15", "S16", "S18")
)
counts <- countAlleles(genotypes, groups)

test_that("countsToProps returns proportions", {
    props <- countsToProps(counts)
    rowsums <- lapply(props, function(x){
        x <- mcols(x)
        x <- as.matrix(x)
        sums <- rowSums(x)
        unique(sums)
    })
    expect_equal(unique(unlist(rowsums)), 1)
})

test_that("countsToProps retains all ranges when `filter = FALSE`", {
    counts <- filterLoci(counts)
    props <- countsToProps(counts, filter = FALSE)
    expect_equal(granges(counts$group1), granges(props$group1))
})

test_that("countsToProps errors when counts are incorrectly filtered", {
    expect_error(
        countsToProps(counts, filter = FALSE),
        "Detected range\\(s\\) with no counts.+"
    )
})

test_that("countsToProps errors when metadata columns are named incorrectly", {
    ## Need to filter before changing names so we can set filter = FALSE
    ##  for countsToProps(). Otherwise the equivalent name check in filterLoci()
    ##  throws the error instead of countsToProps()
    counts <- filterLoci(counts)
    counts <- lapply(counts, function(x) {
        mc <- mcols(x)
        names(mc) <- c("n_called", "n_missing", "n_1", "n_2", "n_3", "n_0")
        mcols(x) <- mc
        x
    })
    counts <- GRangesList(counts)
    expect_error(
        countsToProps(counts, filter = FALSE),
        'Names of metadata columns must equal c(.+)'
    )
})
