fl <- system.file("extdata", "chr1.vcf.bgz", package="darr")
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

test_that("dar returns expected mcols", {
    dar <- dar(props, contrasts)
    expect_equal(length(mcols(dar[[1]])), 2)
    expect_equal(names(mcols(dar[[1]])), c("dar_origin", "dar_region"))
})

test_that("dar returns GRangesList with same names set in contrasts", {
    expect_equal(names(dar(props, contrasts)), dimnames(contrasts)[[2]])
})

test_that("dar errors when matrix has incorrect structure", {
    contrasts <- t(contrasts)
    expect_error(
        dar(props, contrasts),
        "Levels of `contrasts` must match names of `props`"
    )
})

test_that("dar errors when contrasts are labelled incorrectly", {
    dimnames(contrasts) <- list(Levels = c(), Contrasts = c("group1v2"))
    expect_error(
        dar(props, contrasts),
        "Dimnames of `contrasts` must be labelled"
    )
    dimnames(contrasts) <- list(Levels = c("group1", "group2"), Contrasts = c())
    expect_error(
        dar(props, contrasts),
        "Dimnames of `contrasts` must be labelled"
    )
})

test_that("dar errors when contrasts are defined incorrectly", {
    contrasts <- matrix(
        data = c(1, 0),
        dimnames = list(
            Levels = c("group1", "group2"),
            Contrasts = c("group1v2")
        )
    )
    expect_error(
        dar(props, contrasts),
        "`contrasts` defined incorrectly"
    )
})

test_that("dar errors with incorrect win_loci", {
    expect_error(
        dar(props, contrasts, win_loci = 0),
        "`win_loci` must be an odd integer greater than 0"
    )
    expect_error(
        dar(props, contrasts, win_loci = 10),
        "`win_loci` must be an odd integer greater than 0"
    )
    expect_error(
        dar(props, contrasts, win_loci = 10371),
        "`win_loci` greater than number of ranges"
    )
})

test_that("dar adds win_loci to metadata", {
    win_loci <- 5
    dar <- dar(props, contrasts, win_loci = win_loci)
    winMeta <- metadata(dar[[1]])$win_loci
    expect_equal(win_loci, winMeta)
})
