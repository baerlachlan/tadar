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
    dar <- dar(props, contrasts, region_loci = 5)
    expect_equal(length(mcols(dar[[1]])), 2)
    expect_equal(names(mcols(dar[[1]])), c("dar_origin", "dar_region"))

    dar <- dar(props, contrasts)
    expect_equal(length(mcols(dar[[1]])), 1)
    expect_equal(names(mcols(dar[[1]])), "dar_origin")
})

test_that("dar returns GRangesList with same names set in contrasts", {
    expect_equal(
        names(dar(props, contrasts, region_loci = 5)),
        dimnames(contrasts)[[2]]
    )
})

test_that("dar errors when matrix has incorrect structure", {
    contrasts <- t(contrasts)
    expect_error(
        dar(props, contrasts, region_loci = 5),
        "Levels of `contrasts` must match names of `props`"
    )
})

test_that("dar errors when contrasts are labelled incorrectly", {
    dimnames(contrasts) <- list(Levels = c(), Contrasts = c("group1v2"))
    expect_error(
        dar(props, contrasts, region_loci = 5),
        "Dimnames of `contrasts` must be labelled"
    )
    dimnames(contrasts) <- list(Levels = c("group1", "group2"), Contrasts = c())
    expect_error(
        dar(props, contrasts, region_loci = 5),
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
        dar(props, contrasts, region_loci = 5),
        "`contrasts` defined incorrectly"
    )
})

test_that("dar errors with incorrect region_loci", {
    expect_error(
        dar(props, contrasts, region_loci = 0),
        "`region_loci` must be an odd integer greater than 0"
    )
    expect_error(
        dar(props, contrasts, region_loci = 10),
        "`region_loci` must be an odd integer greater than 0"
    )
    expect_error(
        dar(props, contrasts, region_loci = 10371),
        "`region_loci` greater than number of ranges"
    )
})

test_that("dar errors with incorrect region_loci", {
    expect_error(
        dar(props, contrasts, region_fixed = 0),
        "`region_fixed` must be an integer greater than 0"
    )
    expect_error(
        dar(props, contrasts, region_fixed = -1),
        "`region_fixed` must be an integer greater than 0"
    )
})

test_that("dar adds window info to metadata", {
    region_loci <- 5
    dar <- dar(props, contrasts, region_loci = region_loci)
    region_size <- metadata(dar[[1]])$region_size
    region_type <- metadata(dar[[1]])$region_type
    expect_equal(region_loci, region_size)
    expect_equal(region_type, "elastic")

    region_fixed <- 100
    dar <- dar(props, contrasts, region_fixed = region_fixed)
    region_size <- metadata(dar[[1]])$region_size
    region_type <- metadata(dar[[1]])$region_type
    expect_equal(region_fixed, region_size)
    expect_equal(region_type, "fixed")
})
