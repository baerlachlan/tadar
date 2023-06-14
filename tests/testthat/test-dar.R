fl <- system.file("extdata", "", package="darr")
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
