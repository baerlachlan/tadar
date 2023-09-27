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
dar <- dar(props, contrasts, win_loci = 5)
dar_fixed <- dar(props, contrasts, win_fixed = 101)

test_that("flipRanges returns regions", {
    dar_regions <- flipRanges(dar)
    widths <- width(dar_regions$group1v2)
    expect_true(max(widths) > 1)

    dar_regions <- flipRanges(dar_fixed)
    widths <- width(dar_regions$group1v2)
    expect_true(max(widths) == 101)
})

test_that("flipRanges extends ranges", {
    dar_regions <- flipRanges(dar, extend_edges = TRUE)
    expect_equal(min(start(dar_regions$group1v2)), 1)
    expect_equal(max(end(dar_regions$group1v2)), 59578282)
})

test_that("flipRanges errors when missing metadata", {
    dar <- endoapply(dar, function(x){
        metadata(x) <- list()
        x
    })
    expect_error(
        flipRanges(dar),
        paste0(
            "Required metadata not detected\\. Use `dar\\(\\)` with ",
            "either `win_fixed` or `win_loci` arguments specified ",
            "before `flipRanges\\(\\)`"
        )
    )
})

test_that("flipRanges errors when no seqlengths and extend_edges = TRUE", {
    dar <- endoapply(dar, function(x){
        seqinfo(x) <- Seqinfo(seqnames = "1")
        x
    })
    expect_error(
        flipRanges(dar, extend_edges = TRUE),
        "Cannot extend edges\\."
    )
})

test_that("flipRanges can revert back to output of dar()", {
    dar_regions <- flipRanges(dar)
    expect_identical(flipRanges(dar_regions), dar)
    dar_regions <- flipRanges(dar, extend_edges = TRUE)
    expect_identical(flipRanges(dar_regions), dar)
})
