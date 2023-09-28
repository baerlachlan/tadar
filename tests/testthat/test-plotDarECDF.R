set.seed(230704)
gr <- GRanges(
    paste0(rep(seq(1,25), each = 100), ":", seq(1,100)),
    dar_origin = runif(2500, 0, 1),
    dar_region = runif(2500, 0, 1)
)
plotDarECDF(gr, dar_val = "origin") %>% str()

test_that("plotDarECDF returns expected output", {
    p <- plotDarECDF(gr, dar_val = "origin")
    expect_true(is(p, "gg"))
    expect_true(is(p$layers[[1]]$geom, "GeomStep"))
    expect_equal(length(levels(p$data$line_colour)), 25)
    p <- plotDarECDF(gr, dar_val = "region")
    expect_true(is(p, "gg"))
    expect_true(is(p$layers[[1]]$geom, "GeomStep"))
    expect_equal(length(levels(p$data$line_colour)), 25)
})

test_that("plotDarECDF errors correctly", {
    expect_error(plotDarECDF(gr, dar_val = "notaval"))
    expect_error(plotDarECDF(gr, dar_val = "origin", highlight = "100"))
})

test_that("plotDarECDF highlights correctly", {
    p <- plotDarECDF(gr, dar_val = "origin", highlight = "22")
    expect_equal(length(levels(p$data$line_colour)), 2)
})
