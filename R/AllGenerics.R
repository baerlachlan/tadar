#' @name assignFeatureDar
#' @rdname assignFeatureDar-methods
#' @export
setGeneric(
    "assignFeatureDar",
    function(dar, features, dar_val = c("origin", "region"))
        standardGeneric("assignFeatureDar")
)

#' @name flipRanges
#' @rdname flipRanges-methods
#' @export
setGeneric(
    "flipRanges",
    function(dar, extend_edges = FALSE) standardGeneric("flipRanges")
)

#' @name countAlleles
#' @rdname countAlleles-methods
#' @export
setGeneric(
    "countAlleles",
    function(genotypes, groups) standardGeneric("countAlleles")
)

#' @name countsToProps
#' @rdname countsToProps-methods
#' @export
setGeneric(
    "countsToProps",
    function(counts) standardGeneric("countsToProps")
)

#' @name dar
#' @rdname dar-methods
#' @export
setGeneric(
    "dar",
    function(props, contrasts, win_size = 5) standardGeneric("dar")
)

#' @name filterLoci
#' @rdname filterLoci-methods
#' @export
setGeneric(
    "filterLoci",
    function(counts, filter = n_called > n_missing)
        standardGeneric("filterLoci")
)

#' @name plotDarECDF
#' @rdname plotDarECDF-methods
#' @export
setGeneric(
    "plotDarECDF",
    function(dar, dar_val = c("origin", "region"), highlight = NULL)
        standardGeneric("plotDarECDF")
)

#' @name plotChrDar
#' @rdname plotChrDar-methods
#' @export
setGeneric(
    "plotChrDar",
    function(
        dar, dar_val = c("origin", "region"), chr,
        foi, foi_anno, foi_highlight = TRUE,
        features, features_anno, features_highlight = TRUE,
        title = ""
    ) standardGeneric("plotChrDar")
)

#' @name readGenotypes
#' @rdname readGenotypes-methods
#' @export
setGeneric(
    "readGenotypes",
    function(file, unphase = TRUE, ...) standardGeneric("readGenotypes")
)

#' @name unphaseGT
#' @rdname unphaseGT-methods
#' @export
setGeneric(
    "unphaseGT",
    function(gt) standardGeneric("unphaseGT")
)
