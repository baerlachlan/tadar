#' @name assignFeatureDar
#' @rdname assignFeatureDar-methods
#' @export
setGeneric(
    "assignFeatureDar",
    function(features, dar, ...) standardGeneric("assignFeatureDar")
)

#' @name convertRanges
#' @rdname convertRanges-methods
#' @export
setGeneric(
    "convertRanges",
    function(dar, ...) standardGeneric("convertRanges")
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
    function(counts, ...) standardGeneric("countsToProps")
)

#' @name dar
#' @rdname dar-methods
#' @export
setGeneric(
    "dar",
    function(props, contrasts, ...) standardGeneric("dar")
)

#' @name filterLoci
#' @rdname filterLoci-methods
#' @export
setGeneric(
    "filterLoci",
    function(counts, ...) standardGeneric("filterLoci")
)

#' @name readGenotypes
#' @rdname readGenotypes-methods
#' @export
setGeneric(
    "readGenotypes",
    function(file, ...) standardGeneric("readGenotypes")
)

#' @name unphaseGT
#' @rdname unphaseGT-methods
#' @export
setGeneric(
    "unphaseGT",
    function(gt) standardGeneric("unphaseGT")
)
