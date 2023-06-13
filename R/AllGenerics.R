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

#' @name getWinRanges
#' @rdname getWinRanges-methods
#' @export
setGeneric(
    "getWinRanges",
    function(dar, ...) standardGeneric("getWinRanges")
)

#' @name readGenotypes
#' @rdname readGenotypes-methods
#' @export
setGeneric(
    "readGenotypes",
    function(file, ...) standardGeneric("readGenotypes")
)

#' @name smoothDar
#' @rdname smoothDar-methods
#' @export
setGeneric(
    "smoothDar",
    function(dar, ...) standardGeneric("smoothDar")
)

#' @name unphaseGT
#' @rdname unphaseGT-methods
#' @export
setGeneric(
    "unphaseGT",
    function(gt) standardGeneric("unphaseGT")
)
