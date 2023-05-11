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
