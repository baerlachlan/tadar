#' @title Unphase genotypes
#'
#' @description Remove phasing information from genotype calls
#'
#' @details
#' Phasing information is not required for a simple DAR analysis.
#' Removing this enables easy counting of alleles from genotype calls.
#'
#' @param gt A matrix or data.frame containing only genotype information
#' @param ... Not used
#'
#' @return
#' A \link[base]{data.frame} containing unphased genotype calls
#'
#' @examples
#' library(VariantAnnotation)
#' fl <- system.file("extdata", "chr22.vcf.gz", package="VariantAnnotation")
#' vcf <- readVcf(fl, "hg19")
#' gt <- geno(vcf)$GT
#' unphaseGT(gt)
#'
#' @rdname unphaseGT-methods
#' @aliases unphaseGT
#' @export
setMethod(
    "unphaseGT",
    signature = signature(gt = "matrix"),
    function(gt, ...) {

        rowIds <- row.names(gt)
        gt <- as.data.frame(gt)
        gt <- lapply(gt, .unphase)
        as.data.frame(gt, row.names = rowIds)

    }
)
#' @rdname unphaseGT-methods
#' @aliases unphaseGT
#' @export
setMethod(
    "unphaseGT",
    signature = signature("data.frame"),
    function(gt, ...) {

        rowIds <- row.names(gt)
        gt <- lapply(gt, .unphase)
        as.data.frame(gt, row.names = rowIds)

    }
)
#' @keywords internal
.unphase <- function(x) {
    gsub("\\|", "\\/", x)
}
