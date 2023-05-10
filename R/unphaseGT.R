#' @title Unphase genotypes
#'
#' @description Remove phasing information from genotype calls
#'
#' @details
#' Phasing information is not required for a simple DAR analysis.
#' Removing this enables easy counting of alleles from genotype calls.
#'
#' @param gt A matrix or data.frame containing only genotype information
#'
#' @return
#' A \link[base]{matrix} containing unphased genotype calls
#'
#' @examples
#' library(VariantAnnotation)
#' fl <- system.file("extdata", "chr1.vcf.gz", package="darr")
#' vcf <- readVcf(fl)
#' gt <- geno(vcf)$GT
#' unphaseGT(gt)
#'
#' @rdname unphaseGT-methods
#' @aliases unphaseGT
#' @export
setMethod(
    "unphaseGT",
    signature = signature(gt = "matrix"),
    function(gt) {

        gt <- apply(gt, c(1,2), .unphase)
        gt

    }
)
#' @rdname unphaseGT-methods
#' @aliases unphaseGT
#' @export
setMethod(
    "unphaseGT",
    signature = signature("data.frame"),
    function(gt) {

        gt <- as.matrix(gt)
        gt <- apply(gt, c(1,2), .unphase)
        gt

    }
)
#' @keywords internal
.unphase <- function(x) {
    gsub("\\|", "\\/", x)
}
