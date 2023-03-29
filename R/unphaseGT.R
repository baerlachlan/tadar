#' @title Unphase genotypes
#'
#' @description Remove phasing information from genotype calls
#'
#' @details
#' Phasing information is not required for a simple DAR analysis.
#' Removing this enables easy counting of alleles from genotype calls.
#'
#' @param gt matrix or data.frame containing genotype information
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
#' @importFrom dplyr mutate across
#' @importFrom tidyselect everything
#'
#' @export
unphaseGT <- function(gt) {
    gt <- as.data.frame(gt)
    gt <- dplyr::mutate(gt, dplyr::across(tidyselect::everything(), .unphase))
    gt
}

#' @importFrom stringr str_replace
#' @keywords internal
.unphase <- function(x) {
    str_replace(x, "\\|", "\\/")
}
