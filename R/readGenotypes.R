#' @title Read genotypes from a VCF file
#'
#' @description Extract genotypes from a VCF file into a GRanges object
#'
#' @details
#' Extract genotype information from the GT field
#'
#' @param file The file path of a VCF file containing genotype data
#' @param ... Passed to [VariantAnnotation::readVcf()]
#'
#' @return tmp
#'
#' @examples
#' library(VariantAnnotation)
#' fl <- system.file("extdata", "chr22.vcf.gz", package="VariantAnnotation")
#' readGenotypes(fl, genome = "hg19")
#'
#' @import GenomicRanges
#' @importFrom VariantAnnotation readVcf geno
#' @importFrom MatrixGenerics rowRanges
#' @importFrom S4Vectors 'mcols<-'
#'
#' @rdname readGenotypes-methods
#' @aliases readGenotypes
#'
#' @export
setMethod(
    "readGenotypes",
    signature = signature(file = "character"),
    function(file, ...) {

        vcf <- readVcf(file, ...)
        gt <- geno(vcf)$GT
        gt <- unphaseGT(gt)
        gr <- rowRanges(vcf)
        mcols(gr) <- gt
        gr

    }
)
