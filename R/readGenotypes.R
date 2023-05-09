#' @title Read genotypes from a VCF file
#'
#' @description Extract genotypes from a VCF file into a GRanges object for
#' downstream DAR analysis
#'
#' @details
#' Extract genotypes from a VCF file with the option to remove phasing
#' information for DAR analysis
#'
#' @param file The file path of a VCF file containing genotype data.
#' Alternatively, a \link[Rsamtools]{TabixFile} as described in
#' [VariantAnnotation::readVcf()]
#' @param unphase A `logical` specifying if phasing information should be
#' removed from genotypes. This is recommended if proceeding with DAR analysis
#' @param ... Passed to [VariantAnnotation::readVcf()]
#'
#' @return A `GRanges` object constructed from the CHROM, POS, ID and REF
#' fields of the supplied `VCF` file. Genotype data for each sample present in
#' the `VCF` file is added to the metadata columns
#'
#' @examples
#' library(VariantAnnotation)
#' fl <- system.file("extdata", "chr22.vcf.gz", package="VariantAnnotation")
#' readGenotypes(fl, genome = "hg19")
#'
#' @rdname readGenotypes-methods
#' @aliases readGenotypes
#' @export
setMethod(
    "readGenotypes",
    signature = signature(file = "character"),
    function(file, unphase, ...) {

        .readGenotypes(file, unphase, ...)

    }
)
#' @rdname readGenotypes-methods
#' @aliases readGenotypes
#' @export
setMethod(
    "readGenotypes",
    signature = signature(file = "TabixFile"),
    function(file, unphase, ...) {

        .readGenotypes(file, unphase, ...)

    }
)
#' @import GenomicRanges
#' @importFrom VariantAnnotation readVcf geno
#' @importFrom MatrixGenerics rowRanges
#' @importFrom S4Vectors 'mcols<-'
#' @keywords internal
.readGenotypes <- function(file, unphase, ...) {

    stopifnot(is.logical(unphase))
    vcf <- readVcf(file, ...)
    gt <- geno(vcf)$GT
    if (unphase) gt <- unphaseGT(gt)
    gr <- rowRanges(vcf)
    ## Remove names to reduce object size
    gr <- unname(gr)
    mcols(gr) <- gt
    gr

}
