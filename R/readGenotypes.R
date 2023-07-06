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
#' Alternatively, a `TabixFile` as described in
#' \link[VariantAnnotation]{readVcf}
#' @param unphase A `logical` specifying if phasing information should be
#' removed from genotypes. This is recommended if proceeding with DAR analysis
#' @param ... Passed to \link[VariantAnnotation]{readVcf}
#'
#' @return A `GRanges` object constructed from the CHROM, POS, ID and REF
#' fields of the supplied `VCF` file. Genotype data for each sample present in
#' the `VCF` file is added to the metadata columns
#'
#' @examples
#' fl <- system.file("extdata", "chr1.vcf.bgz", package="darr")
#' readGenotypes(fl)
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
#' @importFrom Rsamtools TabixFile
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
#' @importFrom VariantAnnotation readVcf ScanVcfParam geno vcfGeno
#' @importFrom MatrixGenerics rowRanges
#' @importFrom S4Vectors 'mcols<-'
#' @importFrom GenomeInfoDb 'seqlevels<-' seqlevelsInUse
#' @keywords internal
.readGenotypes <- function(file, unphase, ...) {

    stopifnot(is.logical(unphase))
    dotArgs <- list(...)
    ## Define a minimal ScanVcfParam object for fast loading
    if (!exists("param", where = dotArgs)) dotArgs$param <- ScanVcfParam(
        fixed = NA, info = NA, geno = "GT"
    )
    genoParam <- vcfGeno(dotArgs$param)
    svpChecks <- c(
        ## Check for GT field if user-specified param
        "GT" %in% genoParam,
        ## Or if no geno specified, this will also return GT field
        is.character(genoParam) & length(genoParam) == 0
    )
    stopifnot(any(svpChecks))
    vcf <- readVcf(file, ...)
    gt <- geno(vcf)$GT
    stopifnot(!is.null(gt))
    if (unphase) gt <- unphaseGT(gt)
    gr <- rowRanges(vcf)
    seqlevels(gr) <- seqlevelsInUse(gr)  # Remove seqlevels not present in data
    gr <- unname(gr)  # Remove names to reduce object size
    mcols(gr) <- gt
    gr

}
