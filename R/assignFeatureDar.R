#' @title Assign DAR values to genomic features
#'
#' @description Assign DAR values to genomic features of interest by
#' averaging the DAR values at ranges that overlap the feature range
#'
#' @param dar A GRangesList with ranges representing either single base
#' positions for raw DAR values or windows for smoothed DAR values.
#' DAR values for each range must exist in metadata columns.
#' Windows are recommended to assign the greatest number of features with
#' DAR values.
#' Alternatively, the use of single base positions results in a more precise
#' representation of DAR for overlapping features
#' @param features A GRanges object specifying the features of interest
#' @param darVal A character specifying the whether to use raw or smoothed DAR
#' values for the chosen ranges.
#' Possible options are "raw" and "smooth".
#' Raw DAR values should be chosen for single base positions, while
#' smoothed DAR values should be chosen for windows
#' @param ... Not used
#'
#' @return A GRangesList with ranges representing features of interest that
#' overlap at least one DAR range.
#' Feature metadata columns are retained and an additional column is added
#' for the assigned DAR value.
#' Each element of the list represents a DAR analysis between two sample groups
#'
#' @examples
#' data("chr1_genes")
#' fl <- system.file("extdata", "chr1.vcf.bgz", package="darr")
#' genotypes <- readGenotypes(fl)
#' groups <- list(
#'     group1 = c("S2", "S7", "S9", "S10", "S19", "S20"),
#'     group2 = c("S3", "S6", "S11", "S12", "S15", "S16", "S18")
#' )
#' counts <- countAlleles(genotypes, groups)
#' props <- countsToProps(counts)
#' contrasts <- matrix(
#'     data = c(1, -1),
#'     dimnames = list(
#'         Levels = c("group1", "group2"),
#'         Contrasts = c("group1v2")
#'     )
#' )
#' dar <- dar(props, contrasts)
#' assignFeatureDar(chr1_genes, dar, darVal = "raw")
#'
#' darSmooth <- smoothDar(dar, winSize = 5)
#' darWindows <- getWinRanges(darSmooth, extendEdges = TRUE)
#' assignFeatureDar(chr1_genes, darWindows, darVal = "smooth")
#'
#' @import GenomicRanges
#' @importFrom S4Vectors endoapply from to
#' @rdname assignFeatureDar-methods
#' @aliases assignFeatureDar
#' @export
setMethod(
    "assignFeatureDar",
    signature = signature(dar = "GRangesList", features = "GRanges"),
    function(features, dar, darVal = c("smooth", "raw")) {

        darVal <- match.arg(darVal)
        endoapply(dar, function(x){
            .afdChecks(x, darVal)
            hits <- findOverlaps(features, x)
            queries <- from(hits)
            queries <- unique(queries)
            darMean <- vapply(queries, function(y){
                subjects <- to(hits)[from(hits) == y]
                if (darVal == "smooth") featureDar <- x$dar_smooth[subjects]
                if (darVal == "raw") featureDar <- x$dar[subjects]
                mean(featureDar)
            }, numeric(1))
            features <- features[queries]
            features$dar <- darMean
            features
        })

    }
)

#' @keywords internal#'
#' @importFrom S4Vectors mcols
#' @importFrom BiocGenerics width
.afdChecks <- function(dar, darVal) {
    widths <- width(dar)
    if (darVal == "smooth") {
        if (!"dar_smooth" %in% names(mcols(dar)))
            stop("No smoothed DAR values detected", call. = FALSE)
        if (min(widths) == 1)
            warning(
                "Range(s) detected with width == 1 but darVal = smooth. ",
                "See ?assignGeneDar", call. = FALSE
            )
    }
    if (darVal == "raw") {
        if (!"dar" %in% names(mcols(dar)))
            stop("No DAR values detected", call. = FALSE)
        if (max(widths) > 1)
            warning(
                "Range(s) detected with width > 1 but darVal = raw. ",
                "See ?assignGeneDar", call. = FALSE
            )
    }
}
