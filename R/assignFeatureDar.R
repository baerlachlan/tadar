#' @title Assign DAR values to genomic features
#'
#' @description Assign DAR values to genomic features of interest by
#' averaging the DAR values of ranges that overlap the feature range.
#'
#' @param dar `GRangesList` with DAR values of the associated ranges contained
#' in metadata columns.
#' Ranges that represent DAR regions are recommended to assign the greatest
#' number of features with DAR values.
#' This results in an assigned estimate of DAR in the region surrounding the
#' feature.
#' Alternatively, the use of DAR origin ranges results in an assigned average
#' of DAR solely within the feature.
#' Ranges can be converted between origins and regions with
#' \link{flipRanges}.
#' @param features `GRanges` object specifying the features of interest.
#' @param darVal `character(1)` specifying the whether to use origin or region
#' DAR values for the chosen ranges.
#' Options are "origin" and "region".
#' A warning will be produced if the chosen `darVal` does not match the
#' ranges detected in the object provided to the `dar` argument, as this is
#' likely unintended by the user.
#'
#' @return `GRangesList` with ranges representing features of interest that
#' overlap at least one DAR range.
#' Feature metadata columns are retained and an additional column is added
#' for the assigned DAR value.
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
#' assignFeatureDar(chr1_genes, dar, darVal = "origin")
#'
#' darRegions <- flipRanges(dar, extendEdges = TRUE)
#' assignFeatureDar(chr1_genes, darRegions, darVal = "region")
#'
#' @import GenomicRanges
#' @importFrom S4Vectors endoapply from to
#' @rdname assignFeatureDar-methods
#' @aliases assignFeatureDar
#' @export
setMethod(
    "assignFeatureDar",
    signature = signature(dar = "GRangesList", features = "GRanges"),
    function(features, dar, darVal) {

        darVal <- match.arg(darVal)
        endoapply(dar, function(x){
            .assignFeatureDar_checks(x, darVal)
            hits <- findOverlaps(features, x)
            queries <- from(hits)
            queries <- unique(queries)
            darMean <- vapply(queries, function(y){
                subjects <- to(hits)[from(hits) == y]
                if (darVal == "origin") featureDar <- x$dar_origin[subjects]
                if (darVal == "region") featureDar <- x$dar_region[subjects]
                mean(featureDar)
            }, numeric(1))
            features <- features[queries]
            features$dar <- darMean
            features
        })

    }
)

#' @keywords internal
#' @importFrom S4Vectors mcols
#' @importFrom BiocGenerics width
.assignFeatureDar_checks <- function(dar, darVal) {
    widths <- width(dar)
    if (darVal == "region") {
        if (!"dar_region" %in% names(mcols(dar)))
            stop("No dar_region values detected", call. = FALSE)
        if (min(widths) == 1)
            warning(
                "Range(s) detected with width == 1 but darVal = region. ",
                "See ?assignGeneDar", call. = FALSE
            )
    }
    if (darVal == "origin") {
        if (!"dar_origin" %in% names(mcols(dar)))
            stop("No dar_origin values detected", call. = FALSE)
        if (max(widths) > 1)
            warning(
                "Range(s) detected with width > 1 but darVal = origin. ",
                "See ?assignGeneDar", call. = FALSE
            )
    }
}
