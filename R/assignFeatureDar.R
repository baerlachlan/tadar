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
#' @param dar_val `character(1)` specifying the whether to use origin or region
#' DAR values for the chosen ranges.
#' Options are "origin" and "region".
#' A warning will be produced if the chosen `dar_val` does not match the
#' ranges detected in the object provided to the `dar` argument, as this is
#' likely unintended by the user.
#' @param fill_missing The DAR value to assign features with no overlaps.
#' Defaults to `NA`.
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
#'     group1 = paste0("sample", 1:6),
#'     group2 = paste0("sample", 7:13)
#' )
#' counts <- countAlleles(genotypes, groups)
#' counts_filt <- filterLoci(counts)
#' props <- countsToProps(counts_filt)
#' contrasts <- matrix(
#'     data = c(1, -1),
#'     dimnames = list(
#'         Levels = c("group1", "group2"),
#'         Contrasts = c("group1v2")
#'     )
#' )
#' dar <- dar(props, contrasts, region_loci = 5)
#' assignFeatureDar(dar, chr1_genes, dar_val = "origin")
#'
#' dar_regions <- flipRanges(dar, extend_edges = TRUE)
#' assignFeatureDar(dar_regions, chr1_genes, dar_val = "region")
#'
#' @import GenomicRanges
#' @importFrom S4Vectors endoapply from to
#' @rdname assignFeatureDar-methods
#' @aliases assignFeatureDar
#' @export
setMethod(
    "assignFeatureDar",
    signature = signature(dar = "GRangesList", features = "GRanges"),
    function(dar, features, dar_val, fill_missing) {

        dar_val <- match.arg(dar_val)
        endoapply(dar, function(x){
            .assignFeatureDar_checks(x, dar_val)
            hits <- findOverlaps(features, x)
            queries <- from(hits)
            queries <- unique(queries)
            dar_mean <- vapply(queries, function(y){
                subjects <- to(hits)[from(hits) == y]
                if (dar_val == "origin") featureDar <- x$dar_origin[subjects]
                if (dar_val == "region") featureDar <- x$dar_region[subjects]
                mean(featureDar)
            }, numeric(1))
            features$dar <- fill_missing
            features[queries]$dar <- dar_mean
            features
        })

    }
)

#' @keywords internal
#' @importFrom S4Vectors mcols
#' @importFrom BiocGenerics width
.assignFeatureDar_checks <- function(dar, dar_val) {
    widths <- width(dar)
    if (dar_val == "region") {
        if (!"dar_region" %in% names(mcols(dar)))
            stop("No dar_region values detected", call. = FALSE)
        if (min(widths) == 1)
            warning(
                "Range(s) detected with width == 1 but dar_val = region. ",
                "See ?assignGeneDar", call. = FALSE
            )
    }
    if (dar_val == "origin") {
        if (!"dar_origin" %in% names(mcols(dar)))
            stop("No dar_origin values detected", call. = FALSE)
        if (max(widths) > 1)
            warning(
                "Range(s) detected with width > 1 but dar_val = origin. ",
                "See ?assignGeneDar", call. = FALSE
            )
    }
}
