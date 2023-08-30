#' @title Convert DAR origin ranges to DAR region ranges, or vice versa
#'
#' @description Convert the ranges element associated with origin DAR values
#' to ranges associated with the region DAR values.
#' This function can also be used to revert back to the original object
#' containing origin ranges if desired.
#'
#' @param dar `GRangesList` with ranges representing single nucleotide (origin)
#' positions.
#' @param extendEdges `logical(1)` specifying if region DAR ranges at the edges
#' of each chromosome should be extended to cover the entire chromosome.
#' This argument is only considered when converting from origin DAR ranges to
#' region DAR ranges.
#' Useful for downstream assignment of DAR values to genomic features that
#' exist at the 5' or 3' edges of the chromosome, which would have otherwise
#' been missed.
#'
#' @return `GRangesList` with ranges that represent either DAR regions or
#' DAR origins, depending on the ranges of the input object.
#'
#' @examples
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
#' dar <- dar(props, contrasts, winSize = 5)
#' ## Convert ranges to regions associated with dar_region values
#' convertRanges(dar)
#' ## Extend the outer regions the cover the entire chromosome
#' convertRanges(dar, extendEdges = TRUE)
#'
#' ## Convert back to origin ranges associated with dar_origin values
#' darRegions <- convertRanges(dar)
#' convertRanges(darRegions)
#'
#' @import GenomicRanges
#' @importFrom S4Vectors endoapply metadata
#' @rdname convertRanges-methods
#' @aliases convertRanges
#' @export
setMethod(
    "convertRanges",
    signature = signature(dar = "GRangesList"),
    function(dar, extendEdges) {


        endoapply(dar, function(x){
            if (metadata(x)$rangeType == "origin")
                .switchToRegion(x, extendEdges)
            else if (metadata(x)$rangeType == "region")
                .switchToOrigin(x)
        })

    }
)

#' @keywords internal
#' @importFrom GenomeInfoDb seqnames seqlengths
.extend <- function(regions, dar) {
    chr <- seqnames(dar)
    chr <- unique(chr)
    chr <- as.character(chr)
    seqlen <- seqlengths(dar)[[chr]]
    if (is.na(seqlen))
        stop(
            "Cannot extend edges. Check seqlength for seqname ",
            chr
        )
    start(regions)[1] <- 1
    end(regions)[length(end(regions))] <- seqlen
    regions
}

#' @keywords internal
#' @importFrom IRanges IRanges ranges 'ranges<-'
#' @importFrom S4Vectors endoapply metadata 'metadata<-' 'mcols<-'
#' @importFrom GenomeInfoDb seqnames
.switchToRegion <- function(dar, extendEdges) {

    winSize <- metadata(dar)$winSize
    if (is.null(winSize)) {
        stop(
            "No winSize detected. Use `dar()` before ",
            "`convertRanges()`", call. = FALSE
        )
    }
    removedRanges <- dar[is.na(dar$dar_region),]
    grl <- split(dar, f = seqnames(dar))
    grl <- endoapply(grl, function(x){
        isNA <- is.na(x$dar_region)
        nNA <- sum(isNA)
        n <- NROW(x)
        ## Construct regions while accounting for NA removal
        regions <- IRanges(
            start = start(x)[seq_len(n - nNA)],
            end = end(x)[winSize:n]
        )
        if (extendEdges) regions <- .extend(regions, x)
        x <- x[!isNA,]
        mcols(x)$origin_ranges <- ranges(x)
        ranges(x) <- regions
        x
    })
    gr <- unlist(grl,)
    metadata(gr)$removedRanges <- removedRanges  # Keep for .switchToOrigin()
    metadata(gr)$rangeType <- "region"
    gr

}

#' @keywords internal
#' @importFrom IRanges 'ranges<-'
#' @importFrom S4Vectors metadata 'metadata<-' mcols 'mcols<-'
.switchToOrigin <- function(dar) {

    ranges(dar) <- mcols(dar)$origin_ranges
    mcols(dar) <- mcols(dar)[,c("dar_origin", "dar_region")]
    dar <- c(dar, metadata(dar)$removedRanges)
    dar <- sort(dar)
    metadata(dar) <- metadata(dar)[c("rangeType", "winSize")]
    metadata(dar)$rangeType <- "origin"
    dar

}
