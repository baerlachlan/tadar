#' @title Convert ranges from single base positions to elastic windows,
#' or vice versa
#'
#' @description Convert the ranges element associated with single base position
#' raw DAR values to the elastic sliding window ranges associated with the
#' smoothed DAR values, or vice versa
#'
#' @param dar A GRangesList with ranges representing single base positions
#' and metadata columns containing the smoothed DAR values.
#' @param extendEdges A logical specifying if ranges at the edges of each
#' chromosome should be extended to cover the entire chromosome.
#' This argument is only used when converting to ranges associated with the
#' smooth DAR values (i.e. `to = "smooth`).
#' Useful for downstream assignment of DAR values to genomic features that
#' exist at the 5' or 3' edges of the chromosome.
#' @param ... Not used
#'
#' @return A GRangesList with ranges that represent the elastic sliding windows
#' used to smooth the DAR metric.
#' Each element of the list represents a comparison between two sample groups.
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
#' dar <- dar(props, contrasts)
#' convertRanges(dar)
#'
#' ## Convert back to origin ranges
#' darWindows <- convertRanges(dar)
#' convertRanges(darWindows)
#'
#' @import GenomicRanges
#' @importFrom IRanges IRanges ranges 'ranges<-'
#' @importFrom S4Vectors endoapply metadata 'metadata<-' mcols 'mcols<-'
#' @importFrom GenomeInfoDb seqnames seqinfo seqlengths
#' @rdname convertRanges-methods
#' @aliases convertRanges
#' @export
setMethod(
    "convertRanges",
    signature = signature(dar = "GRangesList"),
    function(dar, extendEdges = FALSE) {


        endoapply(dar, function(x){
            if (metadata(x)$rangeType == "origin")
                .switchToWindows(x, extendEdges)
            else if (metadata(x)$rangeType == "window")
                .switchToOrigins(x)
        })

    }
)

#' @keywords internal
.extend <- function(windows, dar) {
    chr <- seqnames(dar)
    chr <- unique(chr)
    chr <- as.character(chr)
    seqlen <- seqlengths(dar)[[chr]]
    if (is.na(seqlen))
        stop(
            "Cannot extend edges. Check seqlength for seqname ",
            chr
        )
    start(windows)[1] <- 1
    end(windows)[length(end(windows))] <- seqlen
    windows
}

#' @keywords internal
.switchToWindows <- function(dar, extendEdges) {

    winSize <- metadata(dar)$winSize
    if (is.null(winSize)) {
        stop(
            "No winSize detected. Use `dar()` before ",
            "`convertRanges()`", call. = FALSE
        )
    }
    removedRanges <- dar[is.na(dar$dar_smooth),]
    grl <- split(dar, f = seqnames(dar))
    grl <- endoapply(grl, function(x){
        isNA <- is.na(x$dar_smooth)
        nNA <- sum(isNA)
        n <- NROW(x)
        ## Construct windows while accounting for NA removal
        windows <- IRanges(
            start = start(x)[seq_len(n - nNA)],
            end = end(x)[winSize:n]
        )
        if (extendEdges) windows <- .extend(windows, x)
        x <- x[!isNA,]
        mcols(x)$origin <- ranges(x)
        ranges(x) <- windows
        x
    })
    gr <- unlist(grl,)
    metadata(gr)$removedRanges <- removedRanges  # Keep for .switchToOrigins()
    metadata(gr)$rangeType <- "window"
    gr

}

#' @keywords internal
.switchToOrigins <- function(dar) {

    ranges(dar) <- mcols(dar)$origin
    mcols(dar) <- mcols(dar)[,c("dar", "dar_smooth")]
    dar <- c(dar, metadata(dar)$removedRanges)
    dar <- sort(dar)
    metadata(dar) <- metadata(dar)[c("rangeType", "winSize")]
    metadata(dar)$rangeType <- "origin"
    dar

}
