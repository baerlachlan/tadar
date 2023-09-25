#' @title Convert DAR origin ranges to DAR region ranges, or vice versa
#'
#' @description Convert the ranges element associated with origin DAR values
#' to ranges associated with the region DAR values.
#' This function can also be used to revert back to the original object
#' containing origin ranges if desired.
#'
#' @param dar `GRangesList` with ranges representing single nucleotide (origin)
#' positions.
#' @param extend_edges `logical(1)` specifying if region DAR ranges at the edges
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
#' dar <- dar(props, contrasts, win_size = 5)
#' ## Convert ranges to regions associated with dar_region values
#' flipRanges(dar)
#' ## Extend the outer regions the cover the entire chromosome
#' flipRanges(dar, extend_edges = TRUE)
#'
#' ## Convert back to origin ranges associated with dar_origin values
#' darRegions <- flipRanges(dar)
#' flipRanges(darRegions)
#'
#' @import GenomicRanges
#' @importFrom S4Vectors endoapply metadata
#' @rdname flipRanges-methods
#' @aliases flipRanges
#' @export
setMethod(
    "flipRanges",
    signature = signature(dar = "GRangesList"),
    function(dar, extend_edges) {


        endoapply(dar, function(x){
            if (metadata(x)$range_type == "origin")
                .switchToRegion(x, extend_edges)
            else if (metadata(x)$range_type == "region")
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
.switchToRegion <- function(dar, extend_edges) {

    win_size <- metadata(dar)$win_size
    if (is.null(win_size)) {
        stop(
            "No win_size detected. Use `dar()` before ",
            "`flipRanges()`", call. = FALSE
        )
    }
    removed_ranges <- dar[is.na(dar$dar_region),]
    grl <- split(dar, f = seqnames(dar))
    grl <- endoapply(grl, function(x){
        is_NA <- is.na(x$dar_region)
        n_NA <- sum(is_NA)
        n <- NROW(x)
        ## Construct regions while accounting for NA removal
        regions <- IRanges(
            start = start(x)[seq_len(n - n_NA)],
            end = end(x)[win_size:n]
        )
        if (extend_edges) regions <- .extend(regions, x)
        x <- x[!is_NA,]
        mcols(x)$origin_ranges <- ranges(x)
        ranges(x) <- regions
        x
    })
    gr <- unlist(grl,)
    metadata(gr)$removed_ranges <- removed_ranges  # Keep for .switchToOrigin()
    metadata(gr)$range_type <- "region"
    gr

}

#' @keywords internal
#' @importFrom IRanges 'ranges<-'
#' @importFrom S4Vectors metadata 'metadata<-' mcols 'mcols<-'
.switchToOrigin <- function(dar) {

    ranges(dar) <- mcols(dar)$origin_ranges
    mcols(dar) <- mcols(dar)[,c("dar_origin", "dar_region")]
    dar <- c(dar, metadata(dar)$removed_ranges)
    dar <- sort(dar)
    metadata(dar) <- metadata(dar)[c("range_type", "win_size")]
    metadata(dar)$range_type <- "origin"
    dar

}
