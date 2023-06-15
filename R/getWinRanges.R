#' @title Get the elastic sliding window ranges of the smoothed DAR metric
#'
#' @description Replace the ranges associated with single base position raw
#' DAR values to the elastic sliding window ranges associated with the smoothed
#' DAR values
#'
#' @param dar A GRangesList with ranges representing single base positions
#' and metadata columns containing the smoothed DAR values.
#' @param extendEdges A logical specifying if ranges at the edges of each
#' chromosome should be extended to cover the entire chromosome.
#' Useful for assigning DAR values to genomic features that exist
#' at the 5' or 3' edges of the chromosome
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
#' darSmooth <- smoothDar(dar, winSize = 5)
#' getWinRanges(darSmooth)
#'
#' @import GenomicRanges
#' @importFrom IRanges IRanges
#' @importFrom S4Vectors endoapply metadata 'mcols<-' mcols
#' @importFrom GenomeInfoDb seqnames seqinfo seqlengths
#' @rdname getWinRanges-methods
#' @aliases getWinRanges
#' @export
setMethod(
    "getWinRanges",
    signature = signature(dar = "GRangesList"),
    function(dar, extendEdges = FALSE) {

        endoapply(dar, function(x){
            winSize <- metadata(x)$winSize
            if (is.null(winSize)) {
                stop(
                    "No winSize detected. Use `smoothDar()` before ",
                    "`getWinRanges()`", call. = FALSE
                )
            }
            grl <- split(x, f = seqnames(x))
            grl <- endoapply(grl, function(y){
                isNA <- is.na(y$dar_smooth)
                nNA <- sum(isNA)
                # ## Commenting this out while deciding if it's needed
                # ## since using metadata() now as opposed to a winSize arg
                # if (nNA != winSize - 1)
                #     stop(
                #         "`winSize` does not match the window size used to ",
                #         "smooth the DAR metric"
                #     )
                ## Grab edges of windows while accounting for NA removal
                n <- NROW(y)
                start <- start(y)[seq_len(n - nNA)]
                end <- end(y)[winSize:n]
                if (extendEdges) {
                    chr <- seqnames(y)
                    chr <- unique(chr)
                    chr <- as.character(chr)
                    seqlen <- seqlengths(y)[[chr]]
                    if (is.na(seqlen))
                        stop(
                            "Cannot extend edges. Check seqlength for seqname ",
                            chr
                        )
                    start[1] <- 1
                    end[length(end)] <- seqlen
                }
                ## Remove edge ranges with missing smoothed DAR values
                y <- y[!isNA,]
                ## Replace SNP ranges with window ranges
                gr <- GRanges(
                    seqnames = seqnames(y),
                    ranges = IRanges(start = start, end = end),
                    seqinfo = seqinfo(y)
                )
                mcols(gr) <- mcols(y)
                gr
            })
            unlist(grl)
        })

    }
)
