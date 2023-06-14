#' @title Smooth the Differential Allelic Representation metric
#'
#' @description Add a metadata column alongside raw DAR values containing
#' smoothed DAR values that are calculated using an elastic sliding window
#' approach
#'
#' @param dar A GRangesList containing DAR values at each overlapping range
#' between two sample groups.
#' Each element of the list represents a comparison between two sample groups
#' @param winSize integer specifying the number of ranges to include in the
#' elastic sliding window used for smoothing the DAR metric
#' @param ... Not used
#'
#' @return A GRangesList containing smoothed DAR values calculated as a rolling
#' average over a specified number of ranges.
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
#' smoothDar(dar, winSize = 5)
#'
#' @import GenomicRanges
#' @importFrom S4Vectors endoapply 'metadata<-' 'mcols<-' mcols
#' @importFrom GenomeInfoDb seqnames
#' @importFrom stats filter
#' @rdname smoothDar-methods
#' @aliases smoothDar
#' @export
setMethod(
    "smoothDar",
    signature = signature(dar = "GRangesList"),
    function(dar, winSize = 5) {

        endoapply(dar, function(x){
            if (winSize < 1 || winSize %% 2 != 1)
                stop("`winSize` must be an odd integer greater than 0")
            ## Add winSize to metadata for downstream use (e.g. getWinRanges())
            ## so the user doesn't need to specify winSize in two funcs
            ## Alternatively could use the number of NA's to determine winSize
            ## Not sure about this, remember to ask Stevie
            metadata(x)$winSize <- winSize
            grl <- split(x, f = seqnames(x))
            grl <- endoapply(grl, function(y){
                ## Throw a more informative error than filter() in the
                ## rare scenario that this occurs
                if (winSize > NROW(y))
                    stop(
                        "`winSize` greater than number of ranges for seqname ",
                        unique(seqnames(y)),
                        call. = FALSE
                    )
                mcols(y)$dar_smooth <- filter(
                    y$dar,
                    rep(1 / winSize, winSize),
                    sides = 2
                )
                mcols(y)$dar_smooth <- as.numeric(mcols(y)$dar_smooth)
                y
            })
            unlist(grl)
        })

    }
)
