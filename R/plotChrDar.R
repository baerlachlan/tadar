#' @title Plot DAR across a chromosome
#'
#' @description Use `Gviz` to plot the trend in DAR across a chromosomal region
#' with the option to add features of interest as separate tracks.
#'
#' @param dar `GRanges` with DAR values of associated ranges contained
#' in metadata columns.
#' Used to build the DataTrack showing the trend in DAR across the chromosome.
#' If ranges of the input object span regions (i.e. post application of
#' \link{flipRanges}), data points are plotted at the midpoint of the region.
#' @param dar_val `character(1)` specifying the whether to use origin or region
#' DAR values for the chosen ranges.
#' Options are "origin" and "region".
#' The default ("region") represents averaged DAR values across a region and
#' produces a smoother graph.
#' @param chr Optional.
#' `character(1)` specifying the chromosome to subset all `GRanges` objects.
#' Plotting is only possible across a single chromosome and  is therefore
#' required if supplying `GRanges` objects spanning multiple chromosomes.
#' Also controls the track title for the GenomeAxisTrack.
#' @param foi Optional.
#' `GRanges` object of features of interest (foi) to be highlighted and labelled
#' along the GenomeAxisTrack.
#' @param foi_anno Optional.
#' `character(1)` specifying the name of the `mcol` containing labels
#' associated with feature ranges in `foi`.
#' @param foi_highlight `logical(1)` specifying if the positions of `foi`
#' should be overlayed on the DataTrack showing the trend in DAR.
#' Useful for visually inspecting DAR at the location of chosen features.
#' Default is `TRUE`.
#' @param features Optional.
#' `GRanges` object of features to be plotted on a separate AnnotationTrack.
#' @param features_anno Optional.
#' `character(1)` specifying the name of the `mcol` containing labels
#' associated with feature ranges in `features`.
#' @param features_highlight `logical(1)` specifying if the positions of
#' `features` should be overlayed on the DataTrack showing the trend in DAR.
#' Useful for visually inspecting DAR at the location of chosen features.
#' Default is `TRUE`.
#' @param title `character(1)`.
#' The main plot title.
#'
#' @return
#' A Gviz object
#'
#' @examples
#' set.seed(230822)
#' data("chr1_genes")
#' foi <- sample(chr1_genes, 1)
#' features <- sample(chr1_genes, 20)
#' fl <- system.file("extdata", "chr1.vcf.bgz", package="darr")
#' genotypes <- readGenotypes(fl)
#' groups <- list(
#'     group1 = paste0("sample", 1:6),
#'     group2 = paste0("sample", 7:13)
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
#' plotChrDar(
#'     dar = dar$group1v2, dar_val = "region", chr = "1",
#'     foi = foi, foi_anno = "gene_name", foi_highlight = TRUE,
#'     features = features, features_anno = "gene_name",
#'     features_highlight = TRUE,
#'     title = "Example plot of DAR along Chromosome 1"
#' )
#'
#' @import GenomicRanges
#' @importFrom Gviz plotTracks
#' @rdname plotChrDar-methods
#' @aliases plotChrDar
#' @export
setMethod(
    "plotChrDar",
    signature = signature(dar = "GRanges"),
    function(
        dar, dar_val, chr, foi, foi_anno, foi_highlight,
        features, features_anno, features_highlight, title
    ) {

        dar_val <- match.arg(dar_val)
        .plotChrDar_checks(
            dar, chr, foi, foi_anno, foi_highlight,
            features, features_anno, features_highlight
        )
        foi_track <- .foiTrack(chr, foi, foi_anno)
        axis_track <- .axisTrack(chr, foi)
        features_track <- .featuresTrack(chr, features, features_anno)
        dar_track <- .darTrack(
            dar, dar_val, chr, foi, foi_highlight, features, features_highlight
        )
        tracks <- c(foi_track, axis_track, features_track, dar_track)
        plotTracks(
            trackList = tracks, main = title, cex.main = 1, cex.title = 0.6,
            col.title = "black", background.title = "white"
        )

    }
)

#' @keywords internal
#' @importFrom GenomeInfoDb seqnames
#' @importFrom S4Vectors mcols
#' @importFrom Gviz AnnotationTrack
.foiTrack <- function(chr, foi, foi_anno) {

    if (missing(foi)) return(NULL)

    if (!missing(chr))
        foi <- foi[seqnames(foi) == chr]
    ## Assume all ranges are independent if no annotation provided
    group <- seq(1, length(foi))
    anno <- NULL
    if (!missing(foi_anno)){
        foi_anno <- as.character(foi_anno)
        group <- mcols(foi)[[foi_anno]]
        anno <- "group"
    }
    foi_track <- AnnotationTrack(
        range = foi, shape = "box", col= "white",  fill = "white",
        group = group, groupAnnotation = anno, fontcolor.group = "red",
        cex.group = 0.6, size = 0.4, just.group = "below",
        ## Give track name for debugging but don't show on plot
        name = "foi", showTitle = FALSE
    )
    foi_track

}

#' @keywords internal
#' @importFrom Gviz GenomeAxisTrack HighlightTrack
#' @importFrom GenomeInfoDb seqnames
.axisTrack <- function(chr, foi) {

    track_name <- ""  # Defaults to empty string if `chr` not provided
    if (!missing(chr))
        track_name <- paste0("Chr", chr)
    axis_track <- GenomeAxisTrack(
        add53 = TRUE, add35 = TRUE, size = 1.2,
        name = track_name, showTitle = TRUE, rotation.title = 0
    )
    if (missing(foi)) return(axis_track)
    if (!missing(chr))
        foi <- foi[seqnames(foi) == chr]
    axis_track <- HighlightTrack(
        trackList = list(axis_track), range = foi, col = "red", fill = "red",
        name = track_name
    )
    axis_track

}

#' @keywords internal
#' @importFrom GenomeInfoDb seqnames
#' @importFrom S4Vectors mcols
#' @importFrom Gviz AnnotationTrack
.featuresTrack <- function(chr, features, features_anno) {

    if (missing(features)) return(NULL)

    if (!missing(chr))
        features <- features[seqnames(features) == chr]
    ## Assume all ranges are independent if no annotation provided
    group <- seq(1, length(features))
    anno <- NULL
    if (!missing(features_anno)){
        features_anno <- as.character(features_anno)
        group <- mcols(features)[[features_anno]]
        anno <- "group"
    }
    features_track <- AnnotationTrack(
        range = features, name = "Features", shape = "box",  col = "darkgray",
        fill = "darkgray", group = group, groupAnnotation = anno,
        fontcolor.group = "black", cex.group = 0.6, size = 0.4,
        rotation.title = 0
    )
    features_track

}

#' @keywords internal
#' @importFrom GenomeInfoDb seqnames
#' @importFrom Gviz DataTrack HighlightTrack
.darTrack <- function(
        dar, dar_val, chr, foi, foi_highlight, features, features_highlight
) {

    dar <- dar[, paste0("dar_", dar_val)]
    if (!missing(chr))
        dar <- dar[seqnames(dar) == chr]
    dar_track <- DataTrack(
        range = dar, type = "b", name = "DAR", size = 8, window = -1,
        windowSize = 1, cex = 0.4, col = "grey20", col.axis = "black",
        yTicksAt = seq(0, 1, 0.1), ylim = c(0, 1), rotation.title = 90
    )
    ## Assign empty GRanges to allow correct colouring when only one of foi or
    ## features is selected for highlighting
    if (any(missing(foi), !foi_highlight)) foi <- GRanges()
    if (any(missing(features), !features_highlight)) features <- GRanges()
    if (!missing(chr)) {
        foi <- foi[seqnames(foi) == chr]
        features <- features[seqnames(features) == chr]
    }
    ranges <- GRanges()
    if (foi_highlight) ranges <- c(ranges, foi)
    if (features_highlight) ranges <- c(ranges, features)
    if (length(ranges))
        colours <- c(rep("red", length(foi)), rep("#ffe0e0", length(features)))
    dar_track <- HighlightTrack(
        trackList = list(dar_track), range = ranges,
        col = colours, fill = colours
    )
    dar_track

}

#' @keywords internal
.plotChrDar_checks <- function(
        dar, chr, foi, foi_anno, foi_highlight = TRUE,
        features, features_anno, features_highlight = TRUE
) {

    msg <- c()
    msg <- .checkChr(msg, chr)
    msg <- .checkGRanges(msg, dar, chr)
    msg <- .checkGRanges(msg, foi, chr)
    msg <- .checkGRanges(msg, features, chr)
    msg <- .checkConsistentSeqnames(msg, dar, chr, foi, features)
    msg <- .checkAnnotations(msg, foi, foi_anno)
    msg <- .checkAnnotations(msg, features, features_anno)
    msg <- .checkBools(msg, foi_highlight, features_highlight)
    if (!is.null(msg))
        stop(c("\n", msg), call. = FALSE)

}

#' @keywords internal
#' @importFrom methods is
.checkChr <- function(msg, chr) {

    if (missing(chr)) return(msg)

    if (!is(chr, "character"))
        msg <- c(msg, "`chr` must be a character\n")
    msg

}

#' @keywords internal
#' @importFrom GenomeInfoDb seqnames
#' @importFrom methods is
.checkGRanges <- function(msg, gr, chr, arg_name = deparse(substitute(gr))) {

    if (missing(gr)) return(msg)

    if (!is(gr, "GRanges")) {
        msg <- c(msg, paste0("`", arg_name, "` must be a GRanges object\n"))
    } else {
        seqs <- as.character(unique(seqnames(gr)))
        if (missing(chr)) {
            if (length(seqs) != 1)
                msg <- c(
                    msg,
                    paste0(
                        "All ranges of `", arg_name,
                        "` must exist on the same chromosome\n"
                    )
                )
        } else {
            if (!(chr %in% seqs))
                msg <- c(
                    msg,
                    paste0(
                        "Chromosome '", chr, "' not found in ranges of `",
                        arg_name, "`\n"
                    )
                )
        }
    }
    msg

}

#' @keywords internal
#' @importFrom GenomeInfoDb seqnames
.checkConsistentSeqnames <- function(msg, dar, chr, foi, features) {

    ## Only need to check when `chr` missing because .checkGRanges()
    ## already checks that `chr %in% seqs` for each GRanges object
    if (!missing(chr)) return(msg)

    seqs <- c()
    if (length(unique(seqnames(dar))) == 1)  # .checkGRanges will catch if FALSE
        seqs <- c(seqs, as.character(unique(seqnames(dar))))
    if (!missing(foi))
        if (is(foi, "GRanges"))
            if (length(unique(seqnames(foi))) == 1)  # As above
                seqs <- c(seqs, as.character(unique(seqnames(foi))))
    if (!missing(features))
        if (is(features, "GRanges"))
            if (length(unique(seqnames(features))) == 1)  # As above
                seqs <- c(seqs, as.character(unique(seqnames(features))))
    seqs <- unique(seqs)
    if (length(seqs) > 1)
        msg <- c(
            msg,
            paste(
                "All ranges of supplied GRanges objects must exist",
                "on the same chromosome\n"
            )
        )
    msg

}

#' @keywords internal
#' @importFrom S4Vectors mcols
#' @importFrom methods is
.checkAnnotations <- function(
        msg, gr, anno, arg_name = deparse(substitute(gr))
) {

    if (any(missing(gr), missing(anno))) return(msg)

    if (!is(anno, "character")) {
        msg <- c(
            msg,
            paste0(
                "Annotation column for `", arg_name, "` must be a character\n"
            )
        )
    } else {
        if (!(anno %in% names(mcols(gr))))
            msg <- c(
                msg,
                paste0(
                    "Column '", anno, "' not found in mcols of `",
                    arg_name, "`\n"
                )
            )
    }
    msg

}

#' @keywords internal
.checkBools <- function(msg, foi_highlight, features_highlight) {

    if (!is.logical(foi_highlight))
        msg <- c(msg, "`foi_highlight` must be logical\n")
    if (!is.logical(features_highlight))
        msg <- c(msg, "`features_highlight` must be logical\n")
    msg

}
