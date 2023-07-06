#' @title Plot the Empirical Cumulative Distribution Function of DAR
#'
#' @description Plot the ECDF of DAR for each chromosome
#'
#' @param dar `GRanges` object with metadata columns containing the desired
#' DAR values to plot
#' @param darVal `character(1)` specifying the whether to plot dar_origin or
#' dar_region values.
#' Options are "origin" and "region"
#' @param highlight `character(1)` specifying the chromosome to highlight with
#' a different colour
#'
#' @return A ggplot2 object
#'
#' @examples
#' set.seed(230704)
#' gr <- GRanges(
#'     paste0(rep(seq(1,25), each = 100), ":", seq(1,100)),
#'     dar_origin = runif(2500, 0, 1)
#' )
#' ## No highlighting, all chromosomes will be given individual colours
#' plotDarECDF(gr, darVal = "origin") +
#' theme_bw()
#'
#' ## With highlighting
#' plotDarECDF(gr, darVal = "origin", highlight = "25") +
#' scale_colour_manual(values = c("red", "grey")) +
#' theme_bw()
#'
#' @import GenomicRanges ggplot2
#' @importFrom GenomeInfoDb seqnames
#' @rdname plotDarECDF-methods
#' @aliases plotDarECDF
#' @export
setMethod(
    "plotDarECDF",
    signature = signature(dar = "GRanges"),
    function(dar, darVal, highlight) {

        darVal <- match.arg(darVal)
        if (darVal == "origin") plot_x <- dar$dar_origin
        if (darVal == "region") plot_x <- dar$dar_region
        chr <- seqnames(dar)
        chr <- factor(chr)
        if (!is.null(highlight)) {
            stopifnot(highlight %in% chr)
            lvls <- levels(chr)
            chr <- factor(  # Plot highlighted on top
                chr,
                levels = c(lvls[lvls != highlight], lvls[lvls == highlight])
            )
            line_colour <- factor(chr == highlight, levels = c(TRUE, FALSE))
            colour_lab <- paste("Chromosome", highlight)
        } else {
            line_colour <- chr
            colour_lab <- "Chromosome"
        }
        df <- data.frame(
            Chromosome = chr, dar = plot_x, line_colour = line_colour
        )
        df <- df[!is.na(df$dar),]
        p_aes <- aes(x = dar, group = Chromosome, colour = line_colour)
        p <- ggplot(df, p_aes)
        p <- p + stat_ecdf(geom = "step")
        p <- p + labs(x = "DAR", y = "F(DAR)", colour = colour_lab)
        p

    }
)
