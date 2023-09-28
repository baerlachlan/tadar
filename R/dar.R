#' @title Calculate Differential Allelic Representation (DAR)
#'
#' @description Calculate DAR between two sample groups.
#'
#' @details DAR is calculated as the Euclidean distance between the allelic
#' proportions (i.e. proportion of As, Cs, Gs and Ts) of two sample groups at
#' a single nucleotide locus, scaled such that all values range inclusively
#' between 0 and 1.
#' A DAR value of 0 represents identical allelic representation between the
#' two sample groups, while a DAR value of 1 represents complete diversity.
#'
#' @param props `GRangesList` containing a summary of normalised allele counts
#' (i.e. as proportions) at each range.
#' Each element of the list represents a distinct sample group.
#' @param contrasts Contrast `matrix` specifying which sample groups to
#' to calculate DAR between.
#' Each column must represent a single contrast, and rows represent the levels
#' (i.e. sample groups) to be contrasted.
#' The two levels involved with each contrast should be specified with
#' `1` and `-1`.
#' @param region_fixed `integer(1)` specifying the width (in base pairs) of a
#' fixed sliding window used for averaging DAR values within a region, which is
#' centralised around the origin.
#' Must be an integer greater than 1.
#' This argument takes precedence over `region_loci`.
#' @param region_loci `integer(1)` specifying the number of loci to include
#' in an elastic sliding window used for averaging DAR values within a region.
#' Must be an odd integer in order to incorporate the origin locus and an
#' equal number of loci either side.
#' Only used when `region_fixed` is NULL.
#'
#' @return `GRangesList` containing DAR values at each overlapping range
#' between the contrasted sample groups.
#' Two types of DAR values are reported in the metadata columns of each GRanges
#' object:
#'
#' - `dar_origin`: The raw DAR values calculated at single nucleotide positions
#' (the origin) between sample groups.
#' These values represent DAR estimates at a precise locus.
#' - `dar_region`: The mean of raw DAR values in a specified region surrounding
#' the origin.
#' This is optionally returned using either of the `region_fixed` or
#' `region_loci` arguments, which control the strategy and size for
#' establishing regions (more information below).
#' This option exists because eQTLs don't necessarily confer their effects on
#' genes in close proximity.
#' Therefore, DAR estimates that are representative of regions may be more
#' suitable for downstream assignment DAR values to genomic features.
#'
#' Each element of the list represents a single contrast defined in the
#' input contrast matrix.
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
#' dar(props, contrasts, region_loci = 5)
#'
#' @import GenomicRanges
#' @rdname dar-methods
#' @aliases dar
#' @export
setMethod(
    "dar",
    signature = signature(props = "GRangesList", contrasts = "matrix"),
    function(props, contrasts, region_fixed, region_loci) {

        lvls <- dimnames(contrasts)[[1]]
        conts <- dimnames(contrasts)[[2]]
        if (any(is.null(lvls), is.null(conts)))
            stop("Dimnames of `contrasts` must be labelled")
        if (!all(lvls %in% names(props)))
            stop("Levels of `contrasts` must match names of `props`")
        contrasts <- .contrastsAsList(contrasts)
        grl <- .calcDar(props = props, contrasts = contrasts)
        if (!is.null(region_fixed)) {
            grl <- .smoothAcrossFixed(dar = grl, region_fixed = region_fixed)
        } else if (!is.null(region_loci)) {
            grl <- .smoothAcrossLoci(dar = grl, region_loci = region_loci)
        }
        grl

    }
)

#' @keywords internal
.contrastsAsList <- function(contrasts) {

    cols <- seq(ncol(contrasts))
    contrasts_list <- lapply(cols, function(i){
        grp1 <- names(which(contrasts[,i] == 1))
        grp2 <- names(which(contrasts[,i] == -1))
        if (!all(length(grp1) == 1, length(grp2) == 1))
            stop("`contrasts` defined incorrectly")
        c(grp1, grp2)
    })
    names(contrasts_list) <- colnames(contrasts)
    contrasts_list

}

#' @keywords internal
#' @importFrom S4Vectors mcols 'mcols<-' from to 'metadata<-'
#' @importFrom stats dist
#' @importFrom GenomeInfoDb 'seqlevels<-' seqlevelsInUse
.calcDar <- function(props, contrasts) {

    grl <- lapply(contrasts, function(x){
        ## Subset for groups to be contrasted
        props <- props[x]
        ## Subset for ranges present in both groups
        overlaps <- findOverlaps(props[[1]], props[[2]])
        props[[1]] <- props[[1]][from(overlaps),]
        props[[2]] <- props[[2]][to(overlaps),]
        stopifnot(identical(granges(props[[1]]), granges(props[[2]])))
        gr <- granges(props[[1]])
        props <- lapply(props, mcols)
        props <- lapply(props, as.matrix)
        ## Calculate Euclidean distance
        dist <- vapply(seq_along(gr), function(i){
            mat <- rbind(props[[1]][i,], props[[2]][i,])
            dist <- dist(mat, method = "euclidean")
            as.numeric(dist)
        }, numeric(1))
        ## Convert to DAR
        mcols(gr)$dar_origin <- dist / sqrt(2)
        ## Remove seqlevels that may be lost due to no overlap
        seqlevels(gr) <- seqlevelsInUse(gr)
        ## Add range_type to metadata for downstream use
        metadata(gr)$range_type <- "origin"
        gr
    })
    GRangesList(grl)

}

#' @keywords internal
#' @importFrom S4Vectors endoapply 'mcols<-' 'metadata<-' from to
.smoothAcrossFixed <- function(dar, region_fixed) {

    if (region_fixed < 1)
        stop("`region_fixed` must be an integer greater than 0")
    endoapply(dar, function(x){
        ## Add window info to metadata for downstream use
        metadata(x)$region_type <- "fixed"
        metadata(x)$region_size <- region_fixed
        ## Resize ranges to have width = region_fixed
        ## Suppress warnings for out-of-bound ranges because we trim these
        regions <- suppressWarnings(resize(x, region_fixed, fix = "center"))
        regions <- trim(regions)
        ## Find origin ranges that overlap our regions/windows
        hits <- findOverlaps(regions, x)
        queries <- from(hits)
        queries <- unique(queries)
        ## Now average any origin dar values contained within the regions
        dar_region <- vapply(queries, function(y){
            subjects <- to(hits)[from(hits) == y]
            mean(x$dar_origin[subjects])
        }, numeric(1))
        mcols(x)$dar_region <- dar_region
        x
    })

}

#' @keywords internal
#' @importFrom S4Vectors endoapply mcols 'mcols<-' 'metadata<-'
#' @importFrom GenomeInfoDb seqnames
#' @importFrom stats filter
.smoothAcrossLoci <- function(dar, region_loci) {

    if (region_loci < 1 || region_loci %% 2 != 1)
        stop("`region_loci` must be an odd integer greater than 0")
    endoapply(dar, function(x){
        ## Add window info to metadata for downstream use
        metadata(x)$region_type <- "elastic"
        metadata(x)$region_size <- region_loci
        grl <- split(x, f = seqnames(x))
        grl <- endoapply(grl, function(y){
            ## Throw a more informative error than filter() would
            if (region_loci > NROW(y))
                stop(
                    "`region_loci` greater than number of ranges for seqname ",
                    unique(seqnames(y)), call. = FALSE
                )
            mcols(y)$dar_region <- filter(
                y$dar_origin, rep(1 / region_loci, region_loci), sides = 2
            )
            mcols(y)$dar_region <- as.numeric(mcols(y)$dar_region)
            y
        })
        unlist(grl, use.names = FALSE)
    })

}
