#' @title DAR-moderated *p*-values
#'
#' @description Moderate *p*-values from DE testing using assigned DAR values.
#'
#' @param pvals `numeric` of feature-level *p*-values from differential
#' expression testing.
#' @param dar `numeric` of DAR values assigned to corresponding features
#' tested for differential expression.
#' @param slope `numeric(1)` specifying the slope of alpha fit.
#' @param min_dar `numeric(1)` p-values where the DAR is below this value will not be moderated
#'
#' @return `numeric` of DAR-moderated *p*-values of same length as
#' input *p*-values.
#'
#' @examples
#' data("chr1_genes")
#' data("chr1_tt")
#' fl <- system.file("extdata", "chr1.vcf.bgz", package="tadar")
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
#' dar_regions <- flipRanges(dar, extend_edges = TRUE)
#' gene_dar <- assignFeatureDar(dar_regions, chr1_genes, dar_val = "region")
#' chr1_tt <- merge(chr1_tt, mcols(gene_dar$group1v2), sort = FALSE)
#' chr1_tt$darP <- modP(chr1_tt$PValue, chr1_tt$dar)
#'
#' @importFrom stats pbeta
#' @rdname modP-methods
#' @aliases modP
#' @export
setMethod(
    "modP",
    signature = signature(pvals = "numeric", dar = "numeric"),
    function(pvals, dar, slope, min_dar = 0.1) {

        if (length(pvals) != length(dar))
            stop("pvals and dar objects must be of same length")
        alpha <- 1 + (slope * dar)
        alpha <- pmax(alpha, min_dar)
        pbeta(q = pvals, shape1 = alpha, shape2 = 1)

    }
)
