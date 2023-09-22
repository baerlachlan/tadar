#' @title Convert allele counts to proportions
#'
#' @description Normalise allele-level counts across samples by converting to a
#' proportion of total alleles in all samples.
#' Optionally remove ranges using the default filtering criteria described
#' in \link{filterLoci}.
#'
#' @param counts `GRangesList` containing a summary of allele counts at each
#' range.
#' Each element of the list represents a distinct sample group.
#' @param filter `logical(1)` specifying if ranges should be filtered using
#' default criteria.
#' If filtering is not required, or if performed manually with
#' \link{filterLoci}, set this to `FALSE`.
#'
#' @return `GRangesList` containing a summary of normalised allele counts
#' (i.e. as proportions) at each range.
#' Each element of the list represents a distinct sample group.
#'
#' @examples
#' fl <- system.file("extdata", "chr1.vcf.bgz", package="darr")
#' genotypes <- readGenotypes(fl)
#' groups <- list(
#'     group1 = paste0("sample", 1:6),
#'     group2 = paste0("sample", 7:13)
#' )
#' counts <- countAlleles(genotypes, groups)
#' countsToProps(counts)
#'
#' @import GenomicRanges
#' @importFrom S4Vectors mcols 'mcols<-' endoapply
#' @rdname countsToProps-methods
#' @aliases countsToProps
#' @export
setMethod(
    "countsToProps",
    signature = signature(counts = "GRangesList"),
    function(counts, filter) {

        if (filter) counts <- filterLoci(counts)
        endoapply(counts, function(x){
            gr <- granges(x)
            x <- mcols(x)
            if (!filter & min(x$n_called) == 0)
                stop(
                    "Detected range(s) with no counts. ",
                    "Set `filter = TRUE` or filter manually (see ?filterLoci)"
                )
            check_names <- c(
                "n_called", "n_missing", "n_0", "n_1", "n_2", "n_3"
            )
            if (!all(names(x) == check_names))
                stop(
                    'Names of metadata columns must equal c("',
                    paste(check_names, collapse = '", "'),
                    '")'
                )
            x <- as.matrix(x)
            prop_0 <- x[,"n_0"] / x[,"n_called"]
            prop_1 <- x[,"n_1"] / x[,"n_called"]
            prop_2 <- x[,"n_2"] / x[,"n_called"]
            prop_3 <- x[,"n_3"] / x[,"n_called"]
            props <- matrix(
                data = c(prop_0, prop_1, prop_2, prop_3),
                ncol = 4,
                byrow = FALSE,
                dimnames = list(NULL, c("prop_0", "prop_1", "prop_2", "prop_3"))
            )
            mcols(gr) <- props
            gr
        })

    }
)
