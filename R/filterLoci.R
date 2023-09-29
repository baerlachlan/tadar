#' @title Filter loci
#'
#' @description Filter loci based on allele count criteria.
#'
#' @param counts `GRangesList` containing a summary of allele counts at each
#' range.
#' Each element of the list represents a distinct sample group.
#' @param filter A logical expression indicating which rows to keep.
#' Possible values include:
#'
#' - `n_called`
#' The number of total alleles called.
#' - `n_missing`
#' The number of total alleles not reported.
#' - `n_0`, `n_1`, `n_2`, `n_3`
#' The number of ref, alt1, alt2 and alt3 alleles respectively.
#'
#' All values represent the sum of counts across all samples within the group.
#' Defaults to return loci where the number of samples containing allele
#' information is greater than number samples with missing information.
#'
#' @return `GRangesList` containing a summary of allele counts at each range
#' passing the filter criteria.
#' Each element of the list represents a distinct sample group.
#'
#' @examples
#' fl <- system.file("extdata", "chr1.vcf.bgz", package="tadar")
#' genotypes <- readGenotypes(fl)
#' groups <- list(
#'     group1 = paste0("sample", 1:6),
#'     group2 = paste0("sample", 7:13)
#' )
#' counts <- countAlleles(genotypes, groups)
#' filterLoci(counts)
#'
#' @import GenomicRanges
#' @importFrom S4Vectors mcols endoapply
#' @importFrom rlang enquo eval_tidy
#' @importFrom GenomeInfoDb 'seqlevels<-' seqlevelsInUse
#' @rdname filterLoci-methods
#' @aliases filterLoci
#' @export
setMethod(
    "filterLoci",
    signature = signature(counts = "GRangesList"),
    function(counts, filter) {

        ## Defuse filter expression
        filter <- enquo(filter)
        endoapply(counts, function(x){
            check_names <- c(
                "n_called", "n_missing", "n_0", "n_1", "n_2", "n_3"
            )
            if (!all(names(mcols(x)) == check_names))
                stop(
                    'Names of metadata columns must equal c("',
                    paste(check_names, collapse = '", "'),
                    '")'
                )
            df <- as.data.frame(mcols(x))
            ## Evaluate filter expression using data mask
            keep <- eval_tidy(filter, data = df)
            checks <- c(is.logical(keep), length(keep) == nrow(df))
            if (!all(checks))
                stop(
                    "`filter` expression must return a logical vector",
                    " of correct length"
                )
            gr <- x[keep,]
            seqlevels(gr) <- seqlevelsInUse(gr)  # In case lost with filter
            gr
        })

    }
)
