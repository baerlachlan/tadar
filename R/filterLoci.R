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
#' fl <- system.file("extdata", "chr1.vcf.bgz", package="darr")
#' genotypes <- readGenotypes(fl)
#' groups <- list(
#'     group1 = c("S2", "S7", "S9", "S10", "S19", "S20"),
#'     group2 = c("S3", "S6", "S11", "S12", "S15", "S16", "S18")
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
            checkNames <- c("n_called", "n_missing", "n_0", "n_1", "n_2", "n_3")
            if (!all(names(mcols(x)) == checkNames))
                stop(
                    'Names of metadata columns must equal c("',
                    paste(checkNames, collapse = '", "'),
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
