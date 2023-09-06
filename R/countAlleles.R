#' @title Count alleles within each experimental group
#'
#' @description Summarise the alleles from genotype calls at each single
#' nucleotide locus within each sample group.
#'
#' @param genotypes `GRanges` object with metadata columns containing genotype
#' information for all samples.
#' @param groups Named `list` specifying the sample grouping structure, where
#' each element contains a character vector of sample names.
#'
#' @return `GRangesList` containing a summary of allele counts at each range.
#' Each element of the list represents a distinct sample group.
#'
#' @examples
#' fl <- system.file("extdata", "chr1.vcf.bgz", package="darr")
#' genotypes <- readGenotypes(fl)
#' groups <- list(
#'     group1 = paste0("sample", 1:6),
#'     group2 = paste0("sample", 7:13)
#' )
#' countAlleles(genotypes, groups)
#'
#' @import GenomicRanges
#' @importFrom S4Vectors mcols 'mcols<-'
#' @rdname countAlleles-methods
#' @aliases countAlleles
#' @export
setMethod(
    "countAlleles",
    signature = signature(genotypes = "GRanges", groups = "list"),
    function(genotypes, groups) {

        gr <- granges(genotypes)
        genotypes <- as.matrix(mcols(genotypes))
        stopifnot(is.character(genotypes))
        grl <- lapply(groups, function(x) {
            stopifnot(is.character(x))
            gtInGroup <- genotypes[,x]
            ## Split so we can lapply()
            gtInGroup <- split(gtInGroup, seq(nrow(gtInGroup)))
            counts <- lapply(gtInGroup, function(y) {
                alleles <- strsplit(y, "/")
                alleles <- unlist(alleles)
                n_called <- sum(alleles != ".")
                n_missing <- sum(alleles == ".")
                n_0 <- sum(alleles == "0")
                n_1 <- sum(alleles == "1")
                n_2 <- sum(alleles == "2")
                n_3 <- sum(alleles == "3")
                c(n_called, n_missing, n_0, n_1, n_2, n_3)
            })
            counts <- do.call(rbind, counts)
            dimnames(counts) <- list(
                NULL,
                c("n_called", "n_missing", "n_0", "n_1", "n_2", "n_3")
            )
            mcols(gr) <- counts
            gr
        })
        GRangesList(grl)

    }
)
