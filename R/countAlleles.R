#' @title Count alleles within each experimental group
#'
#' @description Summarise the alleles from genotype calls at each single
#' nucleotide locus within each sample group
#'
#' @param x A GRanges object with metadata columns containing genotype
#' information for all samples
#' @param groups A named list specifying the sample grouping structure, where
#' each element contains a character vector of sample names
#'
#' @return A GRangesList containing a summary of allele counts at each range
#' within each sample group
#'
#' @examples
#' fl <- system.file("extdata", "chr1.vcf.gz", package="darr")
#' genotypes <- readGenotypes(fl)
#' sampleGroups <- list(
#'     group1 = c("S2", "S7", "S9", "S10", "S19", "S20"),
#'     group2 = c("S3", "S6", "S11", "S12", "S15", "S16", "S18")
#' )
#' alleleCounts <- countAlleles(genotypes, sampleGroups)
#' alleleCounts
#'
#' @import GenomicRanges
#' @importFrom S4Vectors 'mcols<-'
#' @rdname countAlleles-methods
#' @aliases countAlleles
#' @export
setMethod(
    "countAlleles",
    signature = signature(x = "GRanges", groups = "list"),
    function(x, groups) {

        gr <- granges(x)
        gt <- as.matrix(mcols(x))
        stopifnot(is.character(gt))
        grList <- lapply(groups, function(group) {
            stopifnot(is.character(group))
            gtInGroup <- gt[,group]
            ## Split so we can lapply()
            gtInGroup <- split(gtInGroup, seq(NROW(gtInGroup)))
            countsByLoc <- lapply(gtInGroup, function(genotypes) {
                alleles <- strsplit(genotypes, "/")
                alleles <- unlist(alleles)
                n_called <- sum(alleles != ".")
                n_nocall <- sum(alleles == ".")
                n_0 <- sum(alleles == "0")
                n_1 <- sum(alleles == "1")
                n_2 <- sum(alleles == "2")
                n_3 <- sum(alleles == "3")
                matrix(
                    data = c(n_called, n_nocall, n_0, n_1, n_2, n_3),
                    nrow = 1,
                    dimnames = list(
                        NULL,
                        c("n_called", "n_nocall", "n_0", "n_1", "n_2", "n_3")
                    )
                )
            })
            counts <- do.call(rbind, countsByLoc)
            mcols(gr) <- counts
            gr
        })
        GRangesList(grList)

    }
)
