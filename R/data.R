#' @title Genomic feature example data
#'
#' @description Gene features for example usage. Generation of this data is
#' documented in `system.file("data-raw/chr1_genes.R", package = "tadar")`.
#'
#' @return GRanges object with 1456 ranges and 2 metadata columns.
#' \describe{
#'   \item{chr1_genes}{Ranges represent gene features for chromosome 1 of
#'   zebrafish GRCz11 genome.}
#' }
#'
#' @source \url{https://www.ensembl.org}
#' @name chr1_genes
#' @rdname chr1_genes
#' @usage data(chr1_genes)
"chr1_genes"

#' @title Differential expression example data
#'
#' @description *p*-values from differential expressiong testing for
#'  example usage.
#'
#' @return A 716 x 5 tibble object.
#'
#' @name chr1_tt
#' @rdname chr1_tt
#' @usage data(chr1_tt)
"chr1_tt"
