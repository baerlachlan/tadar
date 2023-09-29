#' @title
#' tadar: A package for Differential Allelic Representation (DAR) analysis
#'
#' @description
#' This package enables DAR analysis by providing functions that address
#' discrete steps of the analysis workflow.
#'
#' DAR analysis is intended to be performed using functions in the following
#' order:
#'
#' 1. [readGenotypes()] parses a multi-sample VCF file and returns a `GRanges`
#' object containing only the data that is required for DAR analysis.
#' 2. [countAlleles()] summarises the alleles from genotype data at each range
#' for each sample group.
#' 3. [filterLoci()] removes ranges that do not match a specified criterion.
#' 4. [countsToProps()] normalises the allele counts to account for missing
#' data and sample groups of different sizes.
#' 5. [dar()] calculates the DAR between two sample groups.
#' 6. [flipRanges()] is an optional step that enables the conversion of ranges
#' output by [dar()] from origins to regions, or vice versa.
#' 7. [assignFeatureDar()] assigns DAR values to features of interest.
#'
#' `tadar` also provides visualisation functions that allow quick inspection
#' of DAR within the dataset:
#'
#' - [plotChrDar()] produces a `Gviz` plot that displays the trend in DAR across
#' a chromosome.
#' - [plotDarECDF()] produces a `ggplot2` figure comparing DAR between
#' chromosomes.
#'
#' @author
#' Lachlan Baer, Stevie Pederson
#'
#' @docType package
#' @name tadar-package
NULL

## Supress R CMD check note
n_called <- n_missing <- NULL
Chromosome <- NULL
