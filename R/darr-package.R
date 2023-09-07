#' @title
#' darr: A package for Differential Allelic Representation (DAR) analysis
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
#' 3. [countsToProps()] normalises the allele counts and filters ranges based
#' on default criteria.
#'     - Optionally, the user can specify their own custom filter with
#'     [filterLoci()]
#' 4. [dar()] calculates the DAR between two sample groups.
#' 5. [flipRanges()] is an optional step that enables the conversion of ranges
#' output by [dar()] from origins to regions, or vice versa.
#' 6. [assignFeatureDar()] assigns DAR values to features of interest.
#'
#' `darr` also provides visualisation functions that allow quick inspection
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
#' @name darr-package
NULL

## Supress R CMD check note
n_called <- n_missing <- NULL
Chromosome <- NULL
