# `chr1.vcf.bgz` data

The data in `chr1.vcf.bgz` is a subset of variant calls resulting from processed RNA-seq data.
`chr1.vcf.bgz.tbi` is the accompanying index file.

The raw data can be found in `FASTQ` format under GEO accession [GSE244310](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE244310). 

Further information on the biology and generation of the RNA-seq data currently exists as a pre-print: [Zebrafish models of Mucopolysaccharidosis types IIIA, B, & C show hyperactivity and changes in oligodendrocyte state](https://www.biorxiv.org/content/10.1101/2023.08.02.550904v1).

Raw `FASTQ` files were processed using the Snakemake workflow located at [this GitHub repository](https://github.com/karissa-b/2022_MPSIII_3mBrainRNAseq/tree/main/code/analysis-variants_AC), resulting in a multi-sample `VCF` file.
The `VCF` does not exist on the repository due to its large size, however can be recreated by applying the Snakemake workflow on the raw data located at the GEO accession above.

The `VCF` was then modified and subset to a single chromosome (Chromosome 1) to reduce size using the following R script. Note that the file paths in this script will need to be updated by the user.

```r
library(tidyverse)
library(VariantAnnotation)
library(Rsamtools)

## These sample names my need to be changed if using the filenames directly
## from GEO
## Refer to the metadata sheet on GEO for associated filenames
samples <- c(
    "22-01840_S2", "22-01845_S7", "22-01847_S9", "22-01848_S10",
    "22-01857_S19", "22-01858_S20", "22-01841_S3", "22-01844_S6",
    "22-01849_S11", "22-01850_S12", "22-01853_S15", "22-01854_S16",
    "22-01856_S18"
)
## Update the following path to the VCF produced by the Snakemake workflow
vcf_path <- file.path("/path/to/original/VCF.vcf.gz")
svp <- ScanVcfParam(
    fixed = c("ALT", "QUAL", "FILTER"),
    info = NA,
    geno = "GT",
    samples = samples
)
## Read VCF, subset, and remove unnecessary data
vcf <- readVcf(vcf_path, param = svp, genome = "GRCz11")
vcf <- vcf[seqnames(rowRanges(vcf)) == "1"]
seqlevels(vcf) <- seqlevelsInUse(vcf)
dimnames(vcf) <- list(
    NULL,
    paste0("sample", 1:13) # Rename samples
)
header(vcf)@header$contig <- NULL
header(vcf)@header$GATKCommandLine <- NULL
header(vcf)@header$reference <- NULL
header(vcf)@header$source <- NULL
header(vcf)@header$source.1 <- NULL
header(vcf)@header$INFO <- NULL
header(vcf)@header$FORMAT <- dplyr::filter(
    as.data.frame(header(vcf)@header$FORMAT),
    Description == "Genotype"
)
## Path to save the subsetted VCF
save_path <- "/path/to/modified/VCF.vcf"
writeVcf(vcf, save_path, index = FALSE)
## The saved VCF has the assembly quoted in its header (i.e. assembly="GRCz11")
## This will be parsed incorrectly with VariantAnnotation::readVCF()
## The following is a workaround to correct this
lines <- readLines(save_path)
ind <- which(lines == "##contig=<ID=1,length=59578282,assembly=\"GRCz11\">")
lines[ind] <- gsub('\\"', "", lines[ind])
writeLines(lines, con = save_path)
## Now compress and index
bgzip(save_path)
file.remove(save_path)
indexVcf(paste0(save_path, ".bgz"))
```

