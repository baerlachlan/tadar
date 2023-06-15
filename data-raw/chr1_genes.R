library(AnnotationHub)
library(stringr)
ens_species <- "Danio rerio"
ens_release <- "101"
ah <- AnnotationHub()
ah <- subset(ah, species == ens_species)
ah <- subset(ah, rdataclass == "EnsDb")
ahId <- ah$ah_id[str_detect(ah$title, ens_release)]
ensDb <- ah[[ahId]]
chr1_genes <- genes(ensDb, filter = SeqNameFilter("1"))
names(chr1_genes) <- NULL
chr1_genes <- chr1_genes[,c("gene_id", "gene_name")]
save(chr1_genes, file = "data/chr1_genes.rda", compress = TRUE)
