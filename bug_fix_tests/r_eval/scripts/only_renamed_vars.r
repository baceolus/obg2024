# Identify genes where significantly more variants are observed than expected
# Output gene with most excess variants as measured by log likelihood.
## Cluster Analysis steps
# * Classify nonsense and frameshift mutations
# * Filter based on CADD > 20, [Combined Annotation Dependent Depletion (CADD)](https://cadd.gs.washington.edu/) method
# * Identify likely deleterious variant clusters in patient genes
# * Generate statistics on gene clustering

gemini_out <- read.table("r_eval/gemini_out.txt",sep = "\t", header = TRUE)
gemini_out$cadd_scaled <- as.numeric(as.character(gemini_out$cadd_scaled))
gemini_out$cadd_scaled[is.na(gemini_out$cadd_scaled)] <- 0
gemini_filtered <- gemini_out[(gemini_out$biotype=="protein_coding")& (gemini_out$cadd_scaled>=20),]
CDS_size <- read.table("r_eval/ENS75.txt",sep = "\t", header = TRUE)
gene_variants<-aggregate(x = gemini_filtered$gene, by = list(gene = gemini_filtered$gene), FUN = length)
colnames(gene_variants) <- c("gene","variant_count")
variants_per_CDS<-rbind.data.frame(gene_variants[gene_variants$gene %in% CDS_size$gene,], data.frame(gene = setdiff(CDS_size$gene,gene_variants$gene), variant_count = 0))                   
sorted_genes <- variants_per_CDS[order(as.character(variants_per_CDS$gene)),]
variants_per_CDS_withsizes <- data.frame(sorted_genes,length_proportion=CDS_size$length_proportion)

Cluster.Test <- function(NbrInTarget,TargetProbability){
N <- length(NbrInTarget)
TotalHits <- sum(NbrInTarget)
ExpectedHits <- TotalHits*TargetProbability
# N*ppois(NbrInTarget,ExpectedHits,lower=F)
1-exp(ppois(NbrInTarget,ExpectedHits,lower=T,log=T)*N)
}

combineddata <- data.frame(variants_per_CDS_withsizes, cluster.test.result = Cluster.Test(variants_per_CDS_withsizes$variant_count,variants_per_CDS_withsizes$length_proportion))
combineddata <- combineddata[order(combineddata$cluster.test.result),]

output_file <- commandArgs(trailingOnly = TRUE)[1]
gene_value <- combineddata[1, "gene"]
writeLines(paste0('"', gene_value, '"'), output_file)