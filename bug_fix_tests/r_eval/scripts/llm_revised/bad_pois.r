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
variants_per_gene<-aggregate(x = gemini_filtered$gene, by = list(gene = gemini_filtered$gene), FUN = length)
colnames(variants_per_gene) <- c("gene","variant_count")
variants_per_CDS<-rbind.data.frame(variants_per_gene[variants_per_gene$gene %in% CDS_size$gene,], data.frame(gene = setdiff(CDS_size$gene,variants_per_gene$gene), variant_count = 0))                   
sorted <- variants_per_CDS[order(as.character(variants_per_CDS$gene)),]
CDS_size_sorted <- CDS_size[order(as.character(CDS_size$gene)),]
variants_per_CDS_withsizes <- data.frame(sorted,length_proportion=CDS_size_sorted$length_proportion)

Cluster.Test <- function(NbrInTarget,TargetProbability){
N <- length(NbrInTarget)
TotalHits <- sum(NbrInTarget)
ExpectedHits <- TotalHits*TargetProbability
exp(ppois(NbrInTarget,ExpectedHits,lower=T,log=T)*N)
}

newdata <- data.frame(variants_per_CDS_withsizes, cluster.test.result = Cluster.Test(variants_per_CDS_withsizes$variant_count,variants_per_CDS_withsizes$length_proportion))
newdata <- newdata[order(newdata$cluster.test.result),]

output_file <- commandArgs(trailingOnly = TRUE)[1]
gene_value <- newdata[1, "gene"]
writeLines(paste0('"', gene_value, '"'), output_file)