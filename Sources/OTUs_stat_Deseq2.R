
#Deseq difference between OTUs (phylum and families) AND Regions

library("DESeq2")
library(ggplot2)
library(FactoMineR)
family_abundance <- read.csv('~/work/Monkey/Deseq/Deseq_family.tsv', header = TRUE, sep = "\t", row.names=1)
head(family_abundance)
phylum_abundance <- read.csv('~/work/Monkey/Deseq/Deseq_phylum.tsv', header = TRUE, sep = "\t", row.names=1)
head(phylum_abundance)

metadata <- read.csv('~/work/Monkey/Deseq/Deseq_metadata.tsv', header = TRUE, sep = "\t", row.names=1)
head(metadata)

dds_family <- DESeqDataSetFromMatrix(countData=round(family_abundance), 
                                     colData=metadata, 
                                     design=~Region)

dds_phylum <- DESeqDataSetFromMatrix(countData=round(phylum_abundance), 
                                     colData=metadata, 
                                     design=~Region)

dds_family <- DESeq(dds_family)
dds_phylum <- DESeq(dds_phylum)

res_family <- results(dds_family)
res_phylum <- results(dds_phylum)

#Sort summary list by p-value
p_family <- res_family[order(res_family$padj),]
head(p_family)

p_phylum <- res_phylum[order(res_phylum$padj),]
head(p_phylum)

