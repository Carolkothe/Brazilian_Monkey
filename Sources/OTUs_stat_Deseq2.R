
#Deseq difference between OTUs (phylum and families) AND Regions

library("DESeq2")
library(ggplot2)
library(FactoMineR)

# working directory
wd <- list()

# The path for the data in my working directory
wd$data   <- "C:/Users/myUser/data/"

family_abundance <- read.csv(paste(wd$data,'Deseq_family.tsv', sep=""), header = TRUE, sep = "\t", row.names=1)
head(family_abundance)
phylum_abundance <- read.csv(paste(wd$data, 'Deseq_phylum.tsv', sep=""), header = TRUE, sep = "\t", row.names=1)
head(phylum_abundance)

metadata <- read.csv(paste(wd$data, 'Deseq_metadata.tsv', sep=""), header = TRUE, sep = "\t", row.names=1)
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

