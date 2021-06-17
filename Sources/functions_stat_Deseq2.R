#R script for Deseq2 metabolic functions – Monkey

library("DESeq2")
library(ggplot2)
library(FactoMineR)
monkey_countData <- read.csv('~/work/Monkey/Deseq/Deseq_monkey.csv', header = TRUE, sep = "\t", row.names=1)
head(monkey_countData)

metadata <- read.csv('~/work/Monkey/Deseq/Deseq_metadata.tsv', header = TRUE, sep = "\t", row.names=1)
head(metadata)


dds <- DESeqDataSetFromMatrix(countData=round(monkey_countData), 
                              colData=metadata, 
                              design=~Region)
dds
dds <- DESeq(dds)
res <- results(dds)
head(results(dds, tidy=TRUE))
summary(res) #summary of results
#Sort summary list by p-value
res <- res[order(res$padj),]
head(res)

