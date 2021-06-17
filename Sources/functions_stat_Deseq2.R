#R script for Deseq2 metabolic functions – Monkey

library("DESeq2")
library(ggplot2)
library(FactoMineR)

# working directory
wd <- list()

# The path for the data in my working directory (must be changed according to your environment)
wd$data   <- "C:/Users/myUser/data/"


monkey_countData <- read.csv(paste(wd$data, 'Deseq_monkey.csv', sep=""), header = TRUE, sep = "\t", row.names=1)
head(monkey_countData)

metadata <- read.csv(paste(wd$data, 'Deseq_metadata.tsv', sep=""), header = TRUE, sep = "\t", row.names=1)
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

