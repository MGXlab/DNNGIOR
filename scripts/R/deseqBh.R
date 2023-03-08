library(DESeq2)
setwd("C:/Users/danie/OneDrive/Documentos/GitHub/hungerGamesModel/files/strainSummaries/bh/genes")

library("tximport")
library("readr")
library("tximportData")


samples <- read.table("bh_samples_t14vst32.txt", header=TRUE, row.names = 1)
samples$condition <- factor(samples$condition)

cts <- as.matrix(read.csv('bh_counts_t14vst32.txt',sep="\t",row.names="geneID"))

all(rownames(samples) == colnames(cts))

dds <- DESeqDataSetFromMatrix(countData = round(cts),
                              colData = samples,
                              design = ~ condition)


dds <- DESeq(dds)
res <- results(dds, ,alpha = 0.05, lfcThreshold = 1.5)
dds <- estimateSizeFactors(dds)
counts.norm <- counts(dds, normalized = T)

write.table(as.data.frame(res), 
            file="t14vst32_deseq.txt", sep='\t')



samples <- read.table("bh_samples_t14vst72.txt", header=TRUE, row.names = 1)
samples$condition <- factor(samples$condition)

cts <- as.matrix(read.csv('bh_counts_t14vst72.txt',sep="\t",row.names="geneID"))

all(rownames(samples) == colnames(cts))

dds <- DESeqDataSetFromMatrix(countData = round(cts),
                              colData = samples,
                              design = ~ condition)


dds <- DESeq(dds)
res <- results(dds,alpha = 0.05, lfcThreshold = 1.5)
dds <- estimateSizeFactors(dds)
counts.norm <- counts(dds, normalized = T)

write.table(as.data.frame(res), 
            file="t14vst72_deseq.txt", sep='\t')


samples <- read.table("bh_samples_t32vst72.txt", header=TRUE, row.names = 1)
samples$condition <- factor(samples$condition)

cts <- as.matrix(read.csv('bh_counts_t32vst72.txt',sep="\t",row.names="geneID"))

all(rownames(samples) == colnames(cts))

dds <- DESeqDataSetFromMatrix(countData = round(cts),
                              colData = samples,
                              design = ~ condition)


dds <- DESeq(dds)
res <- results(dds,alpha = 0.05, lfcThreshold = 1.5)
dds <- estimateSizeFactors(dds)
counts.norm <- counts(dds, normalized = T)

write.table(as.data.frame(res), 
            file="t32vst72_deseq.txt", sep='\t')

