# Set the working directory
setwd("C:/Users/dimit/Desktop/ergasia3_bioinformatics/data/workspace/htseq_output/")

# Packages
library("DESeq2")
library("ggplot2")

# Upload our data from the NGS analysis
origData <- read.csv("SRR10045016-17-18-19-20-21_counts.csv", sep = "\t", header = FALSE)
head(origData)


countdata <- origData[,2:7]
head(countdata)
sampleName <- c("rep1","rep2","rep3","_rep1","_rep2","_rep3")
sampleName <- colnames(countdata)
head(sampleName)
condition=as.factor(c("H","H","H","D","D","D")) ## set the header
coldata<-data.frame(sampleName,condition)

head(coldata)
ddsTable <- DESeqDataSetFromMatrix(countData = countdata, colData = coldata,
                                   design = ~ condition )
??DESeqDataSetFromMatrix

#########.....#########
dds <- DESeq(ddsTable)
dds
res <- results(dds)
res
summary(res)

# order by BaseMean
results_ordered <- res[order(res$baseMean, decreasing = TRUE), ]
results_ordered_DF <- data.frame(results_ordered)
head(results_ordered_DF)

# order by log2FoldChange
loG2 <- res[order(res$log2FoldChange, decreasing = TRUE), ]
loG2_DF <- data.frame(loG2)
head(loG2)

# order by padj
padj_results <- res[order(res$padj, decreasing = TRUE), ]
padj_DF <- data.frame(padj_results)
head(padj_DF)

### We are setting our criteria ###

# 1st criteria : padj < 0.01
padj.cutoff <- 0.01
padj_001 <- res[order(res$padj < padj.cutoff), ]
padj_001


# 2nd criteria: log2FoldChange > 0.58
# > 0.58 means 1.5 based on log2
log2.cutoff <- 0.58
log2_058 <- res[order(res$log2FoldChange > log2.cutoff), ]
log2_058


# Now we will concatenate our criteria creating a absolute value
intro_res <- res[which(res$padj < padj.cutoff), ]
# abs() method is used to get the absolute value that is positive value doesn't change and negative value converted into positive value.
res_final <- intro_res[which(abs(intro_res$log2FoldChange) > log2.cutoff), ]
res_final

deseq2ResDF <- as.data.frame(res)
# Examine this data frame
head(deseq2ResDF)

DESeq2::plotMA(res,ylim=c(-2,2))
##plot.default(deseq2ResDF$baseMean,deseq2ResDF$log2FoldChange,xlim=c(10,10000),ylim=c(-1,1),log = "x") ##same as plotMA(res)
DESeq2::plotMA(res_final, ylim=c(-2,2))

idx <- identify(res$baseMean, res$log2FoldChange)
rownames(res)[idx]
idx

rld <- rlog(dds)
plotPCA(rld)

# Differential Expressed Genes
library("ggplot2")

qwe <- results(dds)
names(qwe_DF)
labels(qwe)
qwe_DF <- as.data.frame(qwe)

# add a column of NAs
qwe_DF$diffexpressed <- "NO"

# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
qwe_DF$diffexpressed[qwe_DF$log2FoldChange > 0.6 & qwe_DF$pvalue < 0.05] <- "UP"

# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
qwe_DF$diffexpressed[qwe_DF$log2FoldChange < -0.6 & qwe_DF$pvalue < 0.05] <- "DOWN"

# Now write down the name of genes beside the points...
# Create a new column "delabel" to qwe_DF, that will contain the name of genes differentially expressed (NA in case they are not)
qwe_DF$delabel <- NA
qwe_DF$delabel[qwe_DF$diffexpressed != "NO"] <- rownames(qwe_DF)[qwe_DF$diffexpressed != "NO"]

# Visualization with ggplot package
ggplot(data=qwe_DF, aes(x=log2FoldChange, y=-log10(pvalue), col=diffexpressed, label= delabel)) + 
  geom_point() +
  theme_minimal() +
  geom_text(aes(label=delabel), nudge_y = 0.5)

label_df <- subset(qwe_DF, rownames(qwe_DF) == "2060")
label_df
