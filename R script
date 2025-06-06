# Install necessary packages (uncomment if needed)
# install.packages(c("RColorBrewer", "lattice", "ggplot2", "limma"))
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("DESeq2")
# BiocManager::install("genefilter")
# BiocManager::install("ashr")
# BiocManager::install("vsn")

# Set working directory - replace with your working directory path
setwd('YOUR_WORKING_DIRECTORY_PATH')

# Load required libraries
library(DESeq2)
library(ggplot2)
library(ashr)
library(limma)
library(vsn)

# FUNCTIONS

# Function to plot PCA with sample names
plotPCAWithSampleNames = function(x, targets=targets, intgroup=colnames(targets)[1], ntop=500)
{
  library(RColorBrewer)
  library(genefilter)
  library(lattice)
  
  # Calculate row variances for feature selection
  rv = rowVars(x)
  select = order(rv, decreasing=TRUE)[seq_len(min(ntop, length(rv)))]
  
  # Perform PCA on selected features
  pca = prcomp(t(x[select,]))
  
  # Proportion of variance explained by PCs
  variance = pca$sdev^2 / sum(pca$sdev^2)
  variance = round(variance, 3) * 100
  
  # Sample names
  names = colnames(x)
  
  # Factor grouping for plotting
  fac = factor(apply(as.data.frame(targets[, intgroup, drop=FALSE]), 1, paste, collapse=" : "))
  
  # Colors depending on number of groups
  if( nlevels(fac) >= 10 )
    colors = rainbow(nlevels(fac))
  else if( nlevels(fac) >= 3 )
    colors = brewer.pal(nlevels(fac), "Set1")
  else
    colors = c( "dodgerblue3", "firebrick3" )
  
  # Plot PCA using lattice
  xyplot(
    PC2 ~ PC1, groups=fac, data=as.data.frame(pca$x), pch=16, cex=1.5,
    aspect = "fill",
    col = colors,
    xlab = list(paste("PC1 (", variance[1], "%)", sep=""), cex=0.8),
    ylab = list(paste("PC2 (", variance[2], "%)", sep=""), cex=0.8),
    panel = function(x, y, ...) {
      panel.xyplot(x, y, ...);
      ltext(x=x, y=y, labels=names, pos=1, offset=0.8, cex=0.7)
    },
    main = draw.key(
      key = list(
        rect = list(col = colors),
        text = list(levels(fac)),
        rep = FALSE
      )
    )
  )
}

# Load count data and design table - replace filenames as needed
gene_count_data <- read.csv("gene_count_matrix.csv", row.names = 1)
transcript_count_data <- read.csv("transcript_count_matrix.csv", row.names = 1)
colData <- read.csv("Design.Group.Batch.csv", row.names = 1)

# DESeq2 dataset creation
dds_gene <- DESeqDataSetFromMatrix(countData = gene_count_data, colData = colData, design = ~ group)
dds_trans <- DESeqDataSetFromMatrix(countData = transcript_count_data, colData = colData, design = ~ group)

# Prefiltering: keep genes/transcripts with at least 10 reads in at least 3 samples
filter_gene <- rowSums(counts(dds_gene) >= 10) >= 3
dds_gene <- dds_gene[filter_gene, ]

filter_trans <- rowSums(counts(dds_trans) >= 10) >= 3
dds_trans <- dds_trans[filter_trans, ]

# Run DESeq2
dds_gene <- DESeq(dds_gene)
dds_trans <- DESeq(dds_trans)

# Generate dispersion plots
png("dispersion_plot_gene.png", width = 1000, height = 800)
plotDispEsts(dds_gene)
dev.off()

png("dispersion_plot_transcript.png", width = 1000, height = 800)
plotDispEsts(dds_trans)
dev.off()

# rlog transformation
rld_gene <- rlog(dds_gene, blind = TRUE)
rld_trans <- rlog(dds_trans, blind = TRUE)

# Log2 transform raw and normalized counts matrices
lgc.raw_gene <- log2(counts(dds_gene, normalized = FALSE) + 1)
lgc.norm_gene <- log2(counts(dds_gene, normalized = TRUE) + 1)

lgc.raw_trans <- log2(counts(dds_trans, normalized = FALSE) + 1)
lgc.norm_trans <- log2(counts(dds_trans, normalized = TRUE) + 1)

# PCA plots (basic)
png("pca_gene_plot.png", width = 800, height = 600)
plotPCA(rld_gene, intgroup = c("group"))
dev.off()

png("pca_transcript_plot.png", width = 800, height = 600)
plotPCA(rld_trans, intgroup = c("group"))
dev.off()

# PCA plots with sample names
png("pca_gene_plot_with_sample_names.png", width = 600, height = 500)
plotPCAWithSampleNames(assay(rld_gene), targets = colData(dds_gene), intgroup = 'group')
dev.off()

png("pca_transcript_plot_with_sample_names.png", width = 800, height = 600)
plotPCAWithSampleNames(assay(rld_trans), targets = colData(dds_trans), intgroup = 'group')
dev.off()

# Mean-SD plots
png("SD_vs_Mean_log2.png", width = 800, height = 800)
meanSdPlot(lgc.norm_gene)
dev.off()

png("SD_vs_Mean_rlog.png", width = 800, height = 800)
meanSdPlot(assay(rld_gene))
dev.off()

# Differential expression analysis for contrasts B vs A and C vs A with LFC=0
res_BvsA <- results(dds_gene, contrast = c("group", "B", "A"))
res_CvsA <- results(dds_gene, contrast = c("group", "C", "A"))

# Differential expression with LFC threshold 1
res_BvsA_LFC1 <- results(dds_gene, contrast = c("group", "B", "A"), lfcThreshold = 1)
res_CvsA_LFC1 <- results(dds_gene, contrast = c("group", "C", "A"), lfcThreshold = 1)

# MA plots for LFC=0 and LFC=1
png("MA_plots_LFC0.png", width = 1000, height = 800)
par(mfrow = c(1, 2))
DESeq2::plotMA(res_BvsA, alpha = 0.05)
abline(h = c(-2, 2), col = "black")
title("MA Plot: B vs A (LFC=0)")
DESeq2::plotMA(res_CvsA, alpha = 0.05)
abline(h = c(-2, 2), col = "black")
title("MA Plot: C vs A (LFC=0)")
dev.off()

png("MA_plots_LFC1.png", width = 1000, height = 800)
par(mfrow = c(1, 2))
DESeq2::plotMA(res_BvsA_LFC1, alpha = 0.05)
abline(h = c(-2, 2), col = "black")
title("MA Plot: B vs A (LFC≥1)")
DESeq2::plotMA(res_CvsA_LFC1, alpha = 0.05)
abline(h = c(-2, 2), col = "black")
title("MA Plot: C vs A (LFC≥1)")
dev.off()

# Shrink log2 fold changes with ashr method
res_BvsA_shrunk <- lfcShrink(dds_gene, coef = "group_B_vs_A", type = "ashr")
res_CvsA_shrunk <- lfcShrink(dds_gene, coef = "group_C_vs_A", type = "ashr")

# Subset for LFC ≥ 1
res_BvsA_shrunk_LFC1 <- res_BvsA_shrunk[abs(res_BvsA_shrunk$log2FoldChange) >= 1, ]
res_CvsA_shrunk_LFC1 <- res_CvsA_shrunk[abs(res_CvsA_shrunk$log2FoldChange) >= 1, ]

# Apply adjusted p-value cutoff <= 0.05
res_BvsA_shrunk_LFC1_padj <- res_BvsA_shrunk_LFC1[!is.na(res_BvsA_shrunk_LFC1$padj) & res_BvsA_shrunk_LFC1$padj <= 0.05, ]
res_CvsA_shrunk_LFC1_padj <- res_CvsA_shrunk_LFC1[!is.na(res_CvsA_shrunk_LFC1$padj) & res_CvsA_shrunk_LFC1$padj <= 0.05, ]

# Sort by log2FoldChange descending
res_BvsA_shrunk_LFC1_sorted <- res_BvsA_shrunk_LFC1_padj[order(res_BvsA_shrunk_LFC1_padj$log2FoldChange, decreasing = TRUE), ]
res_CvsA_shrunk_LFC1_sorted <- res_CvsA_shrunk_LFC1_padj[order(res_CvsA_shrunk_LFC1_padj$log2FoldChange, decreasing = TRUE), ]

# Save significant genes
write.csv(res_BvsA_shrunk_LFC1_sorted, file = "Significant_Genes_BvsA.csv", quote = FALSE)
write.csv(res_CvsA_shrunk_LFC1_sorted, file = "Significant_Genes_CvsA.csv", quote = FALSE)

# MA plots for shrunken LFC=0 and LFC≥1
png("MA_plot_BvsA_CvsA_shrunk_LFC0.png", width = 1000, height = 800)
par(mfrow = c(1, 2))
DESeq2::plotMA(res_BvsA_shrunk, alpha = 0.05)
abline(h = c(-2, 2), col = "black")
title("MA Plot: Shrunken B vs A (LFC=0)")
DESeq2::plotMA(res_CvsA_shrunk, alpha = 0.05)
abline(h = c(-2, 2), col = "black")
title("MA Plot: Shrunken C vs A (LFC=0)")
dev.off()

png("MA_plot_BvsA_CvsA_shrunk_LFC1.png", width = 1000, height = 800)
par(mfrow = c(1, 2))
DESeq2::plotMA(res_BvsA_shrunk_LFC1, alpha = 0.05)
abline(h = c(-2, 2), col = "black")
title("MA Plot: Shrunken B vs A (LFC≥1)")
DESeq2::plotMA(res_CvsA_shrunk_LFC1, alpha = 0.05)
abline(h = c(-2, 2), col = "black")
title("MA Plot: Shrunken C vs A (LFC≥1)")
dev.off()

# Remove batch effect using limma
mydesign <- model.matrix(design(dds_gene), colData(dds_gene))
b.corrected <- limma::removeBatchEffect(assay(rld_gene), batch = colData(dds_gene)$batch, design = mydesign)

# PCA plot of batch-corrected data with sample names
png("PCA_plot_batchcorrected_with_sample_names.png", width = 800, height = 600)
plotPCAWithSampleNames(b.corrected, targets = colData, intgroup = 'group')
dev.off()

# Save batch-corrected data
write.csv(b.corrected, file = "BatchCorrected_Rlog.csv", quote = FALSE)
