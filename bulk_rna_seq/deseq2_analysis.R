
# Command line args
args <- commandArgs(trailingOnly = TRUE)
output_dir <- args[1]
kallisto_dir <- args[2]
gtf_file <- args[3]

# Make output dir
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Load libraries
library(DESeq2)
library(tximport)
library(ggplot2)
library(pheatmap)
library(dplyr)
library(tibble)
library(GenomicFeatures)

# 1. Paths and metadata
data_path <- kallisto_dir
sample_dirs <- list.dirs(path = data_path, recursive = FALSE)
sample_dirs <- sample_dirs[grep("SRR", sample_dirs)]

sample_names <- basename(sample_dirs)
files <- file.path(sample_dirs, "abundance.tsv")
names(files) <- sample_names

# 2. Generate table with estimated counts
sample_table <- data.frame(
  sample = sample_names,
  group = factor(rep(c("Hpf192", "Hpf48", "Hpf24"), each = 2)),
  replicate = rep(1:2, 3)
)

rownames(sample_table) <- sample_names

# 3. Tximport (kallisto to DESeq2)
# Generate tx2gene from GTF file
txdb <- makeTxDbFromGFF(gtf_file)
tx2gene <- transcripts(txdb, columns = c("tx_name", "gene_id"))

tx2gene_df <- as.data.frame(tx2gene)[, c("tx_name", "gene_id")]
tx2gene_df$gene_id <- as.character(tx2gene_df$gene_id)
txi <- tximport(files, type = "kallisto", tx2gene = tx2gene_df)

# 4. Make DESeqDataSet object
dds <- DESeqDataSetFromTximport(txi,
                                colData = sample_table,
                                design = ~ group)

# 5. Filtering
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

# 6. DE analysis
dds <- DESeq(dds)

# 7. PCA analysis
vsd <- vst(dds, blind = FALSE)

# PCA plot
pca_data <- plotPCA(vsd, intgroup = "group", returnData = TRUE)
percentVar <- round(100 * attr(pca_data, "percentVar"))

pca_plot <- ggplot(pca_data, aes(PC1, PC2, color = group, label = name)) +
  geom_point(size = 3) +
  geom_text(hjust = 0, vjust = 0, nudge_x = 1) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed() +
  ggtitle("PCA plot") +
  theme_minimal()

ggsave(file.path(output_dir, "pca_plot.png"), pca_plot, width = 8, height = 6, dpi = 300)

# 8. DE
# 24Hpf vs 48Hpf
res_24Hpf_vs_48Hpf <- results(dds, contrast = c("group", "Hpf24", "Hpf48"))
res_24Hpf_vs_48Hpf <- res_24Hpf_vs_48Hpf[order(res_24Hpf_vs_48Hpf$padj), ]

# 24Hpf vs 192Hpf
res_24Hpf_vs_192Hpf <- results(dds, contrast = c("group", "Hpf24", "Hpf192"))
res_24Hpf_vs_192Hpf <- res_24Hpf_vs_192Hpf[order(res_24Hpf_vs_192Hpf$padj), ]

# 9. Volcano plots
# 24Hpf vs 48Hpf
res_24_48_df <- as.data.frame(res_24Hpf_vs_48Hpf)
res_24_48_df$significant <- res_24_48_df$padj < 0.05 & abs(res_24_48_df$log2FoldChange) > 1

volcano_24_48 <- ggplot(res_24_48_df, aes(x = log2FoldChange, y = -log10(pvalue), color = significant)) +
  geom_point(alpha = 0.6) +
  scale_color_manual(values = c("gray", "red")) +
  theme_minimal() +
  ggtitle("Volcano plot: 24Hpf vs 48Hpf")

ggsave(file.path(output_dir, "volcano_24Hpf_vs_48Hpf.png"), volcano_24_48, width = 8, height = 6, dpi = 300)

# 24Hpf vs 192Hpf
res_24_192_df <- as.data.frame(res_24Hpf_vs_192Hpf)
res_24_192_df$significant <- res_24_192_df$padj < 0.05 & abs(res_24_192_df$log2FoldChange) > 1

volcano_24_192 <- ggplot(res_24_192_df, aes(x = log2FoldChange, y = -log10(pvalue), color = significant)) +
  geom_point(alpha = 0.6) +
  scale_color_manual(values = c("gray", "red")) +
  theme_minimal() +
  ggtitle("Volcano plot: 24Hpf vs 192Hpf")

ggsave(file.path(output_dir, "volcano_24Hpf_vs_192Hpf.png"), volcano_24_192, width = 8, height = 6, dpi = 300)

# 10. Save results
# Write tables
write.csv(as.data.frame(res_24Hpf_vs_48Hpf), 
          file = file.path(output_dir, "DESeq2_results_24Hpf_vs_48Hpf.csv"))
write.csv(as.data.frame(res_24Hpf_vs_192Hpf), 
          file = file.path(output_dir, "DESeq2_results_24Hpf_vs_192Hpf.csv"))

# Save RDS 
saveRDS(dds, file = file.path(output_dir, "dds_object.rds"))
saveRDS(vsd, file = file.path(output_dir, "vsd_object.rds"))

cat("DESeq2 analysis completed successfully!\n")
cat("Results saved to:", output_dir, "\n")
cat("Significant genes (padj < 0.05 & |log2FC| > 1):\n")
cat("24Hpf vs 48Hpf:", sum(res_24Hpf_vs_48Hpf$padj < 0.05 & abs(res_24Hpf_vs_48Hpf$log2FoldChange) > 1, na.rm = TRUE), "\n")
cat("24Hpf vs 192Hpf:", sum(res_24Hpf_vs_192Hpf$padj < 0.05 & abs(res_24Hpf_vs_192Hpf$log2FoldChange) > 1, na.rm = TRUE), "\n")
