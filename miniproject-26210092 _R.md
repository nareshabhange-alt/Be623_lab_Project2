
setwd("~/be623_labpro")
# loading library
library(DESeq2)
library(DESeq2)
library(ggplot2)
library(pheatmap)
# Creating files
dir.create("plots", showWarnings = FALSE)
dir.create("results", showWarnings = FALSE)

samples <- read.csv("airway_samples.csv",  row.names = FALSE)
counts <- read.csv("airway_counts.csv", row.names = 1)

dim(counts)
str(samples)

# creating summary for p_bar
samples$lib_size <- colSums(counts)%>%
sam_data <- samples %>%
  group_by(group) %>%
  summarize(
            mean_library_size = mean(lib_size),
            number_of_samples = n(),
            group = "drop")

p_bar <- ggplot(samples, aes(x = run, y = lib_size, fill = group)) +
  geom_col() +
  labs(
    title = "Library size per sample",
    x = "RNA-seq sample (run)",
    y = "Total reads") +
 theme_minimal()
# saving bar
ggsave("plots/library_sizes.png", p_bar, width = 6, height = 4, dpi = 300)

# Building the DESeq2 and checking the sample
rownames(samples) <- samples$run
dds <- DESeqDataSetFromMatrix(
  countData = as.matrix(counts),
  colData = samples,
  design = ~ group
)
levels(dds$group)
# Removing the  low count gene
keep <- rowSums(counts(dds) >= 10) >= 2
dds <- dds[keep, ]
cat("Genes remaining after removing :", nrow(dds), "\n")

vsd <- vst(dds, blind = TRUE)
pca_data <- plotPCA(vsd, intgroup = "group", returnData = TRUE)
percent_var <- round(100 * attr(pca_data, "percentVar"), 1)
p_pca <- ggplot(pca_data, aes(PC1, PC2, color = group)) +
  geom_point(size = 4) +
  xlab(paste0("PC1: ", percent_var[1], "% variance")) +
  ylab(paste0("PC2: ", percent_var[2], "% variance")) +
  ggtitle("PCA of airway samples") +
  theme_bw()
# Saving the plot
ggsave("plots/PCA.png", p_bar, width = 6, height = 4, dpi = 300)

dds <- DESeq(dds)
write.csv(counts(dds, normalized = TRUE), "results/normalized_counts.csv", quote = FALSE)

trt <- "dmtrt"
untrt <- "utrt"
comparison <- paste0(trt, "_vs_", untrt)
res <- results(dds, contrast = c("group", "trt", "untrt"))
head (res)

summary(res)

res <- res[order(res$padj), ]
res_df <- as.data.frame(res)
res_df$gene <- rownames(res_df)

res_df$status <- "gene colomns"
res_df$status[res_df$padj < 0.05 & res_df$log2FoldChange >= 1] <- "Up"
res_df$status[res_df$padj < 0.05 & res_df$log2FoldChange >= 1] <- "Down"
table(res_df$status)

volcano_data <- res_df[!is.na(res_df$padj), ]
p_volcano <- ggplot(volcano_data,
                    aes(x = log2FoldChange, y = -log10(padj), color = status)) +
  geom_point(alpha = 0.6) +
  scale_color_manual(values = c("Up" = "firebrick",
                                "Down" = "steelblue",
                                "Not significant" = "grey70")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  labs(title = comparison,
       x = "log2 fold change",
       y = "-log10 adjusted p-value") +
  theme_bw()
p_volcano 
# saving volcano plot
ggsave(paste0("plots/volcano_", comparison, ".png"), p_volcano,
       width = 6, height = 5, dpi = 300)

res_df$ENSG00000120129




