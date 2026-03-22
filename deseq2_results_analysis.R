### DEG analysis using DESEQ2
library(DESeq2)
library(dplyr)
library(ggplot2)
library(pheatmap)


## For the dataset, the BRCA subtypes were selected as the key group

## The contrasts - 'Normal' (ref group) and other subtypes as test groups

# tcga_raw <- readRDS(file = './deseq_demo/tcga_dds_raw.RDS')
# str(tcga_raw)

# tcga_dds <- DESeq(tcga_raw, parallel=TRUE)

tcga_dds <- readRDS(file = './deseq_demo/tcga_dds_wald.RDS')
resultsNames(tcga_dds)

# extract results for one comparison
deg_res <- results(tcga_dds, contrast = c('paper_BRCA_Subtype_PAM50', 
                                          'Her2', 'Normal'))
# get a summarised view
summary(deg_res)

# get degs
logFC_cutoff <- 2
adjp_cutoff <- 0.001

deg_tops <- deg_res %>%
  as.data.frame() %>%
  dplyr::filter(abs(log2FoldChange) >= logFC_cutoff & 
           padj < adjp_cutoff)
# rank genes based on logFC and padj
deg_tops$rank_vals <- abs(deg_tops$log2FoldChange) * -log10(deg_tops$padj)

# sort the genes deg based on the ranks
deg_tops_sorted <- deg_tops %>%
  dplyr::arrange(desc(rank_vals))

# get top 100 genes
top100_degs <- deg_tops_sorted %>%
  head(n=100) %>%
  rownames()

# save the top 100 degs to a text file
writeLines(top100_degs, con = './deseq_demo/top100_degs_her2_normal.txt',
           sep = '\n')

#### Visualisations
deg_res_df <- deg_res %>%
  as.data.frame()
# basic volcano plot
ggplot(data =  deg_res_df,
       aes(x = log2FoldChange, y = -log10(padj))) + 
  geom_point()

### add additional info
annotate_genes <- deg_tops_sorted %>%
  head(n=20) %>%
  rownames()

deg_res_df$hgnc_symbol <- rownames(deg_res_df)
deg_res_df <- deg_res_df %>%
  dplyr::mutate(
    regulation = dplyr::case_when(
      log2FoldChange > 2 & padj < 0.001 ~ "Upregulated",
      log2FoldChange < -2 & padj < 0.001 ~ "Downregulated",
      TRUE ~ "Not significant"
    ),
    label = ifelse(hgnc_symbol %in% annotate_genes, hgnc_symbol, NA)
  )

library(ggrepel) # for pretty labels

ggplot(deg_res_df,
       aes(x = log2FoldChange, y = -log10(padj), color = regulation)) + 
  geom_point(alpha = 0.7) + 
  theme_classic() + 
  
  # Threshold lines
  geom_hline(yintercept = -log10(0.001), linetype = 'dashed', color = 'red') + 
  geom_vline(xintercept = c(-2, 2), linetype = 'dashed', color = 'blue') + 
  
  # Labels (repelled to avoid overlap)
  geom_point(aes(color = regulation)) +
  geom_text_repel(aes(label = label, color = regulation), show.legend = FALSE) +
  
  # Custom colors + legend
  scale_color_manual(
    values = c(
      "Upregulated" = "green",
      "Downregulated" = "red",
      "Not significant" = "grey"
    )
  ) +
  
  labs(color = "DEG type")


## heatmap
deg_tops$hgnc_symbol <- rownames(deg_tops)
deg_tops <- deg_tops %>%
  dplyr::mutate(
    regulation = dplyr::case_when(
      log2FoldChange > 2 & padj < 0.001 ~ "Upregulated",
      log2FoldChange < -2 & padj < 0.001 ~ "Downregulated",
      TRUE ~ "Not significant"
    ),
    label = ifelse(hgnc_symbol %in% annotate_genes, hgnc_symbol, NA)
  )

top5_up <- deg_tops %>%
  filter(regulation == 'Upregulated') %>%
  dplyr::arrange(desc(rank_vals)) %>%
  head(n=5) %>%
  pull(hgnc_symbol)

top5_down <- deg_tops %>%
  filter(regulation == 'Downregulated') %>%
  dplyr::arrange(desc(rank_vals)) %>%
  head(n=5) %>%
  pull(hgnc_symbol)

top10_degs <- c(top5_up, top5_down)

# get normal and her2 samples
samples_interest_meta <- colData(tcga_dds) %>%
  as.data.frame() %>%
  filter(paper_BRCA_Subtype_PAM50 == 'Normal' | paper_BRCA_Subtype_PAM50 == 'Her2') 
# get normalised counts
# 1. Log transform (very important)
norm_counts <- counts(tcga_dds, normalized = TRUE) 
deg_norms <- norm_counts[top10_degs, rownames(samples_interest_meta)]
deg_norms_log <- log2(deg_norms + 1)

# 2. Optional: Z-score per gene (better visualization)
deg_norms_scaled <- t(scale(t(deg_norms_log)))

# 3. Annotation (example: phenotype)
annotation_col <- data.frame(
  condition = samples_interest_meta$paper_BRCA_Subtype_PAM50,
  race = samples_interest_meta$race
)
rownames(annotation_col) <- colnames(deg_norms_scaled)
annotation_col <- annotation_col %>%
  dplyr::arrange(condition)
# 4. Color palette
heat_colors <- colorRampPalette(c("green", "white", "red"))(50)

# 5. Plot
pheatmap(
  deg_norms_scaled[, rownames(annotation_col)],
  color = heat_colors,breaks = seq(-3, 3, length.out = 50),
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  show_rownames = TRUE,
  show_colnames = FALSE,
  fontsize_row = 8,
  annotation_col = annotation_col,
  border_color = NA,
  main = "Top DEGs Heatmap"
)
