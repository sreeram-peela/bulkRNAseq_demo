### GSEA using ranked list of genes
library(dplyr)
library(clusterProfiler)
library(ggplot2)
deg_res <- readRDS(file = './deseq_demo/deg_results_Her2_normal.RDS')
head(deg_res)

# get ranked list
deg_res$rank_val <- abs(deg_res$log2FoldChange) * -log10(deg_res$padj)
deg_res$gene <- rownames(deg_res)
library(tibble)
lfc_vec <- deg_res %>%
  arrange(desc(log2FoldChange)) %>%
  select(gene, log2FoldChange) %>%
  filter(!is.na(log2FoldChange)) %>%
  deframe()

head(lfc_vec)


# get the gene lists
library(msigdbr)
hs_msigdb_df <- msigdbr(species = "Homo sapiens")
wikipath_gs <- hs_msigdb_df %>%
  filter(gs_subcollection == 'CP:WIKIPATHWAYS')
term_gene_df <- wikipath_gs %>%
  select(c('gs_name', 'gene_symbol'))
colnames(term_gene_df) <- c("term", "gene")

# run gsea
# Set the seed so our results are reproducible:
set.seed(1987)
gsea_results <- GSEA(
  geneList = lfc_vec, # Ordered ranked gene list
  minGSSize = 25, # Minimum gene set size
  maxGSSize = 500, # Maximum gene set set
  pvalueCutoff = 0.05, # p-value cutoff
  eps = 0, # Boundary for calculating the p value
  seed = TRUE, # Set seed to make results reproducible
  pAdjustMethod = "BH", # Benjamini-Hochberg correction
  TERM2GENE = term_gene_df
)
head(gsea_results@result)
gsea_result_df <- data.frame(gsea_results@result)
gsea_result_df %>%
  # This returns the 3 rows with the largest NES values
  dplyr::slice_max(NES, n = 3)
most_positive_nes_plot <- enrichplot::gseaplot(
  gsea_results,
  geneSetID = "WP_RETINOBLASTOMA_GENE_IN_CANCER",
  title = "WP: RETINOBLASTOMA and CANCER",
  color.line = "#0d76ff"
)

# Plot the GSEA results per term
enrichplot::gseaplot(
  gsea_results,
  geneSetID = "WP_EMBRYONIC_STEM_CELL_PLURIPOTENCY_PATHWAYS",
  title = "WP: EMBRYONIC STEM CELL PLURIPOTENCY PATHWAYS",
  color.line = "#0d76ff"
)
 
