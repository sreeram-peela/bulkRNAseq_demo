## Over-representation analysis using ClusterProfiler
BiocManager::install("clusterProfiler", update = FALSE)
BiocManager::install("ggupset", update = FALSE)
BiocManager::install("msigdbr", update = FALSE)
BiocManager::install("org.Hs.eg.db", update = FALSE)

# load the libs
library(clusterProfiler)
library(DESeq2)
library(ggupset)
library(msigdbr)
library(dplyr)
library(ggplot2)

# get the deg computed object
tcga_dds <- readRDS(file = './deseq_demo/tcga_dds_wald.RDS')
resultsNames(tcga_dds)

# # extract results for one comparison
# deg_res <- results(tcga_dds, contrast = c('paper_BRCA_Subtype_PAM50', 
#                                           'Her2', 'Normal'))
# deg_res <- deg_res %>%
#   as.data.frame()
# 
# head(deg_res)
# saveRDS(object = deg_res, file = './deseq_demo/deg_results_Her2_normal.RDS')

# get the human msigdb
hs_msigdb_df <- msigdbr(species = "Homo sapiens")
table(hs_msigdb_df$gs_subcollection) # all subcollections
# get a smaller gene set of wikipathways (test)
wikipath_gs <- hs_msigdb_df %>%
  filter(gs_subcollection == 'CP:WIKIPATHWAYS')
# head(wikipath_gs)

# create a df with term-gene 
term_gene_df <- wikipath_gs %>%
  select(c('gs_name', 'gene_symbol'))
colnames(term_gene_df) <- c("term", "gene")
table(term_gene_df$term) # counts of genes in all terms

# get a list of genes (N=130)
key_genes <- deg_res %>%
  filter(abs(log2FoldChange) > 4 & padj < 0.0001) %>%
  row.names()

# run enrichement OVA
wiki_ova <- enricher(
  gene = key_genes, # A vector of your genes of interest
  pvalueCutoff = 0.3, # Can choose a FDR cutoff
  pAdjustMethod = "BH", # Method to be used for multiple testing correction
  universe = rownames(deg_res),
  TERM2GENE = term_gene_df
)

# probe the results
wiki_result_df <- data.frame(wiki_ova@result)
wiki_result_df %>%
  dplyr::filter(p.adjust < 0.2)
# plot the results
enrichplot::dotplot(wiki_ova) # plot the result

#### Probe for other gene sets as exercises
