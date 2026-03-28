### GSVA enrichment
# load the datasets if not available in the env

BiocManager::install("GSVA", update = FALSE)
BiocManager::install("qusage", update = FALSE)

## load the required libs
library(DESeq2)
library(dplyr)
library(ggplot2)
library(pheatmap)
library(GSVA)

# load the deseq object and perform VST
dds <- readRDS(file = './deseq_demo/tcga_dds_wald.RDS')
dds_norm <- vst(dds, blind = TRUE)
vst_mat <- assay(dds_norm) %>%
  as.matrix() 

# get the gene set (hallmark genes)
hallmark_gene_sets <- msigdbr::msigdbr(
  species = "Homo sapiens", # Can change this to what species you need
  category = "H" # Only hallmark gene sets
)
hallmarks_list <- split(
  hallmark_gene_sets$gene_symbol, # The genes we want split into pathways
  hallmark_gene_sets$gs_name # The pathways made as the higher levels of the list
)

gsva_par <- gsvaParam(
  exprData = vst_mat,
  geneSets = hallmarks_list,
  kcdf = "Gaussian",
  minSize = 15,
  maxSize = 500
)
# run gsva
gsva_results <- gsva(gsva_par)
# get the results
gsva_results_mat <- gsva_results %>%
  as.data.frame() %>%
  as.matrix()
gsva_sub <- gsva_results_mat[1:10,]
# plot the heatmap
pheatmap(mat=gsva_sub, cluster_rows = TRUE, cluster_cols = TRUE,
         show_colnames = FALSE)
