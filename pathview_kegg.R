### pathview to annotate KEGG pathways
BiocManager::install("pathview")

library(pathview)
library(ggplot2)
library(dplyr)

# load the deg results
deg_res <- read.csv('./deseq_demo/deg_res_her2_normal.csv', 
                    header = TRUE)


# get entrez ids for using with KEGG
gene_data <- deg_res$log2FoldChange
names(gene_data) <- deg_res$gene


## hs db
library(org.Hs.eg.db)
dge_mapped_df <- data.frame(
  gene_symbol = mapIds(
    # Replace with annotation package for the organism relevant to your data
    org.Hs.eg.db,
    keys = deg_res$gene,
    # Replace with the type of gene identifiers in your data
    keytype = "SYMBOL",
    # Replace with the type of gene identifiers you would like to map to
    column = "ENTREZID",
    # This will keep only the first mapped value for each Ensembl ID
    multiVals = "first"
  )
)
  
deg_res$symbol <- deg_res$X
  
  dge_mapped_df$symbol <- rownames(dge_mapped_df)
  colnames(dge_mapped_df) <- c('ENTREZ', 'symbol')
  deg_res_merged <- merge(x = deg_res, y = dge_mapped_df, by.x = 'symbol', by.y = 'symbol')

gene_data <- deg_res_merged$log2FoldChange
names(gene_data) <- deg_res_merged$ENTREZ
# hsa05224 breast cancer KEGG pathway
pathview(gene.data = gene_data, pathway.id = "05224",
         species = "hsa", out.suffix = "logFC_her2")
