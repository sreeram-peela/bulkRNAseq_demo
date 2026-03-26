## re-annotate gene names and ids
library(org.Hs.eg.db)
library(dplyr)

deg_res <- readRDS('./deseq_demo/deg_results_Her2_normal.RDS')
head(deg_res)

# type of ids available
keytypes(org.Hs.eg.db)

# create a new df
dge_mapped_df <- data.frame(
  gene_symbol = mapIds(
    # Replace with annotation package for the organism relevant to your data
    org.Hs.eg.db,
    keys = rownames(deg_res),
    # Replace with the type of gene identifiers in your data
    keytype = "SYMBOL",
    # Replace with the type of gene identifiers you would like to map to
    column = "ENTREZID",
    # This will keep only the first mapped value for each Ensembl ID
    multiVals = "first"
  )
)

deg_res$symbol <- rownames(deg_res)

dge_mapped_df$symbol <- rownames(dge_mapped_df)
colnames(dge_mapped_df) <- c('ENTREZ', 'symbol')
deg_res_merged <- merge(x = deg_res, y = dge_mapped_df, by.x = 'symbol', by.y = 'symbol')
