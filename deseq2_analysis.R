#### DESEQ2 DEG analysis
library(cBioPortalData)
library(SummarizedExperiment)
library(dplyr)
# get the data
cbio <- cBioPortal()
# studies <- getStudies(cbio, buildReport = TRUE)
# head(studies)

# download studies as summarised experiment class
#### Study: BRCA TCGA pancancer
tcga_pan <- cBioDataPack(cancer_study_id = "brca_tcga_pan_can_atlas_2018",
                     use_cache = TRUE,
                     names.field = c("Hugo_Symbol", "Entrez_Gene_Id", "Gene"),
                     cleanup = TRUE, ask = FALSE)
tcga_pan_exp <- tcga_pan[['mrna_seq_v2_rsem']] # get the expression data

metadata_tcga <- colData(tcga_pan) %>%
  as.data.frame() %>%
  dplyr::filter(SAMPLE_ID %in% colnames(tcga_pan_exp))
# R expects a Bioconductor DataFrame
colData(tcga_pan_exp) <- S4Vectors::DataFrame(metadata_tcga)

rm(tcga_pan) # remove local vars to reduce mem usage

# DEG analysis using DESEQ2
library(DESeq2)

## create deseq2 object from the above object
counts_mat <- assay(tcga_pan_exp)
counts_mat <- round(counts_mat, digits = 0) # round off to integers
## Although RSEM values can be imported using tximport
# metadata cleaning
table(metadata_tcga$RACE)
table(metadata_tcga$SUBTYPE)
table(metadata_tcga$TISSUE_SOURCE_SITE)
table(metadata_tcga$OS_STATUS)

metadata_tcga_clean <- metadata_tcga %>%
  filter(
    !is.na(RACE),
    !is.na(SUBTYPE),
    !is.na(TISSUE_SOURCE_SITE),
    !is.na(OS_STATUS)
  )

length(intersect(metadata_tcga_clean$PATIENT_ID, colnames(counts_mat)))

counts_mat_filtered <- counts_mat[, metadata_tcga_clean$SAMPLE_ID]
rownames(metadata_tcga_clean) <- metadata_tcga_clean$SAMPLE_ID # align the data
all(colnames(counts_mat_filtered) == rownames(metadata_tcga_clean)) # TRUE
# 897 samples in total

# construct deseq object
library(DESeq2)
tcga_dds <- DESeqDataSetFromMatrix(countData = counts_mat_filtered,
                                   colData = metadata_tcga_clean,
                                   design = ~ RACE + TISSUE_SOURCE_SITE + 
                                     OS_STATUS + SUBTYPE) 
# delete all the previous objects and free mem using rm() and gc()
str(tcga_dds)

# enable parallelisation
library(BiocParallel)
register(SnowParam(workers = 6)) # MultiCoreParam() for linux

# single command: DESeq(tcga_dds) runs all the steps in a reproducible fashion
# We plan to do step-wise analyis
# steps: estimateSizeFactors (library size artefacts) → 
#  estimateDispersions (gene variances) → nbinomWaldTest (model fitting and testing)

# preprocessing deseq data
## remove low count genes
smallestGroupSize <- 3
keep <- rowSums(counts(tcga_dds) >= 10) >= smallestGroupSize
tcga_dds_filtered <- tcga_dds[keep,] # 2k genes removed
# check for batch effects
# common areas: tissue source site
# VST is the recommended normalisation method
# assign to a new object to preserve original data
vsd <- vst(object = tcga_dds_filtered, blind = FALSE)
plotPCA(vsd, intgroup = "TISSUE_SOURCE_SITE")
plotPCA(vsd, intgroup = "RACE")
plotPCA(vsd, intgroup = "OS_STATUS")
plotPCA(vsd, intgroup = "SUBTYPE") # more distinction

# by including the first 3, we are telling DESeq2 to adjust for differences
# before inferring deg for subtypes (our 1st goal) and then for survival
# (2nd goal)
dds <- estimateSizeFactors(tcga_dds_filtered)
dds <- estimateDispersions(dds)
dds <- nbinomWaldTest(dds)

# Extract results
# get the results for normal vs lumA
deg_lumA <- results(dds, contrast = c('SUBTYPE', 'BRCA_LumA', 'BRCA_Normal'))



