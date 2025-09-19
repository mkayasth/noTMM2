library(uwot)
library(cluster)

## cluster: with cell of origin genes.
gene_sets <- list(
  SCP = c("PLP1", "MPZ", "CDH19", "ERBB3", "ERBB4"),
  bridge = c("ERBB4", "CDH9", "CTTNBP2", "ASCL1"),
  chromaffin = c("DBH", "TH", "CHGA", "DDC"),
  late_chromaffin = c("TH", "CHGA", "PNMT"),
  cycling_neuroblast = c("MKI67", "TOP2A", "ALK", "STMN2", "ISL1"),
  neuroblast = c("NEFM", "GAP43", "STMN2", "ISL1", "ALK"),
  late_neuroblast = c("SYN3", "IL7", "GAP43", "STMN2", "ISL1")
)



source("DataCleaning.R")

gsvapar <- gsvaParam(as.matrix(Expression), gene_sets, kcdf = "Gaussian")
gsva_result <- gsva(gsvapar)
GSVA_df <- as.data.frame(gsva_result)
GSVA_df <- t(GSVA_df)

# scaling.
GSVA_df_scaled <- scale(GSVA_df)

###

# umap and kmeans.
umap_res <- umap(GSVA_df_scaled)
set.seed(123)
km_res <- kmeans(umap_res, centers = 7, nstart = 50)


###
  
# Labeling cluster by dominant cell of origin.
cluster_means <- aggregate(GSVA_df, by = list(cluster = km_res$cluster), FUN = mean)


clusters <- km_res$cluster

# Pick a color palette (1 color per cluster)
colors <- rainbow(length(unique(clusters)))


plot(umap_res, col = colors[clusters], pch = 19,
     xlab = "UMAP 1", ylab = "UMAP 2",
     main = "K-means clusters on GSVA scores")
legend("topright", 
       legend = paste("Cluster", sort(unique(clusters))), 
       col = colors, pch = 19)

###################################################################################

## cluster: with ADR/MES genes.
gene_sets_phenotype <- list(
  ADR = c("DLK1", "DBH", "PHOX2A", "PHOX2B", "GATA2", "GATA3"),
  MES = c("VIM", "FN1", "YAP1", "SNAI2", "PRRX1", "WWTR1")
)

gene_sets_phenotype <- list(
  ADR = c("KLF7", "GATA3", "HAND2", "PHOX2A", "ISL1", "HAND1",
          "PHOX2B", "TFAP2B", "GATA2", "SATB1", "SIX3", "EYA1",
          "SOX11", "DACH1", "ASCL1", "HEY1", "KLF13", "PBX3"),
  MES = c("VIM", "FN1", "MEOX2", "ID1", "EGR3", "AEBP1", "CBFB", "IRF3", "IRF2", "IRF1", 
          "TBX18", "MAFF", "RUNX2", "ZFP36L1", "NR3C1", "BHLHE41", "GLIS3", "RUNX1", "FOSL1",
          "FOSL2", "ELK4", "IFI16", "SIX4", "FLI1", "MAML2", "SMAD3", "DCAF6", "WWTR1", "SOX9", 
          "MEF2D", "ZNF217", "PRRX1", "CREG1", "NOTCH2", "SIX1", "MEOX1"))


source("DataCleaning.R")

gsvapar <- gsvaParam(as.matrix(Expression), gene_sets_phenotype, kcdf = "Gaussian")
gsva_result <- gsva(gsvapar)
GSVA_df <- as.data.frame(gsva_result)
GSVA_df <- t(GSVA_df)

# scaling.
GSVA_df_scaled <- scale(GSVA_df)

###

# umap and kmeans.
umap_res <- umap(GSVA_df_scaled)
set.seed(123)
km_res <- kmeans(umap_res, centers = 4, nstart = 50)


###

# Labeling cluster by dominant cell of origin.
cluster_means <- aggregate(GSVA_df, by = list(cluster = km_res$cluster), FUN = mean)


clusters <- km_res$cluster


colors <- rainbow(length(unique(clusters)))


plot(umap_res, col = colors[clusters], pch = 19,
     xlab = "UMAP 1", ylab = "UMAP 2",
     main = "K-means clusters on GSVA scores")
legend("topright", 
       legend = paste("Cluster", sort(unique(clusters))), 
       col = colors, pch = 19)

