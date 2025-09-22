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

# color palette for cluster.
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


source("DataCleaning.R")

gsvapar <- gsvaParam(as.matrix(Expression), gene_sets_phenotype, kcdf = "Gaussian")
gsva_result <- gsva(gsvapar)
GSVA_df <- as.data.frame(gsva_result)
GSVA_df <- t(GSVA_df)

# scaling.
GSVA_df_scaled <- scale(GSVA_df)

###

# umap and kmeans.
set.seed(123)
umap_res2 <- umap(GSVA_df_scaled)
km_res2 <- kmeans(umap_res2, centers = 4, nstart = 50)


###

# Labeling cluster by dominant cell of origin.
cluster_means2 <- aggregate(GSVA_df, by = list(cluster = km_res2$cluster), FUN = mean)


clusters2 <- km_res2$cluster


colors <- rainbow(length(unique(clusters2)))


plot(umap_res2, col = colors[clusters2], pch = 19,
     xlab = "UMAP 1", ylab = "UMAP 2",
     main = "K-means clusters on GSVA scores")
legend("topright", 
       legend = paste("Cluster", sort(unique(clusters2))), 
       col = colors, pch = 19)

