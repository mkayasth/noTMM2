ssource("NO_TMM/DataCleaning.R")

### Adrenergic markers: "PHOX2A","PHOX2B", "HAND1", "HAND2", "GATA2", "GATA3", "ISL1", "TBX2", "ASCL1", "TH", "DBH", "DLL3", "ATOH8"
### Mesenchymal markers: "CD44","VIM", "FN1", "TWIST2", "SNAI2", "PRRX1", "NOTCH3"

library(tidyverse)
library(uwot)
library(cluster)
library(cmapR)
library(GSVA)

ExpressionFpkm <- readRDS("TARGET_GEdata_062024.RDS")


# Data filtering: only including sample IDs from metadata present in the fpkm dataset & vice-versa.
metadata <- metadata %>%
  filter(SampleID %in% colnames(ExpressionFpkm))
ExpressionFpkm <- ExpressionFpkm[, colnames(ExpressionFpkm) %in% metadata$SampleID, drop = FALSE]

metadata <- metadata %>%
  arrange(TMM)

ExpressionFpkm <- ExpressionFpkm[, match(metadata$SampleID, colnames(ExpressionFpkm))]

# ranking the genes in ExpressionFpkm.
ranked_TARGET_NBL <- apply(ExpressionFpkm, 2, function(x) rank(x, ties.method = "average"))

## only looking at the markers genes for adrenergic phenotypes first.
ranked_TARGET_NBL_adrenergic <- ranked_TARGET_NBL[c("PHOX2A","PHOX2B", "HAND1", "HAND2", "GATA2", "GATA3", 
                                                    "ISL1", "TBX2", "ASCL1", "TH", "DBH", "DLL3", "ATOH8"), ]

##########################################################################################

#### Now making the clusters only based on these markers genes.

# Scale and transpose first.
pca_scale <- scale(t(ranked_TARGET_NBL_adrenergic))

# Run PCA
pca_result <- prcomp(pca_scale, center = TRUE, scale. = TRUE)

# plotting the visualize the variance.
pca_var <- pca_result$sdev^2
pca_var_explained <- pca_var / sum(pca_var)


# Make full scree plot to identify plateau.
plot(pca_var_explained * 100, type = "b", pch = 19,
     xlab = "Principal Component",
     ylab = "Variance Explained (%)",
     main = "Scree Plot: Full PCA")

# From the elbow plot, graph seems to plateau at PC:1-4.
pc_scores <- pca_result$x[, 1:4]

##### Clustering the whole expression data.

# k-means clustering.
set.seed(123)
umap_res <- umap(pc_scores)

km_res <- kmeans(umap_res, centers = 2)


plot(umap_res, col = km_res$cluster, pch = 19,
     xlab = "UMAP 1", ylab = "UMAP 2", main = "k-means Clusters on UMAP")



clusters_tmm <- as.data.frame(km_res$cluster)
clusters_tmm$SampleID <- rownames(clusters_tmm)
rownames(clusters_tmm) <- NULL

# coloring by TMM to see if there is any association between TMM and adrenergic/mesenchymal.
clusters_tmm <- left_join(clusters_tmm, metadata[, c("SampleID", "TMM", "TMM_Case", "MYCN.status", "COG.Risk.Group")], by = "SampleID")


tmm_colors <- c("TMM" = "red", "NO_TMM" = "green")
plot(umap_res, col = tmm_colors[metadata$TMM_Case], pch = 19,
     xlab = "UMAP 1", ylab = "UMAP 2", main = "k-means Clusters on UMAP")

######################################################################################

### Now, making the cluster for mesenchymal markers.

ranked_TARGET_NBL_mesenchymal <- ranked_TARGET_NBL[c("CD44","VIM", "FN1", "TWIST2", "SNAI2", "PRRX1", 
                                                    "NOTCH3"), ]

# Scale and transpose first.
pca_scale <- scale(t(ranked_TARGET_NBL_mesenchymal))

# Run PCA
pca_result <- prcomp(pca_scale, center = TRUE, scale. = TRUE)

# plotting the visualize the variance.
pca_var <- pca_result$sdev^2
pca_var_explained <- pca_var / sum(pca_var)


# Make full scree plot to identify plateau.
plot(pca_var_explained * 100, type = "b", pch = 19,
     xlab = "Principal Component",
     ylab = "Variance Explained (%)",
     main = "Scree Plot: Full PCA")

# From the elbow plot, graph seems to plateau at PC:1-4.
pc_scores <- pca_result$x[, 1:4]

##### Clustering the whole expression data.

# k-means clustering.
set.seed(123)
umap_res <- umap(pc_scores)

km_res2 <- kmeans(umap_res, centers = 2)


plot(umap_res, col = km_res2$cluster, pch = 19,
     xlab = "UMAP 1", ylab = "UMAP 2", main = "k-means Clusters on UMAP")



clusters_tmm2 <- as.data.frame(km_res2$cluster)
clusters_tmm2$SampleID <- rownames(clusters_tmm2)
rownames(clusters_tmm2) <- NULL

# coloring by TMM to see if there is any association between TMM and adrenergic/mesenchymal.
clusters_tmm2 <- left_join(clusters_tmm2, metadata[, c("SampleID", "TMM", "TMM_Case", "MYCN.status", "COG.Risk.Group")], by = "SampleID")


tmm_colors <- c("TMM" = "red", "NO_TMM" = "green")
plot(umap_res, col = tmm_colors[metadata$TMM_Case], pch = 19,
     xlab = "UMAP 1", ylab = "UMAP 2", main = "k-means Clusters on UMAP")


clusters_tmm2 <- left_join(clusters_tmm2,  clusters_tmm[, c("SampleID", "km_res$cluster")], by = "SampleID")

clusters_tmm2 <- clusters_tmm2 %>%
  arrange(`km_res2$cluster`)



#########################################################################################

# combining both markers together.
ranked_TARGET_NBL_states <- ranked_TARGET_NBL[c("CD44","VIM", "FN1", "TWIST2", "SNAI2", "PRRX1", 
                                                     "NOTCH3", "PHOX2A","PHOX2B", "HAND1", "HAND2", "GATA2", "GATA3", 
                                                "ISL1", "TBX2", "ASCL1", "TH", "DBH", "DLL3", "ATOH8"), ]

#### Now making the clusters only based on these markers genes.

# Scale and transpose first.
pca_scale <- scale(t(ranked_TARGET_NBL_states))

# Run PCA
pca_result <- prcomp(pca_scale, center = TRUE, scale. = TRUE)

# plotting the visualize the variance.
pca_var <- pca_result$sdev^2
pca_var_explained <- pca_var / sum(pca_var)


# Make full scree plot to identify plateau.
plot(pca_var_explained * 100, type = "b", pch = 19,
     xlab = "Principal Component",
     ylab = "Variance Explained (%)",
     main = "Scree Plot: Full PCA")

# From the elbow plot, graph seems to plateau at PC:1-9.
pc_scores <- pca_result$x[, 1:4]

##### Clustering the whole expression data.

# k-means clustering.
set.seed(123)
umap_res <- umap(pc_scores)

km_res3 <- kmeans(umap_res, centers = 2)


plot(umap_res, col = km_res3$cluster, pch = 19,
     xlab = "UMAP 1", ylab = "UMAP 2", main = "k-means Clusters on UMAP")



clusters_tmm3 <- as.data.frame(km_res3$cluster)
clusters_tmm3$SampleID <- rownames(clusters_tmm3)
rownames(clusters_tmm) <- NULL

# coloring by TMM to see if there is any association between TMM and adrenergic/mesenchymal.
clusters_tmm3 <- left_join(clusters_tmm3, metadata[, c("SampleID", "TMM", "TMM_Case", "MYCN.status", "COG.Risk.Group")], by = "SampleID")
clusters_tmm3 <- left_join(clusters_tmm3, clusters_tmm2[, c("SampleID", "km_res$cluster", "km_res2$cluster")], by = "SampleID")




tmm_colors <- c("TMM" = "red", "NO_TMM" = "green")
plot(umap_res, col = tmm_colors[metadata$TMM_Case], pch = 19,
     xlab = "UMAP 1", ylab = "UMAP 2", main = "k-means Clusters on UMAP")


########################################################################################


### using GSVA on the log counts using these markers.
gene_adr <- list("ADR Genes" = c("PHOX2A","PHOX2B", "HAND1", "HAND2", "GATA2", "GATA3", 
                                 "ISL1", "TBX2", "ASCL1", "TH", "DBH", "DLL3", "ATOH8"))
gene_mes <- list("MES Genes" = c("CD44","VIM", "FN1", "TWIST2", "SNAI2", "PRRX1", "NOTCH3"))

gene_combined <- list(
  "ADR Genes" = c("PHOX2A","PHOX2B", "HAND1", "HAND2", "GATA2", "GATA3", 
                  "ISL1", "TBX2", "ASCL1", "TH", "DBH", "DLL3", "ATOH8"),
  "MES Genes" = c("CD44","VIM", "FN1", "TWIST2", "SNAI2", "PRRX1", "NOTCH3"))

# GSVA for adrenergic marker.
Expression <- as.matrix(Expression)
adr_gsva <- gsvaParam(Expression, gene_adr, kcdf = "Gaussian")
GSVA_result_adr <- gsva(adr_gsva)
GSVA_df_adr <- as.data.frame(GSVA_result_adr)
GSVA_long_adr <- pivot_longer(GSVA_df_adr, cols = everything(), names_to = "SampleID", values_to = "GSVA_Score")



# GSVA for mesenchymal marker.
mes_gsva <- gsvaParam(Expression, gene_mes, kcdf = "Gaussian")
GSVA_result_mes <- gsva(mes_gsva)
GSVA_df_mes <- as.data.frame(GSVA_result_mes)
GSVA_long_mes <- pivot_longer(GSVA_df_mes, cols = everything(), names_to = "SampleID", values_to = "GSVA_Score")


# putting both together.
GSVA_long_state <- left_join(GSVA_long_adr, GSVA_long_mes, by = "SampleID")
colnames(GSVA_long_state) <- c("SampleID", "GSVAadr", "GSVAmes")


GSVA_long_state <- pivot_longer(
  GSVA_long_state,
  cols = c(GSVAadr, GSVAmes),
  names_to = "ScoreType",
  values_to = "Score"
)


# making bar graphs with GSVA scores for both markers.

# Open PDF device
pdf("GSVA_states.pdf", width = 3, height = 3)

# Looping through each sample and create a plot per page
for (sample in unique(GSVA_long_state$SampleID)) {
  p <- ggplot(subset(GSVA_long_state, SampleID == sample),
              aes(x = ScoreType, y = Score, fill = ScoreType)) +
    geom_bar(stat = "identity", width = 0.5) +
    labs(title = sample, x = NULL, y = "GSVA Score") +
    theme_classic() +
    theme(legend.position = "none")
  
  print(p)
}


# Close PDF device
dev.off()


