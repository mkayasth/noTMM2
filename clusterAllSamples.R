# run DGE-Limma first ::)
source("NO_TMM/DataCleaning.R")


## in fpkm data, rank the genes per sample. Cluster them based on NO_TMM differential gene expression and see if the samples cluster out??

library(tidyverse)
library(uwot)
library(cluster)
library(cmapR)

ExpressionFpkm <- readRDS("TARGET_GEdata_062024.RDS")
metadata <- read.table(file = 'Metadata_TARGETFinal_08012024.txt', header = TRUE, sep = '\t')

# Data filtering: only including sample IDs from metadata present in the fpkm dataset & vice-versa.
metadata <- metadata %>%
  filter(SampleID %in% colnames(ExpressionFpkm))
ExpressionFpkm <- ExpressionFpkm[, colnames(ExpressionFpkm) %in% metadata$SampleID, drop = FALSE]

metadata <- metadata %>%
  arrange(TMM)

ExpressionFpkm <- ExpressionFpkm[, match(metadata$SampleID, colnames(ExpressionFpkm))]

# ranking the genes in ExpressionFpkm.
ranked_TARGET_NBL <- apply(ExpressionFpkm, 2, function(x) rank(x, ties.method = "average"))

## only looking at the DGE genes -- common in NO_TMM vs. ALT and NO_TMM vs. Telomerase.
ranked_TARGET_NBL <- ranked_TARGET_NBL[candidate_genes, ]


##########################################################################################

#### Now making the clusters only based on these candidate genes.

# Scale and transpose first.
pca_scale <- scale(t(ranked_TARGET_NBL))

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

# From the elbow plot, graph seems to plateau at PC:1-10.
pc_scores <- pca_result$x[, 1:10]

##### Clustering the whole expression data with batch effects.

# k-means clustering.
set.seed(123)
umap_res <- umap(pc_scores)

km_res <- kmeans(umap_res, centers = 2)

# Plotting with clusters.
tmm_colors <- c("TMM" = "red", "NO_TMM" = "green")

plot(umap_res, col = tmm_colors[metadata$TMM_Case], pch = 19,
     xlab = "UMAP 1", ylab = "UMAP 2", main = "k-means Clusters on UMAP")



clusters_tmm <- as.data.frame(km_res$cluster)
clusters_tmm$SampleID <- rownames(clusters_tmm)
rownames(clusters_tmm) <- NULL

clusters_tmm <- left_join(clusters_tmm, metadata[, c("SampleID", "TMM", "TMM_Case", "MYCN.status", "COG.Risk.Group")], by = "SampleID")


# Labeling NO_TMM high-risk samples.

plot(umap_res, col = tmm_colors[metadata$TMM_Case], pch = 19,
     xlab = "UMAP 1", ylab = "UMAP 2", main = "k-means Clusters on UMAP")


noTMM_idx <- which(metadata$TMM == "NO_TMM" & metadata$COG.Risk.Group == "High Risk")

text(umap_res[noTMM_idx, 1], umap_res[noTMM_idx, 2], 
     labels = metadata$SampleID[noTMM_idx], 
     pos = 3, cex = 0.7, col = "black")

######################################################################################


#### trying clustering with Ackerman data.
gct_file<- parse_gctx("Neuroblastoma_208Samples.gct")
ackerman_NB <- gct_file@mat

# Loading metadata.
ackerman_metadata <- read.table(file = 'NBL_Ackerman_CompleteMeta.txt', header = TRUE, sep = '\t')

# only including SampleID in microarray data present in metadata.
ackerman_NB <- ackerman_NB[, colnames(ackerman_NB) %in% ackerman_metadata$SampleID]
ackerman_metadata <- ackerman_metadata[ackerman_metadata$SampleID %in% colnames(ackerman_NB), ]

ackerman_metadata <- ackerman_metadata %>%
  arrange(TMM_Case)

ackerman_NB <- ackerman_NB[, match(ackerman_metadata$SampleID, colnames(ackerman_NB))]

# ranking the genes.
ranked_ackerman_NB <- apply(ackerman_NB, 2, function(x) rank(x, ties.method = "average"))

## only looking at the DGE genes -- common in NO_TMM vs. ALT and NO_TMM vs. Telomerase.
candidate_genes2 <- candidate_genes[candidate_genes %in% rownames(ranked_ackerman_NB)]
ranked_ackerman_NB <- ranked_ackerman_NB[candidate_genes2, ]

#### Now making the clusters only based on these candidate genes.

# Scale and transpose first.
pca_scale <- scale(t(ranked_ackerman_NB))

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

# From the elbow plot, graph seems to plateau at PC:1-10.
pc_scores <- pca_result$x[, 1:5]

##### Clustering the whole expression data.

# k-means clustering.
set.seed(123)
umap_res <- umap(pc_scores)

km_res <- kmeans(umap_res, centers = 2)

# Plotting with clusters.
tmm_colors <- c("TMM" = "red", "NO_TMM" = "green")

plot(umap_res, col = tmm_colors[ackerman_metadata$TMM_Case], pch = 19,
     xlab = "UMAP 1", ylab = "UMAP 2", main = "k-means Clusters on UMAP")



clusters_tmm <- as.matrix(km_res$cluster)
clusters_tmm <- as.data.frame(clusters_tmm)

clusters_tmm$SampleID <- rownames(clusters_tmm)
rownames(clusters_tmm) <- NULL

clusters_tmm <- left_join(clusters_tmm, ackerman_metadata[, c("SampleID", "TMM_Category", "TMM_Case", "MYCNStatus", "Risk")], by = "SampleID")


# Labeling NO_TMM high-risk samples.

plot(umap_res, col = tmm_colors[ackerman_metadata$TMM_Case], pch = 19,
     xlab = "UMAP 1", ylab = "UMAP 2", main = "k-means Clusters on UMAP")


noTMM_idx <- which(ackerman_metadata$TMM_Case == "NO_TMM" & ackerman_metadata$Risk == "YES")

text(umap_res[noTMM_idx, 1], umap_res[noTMM_idx, 2], 
     labels = metadata$SampleID[noTMM_idx], 
     pos = 3, cex = 0.7, col = "black")



#####################################################################################
# using candidate_genes_minus to hopefully separate TMM and NO_TMM.

ExpressionFpkm <- readRDS("TARGET_GEdata_062024.RDS")
metadata <- read.table(file = 'Metadata_TARGETFinal_08012024.txt', header = TRUE, sep = '\t')

# Data filtering: only including sample IDs from metadata present in the fpkm dataset & vice-versa.
metadata <- metadata %>%
  filter(SampleID %in% colnames(ExpressionFpkm))
ExpressionFpkm <- ExpressionFpkm[, colnames(ExpressionFpkm) %in% metadata$SampleID, drop = FALSE]

metadata <- metadata %>%
  arrange(TMM)

ExpressionFpkm <- ExpressionFpkm[, match(metadata$SampleID, colnames(ExpressionFpkm))]

# ranking the genes in ExpressionFpkm.
ranked_TARGET_NBL <- apply(ExpressionFpkm, 2, function(x) rank(x, ties.method = "average"))

## only looking at the DGE genes -- common in NO_TMM vs. ALT and NO_TMM vs. Telomerase.
ranked_TARGET_NBL <- ranked_TARGET_NBL[candidate_genes_minus, ]


##########################################################################################

#### Now making the clusters only based on these candidate genes.

# Scale and transpose first.
pca_scale <- scale(t(ranked_TARGET_NBL))

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

# From the elbow plot, graph seems to plateau at PC:1-5.
pc_scores <- pca_result$x[, 1:5]

##### Clustering the whole expression data with batch effects.

# k-means clustering.
set.seed(123)
umap_res <- umap(pc_scores)

km_res <- kmeans(umap_res, centers = 2)

# Plotting with clusters.
tmm_colors <- c("TMM" = "red", "NO_TMM" = "green")

plot(umap_res, col = tmm_colors[metadata$TMM_Case], pch = 19,
     xlab = "UMAP 1", ylab = "UMAP 2", main = "k-means Clusters on UMAP")



clusters_tmm <- as.data.frame(km_res$cluster)
clusters_tmm$SampleID <- rownames(clusters_tmm)
rownames(clusters_tmm) <- NULL

clusters_tmm <- left_join(clusters_tmm, metadata[, c("SampleID", "TMM", "TMM_Case", "MYCN.status", "COG.Risk.Group")], by = "SampleID")


# Labeling NO_TMM high-risk samples.

plot(umap_res, col = tmm_colors[metadata$TMM_Case], pch = 19,
     xlab = "UMAP 1", ylab = "UMAP 2", main = "k-means Clusters on UMAP")


noTMM_idx <- which(metadata$TMM == "NO_TMM" & metadata$COG.Risk.Group == "High Risk")

text(umap_res[noTMM_idx, 1], umap_res[noTMM_idx, 2], 
     labels = metadata$SampleID[noTMM_idx], 
     pos = 3, cex = 0.7, col = "black")



######################################################################################

## using DGE between High Risk and Low Risk, separating NO_TMM.

ExpressionFpkm <- readRDS("TARGET_GEdata_062024.RDS")
metadata <- read.table(file = 'Metadata_TARGETFinal_08012024.txt', header = TRUE, sep = '\t')

# Data filtering: only including sample IDs from metadata present in the fpkm dataset & vice-versa.
metadata <- metadata %>%
  filter(SampleID %in% colnames(ExpressionFpkm))

metadata <- metadata %>%
  filter(TMM == "NO_TMM")

metadata <- metadata %>%
  arrange(COG.Risk.Group)

ExpressionFpkm <- ExpressionFpkm[, colnames(ExpressionFpkm) %in% metadata$SampleID, drop = FALSE]

ExpressionFpkm <- ExpressionFpkm[, match(metadata$SampleID, colnames(ExpressionFpkm))]

# ranking the genes in ExpressionFpkm.
ranked_TARGET_NBL <- apply(ExpressionFpkm, 2, function(x) rank(x, ties.method = "average"))

## only looking at a$Genes.
ranked_TARGET_NBL <- ranked_TARGET_NBL[regression_results_sig_highRisk$Gene, ]



#### Now making the clusters only based on these candidate genes.

# Scale and transpose first.
pca_scale <- scale(t(ranked_TARGET_NBL))

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

# From the elbow plot, graph seems to plateau at PC:1-10.
pc_scores <- pca_result$x[, 1:5]

##### Clustering the whole expression data with batch effects.

# k-means clustering.
set.seed(123)
umap_res <- umap(pc_scores)

km_res <- kmeans(umap_res, centers = 2)

# Plotting with clusters.
tmm_colors <- c("High Risk" = "red", "Low Risk" = "green", "Intermediate Risk" = "green")

plot(umap_res, col = tmm_colors[metadata$COG.Risk.Group], pch = 19,
     xlab = "UMAP 1", ylab = "UMAP 2", main = "k-means Clusters on UMAP")



clusters_tmm <- as.data.frame(km_res$cluster)
clusters_tmm$SampleID <- rownames(clusters_tmm)
rownames(clusters_tmm) <- NULL

clusters_tmm <- left_join(clusters_tmm, metadata[, c("SampleID", "TMM", "TMM_Case", "MYCN.status", "COG.Risk.Group")], by = "SampleID")


# Labeling NO_TMM high-risk samples.

plot(umap_res, col = tmm_colors[metadata$COG.Risk.Group], pch = 19,
     xlab = "UMAP 1", ylab = "UMAP 2", main = "k-means Clusters on UMAP")


noTMM_idx <-  #samples not following the trend.
idx <- match(noTMM_idx, metadata$SampleID)

text(umap_res[idx, 1], umap_res[idx, 2], 
     labels = metadata$SampleID[idx], 
     pos = 3, cex = 0.7, col = "black")

### not the most distinct cluster.
