# run DGE-Limma first ::)
source("DataCleaning.R")


## in fpkm data, rank the genes per sample. Cluster them based on NO_TMM differential gene expression and see if the samples cluster out??

library(tidyverse)
library(uwot)
library(cluster)
library(cmapR)

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

source("DGE-Limma.R") # for candidate_genes.

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

# From the elbow plot, graph seems to plateau at PC:1-3.
pc_scores <- pca_result$x[, 1:3]

##### Clustering the whole expression data.

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

# From the elbow plot, graph seems to plateau at PC:1-5.
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


######################################################################################


##########################################################################################

## clustering based on the best signature from AUC-GSVA.

source("DataCleaning.R")
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

gene_list_together <- c("CPNE8", "PGM2L1", "LIFR", "CNR1", "HECW2", "CPNE3", "HOXC9", "SNX16", "IGSF10",
                                           "PRR7", "IGLV6-57", "SAC3D1", "CCDC86", "DDN") # from aucGSVA: genes most upregulated and downregulated in NO_TMM.


## only looking at the DGE genes -- common in NO_TMM vs. ALT and NO_TMM vs. Telomerase.
ranked_TARGET_NBL <- ranked_TARGET_NBL[gene_list_together, ]

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

# From the elbow plot, graph seems to plateau at PC:1-4.
pc_scores <- pca_result$x[, 1:4]

##### Clustering the whole expression data.

# k-means clustering.
set.seed(123)
umap_res <- umap(pc_scores)

km_res_sig <- kmeans(umap_res, centers = 2)

# Plotting with clusters.
tmm_colors <- c("Telomerase" = "red", "ALT" = "blue", "NO_TMM" = "green")

plot(umap_res, col = tmm_colors[metadata$TMM], pch = 19,
     xlab = "UMAP 1", ylab = "UMAP 2", main = "k-means Clusters on UMAP")


noTMM_idx <- which(metadata$TMM == "NO_TMM" & metadata$COG.Risk.Group == "High Risk")

text(umap_res[noTMM_idx, 1], umap_res[noTMM_idx, 2], 
     labels = metadata$SampleID[noTMM_idx], 
     pos = 3, cex = 0.7, col = "black")

#######################################################################################

# clustering on the expression log counts data.
source("DataCleaning.R")


gene_list <- c("CPNE8", "PGM2L1", "CNR1", "LIFR", "HOXC9", "SNX16", "HECW2", "ALDH3A2",
               "THSD7A", "CPNE3", "IGSF10")
gene_list_ <- c("PRR7", "IGLV6-57", "SAC3D1", "CCDC86", "DDN")
gene_list_together <- c("CPNE8", "PGM2L1", "LIFR", "CNR1", "HECW2", "CPNE3", "HOXC9", "SNX16", "IGSF10",
                        "PRR7", "IGLV6-57", "SAC3D1", "CCDC86", "DDN")

## only looking at the DGE genes -- common in NO_TMM vs. ALT and NO_TMM vs. Telomerase.
TARGET_NBL <- Expression[gene_list_together, ]


##########################################################################################

#### Now making the clusters only based on these candidate genes.

# Scale and transpose first.
pca_scale <- scale(t(TARGET_NBL))

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

# Plotting with clusters.
tmm_colors <- c("TMM" = "red", "NO_TMM" = "green")

plot(umap_res, col = tmm_colors[metadata$TMM_Case], pch = 19,
     xlab = "UMAP 1", ylab = "UMAP 2", main = "k-means Clusters on UMAP")


# clustering based on log expression  does not really work...

########################################################################################

# clustering on Ackerman samples using the signature GSVA list.

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
ranked_ackerman_NB <- ranked_ackerman_NB[gene_list_together, ]

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

# From the elbow plot, graph seems to plateau at PC:1-4.
pc_scores <- pca_result$x[, 1:4]

##### Clustering the whole expression data.

# k-means clustering.
set.seed(123)
umap_res <- umap(pc_scores)

km_res_sig2 <- kmeans(umap_res, centers = 2)

# Plotting with clusters.
tmm_colors <- c("Telomerase" = "red", "NO_TMM" = "green", "ALT" = "blue")

plot(umap_res, col = tmm_colors[ackerman_metadata$TMM_Category], pch = 19,
     xlab = "UMAP 1", ylab = "UMAP 2", main = "k-means Clusters on UMAP")

noTMM_idx <- which(ackerman_metadata$TMM_Case == "NO_TMM" & ackerman_metadata$Risk == "YES")

text(umap_res[noTMM_idx, 1], umap_res[noTMM_idx, 2], 
     labels = metadata$SampleID[noTMM_idx], 
     pos = 3, cex = 0.7, col = "black")



###################################################################################

# clustering samples based on MES & ADR markers -- divide samples hopefully into MES, ADR and intermediate samples.
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

# #mes_genes <- c("MEOX1", "WWTR1", "VIM", "CD44", "CBFB", "FOSL2", "MEOX2", "GLIS3", "TBX18", "NR3C1",
#                           "PRRX1", "MEF2D", "BHLHE41", "RUNX2", "IRF1", "NOTCH2", "YAP1", "CREG1", "DCAF6", "FLI1",
#                           "RUNX1", "IRF2", "JUN", "MAML2", "ZFP36L1")

mes_genes <- c("FN1", "VIM", "SNAI2", "PRRX1", "YAP1", "WWTR1", "DLK1", "DBH", "PHOX2A", "PHOX2B", "GATA2", "GATA3")

## only looking at the MES markers correlated with EXTEND scores.
ranked_TARGET_NBL <- ranked_TARGET_NBL[mes_genes, ]


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

##### Clustering the whole expression data.

# k-means clustering.
set.seed(123)
umap_res <- umap(pc_scores)

km_res <- kmeans(umap_res, centers = 3)

colors <- rainbow(length(unique(km_res$cluster)))


# Plot UMAP with colors
plot(umap_res, col = km_res$cluster, pch = 19,
     xlab = "UMAP 1", ylab = "UMAP 2",
     main = "k-means Clusters on UMAP")
legend("topright", 
       legend = paste("Cluster", sort(unique(km_res$cluster))), 
       col = colors, pch = 19)












