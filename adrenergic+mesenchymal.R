source("DataCleaning.R")

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
ranked_TARGET_NBL_adrenergic <- ranked_TARGET_NBL[c("PHOX2B", "HAND2", "GATA3", 
                                                    "ISL1", "ASCL1", "TH", "DBH"), ]

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

ranked_TARGET_NBL_mesenchymal <- ranked_TARGET_NBL[c("VIM", "FN1", "TWIST1", "SNAI2", "PRRX1", 
                                                    "ALDH1A3"), ]

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

# From the elbow plot, graph seems to plateau at PC:1-3.
pc_scores <- pca_result$x[, 1:3]

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

# ranked_TARGET_NBL_states <- ranked_TARGET_NBL[c("CD44","VIM", "FN1", "TWIST2", "SNAI2", "PRRX1", 
#                                                      "NOTCH3", "PHOX2A","PHOX2B", "HAND1", "HAND2", "GATA2", "GATA3", 
#                                                 "ISL1", "TBX2", "ASCL1", "TH", "DBH", "DLL3", "ATOH8"), ]

ranked_TARGET_NBL_states <- ranked_TARGET_NBL[c("PHOX2B", "HAND2", "GATA3", 
                                                "ISL1", "ASCL1", "TH", "DBH", "VIM", "FN1", "TWIST1", "SNAI2", "PRRX1", 
                                                "ALDH1A3"), ]

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

# From the elbow plot, graph seems to plateau at PC:1-5.
pc_scores <- pca_result$x[, 1:5]

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
gene_adr <- list("ADR Genes" = c("PHOX2B", "HAND2", "GATA3", 
                                 "ISL1", "ASCL1", "TH", "DBH"))
gene_mes <- list("MES Genes" = c("VIM", "FN1", "TWIST1", "SNAI2", "PRRX1", 
                                 "ALDH1A3"))

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

# fixing the bad column name in clusters_tmm3.
names(clusters_tmm3)[names(clusters_tmm3) == "km_res3$cluster"] <- "km_res3_cluster"

# wide GSVA
gsva_wide <- GSVA_long_state %>%
  pivot_wider(names_from = ScoreType, values_from = Score)

# computing diffs
gsva_wide$mes_adr <- gsva_wide$GSVAmes - gsva_wide$GSVAadr
gsva_wide$adr_mes <- gsva_wide$GSVAadr - gsva_wide$GSVAmes

# join and filter
tmp <- merge(clusters_tmm3, gsva_wide, by = "SampleID", all.x = TRUE)

results <- tmp[
  (tmp$mes_adr >= 0.3 & tmp$km_res3_cluster == 1) |
    (tmp$adr_mes >= 0.3 & tmp$km_res3_cluster == 2),
]

## plotting the clusters back to see if better GSVA results mean anything to the cluster.
gsva_colors <- case_when(
  clusters_tmm3$SampleID %in% results$SampleID[results$mes_adr > 0] ~ "darkblue",
  clusters_tmm3$SampleID %in% results$SampleID[results$adr_mes > 0] ~ "lightblue",
  TRUE ~ "red"
)

plot(umap_res, col = gsva_colors, pch = 19,
     xlab = "UMAP 1", ylab = "UMAP 2", main = "k-means Clusters on UMAP")

## looking at only NO_TMM/risk group samples.
umap_res_notmm <- umap_res[rownames(umap_res)  %in% metadata$SampleID[metadata$TMM_Case == "NO_TMM"], ]
umap_res_gsva <- umap_res[rownames(umap_res)  %in% results$SampleID, ]


notmm_colors <- c("High Risk" = "red", "Intermediate Risk" = "green", "Low Risk" = "green")

plot(umap_res, col =  notmm_colors[metadata$COG.Risk.Group], pch = 19,
     xlab = "UMAP 1", ylab = "UMAP 2", main = "k-means Clusters on UMAP")

plot(umap_res_notmm, col =  notmm_colors[metadata$COG.Risk.Group], pch = 19,
     xlab = "UMAP 1", ylab = "UMAP 2", main = "k-means Clusters on UMAP")


plot(umap_res_gsva, col =  notmm_colors[metadata$COG.Risk.Group], pch = 19,
     xlab = "UMAP 1", ylab = "UMAP 2", main = "k-means Clusters on UMAP")



#######################################################################################

###################################################################################################

# heatmap based on adrenergic and mesenchymal elements from single cell.

library(ComplexHeatmap)
library(grid)

source("DataCleaning.R")

metadata <- metadata %>%
  arrange(TMM_Case, TMM)

Expression <- Expression[, match(metadata$SampleID, colnames(Expression))]

# common signatures from single & bulk cell studies.
heatmap_genes <- c("PHOX2A", "PHOX2B", "HAND2", "ISL1", "GATA3", "ASCL1", "TBX2", "DBH", "TH",
                  "PRRX1", "RUNX1", "FOSL1", "JUN", "YAP1", "WWTR1", "MEOX1", "MEOX2", "SNAI2", "CD44", "FN1", "VIM")

  
  
  
expression_matrix <- t(scale(t(Expression[rownames(Expression) %in% heatmap_genes, ,drop = FALSE])))

# adding gsva as annotation.

# building split vector aligned to rownames of expression matrix.
grp_map <- setNames(c(rep("Adrenergic", 9), rep("Mesenchymal", 21)), heatmap_genes)
row_groups <- factor(grp_map[rownames(expression_matrix)], levels = c("Adrenergic","Mesenchymal"))

# Running GSVA.

gene_set_list_adrenergic <-list(ADRN = c("PHOX2A", "PHOX2B", "HAND2", "ISL1", "GATA3", "ASCL1", "TBX2", "DBH", "TH"))
gene_set_list_mesenchymal <- list(MES = c("PRRX1", "RUNX1", "FOSL1", "JUN", "YAP1", "WWTR1", "MEOX1", "MEOX2", "SNAI2", "CD44", "FN1", "VIM"))


gsvapar <- gsvaParam(as.matrix(Expression), gene_set_list_adrenergic, kcdf = "Gaussian")
gsva_result <- gsva(gsvapar)

gsvapar2 <- gsvaParam(as.matrix(Expression), gene_set_list_mesenchymal, kcdf = "Gaussian")
gsva_result2 <- gsva(gsvapar2)


adr <- as.numeric(gsva_result[1, ])
names(adr) <- colnames(gsva_result)

mes <- as.numeric(gsva_result2[1, ])
names(mes) <- colnames(gsva_result2)


stopifnot(all(colnames(Expression) %in% names(adr)))
stopifnot(all(colnames(Expression) %in% names(mes)))

# computing clustering order first(no annotation yet).
ht_pre <- Heatmap(
  expression_matrix,
  row_split         = row_groups,
  cluster_columns   = TRUE,
  cluster_rows      = FALSE,
  show_column_names = TRUE,
  row_names_gp      = gpar(fontsize = 8, fontface = "bold"),
  column_names_gp   = gpar(fontsize = 8),
  column_names_rot  = 45
)
ht_pre <- draw(ht_pre)

# column order for annotation.
co <- column_order(ht_pre)

ord_idx <- if (is.list(co)) unlist(co, recursive = TRUE, use.names = FALSE) else as.integer(co)

# sanity check
stopifnot(length(ord_idx) == ncol(expression_matrix))
stopifnot(identical(sort(ord_idx), seq_len(ncol(expression_matrix))))

ord_samples <- colnames(expression_matrix)[ord_idx]

# gsva annotation in the same order now.
adr <- setNames(as.numeric(gsva_result[1, ]),  colnames(gsva_result))
mes <- setNames(as.numeric(gsva_result2[1, ]), colnames(gsva_result2))

# align to the display order
adr_ord <- adr[ord_samples]
mes_ord <- mes[ord_samples]


adr_col <- ifelse(adr_ord > 0, "blue", "red")
mes_col <- ifelse(mes_ord > 0, "blue", "red")

top_anno <- HeatmapAnnotation(
  ADR = anno_barplot(
    adr_ord,
    gp = gpar(fill = adr_col),
    axis_param = list(at = c(-1, 0, 1), labels = c("-1", "0", "1"))
  ),
  MES = anno_barplot(
    mes_ord,
    gp = gpar(fill = mes_col),
    axis_param = list(at = c(-1, 0, 1), labels = c("-1", "0", "1"))
  ),
  which = "column"
)

# final heatmap.
ht_final <- Heatmap(
  expression_matrix[, ord_samples, drop = FALSE],
  row_split         = row_groups,
  cluster_columns   = TRUE,             
  cluster_rows      = FALSE,
  show_column_names = TRUE,
  top_annotation    = top_anno,
  row_names_gp      = gpar(fontsize = 8, fontface = "bold"),
  column_names_gp   = gpar(fontsize = 5),
  column_names_rot  = 45,
  cell_fun = function(j, i, x, y, w, h, fill) {
    grid.rect(x, y, w, h, gp = gpar(col = "black", fill = NA, lwd = 0.25))
  }
)


pdf("adr-mesHeatmap2.pdf", width = 15, height = 7)
draw(ht_final)
dev.off()


###

# heatmap from Thirant single cell markers.

heatmap_genes <- c("KLF7", "GATA3", "HAND2", "PHOX2A", "ISL1", "HAND1", "PHOX2B", "TFAP2B", "GATA2", "SATB1", "SIX3", "EYA1", "SOX11", "DACH1", "ASCL1", "HEY1", "KLF13", "PBX3",
                  "VIM", "FN1", "MEOX2", "ID1", "EGR3", "AEBP1", "CBFB", "IRF3", "IRF2", "IRF1", "TBX18", "MAFF", "RUNX2", "ZFP36L1", "NR3C1", "BHLHE41", "GLIS3", "RUNX1", "FOSL1", "FOSL2", "ELK4", "IFI16", "SIX4", "FLI1", "MAML2", "SMAD3", "DCAF6", "WWTR1", "SOX9", "MEF2D", "ZNF217", "PRRX1", "CREG1", "NOTCH2", "SIX1", "MEOX1")



expression_matrix <- t(scale(t(Expression[rownames(Expression) %in% heatmap_genes, ,drop = FALSE])))

# adding gsva as annotation.

grp_map <- setNames(c(rep("Adrenergic", 18), rep("Mesenchymal", 54)), heatmap_genes)
row_groups <- factor(grp_map[rownames(expression_matrix)], levels = c("Adrenergic","Mesenchymal"))

# Running GSVA.
gene_set_list_adrenergic <-list(ADRN = c("KLF7", "GATA3", "HAND2", "PHOX2A", "ISL1", "HAND1", "PHOX2B", "TFAP2B", "GATA2", "SATB1", "SIX3", "EYA1", "SOX11", "DACH1", "ASCL1", "HEY1", "KLF13", "PBX3"))
gene_set_list_mesenchymal <- list(MES = c("VIM", "FN1", "MEOX2", "ID1", "EGR3", "AEBP1", "CBFB", "IRF3", "IRF2", "IRF1", "TBX18", "MAFF", "RUNX2", "ZFP36L1", "NR3C1", "BHLHE41", "GLIS3", "RUNX1", "FOSL1", "FOSL2", "ELK4", "IFI16", "SIX4", "FLI1", "MAML2", "SMAD3", "DCAF6", "WWTR1", "SOX9", "MEF2D", "ZNF217", "PRRX1", "CREG1", "NOTCH2", "SIX1", "MEOX1"))


gsvapar <- gsvaParam(as.matrix(Expression), gene_set_list_adrenergic, kcdf = "Gaussian")
gsva_result <- gsva(gsvapar)

gsvapar2 <- gsvaParam(as.matrix(Expression), gene_set_list_mesenchymal, kcdf = "Gaussian")
gsva_result2 <- gsva(gsvapar2)


adr <- as.numeric(gsva_result[1, ])
names(adr) <- colnames(gsva_result)

mes <- as.numeric(gsva_result2[1, ])
names(mes) <- colnames(gsva_result2)

stopifnot(all(colnames(Expression) %in% names(adr)))
stopifnot(all(colnames(Expression) %in% names(mes)))

# computing clustering order first(no annotation yet).

ht_pre <- Heatmap(
  expression_matrix,
  row_split         = row_groups,
  cluster_columns   = TRUE,
  cluster_rows      = FALSE,
  show_column_names = TRUE,
  row_names_gp      = gpar(fontsize = 8, fontface = "bold"),
  column_names_gp   = gpar(fontsize = 8),
  column_names_rot  = 45
)
ht_pre <- draw(ht_pre)

# column order for annotation.
co <- column_order(ht_pre)
ord_idx <- if (is.list(co)) unlist(co, recursive = TRUE, use.names = FALSE) else as.integer(co)

# sanity check
stopifnot(length(ord_idx) == ncol(expression_matrix))
stopifnot(identical(sort(ord_idx), seq_len(ncol(expression_matrix))))

ord_samples <- colnames(expression_matrix)[ord_idx]


adr <- setNames(as.numeric(gsva_result[1, ]),  colnames(gsva_result))
mes <- setNames(as.numeric(gsva_result2[1, ]), colnames(gsva_result2))


adr_ord <- adr[ord_samples]
mes_ord <- mes[ord_samples]


adr_col <- ifelse(adr_ord > 0, "blue", "red")
mes_col <- ifelse(mes_ord > 0, "blue", "red")

top_anno <- HeatmapAnnotation(
  ADR = anno_barplot(
    adr_ord,
    gp = gpar(fill = adr_col),
    axis_param = list(at = c(-1, 0, 1), labels = c("-1", "0", "1"))
  ),
  MES = anno_barplot(
    mes_ord,
    gp = gpar(fill = mes_col),
    axis_param = list(at = c(-1, 0, 1), labels = c("-1", "0", "1"))
  ),
  which = "column"
)



# final heatmap.
ht_final <- Heatmap(
  expression_matrix[, ord_samples, drop = FALSE],
  row_split         = row_groups,
  cluster_columns   = TRUE,         
  cluster_rows      = FALSE,
  show_column_names = TRUE,
  top_annotation    = top_anno,
  row_names_gp      = gpar(fontsize = 8, fontface = "bold"),
  column_names_gp   = gpar(fontsize = 5),
  column_names_rot  = 45,
  cell_fun = function(j, i, x, y, w, h, fill) {
    grid.rect(x, y, w, h, gp = gpar(col = "black", fill = NA, lwd = 0.25))
  })


pdf("adr-mesHeatmap.pdf", width = 15, height = 10)
draw(ht_final)
dev.off()

# coloring NO_TMMs.
highlight_samples_NOTMM <- metadata$SampleID[metadata$TMM == "NO_TMM"]
highlight_samples_ALT <- metadata$SampleID[metadata$TMM == "ALT"]
highlight_samples_Telomerase <- metadata$SampleID[metadata$TMM == "Telomerase"]


mat <- expression_matrix[, ord_samples, drop = FALSE]

sample_colors <- setNames(rep("black", ncol(mat)), colnames(mat))
sample_colors[highlight_samples_NOTMM]      <- "black"
sample_colors[highlight_samples_ALT]        <- "blue"
sample_colors[highlight_samples_Telomerase] <- "red"

ht_final <- Heatmap(
  mat,
  name = "expr",
  row_split         = row_groups,
  cluster_columns   = TRUE,
  cluster_rows      = FALSE,
  show_column_names = TRUE,               
  column_names_gp   =  gpar(col = sample_colors[colnames(mat)], fontsize = 5),
  top_annotation    = top_anno,
  row_names_gp      = gpar(fontsize = 8, fontface = "bold"),
  cell_fun = function(j, i, x, y, w, h, fill) {
    grid.rect(x, y, w, h, gp = gpar(col = "black", fill = NA, lwd = 0.25))
  }
)

pdf("adr-mesHeatmap-TMM.pdf", width = 15, height = 10)
draw(ht_final)
dev.off()

# coloring based on risk groups.
highlight_samples_highRisk <- metadata$SampleID[metadata$COG.Risk.Group == "High Risk"]
highlight_samples_intermediate <- metadata$SampleID[metadata$COG.Risk.Group == "Intermediate Risk"]
highlight_samples_lowRisk <- metadata$SampleID[metadata$COG.Risk.Group == "Low Risk"]


mat <- expression_matrix[, ord_samples, drop = FALSE]

sample_colors <- setNames(rep("black", ncol(mat)), colnames(mat))
sample_colors[highlight_samples_highRisk]      <- "red"
sample_colors[highlight_samples_intermediate]        <- "blue"
sample_colors[highlight_samples_lowRisk] <- "black"

ht_final <- Heatmap(
  mat,
  name = "expr",
  row_split         = row_groups,
  cluster_columns   = TRUE,
  cluster_rows      = FALSE,
  show_column_names = TRUE,               
  column_names_gp   =  gpar(col = sample_colors[colnames(mat)], fontsize = 5),
  top_annotation    = top_anno,
  row_names_gp      = gpar(fontsize = 8, fontface = "bold"),
  cell_fun = function(j, i, x, y, w, h, fill) {
    grid.rect(x, y, w, h, gp = gpar(col = "black", fill = NA, lwd = 0.25))
  }
)

pdf("adr-mesHeatmap-RiskGroup.pdf", width = 15, height = 10)
draw(ht_final)
dev.off()

#################################################################################


# heatmap using only the genes with association with TERT/extend. // from playground.R.

mes_genes <- list(MES = c("NR3C1", "CBFB", "GLIS3", "IRF2" , "CD44", "DCAF6", "MAML2", "RUNX2", "WWTR1",  
              "YAP1", "IRF1", "VIM", "ZFP36L1", "BHLHE41",
              "GLIS3",  "CBFB"  , "FLI1", "MEF2D", "RUNX1", "CREG1", "JUN", 
              "NOTCH2", "TBX18", "FOSL2", "MEOX1", "IFI16", "PRRX1" , "MEOX2", "SNAI2",
            "EGR3", "MAFF", "ZNF217", "FN1", "FOSL1", "ELK4", "ID1", "AEBP1"))

expression_matrix <- t(scale(t(Expression[rownames(Expression) %in% unlist(mes_genes), ,drop = FALSE])))


gsvapar <- gsvaParam(as.matrix(Expression), mes_genes, kcdf = "Gaussian")
gsva_result <- gsva(gsvapar)

mes <- as.numeric(gsva_result[1, ])
names(mes) <- colnames(gsva_result)

# computing clustering order first(no annotation yet).

ht_pre <- Heatmap(
  expression_matrix,
  row_order = rownames(expression_matrix),
  cluster_columns   = TRUE,
  cluster_rows      = FALSE,
  show_column_names = TRUE,
  row_names_gp      = gpar(fontsize = 8, fontface = "bold"),
  column_names_gp   = gpar(fontsize = 8),
  column_names_rot  = 45
)
ht_pre <- draw(ht_pre)

# column order for annotation.
co <- column_order(ht_pre)
ord_idx <- if (is.list(co)) unlist(co, recursive = TRUE, use.names = FALSE) else as.integer(co)

# sanity check
stopifnot(length(ord_idx) == ncol(expression_matrix))
stopifnot(identical(sort(ord_idx), seq_len(ncol(expression_matrix))))

ord_samples <- colnames(expression_matrix)[ord_idx]

mes <- setNames(as.numeric(gsva_result[1, ]), colnames(gsva_result))

mes_ord <- mes[ord_samples]



mes_col <- ifelse(mes_ord > 0, "blue", "red")

top_anno <- HeatmapAnnotation(
  
  MES = anno_barplot(
    mes_ord,
    gp = gpar(fill = mes_col),
    axis_param = list(at = c(-1, 0, 1), labels = c("-1", "0", "1"))
  ),
  which = "column"
)

# coloring NO_TMMs.
highlight_samples_NOTMM <- metadata$SampleID[metadata$TMM == "NO_TMM"]
highlight_samples_ALT <- metadata$SampleID[metadata$TMM == "ALT"]
highlight_samples_Telomerase <- metadata$SampleID[metadata$TMM == "Telomerase"]


mat <- expression_matrix[, ord_samples, drop = FALSE]

sample_colors <- setNames(rep("black", ncol(mat)), colnames(mat))
sample_colors[highlight_samples_NOTMM]      <- "black"
sample_colors[highlight_samples_ALT]        <- "blue"
sample_colors[highlight_samples_Telomerase] <- "red"

# final heatmap.
ht_final <- Heatmap(
  expression_matrix[, ord_samples, drop = FALSE],
  row_order = rownames(expression_matrix),
  cluster_columns   = TRUE,         
  cluster_rows      = FALSE,
  show_column_names = TRUE,
  top_annotation    = top_anno,
  row_names_gp      = gpar(fontsize = 8, fontface = "bold"),
  column_names_gp   =  gpar(col = sample_colors[colnames(mat)], fontsize = 5),
  column_names_rot  = 45,
  cell_fun = function(j, i, x, y, w, h, fill) {
    grid.rect(x, y, w, h, gp = gpar(col = "black", fill = NA, lwd = 0.25))
  })


pdf("mesHeatmap.pdf", width = 15, height = 10)
draw(ht_final)
dev.off()# computing clustering order first(no annotation yet).









