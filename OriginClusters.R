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

gene_sets_phenotype <- list(
  ADR = c("KLF7", "GATA3", "HAND2", "PHOX2A", "ISL1", "HAND1", "PHOX2B", "TFAP2B", "GATA2", "SATB1", "SIX3", 
          "EYA1", "SOX11", "DACH1", "ASCL1", "HEY1", "KLF13", "PBX3"),
  MES = c("MEOX1", "WWTR1", "VIM", "CD44", "CBFB", "FOSL2", "MEOX2", "GLIS3", "TBX18", "NR3C1",
                              "PRRX1", "MEF2D", "BHLHE41", "RUNX2", "IRF1", "NOTCH2", "YAP1", "CREG1", "DCAF6", "FLI1",
                               "RUNX1", "IRF2", "JUN", "MAML2", "ZFP36L1")
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


######

phenotype_cluster <- as.data.frame(km_res2$cluster)
colnames(phenotype_cluster) = "Cluster"
phenotype_cluster$SampleID <- rownames(phenotype_cluster)
rownames(phenotype_cluster) = NULL   

phenotype_cluster <- phenotype_cluster %>%
  arrange(Cluster)
phenotype_cluster <- left_join(phenotype_cluster, metadata[, c("SampleID", "TMM", "COG.Risk.Group", "MYCN.status")], by = "SampleID")


table(phenotype_cluster$Cluster, phenotype_cluster$TMM)
table(phenotype_cluster$Cluster, phenotype_cluster$COG.Risk.Group)



