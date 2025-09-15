
### looking at TERT distributions across the phenotype (Expression Levels.)

library(tidyverse)
library(dplyr)
library(ggpubr)
library(survival)
library(survminer)
library(rstatix)


ExpressionTERT <- Expression["TERT", ]


ExpressionTERT<- data.frame(
  Gene = "TERT",
  Expression = as.numeric(ExpressionTERT),
  SampleID = colnames(Expression)
)

ExpressionTERT <- left_join(ExpressionTERT, metadata[, c("SampleID", "TMM_Case", "TMM")], by = "SampleID")

ggplot(ExpressionTERT, aes(x = TMM_Case, y = Expression, fill = TMM_Case, color = TMM)) +
  geom_boxplot(aes(group = TMM_Case), size=0.2, alpha=0.5) +
  geom_point(position = position_jitter(width = .2), size = 4)  + scale_fill_manual(values=c("NO_TMM"="lightblue","TMM"="red")) + 
  scale_color_manual(values=c("ALT"="yellow", "Telomerase"="darkred", "NO_TMM"="darkblue"))+ theme_classic() + 
  labs(x = "Class", y = "Expression Levels") +
  theme(
    axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1),
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 12, face = "bold"), legend.position = "none") +
  stat_compare_means(comparisons = list(c("NO_TMM","TMM")), method= "t.test",
                     method.args = list(alternative ="two.sided"), size = 8, tip.length = 0.01,
                     label.y = 12)

## looking at the difference between survival plot of MYCN-amplified and MYCN-not amplified NO_TMM samples.
ackerman_metadata <- read.table(file = 'NBL_Ackerman_CompleteMeta.txt', header = TRUE, sep = '\t')
survival_metadata <- read_delim("Metadata_TARGETFinal_08012024.txt")

survival_metadata <- survival_metadata[survival_metadata$TMM_Case == "NO_TMM", ]
survival_metadata$Vital.Status <- ifelse(survival_metadata$Vital.Status == "Dead", 1, 0)
fit <- survfit(Surv(Event.Free.Survival.Time.in.Days, Vital.Status) ~ MYCN.status, data = survival_metadata)


# plotting the graph.
ggsurvplot(fit,
           pval = TRUE,
           conf.int = TRUE,
           risk.table = TRUE,
           risk.table.col = "strata",
           linetype = "strata",
           surv.median.line = "hv",
           ncensor.plot = TRUE,
           ggtheme = theme_bw())


###########################################################################################

##### making cls and gct file for running in GSEA for pathway analysis.

# run DataCleaning.R first.

########### First - NO_TMM vs. ALT.

### gct file from Expression function.

write_gct <- function(mat, file, gene_desc = "na") {
  mat <- as.matrix(mat)
  genes <- rownames(mat)
  

  con <- file(file, open = "wt")
  on.exit(close(con), add = TRUE)
  writeLines("#1.2", con)
  writeLines(paste(nrow(mat), ncol(mat), sep = "\t"), con)
  writeLines(paste(c("Name", "Description", colnames(mat)), collapse = "\t"), con)
  
  for (i in seq_len(nrow(mat))) {
    line <- c(genes[i], "na", sprintf("%g", mat[i, ]))
    writeLines(paste(line, collapse = "\t"), con)
  }
}

metadata_ALT <- metadata_ALT %>%
  arrange(TMM)
ExpressionALT <- ExpressionALT[, match(metadata_ALT$SampleID, colnames(ExpressionALT))]

write_gct(ExpressionALT, "Expression.gct")

### .cls file from metadata.

out_cls   <- "classes.cls"
 
classes <- as.character(metadata_ALT[["TMM"]])

levels <- unique(classes)
k <- length(levels)
class_index <- match(classes, levels) - 1  # 0-based.

con <- file(out_cls, open = "wt")
on.exit(close(con), add = TRUE)
writeLines(sprintf("%d %d 1", ncol(ExpressionALT), k), con)
writeLines(paste("#", paste(levels, collapse = " ")), con)
writeLines(paste(class_index, collapse = " "), con)
close(con)



########### Second - NO_TMM vs. Telomerase.

### gct file from Expression function.

metadata_Telomerase <- metadata_Telomerase %>%
  arrange(TMM)
ExpressionTelomerase <- ExpressionTelomerase[, match(metadata_Telomerase$SampleID, colnames(ExpressionTelomerase))]


write_gct(ExpressionTelomerase, "Expression2.gct")

### .cls file from metadata.

out_cls   <- "classes2.cls"

classes <- as.character(metadata_Telomerase[["TMM"]])

levels <- unique(classes)
k <- length(levels)
class_index <- match(classes, levels) - 1  # 0-based.

con <- file(out_cls, open = "wt")
on.exit(close(con), add = TRUE)
writeLines(sprintf("%d %d 1", ncol(ExpressionTelomerase), k), con)
writeLines(paste("#", paste(levels, collapse = " ")), con)
writeLines(paste(class_index, collapse = " "), con)
close(con)

#######################################################################################

## Looking at the difference between High Risk and Low Risk NO_TMM in metadata.
source("DataCleaning.R")

NO_TMM_metadata <- metadata[metadata$TMM_Case == "NO_TMM", ]
NO_TMM_metadata <- NO_TMM_metadata %>%
  arrange(COG.Risk.Group)

# removing the samples that did not follow the cluster order.
NO_TMM_metadata <- NO_TMM_metadata[!NO_TMM_metadata$SampleID %in% c("TARGET.30.PALCBW.01A", "TARGET.30.PASPER.01A", "TARGET.30.PASJZC.01A", "TARGET.30.PARZIP.01A"), ]


# only selecting few relevant columns from everything -- age, telomere content, c-circles, apbs, tert structural variation, tert fpkm, gender, p53-cdkn2a, alk, ras-mapk, age, survival stats, inss stage, mycn status, ploidy, histology, grade.
NO_TMM_metadata <- NO_TMM_metadata[, colnames(NO_TMM_metadata) %in% c("SampleID", "ATRX.Status", "Telomere.Content", "C.Circles", "APBs",
                                                                      "TERT.SV", "TERT.FPKM", "p53.CDKN2A", "RAS.MAPK", "ALK", "Gender",
                                                                      "Age.at.Diagnosis.in.Days", "Event.Free.Survival.Time.in.Days", "Vital.Status",
                                                                      "Overall.Survival.Time.in.Days", "First.Event", "INSS.Stage", "MYCN.status", "Ploidy", "Histology",
                                                                      "Grade",  "COG.Risk.Group")]


## survival graph for High vs. Intermediate vs. Low Risk NO_TMM.

# OS
NO_TMM_metadata$Vital.Status <- ifelse(NO_TMM_metadata$Vital.Status == "Dead", 1, 0)
fit <- survfit(Surv(Overall.Survival.Time.in.Days, Vital.Status) ~ COG.Risk.Group, data = NO_TMM_metadata)


# plotting the graph.
ggsurvplot(fit,
           pval = TRUE,
           conf.int = TRUE,
           risk.table = TRUE,
           risk.table.col = "strata",
           linetype = "strata",
           surv.median.line = "hv",
           ncensor.plot = TRUE,
           ggtheme = theme_bw())

#EFS.
NO_TMM_metadata <- NO_TMM_metadata %>%
  mutate(
    EFS_event = case_when(
      First.Event %in% c("Relapse", "Progression", "Second Malignant Neoplasm", "Death") ~ 1L,
      First.Event %in% c("Censored") ~ 0L,
      TRUE ~ 0L
    )
  )

fit <- survfit(Surv(Event.Free.Survival.Time.in.Days, EFS_event) ~ COG.Risk.Group,
                   data = NO_TMM_metadata)

ggsurvplot(fit,
  pval = TRUE,
  conf.int = TRUE,
  risk.table = TRUE,
  risk.table.col = "strata",
  surv.median.line = "hv",
  ggtheme = theme_bw())

##########################################################################################

# fisher's test to see relationship between Risk Group and other metadata columns.

NO_TMM_metadata2 <- NO_TMM_metadata[NO_TMM_metadata$Histology != "Unknown", ]
table_var <- table(NO_TMM_metadata2$COG.Risk.Group, NO_TMM_metadata2$Histology)
fisher.test(table_var)

NO_TMM_metadata2 <- NO_TMM_metadata[NO_TMM_metadata$Grade != "Unknown", ]
table_var <- table(NO_TMM_metadata2$COG.Risk.Group, NO_TMM_metadata2$Grade)
fisher.test(table_var)

NO_TMM_metadata2 <- NO_TMM_metadata[NO_TMM_metadata$Ploidy != "Unknown", ]
table_var <- table(NO_TMM_metadata2$COG.Risk.Group, NO_TMM_metadata2$Ploidy)
fisher.test(table_var)

NO_TMM_metadata2 <- NO_TMM_metadata[NO_TMM_metadata$MYCN.status != "Unknown", ]
table_var <- table(NO_TMM_metadata2$COG.Risk.Group, NO_TMM_metadata2$MYCN.status)
fisher.test(table_var)

NO_TMM_metadata2 <- NO_TMM_metadata[NO_TMM_metadata$INSS.Stage != "Unknown", ]
table_var <- table(NO_TMM_metadata2$COG.Risk.Group, NO_TMM_metadata2$INSS.Stage)
fisher.test(table_var)

NO_TMM_metadata2 <- NO_TMM_metadata[NO_TMM_metadata$ALK != "N/A", ]
table_var <- table(NO_TMM_metadata2$COG.Risk.Group, NO_TMM_metadata2$ALK)
fisher.test(table_var)


#############################################################################################

# GSVA score difference between two phenotypes using the best signature.

library(GSVA)
source("DataCleaning.R")

gene_list <- list("Gene List" = c("CPNE8", "PGM2L1", "LIFR", "CNR1", "HECW2", "CPNE3", "HOXC9", "SNX16", "IGSF10"))
gene_list_ <- list("Gene List" = c("PRR7", "IGLV6-57", "SAC3D1", "CCDC86", "DDN"))
gene_list_together <- list("Gene List" = c("CPNE8", "PGM2L1", "LIFR", "CNR1", "HECW2", "CPNE3", "HOXC9", "SNX16", "IGSF10",
                                           "PRR7", "IGLV6-57", "SAC3D1", "CCDC86", "DDN"))


Expression <- as.matrix(Expression)
gsva <- gsvaParam(Expression, gene_list, kcdf = "Gaussian")
GSVA_result <- gsva(gsva)
GSVA_df <- as.data.frame(GSVA_result)
GSVA_long <- pivot_longer(GSVA_df, cols = everything(), names_to = "SampleID", values_to = "GSVA_Score")

GSVA_long <- merge(GSVA_long, metadata[, c("SampleID", "TMM_Case")], by = "SampleID")

gsva <- gsvaParam(Expression, gene_list_, kcdf = "Gaussian")
GSVA_result <- gsva(gsva)
GSVA_df <- as.data.frame(GSVA_result)
GSVA_long_ <- pivot_longer(GSVA_df, cols = everything(), names_to = "SampleID", values_to = "GSVA_Score")

GSVA_long_ <- merge(GSVA_long_, metadata[, c("SampleID", "TMM_Case")], by = "SampleID")


# 3. Boxplot.
ggplot(GSVA_long, aes(x = TMM_Case, y = GSVA_Score, fill = TMM_Case, color = TMM_Case)) +
  geom_boxplot(size = 0.2, alpha = 0.5, outlier.shape = NA) +
  geom_point(position = position_jitter(width = 0.2), size = 3) +
  scale_fill_manual(values = c("TMM" = "lightpink2", 
                               "NO_TMM" = "lightgreen")) +
  scale_color_manual(values = c("TMM"="darkred", 
                                "NO_TMM" = "darkgreen")) +
  theme_classic() +
  labs(x = "Class", y = "gsva Score") +
  theme(
    axis.text.x = element_text(vjust = 1, hjust = 1),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 9, face = "bold"),
    legend.position = "none"
  ) +
  stat_compare_means(comparisons = list(c("TMM","NO_TMM")), method= "t.test",
                     method.args = list(alternative ="two.sided"), size = 4, tip.length = 0.01,
                     label.y = 1.2)

ggplot(GSVA_long_, aes(x = TMM_Case, y = GSVA_Score, fill = TMM_Case, color = TMM_Case)) +
  geom_boxplot(size = 0.2, alpha = 0.5, outlier.shape = NA) +
  geom_point(position = position_jitter(width = 0.2), size = 3) +
  scale_fill_manual(values = c("TMM" = "lightpink2", 
                               "NO_TMM" = "lightgreen")) +
  scale_color_manual(values = c("TMM"="darkred", 
                                "NO_TMM" = "darkgreen")) +
  theme_classic() +
  labs(x = "Class", y = "gsva Score") +
  theme(
    axis.text.x = element_text(vjust = 1, hjust = 1),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 10, face = "bold"),
    legend.position = "none"
  ) +
  stat_compare_means(comparisons = list(c("TMM","NO_TMM")), method= "t.test",
                     method.args = list(alternative ="two.sided"), size = 4, tip.length = 0.01,
                     label.y = 1)



############################################################################################

# t-test for gsva signature genes between NO_TMM clusters.
source('clusterAllSamples.R')

clusters_target <- as.data.frame(km_res_sig$cluster)
colnames(clusters_target) <- c("Cluster")

clusters_target <- clusters_target[rownames(clusters_target) %in% metadata$SampleID[metadata$TMM == "NO_TMM"], , drop = FALSE]
signature_genes <- c("CPNE8", "PGM2L1", "LIFR", "CNR1", "HECW2", "CPNE3", "HOXC9", "SNX16", "IGSF10", "PRR7", "IGLV6-57", "SAC3D1", "CCDC86", "DDN")

# now running t-test for these genes across the two clusters. Need to filter out genes that are statistically different.

##### 2) t-test of the DEGs.
# Initializing an empty data frame for storing t-test results.
dge_gene <- Expression[rownames(Expression) %in% signature_genes, ,drop = FALSE]

t_test_results <- data.frame(
  Gene = character(),
  p_value_t_test = numeric(),
  stringsAsFactors = FALSE
)


# Looping through each gene in dge_gene.
for (i in 1:length(signature_genes)) {
  
  gene_id <- signature_genes[i]
  
  
  # Extracting the expression values for this gene, keeping it as a matrix. We are only selecting significant genes (all samples for this gene).
  gene_Expression <- dge_gene[rownames(dge_gene) == gene_id, , drop = FALSE]
  
  # Creating C1 and C2 group. For the gene we are working with, all NO_TMM sample expression data placed in c1_group and c2_group based on their cluster.
  c1_group <- gene_Expression[, clusters_target == 1, drop = FALSE]
  c2_group <- gene_Expression[, clusters_target == 2, drop = FALSE]
  
  # for a gene, compare c1 and c2 samples.
  t_test <- t.test(c1_group, c2_group)
  p_value_t_test <- t_test$p.value
  
  
  # Storing the results.
  t_test_results <- rbind(t_test_results, data.frame(
    Gene = gene_id,
    p_value_t_test = p_value_t_test
  ))
  
  
  
  # removing un-required intermediates formed while forming the t-test table above.
  rm(gene_id)
  rm(gene_Expression)
  rm(c1_group)
  rm(c2_group)
  rm(p_value_t_test)
  rm(t_test)
}

# filtering to only maintain genes where p-value is greater than 0.1.
t_test_results <- t_test_results[t_test_results$p_value_t_test > 0.1, ]

# remaining 9 genes -- clustering based on these 9 genes.

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

gene_list_together <- t_test_results$Gene


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

# From the elbow plot, graph seems to plateau at PC:1-2.
pc_scores <- pca_result$x[, 1:2]

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


# looking at extend -- how good it is in differentiating telomerase and non-telomerase samples in TARGET and Ackerman.


# first TARGET sample.
source("DataCleaning.R")

extendScores <- RunEXTEND(as.matrix(Expression))

telomeraseScores <- read_delim("TelomeraseScores.txt")

telomeraseScores <- telomeraseScores[, c("SampleID", "NormEXTENDScores")]
telomeraseScores <- as.data.frame(telomeraseScores)
telomeraseScores <- left_join(telomeraseScores, metadata[, c("SampleID", "TMM")], by = "SampleID")

## boxplot showing the difference.
ggplot(telomeraseScores, aes(x = TMM, y = NormEXTENDScores, fill = TMM, color = TMM)) +
  geom_boxplot(size = 0.2, alpha = 0.5, outlier.shape = NA) +
  geom_point(position = position_jitter(width = 0.2), size = 3) +
  scale_fill_manual(values = c("Telomerase" = "lightpink2", 
                               "NO_TMM" = "lightgreen",
                               "ALT" = "blue")) +
  scale_color_manual(values = c("Telomerase"="darkred", 
                                "NO_TMM" = "darkgreen",
                                "ALT" = "darkblue")) +
  theme_classic() +
  labs(x = "Class", y = "Normalized EXTEND Score") +
  theme(
    axis.text.x = element_text(vjust = 1, hjust = 1),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 12, face = "bold"),
    legend.position = "none"
  ) +
  stat_compare_means(comparisons = list(c("Telomerase","NO_TMM")), method= "t.test",
                     method.args = list(alternative ="two.sided"), size = 6, tip.length = 0.01,
                     label.y = 1.2) +
  stat_compare_means(comparisons = list(c("Telomerase","ALT")), method= "t.test",
                     method.args = list(alternative ="two.sided"), size = 6, tip.length = 0.01,
                     label.y = 1.4)


## now, ackerman sample.

gct_file <- parse_gctx("Neuroblastoma_208Samples.gct")
ackerman_NB <- gct_file@mat

# Loading metadata.
ackerman_metadata <- read.table(file = 'NBL_Ackerman_CompleteMeta.txt', header = TRUE, sep = '\t')

# only including SampleID in microarray data present in metadata.
ackerman_NB <- ackerman_NB[, colnames(ackerman_NB) %in% ackerman_metadata$SampleID]
ackerman_metadata <- ackerman_metadata[ackerman_metadata$SampleID %in% colnames(ackerman_NB), ]

ackerman_metadata <- ackerman_metadata %>%
  arrange(TMM_Case)

ackerman_NB <- ackerman_NB[, match(ackerman_metadata$SampleID, colnames(ackerman_NB))]

extendScores <- RunEXTEND(as.matrix(ackerman_NB))

telomeraseScores <- read_delim("TelomeraseScores.txt")

telomeraseScores <- telomeraseScores[, c("SampleID", "NormEXTENDScores")]
telomeraseScores <- as.data.frame(telomeraseScores)
telomeraseScores <- left_join(telomeraseScores, ackerman_metadata[, c("SampleID", "TMM_Category")], by = "SampleID")

## boxplot showing the difference.
ggplot(telomeraseScores, aes(x = TMM_Category, y = NormEXTENDScores, fill = TMM_Category, color = TMM_Category)) +
  geom_boxplot(size = 0.2, alpha = 0.5, outlier.shape = NA) +
  geom_point(position = position_jitter(width = 0.2), size = 3) +
  scale_fill_manual(values = c("Telomerase" = "lightpink2", 
                               "NO_TMM" = "lightgreen",
                               "ALT" = "blue")) +
  scale_color_manual(values = c("Telomerase"="darkred", 
                                "NO_TMM" = "darkgreen",
                                "ALT" = "darkblue")) +
  theme_classic() +
  labs(x = "Class", y = "Normalized EXTEND Score") +
  theme(
    axis.text.x = element_text(vjust = 1, hjust = 1),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 12, face = "bold"),
    legend.position = "none"
  ) +
  stat_compare_means(comparisons = list(c("Telomerase","NO_TMM")), method= "t.test",
                     method.args = list(alternative ="two.sided"), size = 6, tip.length = 0.01,
                     label.y = 1.2) +
  stat_compare_means(comparisons = list(c("Telomerase","ALT")), method= "t.test",
                     method.args = list(alternative ="two.sided"), size = 6, tip.length = 0.01,
                     label.y = 1.4)

#####>>>>>>>>############################################################################
##############################################################################################

# heatmap to see the relationship between cell of origin timeline markers and TMM status.

library(ComplexHeatmap)
library(GSVA)
source("DataCleaning.R")

scp <- c("CDH19", "PLP1", "ERBB3", "MPZ", "ERBB4")
chromaffin <- c("TH", "DBH", "DDC", "CHGA", "PNMT")
neuroblast <- c("NEFM", "GAP43", "STMN2", "ISL1")
early_neuroblast <- c("ALK")
late_neuroblast <- c("SYN3", "IL7", "GAP43", "STMN2", "ISL1")
bridge <- c("ERBB4", "ASCL1", "CDH9", "CTTNBP2")
proliferation <- c("MKI67", "PCNA", "BIRC5", "CEP55", "TOP2A", "CDK1", "CDC20", "CCNB1", "CCNB2", "CCNA2", 
                   "UBE2C", "AURKA", "AURKB", "PLK1", "PTTG1", "NUSAP1")


expression_matrix <- t(scale(t(Expression[rownames(Expression) %in% c( "GAP43", "STMN2", "ISL1", "SYN3", "IL7"), ]))) # chromaffin and chromaffin adjacent genes.

# adding gsva as annotation.

# Running GSVA.

gene_set_list_adrenergic <-list(ADRN = c("PHOX2A", "PHOX2B", "HAND2", "ISL1", "GATA3", "ASCL1", "TBX2", "DBH", "TH"))
gene_set_list_mesenchymal <- list(MES = c("PRRX1", "RUNX1", "FOSL1", "JUN", "YAP1", "WWTR1", "MEOX1", 
                                          "MEOX2", "SNAI2", "CD44", "FN1", "VIM"))

gene_set_list_chromaffin <- list(chromaffin = c("TH", "DBH", "DDC", "CHGA", "PNMT"))
gene_set_list_scp <- list(scp =  c("CDH19", "PLP1", "ERBB3", "MPZ", "ERBB4"))
gene_set_list_neuroblast <- list(neuroblast = c("NEFM", "GAP43", "STMN2", "ISL1"))
gene_set_list_lateNeuroblast <- list(lateNeuroblast = c("SYN3", "IL7", "GAP43", "STMN2", "ISL1"))
gene_set_list_bridge <- list(bridge = c("ERBB4", "ASCL1", "CDH9", "CTTNBP2"))
gene_set_list_proliferation <- list(proliferation = c("MKI67", "PCNA", "BIRC5", "CEP55", "TOP2A", "CDK1", "CDC20", "CCNB1", 
                                                      "CCNB2", "CCNA2", "UBE2C", "AURKA", "AURKB", "PLK1", "PTTG1", "NUSAP1"))


gsvapar <- gsvaParam(as.matrix(Expression), gene_set_list_adrenergic, kcdf = "Gaussian")
gsva_result <- gsva(gsvapar)

gsvapar2 <- gsvaParam(as.matrix(Expression), gene_set_list_mesenchymal, kcdf = "Gaussian")
gsva_result2 <- gsva(gsvapar2)

gsvapar3 <- gsvaParam(as.matrix(Expression), gene_set_list_chromaffin, kcdf = "Gaussian")
gsva_result3 <- gsva(gsvapar3)

gsvapar4 <- gsvaParam(as.matrix(Expression), gene_set_list_scp, kcdf = "Gaussian")
gsva_result4 <- gsva(gsvapar4)

gsvapar5 <- gsvaParam(as.matrix(Expression), gene_set_list_neuroblast, kcdf = "Gaussian")
gsva_result5 <- gsva(gsvapar5)

gsvapar6 <- gsvaParam(as.matrix(Expression), gene_set_list_lateNeuroblast, kcdf = "Gaussian")
gsva_result6 <- gsva(gsvapar6)

gsvapar7 <- gsvaParam(as.matrix(Expression), gene_set_list_bridge, kcdf = "Gaussian")
gsva_result7 <- gsva(gsvapar7)

gsvapar8 <- gsvaParam(as.matrix(Expression), gene_set_list_proliferation, kcdf = "Gaussian")
gsva_result8 <- gsva(gsvapar8)


adr <- as.numeric(gsva_result[1, ])
names(adr) <- colnames(gsva_result)

mes <- as.numeric(gsva_result2[1, ])
names(mes) <- colnames(gsva_result2)

chromaffin <- as.numeric(gsva_result3[1, ])
names(chromaffin) <- colnames(gsva_result3)

scp <- as.numeric(gsva_result4[1, ])
names(scp) <- colnames(gsva_result4)

neuroblast <- as.numeric(gsva_result5[1, ])
names(neuroblast) <- colnames(gsva_result5)

lateNeuroblast <- as.numeric(gsva_result6[1, ])
names(lateNeuroblast) <- colnames(gsva_result6)

bridge <- as.numeric(gsva_result7[1, ])
names(bridge) <- colnames(gsva_result7)

proliferation <- as.numeric(gsva_result8[1, ])
names(proliferation) <- colnames(gsva_result8)



stopifnot(all(colnames(Expression) %in% names(adr)))
stopifnot(all(colnames(Expression) %in% names(mes)))
stopifnot(all(colnames(Expression) %in% names(chromaffin)))
stopifnot(all(colnames(Expression) %in% names(scp)))
stopifnot(all(colnames(Expression) %in% names(neuroblast)))
stopifnot(all(colnames(Expression) %in% names(lateNeuroblast)))
stopifnot(all(colnames(Expression) %in% names(bridge)))
stopifnot(all(colnames(Expression) %in% names(proliferation)))


# computing clustering order first(no annotation yet).
ht_pre <- Heatmap(
  expression_matrix,
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
chromaffin <- setNames(as.numeric(gsva_result3[1, ]), colnames(gsva_result3))
scp <- setNames(as.numeric(gsva_result4[1, ]), colnames(gsva_result4))
neuroblast <- setNames(as.numeric(gsva_result5[1, ]), colnames(gsva_result5))
lateNeuroblast <- setNames(as.numeric(gsva_result6[1, ]), colnames(gsva_result6))
bridge <- setNames(as.numeric(gsva_result7[1, ]), colnames(gsva_result7))
proliferation <- setNames(as.numeric(gsva_result8[1, ]), colnames(gsva_result8))


# align to the display order
adr_ord <- adr[ord_samples]
mes_ord <- mes[ord_samples]
chromaffin_ord <- chromaffin[ord_samples]
scp_ord <- scp[ord_samples]
neuroblast_ord <- neuroblast[ord_samples]
lateNeuroblast_ord <- lateNeuroblast[ord_samples]
bridge_ord <- bridge[ord_samples]
proliferation_ord <- proliferation[ord_samples]


adr_col <- ifelse(adr_ord > 0, "blue", "red")
mes_col <- ifelse(mes_ord > 0, "blue", "red")
chromaffin_col <- ifelse(chromaffin_ord > 0, "blue", "red")
scp_col <- ifelse(scp_ord > 0, "blue", "red")
neuroblast_col <- ifelse(neuroblast_ord > 0, "blue", "red")
lateNeuroblast_col <- ifelse(lateNeuroblast_ord > 0, "blue", "red")
bridge_col <- ifelse(bridge_ord > 0, "blue", "red")
proliferation_col <- ifelse(proliferation_ord > 0, "blue", "red")

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
  Chromaffin = anno_barplot(
    chromaffin_ord,
    gp = gpar(fill = chromaffin_col),
    axis_param = list(at = c(-1, 0, 1), labels = c("-1", "0", "1"))
  ),
  SCP = anno_barplot(
    scp_ord,
    gp = gpar(fill = scp_col),
    axis_param = list(at = c(-1, 0, 1), labels = c("-1", "0", "1"))
  ),
  Neuroblast = anno_barplot(
    neuroblast_ord,
    gp = gpar(fill = neuroblast_col),
    axis_param = list(at = c(-1, 0, 1), labels = c("-1", "0", "1"))
  ),
  LateNeuroblast = anno_barplot(
    lateNeuroblast_ord,
    gp = gpar(fill = lateNeuroblast_col),
    axis_param = list(at = c(-1, 0, 1), labels = c("-1", "0", "1"))
  ),
  Bridge = anno_barplot(
    bridge_ord,
    gp = gpar(fill = bridge_col),
    axis_param = list(at = c(-1, 0, 1), labels = c("-1", "0", "1"))
  ),
   Proliferation = anno_barplot(
    proliferation_ord,
    gp = gpar(fill = proliferation_col),
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
  cluster_columns   = TRUE,             
  cluster_rows      = TRUE,
  show_column_names = TRUE,
  top_annotation    = top_anno,
  row_names_gp      = gpar(fontsize = 8, fontface = "bold"),
  column_names_gp   =  gpar(col = sample_colors[colnames(mat)], fontsize = 5),
  column_names_rot  = 45,
  cell_fun = function(j, i, x, y, w, h, fill) {
    grid.rect(x, y, w, h, gp = gpar(col = "black", fill = NA, lwd = 0.25))
  }
)


pdf("adr-mesHeatmapNeuroblast.pdf", width = 15, height = 7)
draw(ht_final)
dev.off()

######################################################################################################


########### test for:

## NO_TMM: high proliferation vs. low proliferation ~ Low vs. High Risk.
proliferation_ord <- as.data.frame(proliferation_ord)

colnames(proliferation_ord) <- "GSVA_Score"
proliferation_ord$SampleID <- rownames(proliferation_ord)
rownames(proliferation_ord) = NULL

proliferation_ord <- left_join(proliferation_ord, metadata[, c("SampleID", "TMM", "COG.Risk.Group", "MYCN.status")], by = "SampleID")

proliferation_ord_noTMM <- proliferation_ord[proliferation_ord$TMM == "NO_TMM", ]


# t-test doing difference between GSVA Score of high, intermediate and low risk in NO_TMM.
ggplot(proliferation_ord_noTMM, aes(x = COG.Risk.Group, y = GSVA_Score, fill = COG.Risk.Group, color = COG.Risk.Group)) +
  geom_boxplot(size = 0.2, alpha = 0.5, outlier.shape = NA) +
  geom_point(position = position_jitter(width = 0.2), size = 3) +
  scale_fill_manual(values = c("High Risk" = "lightpink2", 
                               "Low Risk" = "lightgreen",
                               "Intermediate Risk" = "blue")) +
  scale_color_manual(values = c("High Risk"="darkred", 
                                "Low Risk" = "darkgreen",
                                "Intermediate Risk" = "darkblue")) +
  theme_classic() +
  labs(x = "Risk Group", y = "Proliferation Marker GSVA Score") +
  theme(
    axis.text.x = element_text(vjust = 1, hjust = 1),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 12, face = "bold"),
    legend.position = "none"
  ) +
  stat_compare_means(comparisons = list(c("High Risk","Low Risk")), method= "t.test",
                     method.args = list(alternative ="two.sided"), size = 6, tip.length = 0.01,
                     label.y = 1.2) + 
  stat_compare_means(comparisons = list(c("High Risk","Intermediate Risk")), method= "t.test",
                     method.args = list(alternative ="two.sided"), size = 6, tip.length = 0.01,
                     label.y = 1.4) +
  stat_compare_means(comparisons = list(c("Low Risk","Intermediate Risk")), method= "t.test",
                     method.args = list(alternative ="two.sided"), size = 6, tip.length = 0.01,
                     label.y = 1.1)

# t-test doing difference between GSVA Score of Telomerase, ALT and NO_TMM.
ggplot(proliferation_ord, aes(x = TMM, y = GSVA_Score, fill = TMM, color = TMM)) +
  geom_boxplot(size = 0.2, alpha = 0.5, outlier.shape = NA) +
  geom_point(position = position_jitter(width = 0.2), size = 3) +
  scale_fill_manual(values = c("Telomerase" = "lightpink2", 
                               "NO_TMM" = "lightgreen",
                               "ALT" = "blue")) +
  scale_color_manual(values = c("Telomerase"="darkred", 
                                "NO_TMM" = "darkgreen",
                                "ALT" = "darkblue")) +
  theme_classic() +
  labs(x = "Class", y = "Proliferation Marker GSVA Score") +
  theme(
    axis.text.x = element_text(vjust = 1, hjust = 1),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 12, face = "bold"),
    legend.position = "none"
  ) +
  stat_compare_means(comparisons = list(c("ALT","NO_TMM")), method= "t.test",
                     method.args = list(alternative ="two.sided"), size = 6, tip.length = 0.01,
                     label.y = 1.2) +
  stat_compare_means(comparisons = list(c("Telomerase","NO_TMM")), method= "t.test",
                     method.args = list(alternative ="two.sided"), size = 6, tip.length = 0.01,
                     label.y = 1.4)


#######



## Late neuroblast ~ NO_TMM vs. others?
lateNeuroblast_ord <- as.data.frame(lateNeuroblast_ord)

colnames(lateNeuroblast_ord) <- "GSVA_Score"
lateNeuroblast_ord$SampleID <- rownames(lateNeuroblast_ord)
rownames(lateNeuroblast_ord) = NULL

lateNeuroblast_ord <- left_join(lateNeuroblast_ord, metadata[, c("SampleID", "TMM", "COG.Risk.Group", "MYCN.status")], by = "SampleID")

lateNeuroblast_ord_noTMM <- lateNeuroblast_ord[lateNeuroblast_ord$TMM == "NO_TMM", ]

# t-test doing difference between GSVA Score of high, intermediate and low risk in NO_TMM.
ggplot(lateNeuroblast_ord_noTMM, aes(x = COG.Risk.Group, y = GSVA_Score, fill = COG.Risk.Group, color = COG.Risk.Group)) +
  geom_boxplot(size = 0.2, alpha = 0.5, outlier.shape = NA) +
  geom_point(position = position_jitter(width = 0.2), size = 3) +
  scale_fill_manual(values = c("High Risk" = "lightpink2", 
                               "Low Risk" = "lightgreen",
                               "Intermediate Risk" = "blue")) +
  scale_color_manual(values = c("High Risk"="darkred", 
                                "Low Risk" = "darkgreen",
                                "Intermediate Risk" = "darkblue")) +
  theme_classic() +
  labs(x = "Risk Group", y = "Late Neuroblasts Markers GSVA Score") +
  theme(
    axis.text.x = element_text(vjust = 1, hjust = 1),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 12, face = "bold"),
    legend.position = "none"
  ) +
  stat_compare_means(comparisons = list(c("High Risk","Low Risk")), method= "t.test",
                     method.args = list(alternative ="two.sided"), size = 6, tip.length = 0.01,
                     label.y = 1.2) +
  stat_compare_means(comparisons = list(c("High Risk","Intermediate Risk")), method= "t.test",
                     method.args = list(alternative ="two.sided"), size = 6, tip.length = 0.01,
                     label.y = 1.4)


# t-test doing difference between GSVA Score of Telomerase, ALT and NO_TMM.
ggplot(lateNeuroblast_ord, aes(x = TMM, y = GSVA_Score)) +
  geom_boxplot(aes(fill = TMM, group = TMM), size = 0.2, alpha = 0.5, outlier.shape = NA) +
  geom_point(aes(color = TMM), position = position_jitter(width = 0.2), size = 3) +
  scale_fill_manual(values = c("Telomerase" = "lightpink2", 
                               "NO_TMM" = "lightgreen",
                               "ALT" = "blue")) +
  scale_color_manual(values = c("Telomerase"="darkred", 
                                "ALT" = "darkblue",
                                "NO_TMM" = "darkgreen")) +
  theme_classic() +
  labs(x = "Class", y = "Late Neuroblasts Markers GSVA Score") +
  theme(
    axis.text.x = element_text(vjust = 1, hjust = 1),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 12, face = "bold"),
    legend.position = "none"
  ) +
  stat_compare_means(comparisons = list(c("ALT","NO_TMM")), method= "t.test",
                     method.args = list(alternative ="two.sided"), size = 6, tip.length = 0.01,
                     label.y = 1.2) +
  stat_compare_means(comparisons = list(c("Telomerase","NO_TMM")), method= "t.test",
                     method.args = list(alternative ="two.sided"), size = 6, tip.length = 0.01,
                     label.y = 1.4)


# t-test doing difference between GSVA Score of mycn amplified, non-amplified.
lateNeuroblast_ord <- as.data.frame(lateNeuroblast_ord)
lateNeuroblast_ord <- lateNeuroblast_ord[lateNeuroblast_ord$MYCN.status %in% c("Amplified", "Not Amplified"), ] 

ggplot(lateNeuroblast_ord, aes(x = MYCN.status, y = GSVA_Score, fill = MYCN.status, color = MYCN.status)) +
  geom_boxplot(size = 0.2, alpha = 0.5, outlier.shape = NA) +
  geom_point(position = position_jitter(width = 0.2), size = 3) +
  scale_fill_manual(values = c("Amplified" = "lightpink2", 
                               "Not Amplified" = "lightgreen")) +
  scale_color_manual(values = c("Amplified"="darkred", 
                                "Not Amplified" = "darkgreen")) +
  theme_classic() +
  labs(x = "Risk Group", y = "Late Neuroblasts Markers GSVA Score") +
  theme(
    axis.text.x = element_text(vjust = 1, hjust = 1),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 12, face = "bold"),
    legend.position = "none"
  ) +
  stat_compare_means(comparisons = list(c("Amplified", "Not Amplified")), method= "t.test",
                     method.args = list(alternative = "two.sided"), size = 6, tip.length = 0.01,
                     label.y = 1.2)


# t-test doing difference between GSVA Score of high, intermediate and low risk in ALL samples.
ggplot(lateNeuroblast_ord, aes(x = COG.Risk.Group, y = GSVA_Score, fill = COG.Risk.Group, color = COG.Risk.Group)) +
  geom_boxplot(size = 0.2, alpha = 0.5, outlier.shape = NA) +
  geom_point(position = position_jitter(width = 0.2), size = 3) +
  scale_fill_manual(values = c("High Risk" = "lightpink2", 
                               "Low Risk" = "lightgreen",
                               "Intermediate Risk" = "blue")) +
  scale_color_manual(values = c("High Risk"="darkred", 
                                "Low Risk" = "darkgreen",
                                "Intermediate Risk" = "darkblue")) +
  theme_classic() +
  labs(x = "Risk Group", y = "Late Neuroblasts Markers GSVA Score") +
  theme(
    axis.text.x = element_text(vjust = 1, hjust = 1),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 12, face = "bold"),
    legend.position = "none"
  ) +
  stat_compare_means(comparisons = list(c("High Risk","Low Risk")), method= "t.test",
                     method.args = list(alternative ="two.sided"), size = 6, tip.length = 0.01,
                     label.y = 1.2)

## 


########################################################################################################
##### Now, correlation studies.

# GSVA NO_TMM vs. Proliferation Rank.
source("DataCleaning.R")

gene_list_noTMM <- list("NO_TMM Genes" = c("CPNE8", "PGM2L1", "LIFR", "CNR1", "HECW2", "CPNE3", "HOXC9", "SNX16", "IGSF10"))
gene_list_proliferation <- c("MKI67", "PCNA", "BIRC5", "CEP55", "TOP2A", "CDK1", "CDC20", "CCNB1", 
                                                    "CCNB2", "CCNA2", "UBE2C", "AURKA", "AURKB", "PLK1", "PTTG1", "NUSAP1")
# ranks of proliferation gene.
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

gene_list_proliferation <- intersect(gene_list_proliferation, rownames(ranked_TARGET_NBL))
gene_proliferation <- t(as.matrix(ranked_TARGET_NBL[gene_list_proliferation, , drop = FALSE]))

gsva <- gsvaParam(as.matrix(Expression), gene_list_noTMM, kcdf = "Gaussian")
GSVA_result <- gsva(gsva)
GSVA_df <- as.data.frame(GSVA_result)
GSVA_long_noTMM <- pivot_longer(GSVA_df, cols = everything(), names_to = "SampleID", values_to = "GSVA_Score")


GSVA_long_noTMM <- GSVA_long_noTMM[match(rownames(gene_proliferation),
                                         GSVA_long_noTMM$SampleID), , drop = FALSE]




y <- GSVA_long_noTMM$GSVA_Score


cor_values <- numeric(ncol(gene_proliferation))
p_values   <- numeric(ncol(gene_proliferation))


for (i in seq_len(ncol(gene_proliferation))) {
  a <- as.numeric(gene_proliferation[, i])
  b <- complete.cases(a, y)
  correlation <- cor.test(a[b], y[b], method = "spearman", exact = FALSE)
  cor_values[i] <- unname(correlation$estimate)
  p_values[i]   <- correlation$p.value
}

cor_results <- data.frame(
  Gene = colnames(gene_proliferation),
  estimate = cor_values,
  p_value = p_values,
  fdr = p.adjust(p_values, "BH"),
  stringsAsFactors = FALSE
)[order(p.adjust(p_values, "BH"), -abs(cor_values)), ]

noTMM_proliferation_sig <- cor_results[cor_results$fdr < 0.01, ]

### GSVA Proliferation vs. extend
source("DataCleaning.R")


gene_list_proliferation <- c("MKI67", "PCNA", "BIRC5", "CEP55", "TOP2A", "CDK1", "CDC20", "CCNB1", 
                             "CCNB2", "CCNA2", "UBE2C", "AURKA", "AURKB", "PLK1", "PTTG1", "NUSAP1")


gene_list_proliferation <- intersect(gene_list_proliferation, rownames(Expression))
gene_proliferation <- t(as.matrix(Expression[gene_list_proliferation, , drop = FALSE]))

extendScores <- RunEXTEND(as.matrix(Expression))

telomeraseScores <- read_delim("TelomeraseScores.txt")

telomeraseScores <- telomeraseScores[, c("SampleID", "NormEXTENDScores")]
telomeraseScores <- as.data.frame(telomeraseScores)
telomeraseScores <- left_join(telomeraseScores, metadata[, c("SampleID", "TMM")], by = "SampleID")
rownames(telomeraseScores) <- telomeraseScores$SampleID
telomeraseScores$SampleID = NULL

telomeraseScores[match(rownames(gene_proliferation), rownames(telomeraseScores)), ]

y <- telomeraseScores$NormEXTENDScores


cor_values <- numeric(ncol(gene_proliferation))
p_values   <- numeric(ncol(gene_proliferation))



for (i in seq_len(ncol(gene_proliferation))) {
  a <- as.numeric(gene_proliferation[, i])
  b <- complete.cases(a, y)
  correlation <- cor.test(a[b], y[b], method = "spearman", exact = FALSE)
  cor_values[i] <- unname(correlation$estimate)
  p_values[i]   <- correlation$p.value
}

cor_results <- data.frame(
  Gene = colnames(gene_proliferation),
  estimate = cor_values,
  p_value = p_values,
  fdr = p.adjust(p_values, "BH"),
  stringsAsFactors = FALSE
)[order(p.adjust(p_values, "BH"), -abs(cor_values)), ]

extend_proliferation_sig <- cor_results[cor_results$fdr < 0.01, ]

############################

##### Proliferation vs. GSVA ADR.

source("DataCleaning.R")


adrenergic_markers <- list(ADR = c("PHOX2A", "PHOX2B", "HAND2", "ISL1", "GATA3", "ASCL1", "TBX2", "DBH", "TH"))




# GSVA for adrenergic marker.
Expression <- as.matrix(Expression)
adr_gsva <- gsvaParam(Expression, adrenergic_markers, kcdf = "Gaussian")
GSVA_result_adr <- gsva(adr_gsva)
GSVA_df_adr <- as.data.frame(GSVA_result_adr)
GSVA_long_adr <- pivot_longer(GSVA_df_adr, cols = everything(), names_to = "SampleID", values_to = "GSVA_Score")

GSVA_long_adr <- GSVA_long_adr[match(rownames(gene_proliferation), GSVA_long_adr$SampleID), ]

y <- GSVA_long_adr$GSVA_Score

cor_values <- numeric(ncol(gene_proliferation))
p_values   <- numeric(ncol(gene_proliferation))



for (i in seq_len(ncol(gene_proliferation))) {
  a <- as.numeric(gene_proliferation[, i])
  b <- complete.cases(a, y)
  correlation <- cor.test(a[b], y[b], method = "spearman", exact = FALSE)
  cor_values[i] <- unname(correlation$estimate)
  p_values[i]   <- correlation$p.value
}

cor_results <- data.frame(
  Gene = colnames(gene_proliferation),
  estimate = cor_values,
  p_value = p_values,
  fdr = p.adjust(p_values, "BH"),
  stringsAsFactors = FALSE
)[order(p.adjust(p_values, "BH"), -abs(cor_values)), ]

adr_proliferation_sig <- cor_results[cor_results$fdr < 0.01, ]


##### Proliferation vs. GSVA MES.

source("DataCleaning.R")

mesenchymal_markers <- list(MES = c("PRRX1", "RUNX1", "FOSL1", "JUN", "YAP1", "WWTR1", "MEOX1", "MEOX2",
                                    "SNAI2", "CD44", "FN1", "VIM"))




# GSVA for mesenchymal marker.
Expression <- as.matrix(Expression)
mes_gsva <- gsvaParam(Expression, mesenchymal_markers, kcdf = "Gaussian")
GSVA_result_mes <- gsva(mes_gsva)
GSVA_df_mes <- as.data.frame(GSVA_result_mes)
GSVA_long_mes <- pivot_longer(GSVA_df_mes, cols = everything(), names_to = "SampleID", values_to = "GSVA_Score")

GSVA_long_mes <- GSVA_long_mes[match(rownames(gene_proliferation), GSVA_long_mes$SampleID), ]

y <- GSVA_long_mes$GSVA_Score

cor_values <- numeric(ncol(gene_proliferation))
p_values   <- numeric(ncol(gene_proliferation))



for (i in seq_len(ncol(gene_proliferation))) {
  a <- as.numeric(gene_proliferation[, i])
  b <- complete.cases(a, y)
  correlation <- cor.test(a[b], y[b], method = "spearman", exact = FALSE)
  cor_values[i] <- unname(correlation$estimate)
  p_values[i]   <- correlation$p.value
}

cor_results <- data.frame(
  Gene = colnames(gene_proliferation),
  estimate = cor_values,
  p_value = p_values,
  fdr = p.adjust(p_values, "BH"),
  stringsAsFactors = FALSE
)[order(p.adjust(p_values, "BH"), -abs(cor_values)), ]

mes_proliferation_sig <- cor_results[cor_results$fdr < 0.01, ]


##########

### cell cycle inhibitor vs. NO_TMM GSVA.

cellCycle_inhibitor <- c("FOXO1", "BTG1", "CDKN1B", "CDKN2A", "HES1", "EGR1", "CDKN2C", "CDKN2B", "HOPX", "CDKN1A", 
                         "KLF4", "CDKN1C", "TGFB1", "CDKN2D", "SMAD7", "KLF10")

# ranks of cell cycle inhibitor gene.
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

cellCycle_inhibitor <- intersect(cellCycle_inhibitor, rownames(ranked_TARGET_NBL))
cellCycle_inhibitor_matrix <- t(as.matrix(ranked_TARGET_NBL[cellCycle_inhibitor, , drop = FALSE]))


GSVA_long_noTMM <- GSVA_long_noTMM[match(rownames(cellCycle_inhibitor_matrix), GSVA_long_noTMM$SampleID), ]

y <- GSVA_long_noTMM$GSVA_Score

cor_values <- numeric(ncol(cellCycle_inhibitor_matrix))
p_values   <- numeric(ncol(cellCycle_inhibitor_matrix))



for (i in seq_len(ncol(cellCycle_inhibitor_matrix))) {
  a <- as.numeric(cellCycle_inhibitor_matrix[, i])
  b <- complete.cases(a, y)
  correlation <- cor.test(a[b], y[b], method = "spearman", exact = FALSE)
  cor_values[i] <- unname(correlation$estimate)
  p_values[i]   <- correlation$p.value
}

cor_results <- data.frame(
  Gene = colnames(cellCycle_inhibitor_matrix),
  estimate = cor_values,
  p_value = p_values,
  fdr = p.adjust(p_values, "BH"),
  stringsAsFactors = FALSE
)[order(p.adjust(p_values, "BH"), -abs(cor_values)), ]

noTMM_inhibitor_sig <- cor_results[cor_results$fdr < 0.01, ]




##############

##### ADR GSVA vs. Chromaffin ranks correlation.
source("DataCleaning.R")

adrenergic_markers <- list(ADR = c("PHOX2A", "PHOX2B", "HAND2", "ISL1", "GATA3", "ASCL1", "TBX2", "DBH", "TH"))

chromaffin_genes <- c("DBH", "TH", "CHGA", "DDC", "PNMT")



# GSVA for adrenergic marker.
Expression <- as.matrix(Expression)
adr_gsva <- gsvaParam(Expression, adrenergic_markers, kcdf = "Gaussian")
GSVA_result_adr <- gsva(adr_gsva)
GSVA_df_adr <- as.data.frame(GSVA_result_adr)
GSVA_long_adr <- pivot_longer(GSVA_df_adr, cols = everything(), names_to = "SampleID", values_to = "GSVA_Score")


# ranks of chromaffin gene.
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

chromaffin_genes <- intersect(chromaffin_genes, rownames(ranked_TARGET_NBL))
chromaffin_matrix <- t(as.matrix(ranked_TARGET_NBL[chromaffin_genes, , drop = FALSE]))


GSVA_long_adr <- GSVA_long_adr[match(rownames(chromaffin_matrix), GSVA_long_adr$SampleID), ]

y <- GSVA_long_adr$GSVA_Score

cor_values <- numeric(ncol(chromaffin_matrix))
p_values   <- numeric(ncol(chromaffin_matrix))



for (i in seq_len(ncol(chromaffin_matrix))) {
  a <- as.numeric(chromaffin_matrix[, i])
  b <- complete.cases(a, y)
  correlation <- cor.test(a[b], y[b], method = "spearman", exact = FALSE)
  cor_values[i] <- unname(correlation$estimate)
  p_values[i]   <- correlation$p.value
}

cor_results <- data.frame(
  Gene = colnames(chromaffin_matrix),
  estimate = cor_values,
  p_value = p_values,
  fdr = p.adjust(p_values, "BH"),
  stringsAsFactors = FALSE
)[order(p.adjust(p_values, "BH"), -abs(cor_values)), ]

chromaffin_sig <- cor_results[cor_results$fdr < 0.01, ]


##### Chromaffin vs. NO_TMM GSVA.
source("DataCleaning.R")

gene_list_noTMM <- list("NO_TMM Genes" = c("CPNE8", "PGM2L1", "LIFR", "CNR1", "HECW2", "CPNE3", "HOXC9", "SNX16", "IGSF10"))
chromaffin_genes <- c("DBH", "TH", "CHGA", "DDC", "PNMT")



# GSVA for no_TMM marker.
Expression <- as.matrix(Expression)
noTMM_gsva <- gsvaParam(Expression, gene_list_noTMM, kcdf = "Gaussian")
GSVA_result_noTMM <- gsva(noTMM_gsva)
GSVA_df_noTMM <- as.data.frame(GSVA_result_noTMM)
GSVA_long_noTMM <- pivot_longer(GSVA_df_noTMM, cols = everything(), names_to = "SampleID", values_to = "GSVA_Score")


# ranks of chromaffin gene.
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

chromaffin_genes <- intersect(chromaffin_genes, rownames(ranked_TARGET_NBL))
chromaffin_matrix <- t(as.matrix(ranked_TARGET_NBL[chromaffin_genes, , drop = FALSE]))


GSVA_long_noTMM <- GSVA_long_noTMM[match(rownames(chromaffin_matrix), GSVA_long_noTMM$SampleID), ]

y <- GSVA_long_noTMM$GSVA_Score

cor_values <- numeric(ncol(chromaffin_matrix))
p_values   <- numeric(ncol(chromaffin_matrix))



for (i in seq_len(ncol(chromaffin_matrix))) {
  a <- as.numeric(chromaffin_matrix[, i])
  b <- complete.cases(a, y)
  correlation <- cor.test(a[b], y[b], method = "spearman", exact = FALSE)
  cor_values[i] <- unname(correlation$estimate)
  p_values[i]   <- correlation$p.value
}

cor_results <- data.frame(
  Gene = colnames(chromaffin_matrix),
  estimate = cor_values,
  p_value = p_values,
  fdr = p.adjust(p_values, "BH"),
  stringsAsFactors = FALSE
)[order(p.adjust(p_values, "BH"), -abs(cor_values)), ]

noTMM_chromaffin_sig <- cor_results[cor_results$fdr < 0.01, ]


##### MES GSVA vs. Schwann Cell.
source("DataCleaning.R")

gene_list_mesenchymal <- list(MES = c("PRRX1", "RUNX1", "FOSL1", "JUN", "YAP1", "WWTR1", "MEOX1", 
                                          "MEOX2", "SNAI2", "CD44", "FN1", "VIM"))
scp_genes <- c("CDH19", "PLP1", "ERBB3", "MPZ", "ERBB4")

# GSVA for MES marker.
Expression <- as.matrix(Expression)
mes_gsva <- gsvaParam(Expression, gene_list_mesenchymal, kcdf = "Gaussian")
GSVA_result_mes <- gsva(mes_gsva)
GSVA_df_mes <- as.data.frame(GSVA_result_mes)
GSVA_long_mes <- pivot_longer(GSVA_df_mes, cols = everything(), names_to = "SampleID", values_to = "GSVA_Score")


# ranks of schwann cell markers.
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

scp_genes <- intersect(scp_genes, rownames(ranked_TARGET_NBL))
scp_matrix <- t(as.matrix(ranked_TARGET_NBL[scp_genes, , drop = FALSE]))


GSVA_long_mes <- GSVA_long_mes[match(rownames(scp_matrix), GSVA_long_mes$SampleID), ]

y <- GSVA_long_mes$GSVA_Score

cor_values <- numeric(ncol(scp_matrix))
p_values   <- numeric(ncol(scp_matrix))



for (i in seq_len(ncol(scp_matrix))) {
  a <- as.numeric(scp_matrix[, i])
  b <- complete.cases(a, y)
  correlation <- cor.test(a[b], y[b], method = "spearman", exact = FALSE)
  cor_values[i] <- unname(correlation$estimate)
  p_values[i]   <- correlation$p.value
}

cor_results <- data.frame(
  Gene = colnames(scp_matrix),
  estimate = cor_values,
  p_value = p_values,
  fdr = p.adjust(p_values, "BH"),
  stringsAsFactors = FALSE
)[order(p.adjust(p_values, "BH"), -abs(cor_values)), ]

scp_mes_sig <- cor_results[cor_results$fdr < 0.01, ]


#### Schwann Cell vs. NO_TMM GSVA.
source("DataCleaning.R")

gene_list_noTMM <- list("NO_TMM Genes" = c("CPNE8", "PGM2L1", "LIFR", "CNR1", "HECW2", "CPNE3", "HOXC9", "SNX16", "IGSF10"))

scp_genes <- c("CDH19", "PLP1", "ERBB3", "MPZ", "ERBB4")

# GSVA for NO_TMM marker.
Expression <- as.matrix(Expression)
noTMM_gsva <- gsvaParam(Expression, gene_list_noTMM, kcdf = "Gaussian")
GSVA_result_noTMM <- gsva(noTMM_gsva)
GSVA_df_noTMM <- as.data.frame(GSVA_result_noTMM)
GSVA_long_noTMM <- pivot_longer(GSVA_df_noTMM, cols = everything(), names_to = "SampleID", values_to = "GSVA_Score")


# ranks of schwann cell marker gene.
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

scp_genes <- intersect(scp_genes, rownames(ranked_TARGET_NBL))
scp_matrix <- t(as.matrix(ranked_TARGET_NBL[scp_genes, , drop = FALSE]))


GSVA_long_noTMM <- GSVA_long_noTMM[match(rownames(scp_matrix), GSVA_long_noTMM$SampleID), ]

y <- GSVA_long_noTMM$GSVA_Score

cor_values <- numeric(ncol(scp_matrix))
p_values   <- numeric(ncol(scp_matrix))



for (i in seq_len(ncol(scp_matrix))) {
  a <- as.numeric(scp_matrix[, i])
  b <- complete.cases(a, y)
  correlation <- cor.test(a[b], y[b], method = "spearman", exact = FALSE)
  cor_values[i] <- unname(correlation$estimate)
  p_values[i]   <- correlation$p.value
}

cor_results <- data.frame(
  Gene = colnames(scp_matrix),
  estimate = cor_values,
  p_value = p_values,
  fdr = p.adjust(p_values, "BH"),
  stringsAsFactors = FALSE
)[order(p.adjust(p_values, "BH"), -abs(cor_values)), ]

scp_noTMM_sig <- cor_results[cor_results$fdr < 0.01, ]


### Correlation Late Neuroblast with NO_TMM GSVA.
source("DataCleaning.R")

gene_list_noTMM <- list("NO_TMM Genes" = c("CPNE8", "PGM2L1", "LIFR", "CNR1", "HECW2", "CPNE3", "HOXC9", "SNX16", "IGSF10"))

lateNeuroblast_genes <- c("IL7", "GAP43", "STMN2", "SYN3", "ISL1")

# GSVA for NO_TMM marker.
Expression <- as.matrix(Expression)
noTMM_gsva <- gsvaParam(Expression, gene_list_noTMM, kcdf = "Gaussian")
GSVA_result_noTMM <- gsva(noTMM_gsva)
GSVA_df_noTMM <- as.data.frame(GSVA_result_noTMM)
GSVA_long_noTMM <- pivot_longer(GSVA_df_noTMM, cols = everything(), names_to = "SampleID", values_to = "GSVA_Score")


# ranks of late neuroblast gene.
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

lateNeuroblast_genes <- intersect(lateNeuroblast_genes, rownames(ranked_TARGET_NBL))
lateNeuroblast_matrix <- t(as.matrix(ranked_TARGET_NBL[lateNeuroblast_genes, , drop = FALSE]))


GSVA_long_noTMM <- GSVA_long_noTMM[match(rownames(lateNeuroblast_matrix), GSVA_long_noTMM$SampleID), ]

y <- GSVA_long_noTMM$GSVA_Score

cor_values <- numeric(ncol(lateNeuroblast_matrix))
p_values   <- numeric(ncol(lateNeuroblast_matrix))



for (i in seq_len(ncol(lateNeuroblast_matrix))) {
  a <- as.numeric(lateNeuroblast_matrix[, i])
  b <- complete.cases(a, y)
  correlation <- cor.test(a[b], y[b], method = "spearman", exact = FALSE)
  cor_values[i] <- unname(correlation$estimate)
  p_values[i]   <- correlation$p.value
}

cor_results <- data.frame(
  Gene = colnames(lateNeuroblast_matrix),
  estimate = cor_values,
  p_value = p_values,
  fdr = p.adjust(p_values, "BH"),
  stringsAsFactors = FALSE
)[order(p.adjust(p_values, "BH"), -abs(cor_values)), ]

lateNeuroblast_noTMM_sig <- cor_results[cor_results$fdr < 0.01, ]


### Correlation Late Neuroblast rank with extend SCORE.
source("DataCleaning.R")

lateNeuroblast_genes <- c("IL7", "GAP43", "STMN2", "SYN3", "ISL1")

# extend.
extendScores <- RunEXTEND(as.matrix(Expression))

telomeraseScores <- read_delim("TelomeraseScores.txt")

telomeraseScores <- telomeraseScores[, c("SampleID", "NormEXTENDScores")]
telomeraseScores <- as.data.frame(telomeraseScores)


# ranks of lateNeuroblast gene.
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

lateNeuroblast_genes <- intersect(lateNeuroblast_genes, rownames(ranked_TARGET_NBL))
lateNeuroblast_matrix <- t(as.matrix(ranked_TARGET_NBL[lateNeuroblast_genes, , drop = FALSE]))


telomeraseScores <- telomeraseScores[match(rownames(lateNeuroblast_matrix), telomeraseScores$SampleID), ]

y <- telomeraseScores$NormEXTENDScores

cor_values <- numeric(ncol(lateNeuroblast_matrix))
p_values   <- numeric(ncol(lateNeuroblast_matrix))



for (i in seq_len(ncol(lateNeuroblast_matrix))) {
  a <- as.numeric(lateNeuroblast_matrix[, i])
  b <- complete.cases(a, y)
  correlation <- cor.test(a[b], y[b], method = "spearman", exact = FALSE)
  cor_values[i] <- unname(correlation$estimate)
  p_values[i]   <- correlation$p.value
}

cor_results <- data.frame(
  Gene = colnames(lateNeuroblast_matrix),
  estimate = cor_values,
  p_value = p_values,
  fdr = p.adjust(p_values, "BH"),
  stringsAsFactors = FALSE
)[order(p.adjust(p_values, "BH"), -abs(cor_values)), ]

lateNeuroblast_extend_sig <- cor_results[cor_results$fdr < 0.01, ]

### Correlation Late Neuroblast with MES GSVA. :))))
lateNeuroblast_genes <- c("IL7", "GAP43", "STMN2", "SYN3", "ISL1")
mesenchymal_markers <- list(MES = c("PRRX1", "RUNX1", "FOSL1", "JUN", "YAP1", "WWTR1", "MEOX1", "MEOX2",
                                    "SNAI2", "CD44", "FN1", "VIM"))

# GSVA for mesenchymal marker.
Expression <- as.matrix(Expression)
mes_gsva <- gsvaParam(Expression, mesenchymal_markers, kcdf = "Gaussian")
GSVA_result_mes <- gsva(mes_gsva)
GSVA_df_mes <- as.data.frame(GSVA_result_mes)
GSVA_long_mes <- pivot_longer(GSVA_df_mes, cols = everything(), names_to = "SampleID", values_to = "GSVA_Score")



y <- GSVA_long_mes$GSVA_Score

# ranks of late Neuroblast gene.
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

lateNeuroblast_genes <- intersect(lateNeuroblast_genes, rownames(ranked_TARGET_NBL))
lateNeuroblast_matrix <- t(as.matrix(ranked_TARGET_NBL[lateNeuroblast_genes, , drop = FALSE]))
lateNeuroblast_matrix <- lateNeuroblast_matrix[match(GSVA_long_mes$SampleID, rownames(lateNeuroblast_matrix)), ]


cor_values <- numeric(ncol(lateNeuroblast_matrix))
p_values   <- numeric(ncol(lateNeuroblast_matrix))


for (i in seq_len(ncol(lateNeuroblast_matrix))) {
  a <- as.numeric(lateNeuroblast_matrix[, i])
  b <- complete.cases(a, y)
  correlation <- cor.test(a[b], y[b], method = "spearman", exact = FALSE)
  cor_values[i] <- unname(correlation$estimate)
  p_values[i]   <- correlation$p.value
}

cor_results <- data.frame(
  Gene = colnames(lateNeuroblast_matrix),
  estimate = cor_values,
  p_value = p_values,
  fdr = p.adjust(p_values, "BH"),
  stringsAsFactors = FALSE
)[order(p.adjust(p_values, "BH"), -abs(cor_values)), ]

mes_lateNeuroblast_sig <- cor_results[cor_results$fdr < 0.01, ]



### Correlation Late Neuroblast with ADR GSVA.

adrenergic_markers <- list(ADR = c("PHOX2A", "PHOX2B", "HAND2", "ISL1", "GATA3", "ASCL1", "TBX2", "DBH", "TH"))
lateNeuroblast_genes <- c("IL7", "GAP43", "STMN2", "SYN3", "ISL1")

# GSVA for mesenchymal marker.
Expression <- as.matrix(Expression)
adr_gsva <- gsvaParam(Expression, adrenergic_markers, kcdf = "Gaussian")
GSVA_result_adr <- gsva(adr_gsva)
GSVA_df_adr <- as.data.frame(GSVA_result_adr)
GSVA_long_adr <- pivot_longer(GSVA_df_adr, cols = everything(), names_to = "SampleID", values_to = "GSVA_Score")

y <- GSVA_long_mes$GSVA_Score

# ranks of lasteNeuroblast gene.
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

lateNeuroblast_genes <- intersect(lateNeuroblast_genes, rownames(ranked_TARGET_NBL))
lateNeuroblast_matrix <- t(as.matrix(ranked_TARGET_NBL[lateNeuroblast_genes, , drop = FALSE]))
lateNeuroblast_matrix <- lateNeuroblast_matrix[match(GSVA_long_mes$SampleID, rownames(lateNeuroblast_matrix)), ]


cor_values <- numeric(ncol(lateNeuroblast_matrix))
p_values   <- numeric(ncol(lateNeuroblast_matrix))


for (i in seq_len(ncol(lateNeuroblast_matrix))) {
  a <- as.numeric(lateNeuroblast_matrix[, i])
  b <- complete.cases(a, y)
  correlation <- cor.test(a[b], y[b], method = "spearman", exact = FALSE)
  cor_values[i] <- unname(correlation$estimate)
  p_values[i]   <- correlation$p.value
}

cor_results <- data.frame(
  Gene = colnames(lateNeuroblast_matrix),
  estimate = cor_values,
  p_value = p_values,
  fdr = p.adjust(p_values, "BH"),
  stringsAsFactors = FALSE
)[order(p.adjust(p_values, "BH"), -abs(cor_values)), ]

adr_lateNeuroblast_sig <- cor_results[cor_results$fdr < 0.01, ]




#######


# t-test: cell cycle inhibitors GSVA for different TMM groups.
source("DataCleaning.R")
cell_cycle_inhibitor <- list(Inhibitor = c("FOXO1", "BTG1", "CDKN1B", "CDKN2A", "HES1", "EGR1", "CDKN2C", "CDKN2B", 
                                           "HOPX", "CDKN1A", "KLF4", "CDKN1C", "TGFB1", "CDKN2D", "SMAD7", "KLF10"))
Expression <- as.matrix(Expression)
gsva <- gsvaParam(Expression, cell_cycle_inhibitor, kcdf = "Gaussian")
GSVA_result <- gsva(gsva)
GSVA_df <- as.data.frame(GSVA_result)
GSVA_long <- pivot_longer(GSVA_df, cols = everything(), names_to = "SampleID", values_to = "GSVA_Score")

GSVA_long <- merge(GSVA_long, metadata[, c("SampleID", "TMM", "COG.Risk.Group")], by = "SampleID")

# boxplot: cell cycle inhibitors vs. TMM.
ggplot(GSVA_long, aes(x = TMM, y = GSVA_Score, fill = TMM, color = TMM)) +
  geom_boxplot(size = 0.2, alpha = 0.5, outlier.shape = NA) +
  geom_point(position = position_jitter(width = 0.2), size = 3) +
  scale_fill_manual(values = c("Telomerase" = "lightpink2", 
                               "NO_TMM" = "lightgreen",
                               "ALT" = "blue")) +
  scale_color_manual(values = c("Telomerase"="darkred", 
                                "NO_TMM" = "darkgreen",
                                "ALT" = "darkblue")) +
  theme_classic() +
  labs(x = "Class", y = "Cell Cycle Inhibitor GSVA Score") +
  theme(
    axis.text.x = element_text(vjust = 1, hjust = 1),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 12, face = "bold"),
    legend.position = "none"
  ) +
  stat_compare_means(comparisons = list(c("Telomerase","NO_TMM")), method= "t.test",
                     method.args = list(alternative ="two.sided"), size = 6, tip.length = 0.01,
                     label.y = 1.2) +
  stat_compare_means(comparisons = list(c("NO_TMM","ALT")), method= "t.test",
                     method.args = list(alternative ="two.sided"), size = 6, tip.length = 0.01,
                     label.y = 1.4)

# now, based on risk groups.
ggplot(GSVA_long[GSVA_long$TMM == "NO_TMM", ], aes(x = COG.Risk.Group, y = GSVA_Score, fill = COG.Risk.Group, color = COG.Risk.Group)) +
  geom_boxplot(size = 0.2, alpha = 0.5, outlier.shape = NA) +
  geom_point(position = position_jitter(width = 0.2), size = 3) +
  scale_fill_manual(values = c("High Risk" = "lightpink2", 
                               "Low Risk" = "lightgreen",
                               "Intermediate Risk" = "blue")) +
  scale_color_manual(values = c("High Risk"="darkred", 
                                "Low Risk" = "darkgreen",
                                "Intermediate Risk" = "darkblue")) +
  theme_classic() +
  labs(x = "Class", y = "Cell Cycle Inhibitor GSVA Score") +
  theme(
    axis.text.x = element_text(vjust = 1, hjust = 1),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 12, face = "bold"),
    legend.position = "none"
  ) +
  stat_compare_means(comparisons = list(c("High Risk","Low Risk")), method= "t.test",
                     method.args = list(alternative ="two.sided"), size = 6, tip.length = 0.01,
                     label.y = 1.2) +
  stat_compare_means(comparisons = list(c("Low Risk","Intermediate Risk")), method= "t.test",
                     method.args = list(alternative ="two.sided"), size = 6, tip.length = 0.01,
                     label.y = 1.4)


### t-test: scp gsva for different risk groups and TMM groups.
source("DataCleaning.R")
scp <- list(SCP = c("CDH19", "PLP1", "ERBB3", "MPZ", "ERBB4"))

Expression <- as.matrix(Expression)
gsva <- gsvaParam(Expression, scp, kcdf = "Gaussian")
GSVA_result <- gsva(gsva)
GSVA_df <- as.data.frame(GSVA_result)
GSVA_long <- pivot_longer(GSVA_df, cols = everything(), names_to = "SampleID", values_to = "GSVA_Score")

GSVA_long <- merge(GSVA_long, metadata[, c("SampleID", "TMM", "COG.Risk.Group")], by = "SampleID")

# boxplot: scp gsva vs. TMM.
ggplot(GSVA_long, aes(x = TMM, y = GSVA_Score, fill = TMM, color = TMM)) +
  geom_boxplot(size = 0.2, alpha = 0.5, outlier.shape = NA) +
  geom_point(position = position_jitter(width = 0.2), size = 3) +
  scale_fill_manual(values = c("Telomerase" = "lightpink2", 
                               "NO_TMM" = "lightgreen",
                               "ALT" = "blue")) +
  scale_color_manual(values = c("Telomerase"="darkred", 
                                "NO_TMM" = "darkgreen",
                                "ALT" = "darkblue")) +
  theme_classic() +
  labs(x = "Class", y = "SCP GSVA Score") +
  theme(
    axis.text.x = element_text(vjust = 1, hjust = 1),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 12, face = "bold"),
    legend.position = "none"
  ) +
  stat_compare_means(comparisons = list(c("Telomerase","NO_TMM")), method= "t.test",
                     method.args = list(alternative ="two.sided"), size = 6, tip.length = 0.01,
                     label.y = 1.2) +
  stat_compare_means(comparisons = list(c("NO_TMM","ALT")), method= "t.test",
                     method.args = list(alternative ="two.sided"), size = 6, tip.length = 0.01,
                     label.y = 1.4)

# now, based on risk groups.
ggplot(GSVA_long[GSVA_long$TMM == "NO_TMM", ], aes(x = COG.Risk.Group, y = GSVA_Score, fill = COG.Risk.Group, color = COG.Risk.Group)) +
  geom_boxplot(size = 0.2, alpha = 0.5, outlier.shape = NA) +
  geom_point(position = position_jitter(width = 0.2), size = 3) +
  scale_fill_manual(values = c("High Risk" = "lightpink2", 
                               "Low Risk" = "lightgreen",
                               "Intermediate Risk" = "blue")) +
  scale_color_manual(values = c("High Risk"="darkred", 
                                "Low Risk" = "darkgreen",
                                "Intermediate Risk" = "darkblue")) +
  theme_classic() +
  labs(x = "Class", y = "SCP GSVA Score") +
  theme(
    axis.text.x = element_text(vjust = 1, hjust = 1),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 12, face = "bold"),
    legend.position = "none"
  ) +
  stat_compare_means(comparisons = list(c("High Risk","Low Risk")), method= "t.test",
                     method.args = list(alternative ="two.sided"), size = 6, tip.length = 0.01,
                     label.y = 1.2) +
  stat_compare_means(comparisons = list(c("High Risk","Intermediate Risk")), method= "t.test",
                     method.args = list(alternative ="two.sided"), size = 6, tip.length = 0.01,
                     label.y = 1.4)




### t-test: chromaffin gsva for different risk groups and TMM groups.
source("DataCleaning.R")
chromaffin <- list(Chromaffin = c("TH", "DBH", "DDC", "CHGA", "PNMT"))

Expression <- as.matrix(Expression)
gsva <- gsvaParam(Expression, chromaffin, kcdf = "Gaussian")
GSVA_result <- gsva(gsva)
GSVA_df <- as.data.frame(GSVA_result)
GSVA_long <- pivot_longer(GSVA_df, cols = everything(), names_to = "SampleID", values_to = "GSVA_Score")

GSVA_long <- merge(GSVA_long, metadata[, c("SampleID", "TMM", "COG.Risk.Group")], by = "SampleID")

# boxplot: chromaffin markers GSVA vs. TMM.
ggplot(GSVA_long, aes(x = TMM, y = GSVA_Score, fill = TMM, color = TMM)) +
  geom_boxplot(size = 0.2, alpha = 0.5, outlier.shape = NA) +
  geom_point(position = position_jitter(width = 0.2), size = 3) +
  scale_fill_manual(values = c("Telomerase" = "lightpink2", 
                               "NO_TMM" = "lightgreen",
                               "ALT" = "blue")) +
  scale_color_manual(values = c("Telomerase"="darkred", 
                                "NO_TMM" = "darkgreen",
                                "ALT" = "darkblue")) +
  theme_classic() +
  labs(x = "Class", y = "Chromaffin GSVA Score") +
  theme(
    axis.text.x = element_text(vjust = 1, hjust = 1),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 12, face = "bold"),
    legend.position = "none"
  ) +
  stat_compare_means(comparisons = list(c("Telomerase","NO_TMM")), method= "t.test",
                     method.args = list(alternative ="two.sided"), size = 6, tip.length = 0.01,
                     label.y = 1.2) +
  stat_compare_means(comparisons = list(c("NO_TMM","ALT")), method= "t.test",
                     method.args = list(alternative ="two.sided"), size = 6, tip.length = 0.01,
                     label.y = 1.4)

# now, based on risk groups.
ggplot(GSVA_long, aes(x = COG.Risk.Group, y = GSVA_Score, fill = COG.Risk.Group, color = COG.Risk.Group)) +
  geom_boxplot(size = 0.2, alpha = 0.5, outlier.shape = NA) +
  geom_point(position = position_jitter(width = 0.2), size = 3) +
  scale_fill_manual(values = c("High Risk" = "lightpink2", 
                               "Low Risk" = "lightgreen",
                               "Intermediate Risk" = "blue")) +
  scale_color_manual(values = c("High Risk"="darkred", 
                                "Low Risk" = "darkgreen",
                                "Intermediate Risk" = "darkblue")) +
  theme_classic() +
  labs(x = "Class", y = "Chromaffin GSVA Score") +
  theme(
    axis.text.x = element_text(vjust = 1, hjust = 1),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 12, face = "bold"),
    legend.position = "none"
  ) +
  stat_compare_means(comparisons = list(c("High Risk","Low Risk")), method= "t.test",
                     method.args = list(alternative ="two.sided"), size = 6, tip.length = 0.01,
                     label.y = 1.2) +
  stat_compare_means(comparisons = list(c("High Risk","Intermediate Risk")), method= "t.test",
                     method.args = list(alternative ="two.sided"), size = 6, tip.length = 0.01,
                     label.y = 1.4)


#############>>>>>>###############################




