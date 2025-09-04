
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







