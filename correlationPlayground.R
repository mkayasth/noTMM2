library(tidyverse)
library(dplyr)
library(ggpubr)
library(survival)
library(survminer)
library(rstatix)
library(GSVA)


### TERT association with ADR and MES genes from Thirant paper.

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

# transpose so samples are rows and genes are columns.
t_NBL <- t(ranked_TARGET_NBL)
t_NBL <- as.data.frame(t_NBL)

# TERT expression.
tert_expr <- t_NBL$TERT

# adrenergic and mesenchymal markers from Thirant paper.
mesenchymal_markers <- c("FN1", "VIM", "SNAI2", "PRRX1", "YAP1", "WWTR1")

adrenergic_markers <- c("DLK1", "DBH", "PHOX2A", "PHOX2B", "GATA2", "GATA3")

# first looking at correlation of TERT with adrenergic markers.
adrenergic_ranks <- t_NBL[, colnames(t_NBL) %in% adrenergic_markers, drop = FALSE]
adrenergic_ranks <- adrenergic_ranks[match(rownames(t_NBL), rownames(adrenergic_ranks)), ]

# vector to store results.
cor_values_adr <- numeric(ncol(adrenergic_ranks))
p_values_adr <- numeric(ncol(adrenergic_ranks))

# Looping through each gene for rank-based Spearman correlation with TERT.
for (i in seq_along(adrenergic_ranks)) {
  result <- cor.test(adrenergic_ranks[[i]], tert_expr, method = "spearman")
  cor_values_adr[i] <- result$estimate
  p_values_adr[i] <- result$p.value
}

# result dataframe.
cor_results_adr <- data.frame(
  Gene = colnames(adrenergic_ranks),
  correlation = cor_values_adr,
  pvalue = p_values_adr
)

# p-values for calculating FDR.
cor_results_adr$FDR <- p.adjust(cor_results_adr$pvalue, method = "fdr")
cor_results_adr <- cor_results_adr %>%
  arrange(FDR)

# filtering fdr < 0.01.
tert_correlated_genes_adr <- cor_results_adr %>%
  filter(FDR <= 0.01) %>%
  arrange(desc(abs(correlation)))

#####

# now looking at correlation of TERT with mesenchymal markers.
mesenchymal_ranks <- t_NBL[, colnames(t_NBL) %in% mesenchymal_markers, drop = FALSE]
mesenchymal_ranks <- mesenchymal_ranks[match(rownames(t_NBL), rownames(mesenchymal_ranks)), ]

# vector to store results.
cor_values_mes <- numeric(ncol(mesenchymal_ranks))
p_values_mes <- numeric(ncol(mesenchymal_ranks))

# Looping through each gene for rank-based Spearman correlation with TERT.
for (i in seq_along(mesenchymal_ranks)) {
  result <- cor.test(mesenchymal_ranks[[i]], tert_expr, method = "spearman")
  cor_values_mes[i] <- result$estimate
  p_values_mes[i] <- result$p.value
}

# result dataframe.
cor_results_mes <- data.frame(
  Gene = colnames(mesenchymal_ranks),
  correlation = cor_values_mes,
  pvalue = p_values_mes
)

# p-values for calculating FDR.
cor_results_mes$FDR <- p.adjust(cor_results_mes$pvalue, method = "fdr")
cor_results_mes <-  cor_results_mes %>%
  arrange(FDR)

# filtering fdr < 0.01.
tert_correlated_genes_mes <- cor_results_mes %>%
  filter(FDR <= 0.01) %>%
  arrange(desc(abs(correlation)))

##################################################################################

# seeing correlation between EXTEND scores and MES markers rank.

source("EXTEND/ComponentAndMarkerFunction.r")
source("EXTEND/ComponentOneAndMarkerFunction.r")
source("EXTEND/ComponentTwoAndMarkerFunction.r")
source("EXTEND/InputData.r")
source("EXTEND/IterativeRS.r")
source("EXTEND/IterativeRS.r")
source("EXTEND/MarkerFunction.r")
source("EXTEND/RunEXTEND.r")

source("dataCleaning.R")

Expression_mes <- Expression[rownames(Expression) %in% mesenchymal_markers, ,drop = FALSE]
extendScores <- RunEXTEND(as.matrix(Expression))

telomeraseScores <- read_delim("TelomeraseScores.txt")

telomeraseScores <- telomeraseScores[, c("SampleID", "NormEXTENDScores")]
telomeraseScores <- as.data.frame(telomeraseScores)
telomeraseScores <- left_join(telomeraseScores, metadata[, c("SampleID", "TMM")], by = "SampleID")
rownames(telomeraseScores) <- telomeraseScores$SampleID
telomeraseScores$SampleID = NULL

telomeraseScores <- telomeraseScores[rownames(mesenchymal_ranks), ]

# vector to store results.
cor_values <- numeric(ncol(mesenchymal_ranks))
p_values <- numeric(ncol(mesenchymal_ranks))

# Looping through each gene for rank-based Spearman correlation with TERT.
for (i in seq_along(mesenchymal_ranks)) {
  result <- cor.test(mesenchymal_ranks[[i]], telomeraseScores$NormEXTENDScores, method = "spearman")
  cor_values[i] <- result$estimate
  p_values[i] <- result$p.value
}

# result dataframe.
cor_results <- data.frame(
  Gene = colnames(mesenchymal_ranks),
  correlation = cor_values,
  pvalue = p_values
)

# p-values for calculating FDR.
cor_results$FDR <- p.adjust(cor_results$pvalue, method = "fdr")
cor_results <- cor_results %>%
  arrange(FDR)

# filtering fdr < 0.01.
extend_correlated_genes <- cor_results %>%
  filter(FDR <= 0.01) %>%
  arrange(correlation)

## So, Telomerase is negatively associated with mesenchymal characteristics. what is up with heatmap? ::)


###########################################################################################################################################
# correlation between GSVA of NO_TMM signatures and mesenchymal/adrenergic genes.

# source("aucGSVA.R")
kfold_genes <- list(NO_TMM = c("CPNE8", "PGM2L1", "CNR1", "LIFR", "HOXC9", "SNX16", "HECW2", "ALDH3A2",
                               "THSD7A", "CPNE3", "IGSF10"))
kfold_genes_ <- list(NO_TMM = c("PRR7", "IGLV6-57", "SAC3D1", "DDN", "ZNHIT2", "CCDC86", "TCF15", "HTR6")) # these genes are downregulated in NO_TMM.

# gsvapar <- gsvaParam(as.matrix(Expression), kfold_genes, kcdf = "Gaussian") # unlocking this as necessary.
gsvapar <- gsvaParam(as.matrix(Expression), kfold_genes_, kcdf = "Gaussian")
gsva_result <- gsva(gsvapar)
gsva_result <- as.data.frame(gsva_result)
rownames(gsva_result) <- NULL


gsva_long <- gsva_result %>%
  pivot_longer(cols = everything(),
               names_to = "SampleID",
               values_to = "GSVA_Score"
  )


gsva_long <- left_join(gsva_long, metadata[, c("SampleID", "TMM")], by = "SampleID")

# first looking at correlation of NO_TMM genes with adrenergic markers.

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

# transpose so samples are rows and genes are columns.
t_NBL <- t(ranked_TARGET_NBL)
t_NBL <- as.data.frame(t_NBL)

adrenergic_ranks <- t_NBL[, colnames(t_NBL) %in% adrenergic_markers, drop = FALSE]
gsva_long <- gsva_long[match(rownames(adrenergic_ranks), gsva_long$SampleID), ]
gsva_expr <- gsva_long$GSVA_Score

# vector to store results.
cor_values_adr <- numeric(ncol(adrenergic_ranks))
p_values_adr <- numeric(ncol(adrenergic_ranks))

# Looping through each gene for rank-based Spearman correlation with TERT.
for (i in seq_along(adrenergic_ranks)) {
  result <- cor.test(adrenergic_ranks[[i]], gsva_expr, method = "spearman")
  cor_values_adr[i] <- result$estimate
  p_values_adr[i] <- result$p.value
}

# result dataframe.
cor_results_adr <- data.frame(
  Gene = colnames(adrenergic_ranks),
  correlation = cor_values_adr,
  pvalue = p_values_adr
)

# p-values for calculating FDR.
cor_results_adr$FDR <- p.adjust(cor_results_adr$pvalue, method = "fdr")
cor_results_adr <- cor_results_adr %>%
  arrange(FDR)

# filtering fdr < 0.01.
notmm_correlated_genes_adr <- cor_results_adr %>%
  filter(FDR <= 0.01) %>%
  arrange(desc(abs(correlation)))

###

# now, looking at correlation of NO_TMM genes with mesenchymal markers.
mesenchymal_ranks <- t_NBL[, colnames(t_NBL) %in% mesenchymal_markers, drop = FALSE]
gsva_long <- gsva_long[match(rownames(mesenchymal_ranks), gsva_long$SampleID), ]
gsva_expr <- gsva_long$GSVA_Score

# vector to store results.
cor_values_mes <- numeric(ncol(mesenchymal_ranks))
p_values_mes <- numeric(ncol(mesenchymal_ranks))

# Looping through each gene for rank-based Spearman correlation with TERT.
for (i in seq_along(mesenchymal_ranks)) {
  result <- cor.test(mesenchymal_ranks[[i]], gsva_expr, method = "spearman")
  cor_values_mes[i] <- result$estimate
  p_values_mes[i] <- result$p.value
}

# result dataframe.
cor_results_mes <- data.frame(
  Gene = colnames(mesenchymal_ranks),
  correlation = cor_values_mes,
  pvalue = p_values_mes
)

# p-values for calculating FDR.
cor_results_mes$FDR <- p.adjust(cor_results_mes$pvalue, method = "fdr")
cor_results_mes <- cor_results_mes %>%
  arrange(FDR)

# filtering fdr < 0.01.
notmm_correlated_genes_mes <- cor_results_mes %>%
  filter(FDR <= 0.01) %>%
  arrange(desc(correlation))

########################################################################################

# mes markers and extend scores have negative relationship. Looking if ADR markers have positive relationship with extend scores.

# seeing correlation between EXTEND scores and ADR markers rank.

source("EXTEND/ComponentAndMarkerFunction.r")
source("EXTEND/ComponentOneAndMarkerFunction.r")
source("EXTEND/ComponentTwoAndMarkerFunction.r")
source("EXTEND/InputData.r")
source("EXTEND/IterativeRS.r")
source("EXTEND/IterativeRS.r")
source("EXTEND/MarkerFunction.r")
source("EXTEND/RunEXTEND.r")

ExpressionFpkm <- ExpressionFpkm[, match(metadata$SampleID, colnames(ExpressionFpkm))]

# ranking the genes in ExpressionFpkm.
ranked_TARGET_NBL <- apply(ExpressionFpkm, 2, function(x) rank(x, ties.method = "average"))

# transpose so samples are rows and genes are columns.
t_NBL <- t(ranked_TARGET_NBL)
t_NBL <- as.data.frame(t_NBL)


# looking at correlation of extend score with adrenergic markers -- expecting positive correlation.
adrenergic_ranks <- t_NBL[, colnames(t_NBL) %in% adrenergic_markers, drop = FALSE]
adrenergic_ranks <- adrenergic_ranks[match(rownames(t_NBL), rownames(adrenergic_ranks)), ]

Expression_adr <- Expression[rownames(Expression) %in% adrenergic_markers, ,drop = FALSE]
extendScores <- RunEXTEND(as.matrix(Expression))

telomeraseScores <- read_delim("TelomeraseScores.txt")

telomeraseScores <- telomeraseScores[, c("SampleID", "NormEXTENDScores")]
telomeraseScores <- as.data.frame(telomeraseScores)
telomeraseScores <- left_join(telomeraseScores, metadata[, c("SampleID", "TMM")], by = "SampleID")
rownames(telomeraseScores) <- telomeraseScores$SampleID
telomeraseScores$SampleID = NULL

telomeraseScores <- telomeraseScores[rownames(adrenergic_ranks), ]

# vector to store results.
cor_values <- numeric(ncol(adrenergic_ranks))
p_values <- numeric(ncol(adrenergic_ranks))

# Looping through each gene for rank-based Spearman correlation with extend.
for (i in seq_along(adrenergic_ranks)) {
  result <- cor.test(adrenergic_ranks[[i]], telomeraseScores$NormEXTENDScores, method = "spearman")
  cor_values[i] <- result$estimate
  p_values[i] <- result$p.value
}

# result dataframe.
cor_results <- data.frame(
  Gene = colnames(adrenergic_ranks),
  correlation = cor_values,
  pvalue = p_values
)

# p-values for calculating FDR.
cor_results$FDR <- p.adjust(cor_results$pvalue, method = "fdr")
cor_results <- cor_results %>%
  arrange(FDR)

# filtering fdr < 0.01.
extend_correlated_genes <- cor_results %>%
  filter(FDR <= 0.01) %>%
  arrange(desc(abs(correlation)))


###########################################################################################################
###########################################################################################################

## Looking if the same trend is true for the Ackerman data.

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

# ranking the genes.
ranked_ackerman_NB <- apply(ackerman_NB, 2, function(x) rank(x, ties.method = "average"))

# transpose so samples are rows and genes are columns.
t_NBL <- t(ranked_ackerman_NB)
t_NBL <- as.data.frame(t_NBL)


# now looking at mesenchymal marker ranks.
mesenchymal_ranks <- t_NBL[, colnames(t_NBL) %in% mesenchymal_markers, drop = FALSE]
mesenchymal_ranks <- mesenchymal_ranks[match(rownames(t_NBL), rownames(mesenchymal_ranks)), ]


extendScores <- RunEXTEND(as.matrix(ackerman_NB))

telomeraseScores <- read_delim("TelomeraseScores.txt")

telomeraseScores <- telomeraseScores[, c("SampleID", "NormEXTENDScores")]
telomeraseScores <- as.data.frame(telomeraseScores)
telomeraseScores <- left_join(telomeraseScores, ackerman_metadata[, c("SampleID", "TMM_Case")], by = "SampleID")
telomeraseScores <- telomeraseScores[!duplicated(telomeraseScores), ]
rownames(telomeraseScores) <- telomeraseScores$SampleID

telomeraseScores$SampleID = NULL

telomeraseScores <- telomeraseScores[rownames(mesenchymal_ranks), ]

# vector to store results.
cor_values <- numeric(ncol(mesenchymal_ranks))
p_values <- numeric(ncol(mesenchymal_ranks))

# Looping through each gene for rank-based Spearman correlation with TERT.
for (i in seq_along(mesenchymal_ranks)) {
  result <- cor.test(mesenchymal_ranks[[i]], telomeraseScores$NormEXTENDScores, method = "spearman")
  cor_values[i] <- result$estimate
  p_values[i] <- result$p.value
}

# result dataframe.
cor_results <- data.frame(
  Gene = colnames(mesenchymal_ranks),
  correlation = cor_values,
  pvalue = p_values
)

# p-values for calculating FDR.
cor_results$FDR <- p.adjust(cor_results$pvalue, method = "fdr")
cor_results <- cor_results %>%
  arrange(FDR)

# filtering fdr < 0.01.
extend_correlated_genes_mes <- cor_results %>%
  filter(FDR <= 0.01) %>%
  arrange(desc(abs(correlation)))

###

# Now, looking at EXTEND relationship with ADR.

# adrenergic marker ranks.
adrenergic_ranks <- t_NBL[, colnames(t_NBL) %in% adrenergic_markers, drop = FALSE]
adrenergic_ranks <- adrenergic_ranks[match(rownames(t_NBL), rownames(adrenergic_ranks)), ]


telomeraseScores <- read_delim("TelomeraseScores.txt")

telomeraseScores <- telomeraseScores[, c("SampleID", "NormEXTENDScores")]
telomeraseScores <- as.data.frame(telomeraseScores)
telomeraseScores <- left_join(telomeraseScores, ackerman_metadata[, c("SampleID", "TMM_Case")], by = "SampleID")
telomeraseScores <- telomeraseScores[!duplicated(telomeraseScores), ]
rownames(telomeraseScores) <- telomeraseScores$SampleID

telomeraseScores$SampleID = NULL

telomeraseScores <- telomeraseScores[rownames(adrenergic_ranks), ]

# vector to store results.
cor_values <- numeric(ncol(adrenergic_ranks))
p_values <- numeric(ncol(adrenergic_ranks))

# Looping through each gene for rank-based Spearman correlation with TERT.
for (i in seq_along(adrenergic_ranks)) {
  result <- cor.test(adrenergic_ranks[[i]], telomeraseScores$NormEXTENDScores, method = "spearman")
  cor_values[i] <- result$estimate
  p_values[i] <- result$p.value
}

# result dataframe.
cor_results <- data.frame(
  Gene = colnames(adrenergic_ranks),
  correlation = cor_values,
  pvalue = p_values
)

# p-values for calculating FDR.
cor_results$FDR <- p.adjust(cor_results$pvalue, method = "fdr")
cor_results <- cor_results %>%
  arrange(FDR)

# filtering fdr < 0.01.
extend_correlated_genes_adr <- cor_results %>%
  filter(FDR <= 0.01) %>%
  arrange(desc(abs(correlation)))




#########

# Looking for relationship of MES and ADR markers with TERT in Ackerman dataset.
# first looking at correlation of TERT with adrenergic markers.
t_NBL <- t(ranked_ackerman_NB)
t_NBL <- as.data.frame(t_NBL)

# TERT expression.
tert_expr <- t_NBL$TERT


adrenergic_ranks <- t_NBL[, colnames(t_NBL) %in% adrenergic_markers, drop = FALSE]
adrenergic_ranks <- adrenergic_ranks[match(rownames(t_NBL), rownames(adrenergic_ranks)), ]

# vector to store results.
cor_values_adr <- numeric(ncol(adrenergic_ranks))
p_values_adr <- numeric(ncol(adrenergic_ranks))

# Looping through each gene for rank-based Spearman correlation with TERT.
for (i in seq_along(adrenergic_ranks)) {
  result <- cor.test(adrenergic_ranks[[i]], tert_expr, method = "spearman")
  cor_values_adr[i] <- result$estimate
  p_values_adr[i] <- result$p.value
}

# result dataframe.
cor_results_adr <- data.frame(
  Gene = colnames(adrenergic_ranks),
  correlation = cor_values_adr,
  pvalue = p_values_adr
)

# p-values for calculating FDR.
cor_results_adr$FDR <- p.adjust(cor_results_adr$pvalue, method = "fdr")
cor_results_adr <- cor_results_adr %>%
  arrange(FDR)

# filtering fdr < 0.01.
tert_correlated_genes_adr <- cor_results_adr %>%
  filter(FDR <= 0.01) %>%
  arrange(desc(abs(correlation)))

#####

# now looking at correlation of TERT with mesenchymal markers.
mesenchymal_ranks <- t_NBL[, colnames(t_NBL) %in% mesenchymal_markers, drop = FALSE]
mesenchymal_ranks <- mesenchymal_ranks[match(rownames(t_NBL), rownames(mesenchymal_ranks)), ]

# vector to store results.
cor_values_mes <- numeric(ncol(mesenchymal_ranks))
p_values_mes <- numeric(ncol(mesenchymal_ranks))

# Looping through each gene for rank-based Spearman correlation with TERT.
for (i in seq_along(mesenchymal_ranks)) {
  result <- cor.test(mesenchymal_ranks[[i]], tert_expr, method = "spearman")
  cor_values_mes[i] <- result$estimate
  p_values_mes[i] <- result$p.value
}

# result dataframe.
cor_results_mes <- data.frame(
  Gene = colnames(mesenchymal_ranks),
  correlation = cor_values_mes,
  pvalue = p_values_mes
)

# p-values for calculating FDR.
cor_results_mes$FDR <- p.adjust(cor_results_mes$pvalue, method = "fdr")
cor_results_mes <-  cor_results_mes %>%
  arrange(FDR)

# filtering fdr < 0.01.
tert_correlated_genes_mes <- cor_results_mes %>%
  filter(FDR <= 0.01) %>%
  arrange(correlation)


##################################

########## NO_TMM relationship with MES/ADR phenotype.
kfold_genes <- list(NO_TMM = c("CPNE8", "PGM2L1", "CNR1", "LIFR", "HOXC9", "SNX16", "HECW2", "ALDH3A2",
                               "THSD7A", "CPNE3", "IGSF10"))
kfold_genes <- list(NO_TMM = c("PRR7", "IGLV6-57", "SAC3D1", "DDN", "ZNHIT2", "CCDC86", "TCF15", "HTR6")) #NO_TMM downregulated genes.

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

gsvapar <- gsvaParam(as.matrix(ackerman_NB), kfold_genes, kcdf = "Gaussian")
gsva_result <- gsva(gsvapar)
gsva_result <- as.data.frame(gsva_result)
rownames(gsva_result) <- NULL


gsva_long <- gsva_result %>%
  pivot_longer(cols = everything(),
               names_to = "SampleID",
               values_to = "GSVA_Score"
  )


gsva_long <- left_join(gsva_long, ackerman_metadata[, c("SampleID", "TMM_Category")], by = "SampleID")

# first looking at correlation of NO_TMM genes with adrenergic markers.
# ranking the genes.
ranked_ackerman_NB <- apply(ackerman_NB, 2, function(x) rank(x, ties.method = "average"))

# transpose so samples are rows and genes are columns.
t_NBL <- t(ranked_ackerman_NB)
t_NBL <- as.data.frame(t_NBL)


adrenergic_ranks <- t_NBL[, colnames(t_NBL) %in% adrenergic_markers, drop = FALSE]
adrenergic_ranks <- adrenergic_ranks[match(rownames(t_NBL), rownames(adrenergic_ranks)), ]

gsva_long <- gsva_long[match(rownames(adrenergic_ranks), gsva_long$SampleID), ]
gsva_expr <- gsva_long$GSVA_Score

# vector to store results.
cor_values_adr <- numeric(ncol(adrenergic_ranks))
p_values_adr <- numeric(ncol(adrenergic_ranks))

# Looping through each gene for rank-based Spearman correlation with TERT.
for (i in seq_along(adrenergic_ranks)) {
  result <- cor.test(adrenergic_ranks[[i]], gsva_expr, method = "spearman")
  cor_values_adr[i] <- result$estimate
  p_values_adr[i] <- result$p.value
}

# result dataframe.
cor_results_adr <- data.frame(
  Gene = colnames(adrenergic_ranks),
  correlation = cor_values_adr,
  pvalue = p_values_adr
)

# p-values for calculating FDR.
cor_results_adr$FDR <- p.adjust(cor_results_adr$pvalue, method = "fdr")
cor_results_adr <- cor_results_adr %>%
  arrange(FDR)

# filtering fdr < 0.01.
notmm_correlated_genes_adr <- cor_results_adr %>%
  filter(FDR <= 0.01) %>%
  arrange(desc(abs(correlation)))

#####

# now, looking at correlation of NO_TMM genes with mesenchymal markers.
mesenchymal_ranks <- t_NBL[, colnames(t_NBL) %in% mesenchymal_markers, drop = FALSE]
gsva_long <- gsva_long[match(rownames(mesenchymal_ranks), gsva_long$SampleID), ]
gsva_expr <- gsva_long$GSVA_Score

# vector to store results.
cor_values_mes <- numeric(ncol(mesenchymal_ranks))
p_values_mes <- numeric(ncol(mesenchymal_ranks))

# Looping through each gene for rank-based Spearman correlation with TERT.
for (i in seq_along(mesenchymal_ranks)) {
  result <- cor.test(mesenchymal_ranks[[i]], gsva_expr, method = "spearman")
  cor_values_mes[i] <- result$estimate
  p_values_mes[i] <- result$p.value
}

# result dataframe.
cor_results_mes <- data.frame(
  Gene = colnames(mesenchymal_ranks),
  correlation = cor_values_mes,
  pvalue = p_values_mes
)

# p-values for calculating FDR.
cor_results_mes$FDR <- p.adjust(cor_results_mes$pvalue, method = "fdr")
cor_results_mes <- cor_results_mes %>%
  arrange(FDR)

# filtering fdr < 0.01.
notmm_correlated_genes_mes <- cor_results_mes %>%
  filter(FDR <= 0.01) %>%
  arrange(desc(abs(correlation)))



###################################>>>>####################################################
###################################>

##### Now, correlation studies.

# GSVA NO_TMM vs. Proliferation Rank.
source("DataCleaning.R")

gene_list_noTMM <- list("NO_TMM Genes" = c("CPNE8", "PGM2L1", "LIFR", "CNR1", "HECW2", "CPNE3", 
                                           "HOXC9", "SNX16", "IGSF10", ))
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

###? cell cycle inhibitor vs. NO_TMM GSVA.

cellCycle_inhibitor <- c(
  "CDKN1A","CDKN1B","CDKN1C",        # CIP/KIP: p21, p27, p57
  "CDKN2A","CDKN2B","CDKN2C","CDKN2D", # INK4: p16, p15, p18, p19
  "RB1","RBL1","RBL2",               # RB family
  "TP53","MDM2","MDM4","TP73","TP63", # p53 pathway regulators
  "GADD45A","GADD45B","GADD45G",     # DNA damage-inducible growth arrest
  "BTG1","BTG2","BTG3","BTG4",       # Anti-proliferative BTG/Tob family
  "DUSP1","DUSP4","DUSP5","DUSP6",   # MAPK phosphatase inhibitors
  "PPP1R15A","PPP1R15B",             # Stress-induced cell cycle arrest
  "CDKN3",                           # Cyclin-dependent kinase inhibitor 3
  "TP53BP1","TP53BP2",               # p53 binding proteins
  "CHEK1","CHEK2",                   # Checkpoint kinases
  "ATM","ATR",                       # DDR kinases enforcing arrest
  "RBL1","RBL2",                     # p107, p130
  "ING1","ING2","ING3","ING4","ING5", # Inhibitor of Growth family
  "ZBTB17","ZNF385A",                # p53 pathway repressors
  "NFATC1","NFATC2","NFATC3","NFATC4", # Negative regulators of proliferation
  "KLF4","KLF6","KLF10","KLF11",     # Krüppel-like factors with antiproliferative roles
  "EGR1","EGR2","EGR3",              # Early growth response, growth suppression
  "CD82","CD9","CD81",               # Cell cycle suppressive tetraspanins
  "MIR34A","MIR34B","MIR34C"         # miR-34 family (p53-induced inhibitors)
)


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

adrenergic_markers <- list(ADR = c("DLK1", "DBH", "PHOX2A", "PHOX2B", "GATA2", "GATA3"))

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

gene_list_noTMM <- list("NO_TMM Genes" = c("CPNE8", "PGM2L1", "CNR1", "LIFR", "HOXC9", "SNX16", "HECW2", "ALDH3A2",
                                           "THSD7A", "CPNE3", "IGSF10"))
chromaffin_genes <- c("DBH", "TH", "CHGA", "DDC", "PNMT") # pnmt is a late chromaffin feature.



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

gene_list_mesenchymal <- list(MES = c("FN1", "VIM", "SNAI2", "PRRX1", "YAP1", "WWTR1"))
scp_genes <- c("CDH19", "PLP1", "ERBB3", "MPZ", "ERBB4") #erbb4 not expressed in late scps.

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

gene_list_noTMM <- list("NO_TMM Genes" = c("CPNE8", "PGM2L1", "CNR1", "LIFR", "HOXC9", "SNX16", "HECW2", "ALDH3A2",
                                           "THSD7A", "CPNE3", "IGSF10"))

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

# ### Correlation Late Neuroblast with MES GSVA. :))))
# lateNeuroblast_genes <- c("IL7", "GAP43", "STMN2", "SYN3", "ISL1")
# mesenchymal_markers <- list(MES = c("PRRX1", "RUNX1", "FOSL1", "JUN", "YAP1", "WWTR1", "MEOX1", "MEOX2",
#                                     "SNAI2", "CD44", "FN1", "VIM"))
# 
# # GSVA for mesenchymal marker.
# Expression <- as.matrix(Expression)
# mes_gsva <- gsvaParam(Expression, mesenchymal_markers, kcdf = "Gaussian")
# GSVA_result_mes <- gsva(mes_gsva)
# GSVA_df_mes <- as.data.frame(GSVA_result_mes)
# GSVA_long_mes <- pivot_longer(GSVA_df_mes, cols = everything(), names_to = "SampleID", values_to = "GSVA_Score")
# 
# 
# 
# y <- GSVA_long_mes$GSVA_Score
# 
# # ranks of late Neuroblast gene.
# ExpressionFpkm <- readRDS("TARGET_GEdata_062024.RDS")
# 
# 
# # Data filtering: only including sample IDs from metadata present in the fpkm dataset & vice-versa.
# metadata <- metadata %>%
#   filter(SampleID %in% colnames(ExpressionFpkm))
# ExpressionFpkm <- ExpressionFpkm[, colnames(ExpressionFpkm) %in% metadata$SampleID, drop = FALSE]
# 
# metadata <- metadata %>%
#   arrange(TMM)
# 
# ExpressionFpkm <- ExpressionFpkm[, match(metadata$SampleID, colnames(ExpressionFpkm))]
# 
# # ranking the genes in ExpressionFpkm.
# ranked_TARGET_NBL <- apply(ExpressionFpkm, 2, function(x) rank(x, ties.method = "average"))
# 
# lateNeuroblast_genes <- intersect(lateNeuroblast_genes, rownames(ranked_TARGET_NBL))
# lateNeuroblast_matrix <- t(as.matrix(ranked_TARGET_NBL[lateNeuroblast_genes, , drop = FALSE]))
# lateNeuroblast_matrix <- lateNeuroblast_matrix[match(GSVA_long_mes$SampleID, rownames(lateNeuroblast_matrix)), ]
# 
# 
# cor_values <- numeric(ncol(lateNeuroblast_matrix))
# p_values   <- numeric(ncol(lateNeuroblast_matrix))
# 
# 
# for (i in seq_len(ncol(lateNeuroblast_matrix))) {
#   a <- as.numeric(lateNeuroblast_matrix[, i])
#   b <- complete.cases(a, y)
#   correlation <- cor.test(a[b], y[b], method = "spearman", exact = FALSE)
#   cor_values[i] <- unname(correlation$estimate)
#   p_values[i]   <- correlation$p.value
# }
# 
# cor_results <- data.frame(
#   Gene = colnames(lateNeuroblast_matrix),
#   estimate = cor_values,
#   p_value = p_values,
#   fdr = p.adjust(p_values, "BH"),
#   stringsAsFactors = FALSE
# )[order(p.adjust(p_values, "BH"), -abs(cor_values)), ]
# 
# mes_lateNeuroblast_sig <- cor_results[cor_results$fdr < 0.01, ]



# ### Correlation Late Neuroblast with ADR GSVA.
# 
# adrenergic_markers <- list(ADR = c("PHOX2A", "PHOX2B", "HAND2", "ISL1", "GATA3", "ASCL1", "TBX2", "DBH", "TH"))
# lateNeuroblast_genes <- c("IL7", "GAP43", "STMN2", "SYN3", "ISL1")
# 
# # GSVA for mesenchymal marker.
# Expression <- as.matrix(Expression)
# adr_gsva <- gsvaParam(Expression, adrenergic_markers, kcdf = "Gaussian")
# GSVA_result_adr <- gsva(adr_gsva)
# GSVA_df_adr <- as.data.frame(GSVA_result_adr)
# GSVA_long_adr <- pivot_longer(GSVA_df_adr, cols = everything(), names_to = "SampleID", values_to = "GSVA_Score")
# 
# y <- GSVA_long_mes$GSVA_Score
# 
# # ranks of lasteNeuroblast gene.
# ExpressionFpkm <- readRDS("TARGET_GEdata_062024.RDS")
# 
# 
# # Data filtering: only including sample IDs from metadata present in the fpkm dataset & vice-versa.
# metadata <- metadata %>%
#   filter(SampleID %in% colnames(ExpressionFpkm))
# ExpressionFpkm <- ExpressionFpkm[, colnames(ExpressionFpkm) %in% metadata$SampleID, drop = FALSE]
# 
# metadata <- metadata %>%
#   arrange(TMM)
# 
# ExpressionFpkm <- ExpressionFpkm[, match(metadata$SampleID, colnames(ExpressionFpkm))]
# 
# # ranking the genes in ExpressionFpkm.
# ranked_TARGET_NBL <- apply(ExpressionFpkm, 2, function(x) rank(x, ties.method = "average"))
# 
# lateNeuroblast_genes <- intersect(lateNeuroblast_genes, rownames(ranked_TARGET_NBL))
# lateNeuroblast_matrix <- t(as.matrix(ranked_TARGET_NBL[lateNeuroblast_genes, , drop = FALSE]))
# lateNeuroblast_matrix <- lateNeuroblast_matrix[match(GSVA_long_mes$SampleID, rownames(lateNeuroblast_matrix)), ]
# 
# 
# cor_values <- numeric(ncol(lateNeuroblast_matrix))
# p_values   <- numeric(ncol(lateNeuroblast_matrix))
# 
# 
# for (i in seq_len(ncol(lateNeuroblast_matrix))) {
#   a <- as.numeric(lateNeuroblast_matrix[, i])
#   b <- complete.cases(a, y)
#   correlation <- cor.test(a[b], y[b], method = "spearman", exact = FALSE)
#   cor_values[i] <- unname(correlation$estimate)
#   p_values[i]   <- correlation$p.value
# }
# 
# cor_results <- data.frame(
#   Gene = colnames(lateNeuroblast_matrix),
#   estimate = cor_values,
#   p_value = p_values,
#   fdr = p.adjust(p_values, "BH"),
#   stringsAsFactors = FALSE
# )[order(p.adjust(p_values, "BH"), -abs(cor_values)), ]
# 
# adr_lateNeuroblast_sig <- cor_results[cor_results$fdr < 0.01, ]

### Correlation SCP rank with extend SCORE.
source("DataCleaning.R")

scp_genes <- c("CDH19", "PLP1", "ERBB3", "ERBB4", "MPZ")

# extend.
extendScores <- RunEXTEND(as.matrix(Expression))

telomeraseScores <- read_delim("TelomeraseScores.txt")

telomeraseScores <- telomeraseScores[, c("SampleID", "NormEXTENDScores")]
telomeraseScores <- as.data.frame(telomeraseScores)


# ranks of  gene.
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


telomeraseScores <- telomeraseScores[match(rownames(scp_matrix), telomeraseScores$SampleID), ]

y <- telomeraseScores$NormEXTENDScores

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

scp_extend_sig <- cor_results[cor_results$fdr < 0.01, ]

######
### Correlation early neuroblast rank with extend SCORE.
source("DataCleaning.R")

earlyNeuroblast_genes <- c("MKI67", "TOP2A", "ALK", "ISL1", "STMN2")

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

earlyNeuroblast_genes <- intersect(earlyNeuroblast_genes, rownames(ranked_TARGET_NBL))
earlyNeuroblast_matrix <- t(as.matrix(ranked_TARGET_NBL[earlyNeuroblast_genes, , drop = FALSE]))


telomeraseScores <- telomeraseScores[match(rownames(earlyNeuroblast_matrix), telomeraseScores$SampleID), ]

y <- telomeraseScores$NormEXTENDScores

cor_values <- numeric(ncol(earlyNeuroblast_matrix))
p_values   <- numeric(ncol(earlyNeuroblast_matrix))



for (i in seq_len(ncol(earlyNeuroblast_matrix))) {
  a <- as.numeric(earlyNeuroblast_matrix[, i])
  b <- complete.cases(a, y)
  correlation <- cor.test(a[b], y[b], method = "spearman", exact = FALSE)
  cor_values[i] <- unname(correlation$estimate)
  p_values[i]   <- correlation$p.value
}

cor_results <- data.frame(
  Gene = colnames(earlyNeuroblast_matrix),
  estimate = cor_values,
  p_value = p_values,
  fdr = p.adjust(p_values, "BH"),
  stringsAsFactors = FALSE
)[order(p.adjust(p_values, "BH"), -abs(cor_values)), ]

earlyNeuroblast_extend_sig <- cor_results[cor_results$fdr < 0.01, ]

######
### Correlation chromaffin rank with extend SCORE.
source("DataCleaning.R")

chromaffin_genes <- c("DBH", "TH", "CHGA", "DDC", "PNMT")

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

chromaffin_genes <- intersect(chromaffin_genes, rownames(ranked_TARGET_NBL))
chromaffin_matrix <- t(as.matrix(ranked_TARGET_NBL[chromaffin_genes, , drop = FALSE]))


telomeraseScores <- telomeraseScores[match(rownames(chromaffin_matrix), telomeraseScores$SampleID), ]

y <- telomeraseScores$NormEXTENDScores

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

chromaffin_extend_sig <- cor_results[cor_results$fdr < 0.01, ]


### NO_TMM vs. early Neuroblast.


source("DataCleaning.R")

earlyNeuroblast_genes <- c("MKI67", "TOP2A", "ALK", "STMN2", "ISL1")

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

earlyNeuroblast_genes <- intersect(earlyNeuroblast_genes, rownames(ranked_TARGET_NBL))
earlyNeuroblast_matrix <- t(as.matrix(ranked_TARGET_NBL[earlyNeuroblast_genes, , drop = FALSE]))


GSVA_long_noTMM <- GSVA_long_noTMM[match(rownames(earlyNeuroblast_matrix), GSVA_long_noTMM$SampleID), ]

y <- GSVA_long_noTMM$GSVA_Score

cor_values <- numeric(ncol(earlyNeuroblast_matrix))
p_values   <- numeric(ncol(earlyNeuroblast_matrix))



for (i in seq_len(ncol(earlyNeuroblast_matrix))) {
  a <- as.numeric(earlyNeuroblast_matrix[, i])
  b <- complete.cases(a, y)
  correlation <- cor.test(a[b], y[b], method = "spearman", exact = FALSE)
  cor_values[i] <- unname(correlation$estimate)
  p_values[i]   <- correlation$p.value
}

cor_results <- data.frame(
  Gene = colnames(earlyNeuroblast_matrix),
  estimate = cor_values,
  p_value = p_values,
  fdr = p.adjust(p_values, "BH"),
  stringsAsFactors = FALSE
)[order(p.adjust(p_values, "BH"), -abs(cor_values)), ]

earlyNeuroblast_noTMM_sig <- cor_results[cor_results$fdr < 0.01, ]





 



