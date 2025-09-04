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
adrenergic_markers <- c("TFAP2B", "ISL1", "EYA1", "GATA3", "PHOX2B", "HAND1", "KLF7", "SIX3", "ASCL1", "HAND2", "HEY1",
                        "PHOX2A", "PBX3", "KLF13", "SOX11", "GATA2", "SATB1", "DACH1", "TBX2", "DBH", "ASCL1", "TH")

mesenchymal_markers <- c("VIM", "FN1", "MEOX2", "ID1", "EGR3", "AEBP1", "CBFB", "IRF3", "IRF2", "IRF1", "TBX18", "MAFF", 
                         "RUNX2", "ZFP36L1", "NR3C1", "BHLHE41", "GLIS3", "RUNX1", "FOSL1", "FOSL2", "ELK4", "IFI16", "SIX4", 
                         "FLI1", "MAML2", "SMAD3", "DCAF6", "WWTR1", "SOX9", "MEF2D", "ZNF217", "PRRX1", "CREG1", "NOTCH2", 
                         "SIX1", "MEOX1", "SNAI2", "CD44", "YAP1", "JUN")

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

## So, TERT is negatively associated with mesenchymal characteristics. what is up with heatmap? ::)


###########################################################################################################################################
# correlation between GSVA of NO_TMM signatures and mesenchymal/adrenergic genes.

# source("aucGSVA.R")
kfold_genes <- list(NO_TMM = c("CPNE8", "PGM2L1", "CNR1", "LIFR", "HOXC9", "SNX16", "HECW2", "ALDH3A2",
                               "THSD7A", "CPNE3", "IGSF10"))
gsvapar <- gsvaParam(as.matrix(Expression), kfold_genes, kcdf = "Gaussian")
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

# adrenergic and mesenchymal markers from Thirant paper.
adrenergic_markers <- c("TFAP2B", "ISL1", "EYA1", "GATA3", "PHOX2B", "HAND1", "KLF7", "SIX3", "ASCL1", "HAND2", "HEY1",
                        "PHOX2A", "PBX3", "KLF13", "SOX11", "GATA2", "SATB1", "DACH1", "TBX2", "DBH", "ASCL1", "TH")

mesenchymal_markers <- c("VIM", "FN1", "MEOX2", "ID1", "EGR3", "AEBP1", "CBFB", "IRF3", "IRF2", "IRF1", "TBX18", "MAFF", 
                         "RUNX2", "ZFP36L1", "NR3C1", "BHLHE41", "GLIS3", "RUNX1", "FOSL1", "FOSL2", "ELK4", "IFI16", "SIX4", 
                         "FLI1", "MAML2", "SMAD3", "DCAF6", "WWTR1", "SOX9", "MEF2D", "ZNF217", "PRRX1", "CREG1", "NOTCH2", 
                         "SIX1", "MEOX1", "SNAI2", "CD44", "YAP1", "JUN")

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




# adrenergic and mesenchymal markers from Thirant paper.
adrenergic_markers <- c("TFAP2B", "ISL1", "EYA1", "GATA3", "PHOX2B", "HAND1", "KLF7", "SIX3", "ASCL1", "HAND2", "HEY1",
                        "PHOX2A", "PBX3", "KLF13", "SOX11", "GATA2", "SATB1", "DACH1", "TBX2", "DBH", "ASCL1", "TH")

mesenchymal_markers <- c("VIM", "FN1", "MEOX2", "ID1", "EGR3", "AEBP1", "CBFB", "IRF3", "IRF2", "IRF1", "TBX18", "MAFF", 
                         "RUNX2", "ZFP36L1", "NR3C1", "BHLHE41", "GLIS3", "RUNX1", "FOSL1", "FOSL2", "ELK4", "IFI16", "SIX4", 
                         "FLI1", "MAML2", "SMAD3", "DCAF6", "WWTR1", "SOX9", "MEF2D", "ZNF217", "PRRX1", "CREG1", "NOTCH2", 
                         "SIX1", "MEOX1", "SNAI2", "CD44", "YAP1", "JUN")

# now looking at mesenchymal marker ranks.
mesenchymal_ranks <- t_NBL[, colnames(t_NBL) %in% mesenchymal_markers, drop = FALSE]
mesenchymal_ranks <- mesenchymal_ranks[match(rownames(t_NBL), rownames(mesenchymal_ranks)), ]




source("EXTEND/ComponentAndMarkerFunction.r")
source("EXTEND/ComponentOneAndMarkerFunction.r")
source("EXTEND/ComponentTwoAndMarkerFunction.r")
source("EXTEND/InputData.r")
source("EXTEND/IterativeRS.r")
source("EXTEND/IterativeRS.r")
source("EXTEND/MarkerFunction.r")
source("EXTEND/RunEXTEND.r")

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


# adrenergic and mesenchymal markers from Thirant paper.
adrenergic_markers <- c("TFAP2B", "ISL1", "EYA1", "GATA3", "PHOX2B", "HAND1", "KLF7", "SIX3", "ASCL1", "HAND2", "HEY1",
                        "PHOX2A", "PBX3", "KLF13", "SOX11", "GATA2", "SATB1", "DACH1", "TBX2", "DBH", "ASCL1", "TH")

mesenchymal_markers <- c("VIM", "FN1", "MEOX2", "ID1", "EGR3", "AEBP1", "CBFB", "IRF3", "IRF2", "IRF1", "TBX18", "MAFF", 
                         "RUNX2", "ZFP36L1", "NR3C1", "BHLHE41", "GLIS3", "RUNX1", "FOSL1", "FOSL2", "ELK4", "IFI16", "SIX4", 
                         "FLI1", "MAML2", "SMAD3", "DCAF6", "WWTR1", "SOX9", "MEF2D", "ZNF217", "PRRX1", "CREG1", "NOTCH2", 
                         "SIX1", "MEOX1", "SNAI2", "CD44", "YAP1", "JUN")

# now looking at adrenergic marker ranks.
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






