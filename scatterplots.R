library(ggplot2)
library(ggpubr)
library(dplyr)

## extend tools.
source("EXTEND/ComponentAndMarkerFunction.r")
source("EXTEND/ComponentOneAndMarkerFunction.r")
source("EXTEND/ComponentTwoAndMarkerFunction.r")
source("EXTEND/InputData.r")
source("EXTEND/IterativeRS.r")
source("EXTEND/IterativeRS.r")
source("EXTEND/MarkerFunction.r")
source("EXTEND/RunEXTEND.r")

###

## : scatterplot era.

## scatterplot: Proliferation GSVA vs. EXTEND.

source("DataCleaning.R")
gene_set_list_proliferation <- list(proliferation = c("MKI67", "PCNA", "BIRC5", "CEP55", "TOP2A", "CDK1", "CDC20", "CCNB1", 
                                                      "CCNB2", "CCNA2", "UBE2C", "AURKA", "AURKB", "PLK1", "PTTG1", "NUSAP1"))
gsvapar <- gsvaParam(as.matrix(Expression), gene_set_list_proliferation, kcdf = "Gaussian")
gsva_result <- gsva(gsvapar)
GSVA_df <- as.data.frame(gsva_result)
GSVA_long_proliferation <- pivot_longer(GSVA_df, cols = everything(), names_to = "SampleID", values_to = "GSVA_Score")


extendScores <- RunEXTEND(as.matrix(Expression))
telomeraseScores <- read_delim("TelomeraseScores.txt")
telomeraseScores <- telomeraseScores[, c("SampleID", "NormEXTENDScores")]
telomeraseScores <- as.data.frame(telomeraseScores)
telomeraseScores <- left_join(telomeraseScores, metadata[, c("SampleID", "TMM", "COG.Risk.Group")], by = "SampleID")

telomeraseScores <- telomeraseScores[match(GSVA_long_proliferation$SampleID, telomeraseScores$SampleID), ]

# joining the two results.
GSVA_long_proliferation <- left_join(GSVA_long_proliferation, telomeraseScores, by = "SampleID")


# for scatterplot first computing Pearson correlation value.
cor_val <- cor(GSVA_long_proliferation$GSVA_Score, GSVA_long_proliferation$NormEXTENDScores, method = "pearson")
r_text <- paste("r =", round(cor_val, 3))

ggplot(GSVA_long_proliferation, aes(y = NormEXTENDScores, x = GSVA_Score, colour = TMM)) + 
  scale_color_manual(values = c("ALT" = "darkblue", "Telomerase" = "darkred", "black" = "NO_TMM")) + 
  geom_point(size = 4, alpha = 0.9) + 
  labs(x = "GSVA Signature Score From Proliferation Markers", y = "Normalized EXTEND Scores") + 
  geom_smooth(method = "lm", level = 0.95, se = TRUE, color = "black", fill = "lightgray") + 
  annotate("text", x = 0.25, y = 1.3, label = r_text, size = 7) +
  theme_classic() +
  theme(
    axis.title.x = element_text(size = 24),
    axis.title.y = element_text(size = 24),
    axis.text.x  = element_text(size = 14, face = "bold"),
    axis.text.y  = element_text(size = 14, face = "bold"),
    legend.title = element_text(size = 18, face = "bold"),
    legend.text  = element_text(size = 16),
    legend.position = c(-0.5, 1.3)
  )

# same scatterplot colored with risk group.
ggplot(GSVA_long_proliferation, aes(y = NormEXTENDScores, x = GSVA_Score, colour = COG.Risk.Group)) + 
  scale_color_manual(values = c("Intermediate Risk" = "darkblue", "High Risk" = "darkred", "black" = "Low Risk")) + 
  geom_point(size = 4, alpha = 0.9) + 
  labs(x = "GSVA Signature Score From Proliferation Markers", y = "Normalized EXTEND Scores") + 
  geom_smooth(method = "lm", level = 0.95, se = TRUE, color = "black", fill = "lightgray") + 
  annotate("text", x = 0.25, y = 1.3, label = r_text, size = 7) +
  theme_classic() +
  theme(
    axis.title.x = element_text(size = 24),
    axis.title.y = element_text(size = 24),
    axis.text.x  = element_text(size = 14, face = "bold"),
    axis.text.y  = element_text(size = 14, face = "bold"),
    legend.title = element_text(size = 18, face = "bold"),
    legend.text  = element_text(size = 16),
    legend.position = c(-0.5, 1.3)
  )


## scatterplot: adr gsva vs. proliferation gsva.
source("DataCleaning.R")

gene_set_list_proliferation <- list(proliferation = c("MKI67", "PCNA", "BIRC5", "CEP55", "TOP2A", "CDK1", "CDC20", "CCNB1", 
                                                      "CCNB2", "CCNA2", "UBE2C", "AURKA", "AURKB", "PLK1", "PTTG1", "NUSAP1"))
gsvapar <- gsvaParam(as.matrix(Expression), gene_set_list_proliferation, kcdf = "Gaussian")
gsva_result <- gsva(gsvapar)
GSVA_df <- as.data.frame(gsva_result)
GSVA_long_proliferation <- pivot_longer(GSVA_df, cols = everything(), names_to = "SampleID", values_to = "GSVA_Score")

gene_list_adrenergic <- list("ADR" = c("DLK1", "DBH", "PHOX2A", "PHOX2B", "GATA2", "GATA3"))

gsvapar <- gsvaParam(as.matrix(Expression), gene_list_adrenergic, kcdf = "Gaussian")
gsva_result <- gsva(gsvapar)
GSVA_df <- as.data.frame(gsva_result)
GSVA_long_adr <- pivot_longer(GSVA_df, cols = everything(), names_to = "SampleID", values_to = "GSVA_Score_Adr")

# joining the two results.
GSVA_long_adr <- left_join(GSVA_long_adr, GSVA_long_proliferation, by = "SampleID")
GSVA_long_adr <- left_join(GSVA_long_adr, metadata[, c("SampleID", "COG.Risk.Group", "TMM")], by = "SampleID")

# for scatterplot first computing Pearson correlation value.
cor_val <- cor(GSVA_long_adr$GSVA_Score_Adr, GSVA_long_adr$GSVA_Score, method = "pearson")
r_text <- paste("r =", round(cor_val, 3))

ggplot(GSVA_long_adr, aes(y = GSVA_Score_Adr, x = GSVA_Score, colour = TMM)) + 
  scale_color_manual(values = c("ALT" = "darkblue", "Telomerase" = "darkred")) + 
  geom_point(size = 4, alpha = 0.9) + 
  labs(y = "GSVA Score from ADR Markers", x = "Proliferation GSVA Scores") + 
  geom_smooth(method = "lm", level = 0.95, se = TRUE, color = "black", fill = "lightgray") + 
  annotate("text", x = 0.25, y = 1.3, label = r_text, size = 7) +
  theme_classic() +
  theme(
    axis.title.x = element_text(size = 24),
    axis.title.y = element_text(size = 24),
    axis.text.x  = element_text(size = 14, face = "bold"),
    
    axis.text.y  = element_text(size = 14, face = "bold"),
    legend.title = element_text(size = 18, face = "bold"),
    legend.text  = element_text(size = 16),
    legend.position = c(-0.5, 1.3)
  )


## scatterplot: mes gsva vs. proliferation gsva.

## scatterplot: scp gsva vs. NO_TMM GSVA.
source("DataCleaning.R")
gene_set_list_scp <- list("SCP" = c("PLP1", "MPZ", "CDH19", "ERBB3", "ERBB4"))
gene_list_noTMM <- list("NO_TMM Genes" = c("CPNE8", "PGM2L1", "CNR1", "LIFR", "HOXC9", "SNX16", "HECW2", "ALDH3A2",
                                           "THSD7A", "CPNE3", "IGSF10"))


gsvapar <- gsvaParam(as.matrix(Expression), gene_set_list_scp, kcdf = "Gaussian")
gsva_result <- gsva(gsvapar)
GSVA_df <- as.data.frame(gsva_result)
GSVA_long_scp <- pivot_longer(GSVA_df, cols = everything(), names_to = "SampleID", values_to = "GSVA_Score_scp")


Expression <- as.matrix(Expression)
noTMM_gsva <- gsvaParam(Expression, gene_list_noTMM, kcdf = "Gaussian")
GSVA_result_noTMM <- gsva(noTMM_gsva)
GSVA_df_noTMM <- as.data.frame(GSVA_result_noTMM)
GSVA_long_noTMM <- pivot_longer(GSVA_df_noTMM, cols = everything(), names_to = "SampleID", values_to = "GSVA_Score")


# joining the two results.
GSVA_long_scp <- left_join(GSVA_long_scp, GSVA_long_noTMM, by = "SampleID")
GSVA_long_scp <- left_join(GSVA_long_scp, metadata[, c("SampleID", "COG.Risk.Group", "TMM")], by = "SampleID")


# for scatterplot first computing Pearson correlation value.
cor_val <- cor(GSVA_long_scp$GSVA_Score_scp, GSVA_long_scp$GSVA_Score, method = "pearson")
r_text <- paste("r =", round(cor_val, 3))

ggplot(GSVA_long_scp, aes(y = GSVA_Score_scp, x = GSVA_Score, colour = TMM)) + 
  scale_color_manual(values = c("ALT" = "darkblue", "Telomerase" = "darkred")) + 
  geom_point(size = 4, alpha = 0.9) + 
  labs(x = "GSVA Signature from NO_TMM markers", y = "GSVA Signature from SCP Markers") + 
  geom_smooth(method = "lm", level = 0.95, se = TRUE, color = "black", fill = "lightgray") + 
  annotate("text", x = 0.25, y = 1.3, label = r_text, size = 7) +
  theme_classic() +
  theme(
    axis.title.x = element_text(size = 24),
    axis.title.y = element_text(size = 24),
    axis.text.x  = element_text(size = 14, face = "bold"),
    axis.text.y  = element_text(size = 14, face = "bold"),
    legend.title = element_text(size = 18, face = "bold"),
    legend.text  = element_text(size = 16),
    legend.position = c(-0.5, 1.3)
  )

### Now: risk group.
ggplot(GSVA_long_scp, aes(y = GSVA_Score_scp, x = GSVA_Score, colour = COG.Risk.Group)) + 
  scale_color_manual(values = c("Intermediate Risk" = "darkblue", "High Risk" = "darkred")) + 
  geom_point(size = 4, alpha = 0.9) + 
  labs(x = "GSVA Signature from NO_TMM markers", y = "GSVA Signature from SCP Markers") + 
  geom_smooth(method = "lm", level = 0.95, se = TRUE, color = "black", fill = "lightgray") + 
  annotate("text", x = 0.25, y = 1.3, label = r_text, size = 7) +
  theme_classic() +
  theme(
    axis.title.x = element_text(size = 24),
    axis.title.y = element_text(size = 24),
    axis.text.x  = element_text(size = 14, face = "bold"),
    axis.text.y  = element_text(size = 14, face = "bold"),
    legend.title = element_text(size = 18, face = "bold"),
    legend.text  = element_text(size = 16),
    legend.position = c(-0.5, 1.3)
  )



## scatterplot: late neuroblast gsva vs. NO_TMM GSVA.
source("DataCleaning.R")
late_neuroblast <- list("lateNeuroblast" = c("GAP43", "STMN2", "ISL1", "SYN3", "IL7"))

gsvapar <- gsvaParam(as.matrix(Expression), late_neuroblast, kcdf = "Gaussian")
gsva_result <- gsva(gsvapar)
GSVA_df <- as.data.frame(gsva_result)
GSVA_long_lateNeuroblast <- pivot_longer(GSVA_df, cols = everything(), names_to = "SampleID", values_to = "GSVA_Score_lateNeuroblast")

# joining the two results.
GSVA_long_lateNeuroblast <- left_join(GSVA_long_lateNeuroblast, GSVA_long_noTMM, by = "SampleID")
GSVA_long_lateNeuroblast <- left_join(GSVA_long_lateNeuroblast, metadata[, c("SampleID", "COG.Risk.Group", "TMM")], by = "SampleID")


# for scatterplot first computing Pearson correlation value.
cor_val <- cor(GSVA_long_lateNeuroblast$GSVA_Score_lateNeuroblast, GSVA_long_lateNeuroblast$GSVA_Score, method = "pearson")
r_text <- paste("r =", round(cor_val, 3))

ggplot(GSVA_long_lateNeuroblast, aes(y = GSVA_Score_lateNeuroblast, x = GSVA_Score, colour = TMM)) + 
  scale_color_manual(values = c("ALT" = "darkblue", "Telomerase" = "darkred")) + 
  geom_point(size = 4, alpha = 0.9) + 
  labs(x = "GSVA Signature from NO_TMM markers", y = "GSVA Signature from Late Neuroblast Markers") + 
  geom_smooth(method = "lm", level = 0.95, se = TRUE, color = "black", fill = "lightgray") + 
  annotate("text", x = 0.25, y = 1.3, label = r_text, size = 7) +
  theme_classic() +
  theme(
    axis.title.x = element_text(size = 24),
    axis.title.y = element_text(size = 24),
    axis.text.x  = element_text(size = 14, face = "bold"),
    axis.text.y  = element_text(size = 14, face = "bold"),
    legend.title = element_text(size = 18, face = "bold"),
    legend.text  = element_text(size = 16),
    legend.position = c(-0.5, 1.3)
  )

ggplot(GSVA_long_lateNeuroblast, aes(y = GSVA_Score_lateNeuroblast, x = GSVA_Score, colour = COG.Risk.Group)) + 
  scale_color_manual(values = c("Intermediate Risk" = "darkblue", "High Risk" = "darkred", "black" = "Low Risk")) + 
  geom_point(size = 4, alpha = 0.9) + 
  labs(x = "GSVA Signature from NO_TMM markers", y = "GSVA Signature from Late Neuroblast Markers") + 
  geom_smooth(method = "lm", level = 0.95, se = TRUE, color = "black", fill = "lightgray") + 
  annotate("text", x = 0.25, y = 1.3, label = r_text, size = 7) +
  theme_classic() +
  theme(
    axis.title.x = element_text(size = 24),
    axis.title.y = element_text(size = 24),
    axis.text.x  = element_text(size = 14, face = "bold"),
    axis.text.y  = element_text(size = 14, face = "bold"),
    legend.title = element_text(size = 18, face = "bold"),
    legend.text  = element_text(size = 16),
    legend.position = c(-0.5, 1.3)
  )




## scatterplot: MES vs. SCP.
source("DataCleaning.R")
gene_list_mesenchymal <- list("MES" = c("VIM", "FN1", "YAP1", "SNAI2", "PRRX1", "WWTR1"))
# gene_list_mesenchymal <- list("MES" = c("PRRX1", "RUNX1", "FOSL1", "JUN", "YAP1", "WWTR1", "MEOX1", 
#                                         "MEOX2", "SNAI2", "CD44", "FN1", "VIM"))

gsvapar <- gsvaParam(as.matrix(Expression), gene_list_mesenchymal, kcdf = "Gaussian")
gsva_result <- gsva(gsvapar)
GSVA_df <- as.data.frame(gsva_result)
GSVA_long_mes <- pivot_longer(GSVA_df, cols = everything(), names_to = "SampleID", values_to = "GSVA_Score_Mes")

gene_list_scp <- list("SCP" = c("PLP1", "MPZ", "CDH19", "ERBB3", "ERBB4"))
gsvapar <- gsvaParam(as.matrix(Expression), gene_list_scp, kcdf = "Gaussian")
gsva_result <- gsva(gsvapar)
GSVA_df <- as.data.frame(gsva_result)
GSVA_long_scp <- pivot_longer(GSVA_df, cols = everything(), names_to = "SampleID", values_to = "GSVA_Score_SCP")


# joining the two results.
GSVA_long_mes <- left_join(GSVA_long_mes, GSVA_long_scp, by = "SampleID")
GSVA_long_mes <- left_join(GSVA_long_mes, metadata[, c("SampleID", "COG.Risk.Group", "TMM")], by = "SampleID")


# for scatterplot first computing Pearson correlation value.
cor_val <- cor(GSVA_long_mes$GSVA_Score_Mes, GSVA_long_mes$GSVA_Score_SCP, method = "pearson")
r_text <- paste("r =", round(cor_val, 3))

ggplot(GSVA_long_mes, aes(y = GSVA_Score_Mes, x = GSVA_Score_SCP, colour = TMM)) + 
  scale_color_manual(values = c("ALT" = "darkblue", "Telomerase" = "darkred", "black" = "NO_TMM")) + 
  geom_point(size = 4, alpha = 0.9) + 
  labs(x = "GSVA Signature from SCP markers", y = "GSVA Signature from Mesenchymal Markers") + 
  geom_smooth(method = "lm", level = 0.95, se = TRUE, color = "black", fill = "lightgray") + 
  annotate("text", x = 0.25, y = 1.3, label = r_text, size = 7) +
  theme_classic() +
  theme(
    axis.title.x = element_text(size = 24),
    axis.title.y = element_text(size = 24),
    axis.text.x  = element_text(size = 14, face = "bold"),
    axis.text.y  = element_text(size = 14, face = "bold"),
    legend.title = element_text(size = 18, face = "bold"),
    legend.text  = element_text(size = 16),
    legend.position = c(-0.5, 1.3)
  )

ggplot(GSVA_long_mes, aes(y = GSVA_Score_Mes, x = GSVA_Score_SCP, colour = COG.Risk.Group)) + 
  scale_color_manual(values = c("Intermediate Risk" = "darkblue", "High Risk" = "darkred", "black" = "Low Risk")) + 
  geom_point(size = 4, alpha = 0.9) + 
  labs(x = "GSVA Signature from SCP markers", y = "GSVA Signature from Mesenchymal Markers") + 
  geom_smooth(method = "lm", level = 0.95, se = TRUE, color = "black", fill = "lightgray") + 
  annotate("text", x = 0.25, y = 1.3, label = r_text, size = 7) +
  theme_classic() +
  theme(
    axis.title.x = element_text(size = 24),
    axis.title.y = element_text(size = 24),
    axis.text.x  = element_text(size = 14, face = "bold"),
    axis.text.y  = element_text(size = 14, face = "bold"),
    legend.title = element_text(size = 18, face = "bold"),
    legend.text  = element_text(size = 16),
    legend.position = c(-0.5, 1.3)
  )


#### scatterplot: MES vs. NO_TMM.
gene_list_mesenchymal <- list("MES" = c("FN1", "VIM", "SNAI2", "PRRX1", "YAP1", "WWTR1"))
gene_list_noTMM <- list("NO_TMM Genes" = c("CPNE8", "PGM2L1", "CNR1", "LIFR", "HOXC9", "SNX16", "HECW2", "ALDH3A2",
                                           "THSD7A", "CPNE3", "IGSF10"))


gsvapar <- gsvaParam(as.matrix(Expression), gene_list_mesenchymal, kcdf = "Gaussian")
gsva_result <- gsva(gsvapar)
GSVA_df <- as.data.frame(gsva_result)
GSVA_long_mes <- pivot_longer(GSVA_df, cols = everything(), names_to = "SampleID", values_to = "GSVA_Score_mes")


Expression <- as.matrix(Expression)
noTMM_gsva <- gsvaParam(Expression, gene_list_noTMM, kcdf = "Gaussian")
GSVA_result_noTMM <- gsva(noTMM_gsva)
GSVA_df_noTMM <- as.data.frame(GSVA_result_noTMM)
GSVA_long_noTMM <- pivot_longer(GSVA_df_noTMM, cols = everything(), names_to = "SampleID", values_to = "GSVA_Score")


# joining the two results.
GSVA_long_mes <- left_join(GSVA_long_mes, GSVA_long_noTMM, by = "SampleID")
GSVA_long_mes <- left_join(GSVA_long_mes, metadata[, c("SampleID", "COG.Risk.Group", "TMM")], by = "SampleID")


# for scatterplot first computing Pearson correlation value.
cor_val <- cor(GSVA_long_mes$GSVA_Score_mes, GSVA_long_mes$GSVA_Score, method = "pearson")
r_text <- paste("r =", round(cor_val, 3))

ggplot(GSVA_long_mes, aes(y = GSVA_Score_mes, x = GSVA_Score, colour = TMM)) + 
  scale_color_manual(values = c("ALT" = "darkblue", "Telomerase" = "darkred", "black" = "NO_TMM")) + 
  geom_point(size = 4, alpha = 0.9) + 
  labs(x = "GSVA Signature from NO_TMM markers", y = "GSVA Signature from MES Markers") + 
  geom_smooth(method = "lm", level = 0.95, se = TRUE, color = "black", fill = "lightgray") + 
  annotate("text", x = 0.25, y = 1.3, label = r_text, size = 7) +
  theme_classic() +
  theme(
    axis.title.x = element_text(size = 24),
    axis.title.y = element_text(size = 24),
    axis.text.x  = element_text(size = 14, face = "bold"),
    axis.text.y  = element_text(size = 14, face = "bold"),
    legend.title = element_text(size = 18, face = "bold"),
    legend.text  = element_text(size = 16),
    legend.position = c(-0.5, 1.3)
  )

ggplot(GSVA_long_mes, aes(y = GSVA_Score_mes, x = GSVA_Score, colour = COG.Risk.Group)) + 
  scale_color_manual(values = c("Intermediate Risk" = "darkblue", "High Risk" = "darkred")) + 
  geom_point(size = 4, alpha = 0.9) + 
  labs(x = "GSVA Signature from NO_TMM markers", y = "GSVA Signature from MES Markers") + 
  geom_smooth(method = "lm", level = 0.95, se = TRUE, color = "black", fill = "lightgray") + 
  annotate("text", x = 0.25, y = 1.3, label = r_text, size = 7) +
  theme_classic() +
  theme(
    axis.title.x = element_text(size = 24),
    axis.title.y = element_text(size = 24),
    axis.text.x  = element_text(size = 14, face = "bold"),
    axis.text.y  = element_text(size = 14, face = "bold"),
    legend.title = element_text(size = 18, face = "bold"),
    legend.text  = element_text(size = 16),
    legend.position = c(-0.5, 1.3)
  )


#### scatterplot: early neuroblast vs. NO_TMM.
early_neuroblast <- list("earlyNeuroblast" = c("ALK", "TOP2A", "MKI67", "ISL1", "STMN2"))
gene_list_noTMM <- list("NO_TMM Genes" = c("CPNE8", "PGM2L1", "CNR1", "LIFR", "HOXC9", "SNX16", "HECW2", "ALDH3A2",
                                           "THSD7A", "CPNE3", "IGSF10"))


gsvapar <- gsvaParam(as.matrix(Expression), early_neuroblast, kcdf = "Gaussian")
gsva_result <- gsva(gsvapar)
GSVA_df <- as.data.frame(gsva_result)
GSVA_long_neuroblast <- pivot_longer(GSVA_df, cols = everything(), names_to = "SampleID", 
                                     values_to = "GSVA_Score_neuroblast")


Expression <- as.matrix(Expression)
noTMM_gsva <- gsvaParam(Expression, gene_list_noTMM, kcdf = "Gaussian")
GSVA_result_noTMM <- gsva(noTMM_gsva)
GSVA_df_noTMM <- as.data.frame(GSVA_result_noTMM)
GSVA_long_noTMM <- pivot_longer(GSVA_df_noTMM, cols = everything(), names_to = "SampleID", values_to = "GSVA_Score")


# joining the two results.
GSVA_long_neuroblast <- left_join(GSVA_long_neuroblast, GSVA_long_noTMM, by = "SampleID")
GSVA_long_neuroblast <- left_join(GSVA_long_neuroblast, metadata[, c("SampleID", "COG.Risk.Group", "TMM")], by = "SampleID")


# for scatterplot first computing Pearson correlation value.
cor_val <- cor(GSVA_long_neuroblast$GSVA_Score_neuroblast, GSVA_long_neuroblast$GSVA_Score, method = "pearson")
r_text <- paste("r =", round(cor_val, 3))

ggplot(GSVA_long_neuroblast, aes(y = GSVA_Score_neuroblast, x = GSVA_Score, colour = TMM)) + 
  scale_color_manual(values = c("ALT" = "darkblue", "Telomerase" = "darkred")) + 
  geom_point(size = 4, alpha = 0.9) + 
  labs(x = "GSVA Signature from NO_TMM markers", y = "GSVA Signature from Early Neuroblast Markers") + 
  geom_smooth(method = "lm", level = 0.95, se = TRUE, color = "black", fill = "lightgray") + 
  annotate("text", x = 0.25, y = 1.3, label = r_text, size = 7) +
  theme_classic() +
  theme(
    axis.title.x = element_text(size = 24),
    axis.title.y = element_text(size = 24),
    axis.text.x  = element_text(size = 14, face = "bold"),
    axis.text.y  = element_text(size = 14, face = "bold"),
    legend.title = element_text(size = 18, face = "bold"),
    legend.text  = element_text(size = 16),
    legend.position = c(-0.5, 1.3)
  )

ggplot(GSVA_long_neuroblast, aes(y = GSVA_Score_neuroblast, x = GSVA_Score, colour = COG.Risk.Group)) + 
  scale_color_manual(values = c("Intermediate Risk" = "darkblue", "High Risk" = "darkred")) + 
  geom_point(size = 4, alpha = 0.9) + 
  labs(x = "GSVA Signature from NO_TMM markers", y = "GSVA Signature from Early Neuroblast Markers") + 
  geom_smooth(method = "lm", level = 0.95, se = TRUE, color = "black", fill = "lightgray") + 
  annotate("text", x = 0.25, y = 1.3, label = r_text, size = 7) +
  theme_classic() +
  theme(
    axis.title.x = element_text(size = 24),
    axis.title.y = element_text(size = 24),
    axis.text.x  = element_text(size = 14, face = "bold"),
    axis.text.y  = element_text(size = 14, face = "bold"),
    legend.title = element_text(size = 18, face = "bold"),
    legend.text  = element_text(size = 16),
    legend.position = c(-0.5, 1.3)
  )


#### scatterplot: early neuroblast vs. EXTEND.
source("DataCleaning.R")

early_neuroblast <- list("earlyNeuroblast" = c("ALK", "TOP2A", "MKI67", "ISL1", "STMN2"))

gsvapar <- gsvaParam(as.matrix(Expression), early_neuroblast, kcdf = "Gaussian")
gsva_result <- gsva(gsvapar)
GSVA_df <- as.data.frame(gsva_result)
GSVA_long_neuroblast <- pivot_longer(GSVA_df, cols = everything(), names_to = "SampleID", 
                                     values_to = "GSVA_Score_neuroblast")


extendScores <- RunEXTEND(as.matrix(Expression))
telomeraseScores <- read_delim("TelomeraseScores.txt")
telomeraseScores <- telomeraseScores[, c("SampleID", "NormEXTENDScores")]
telomeraseScores <- as.data.frame(telomeraseScores)
telomeraseScores <- left_join(telomeraseScores, metadata[, c("SampleID", "TMM", "COG.Risk.Group")], by = "SampleID")

# joining the two results.
GSVA_long_neuroblast <- left_join(GSVA_long_neuroblast, telomeraseScores, by = "SampleID")

# for scatterplot first computing Pearson correlation value.
cor_val <- cor(GSVA_long_neuroblast$GSVA_Score_neuroblast, GSVA_long_neuroblast$NormEXTENDScores, method = "pearson")
r_text <- paste("r =", round(cor_val, 3))

ggplot(GSVA_long_neuroblast, aes(y = GSVA_Score_neuroblast, x = NormEXTENDScores, colour = TMM)) + 
  scale_color_manual(values = c("ALT" = "darkblue", "Telomerase" = "darkred")) + 
  geom_point(size = 4, alpha = 0.9) + 
  labs(x = "EXTEND Scores", y = "GSVA Signature from Early Neuroblast Markers") + 
  geom_smooth(method = "lm", level = 0.95, se = TRUE, color = "black", fill = "lightgray") + 
  annotate("text", x = 0.25, y = 1.3, label = r_text, size = 7) +
  theme_classic() +
  theme(
    axis.title.x = element_text(size = 24),
    axis.title.y = element_text(size = 24),
    axis.text.x  = element_text(size = 14, face = "bold"),
    axis.text.y  = element_text(size = 14, face = "bold"),
    legend.title = element_text(size = 18, face = "bold"),
    legend.text  = element_text(size = 16),
    legend.position = c(-0.5, 1.3)
  )

ggplot(GSVA_long_neuroblast, aes(y = GSVA_Score_neuroblast, x = NormEXTENDScores, colour = COG.Risk.Group)) + 
  scale_color_manual(values = c("Intermediate Risk" = "darkblue", "High Risk" = "darkred", "Low Risk" = "gray")) + 
  geom_point(size = 4, alpha = 0.9) + 
  labs(x = "EXTEND Scores", y = "GSVA Signature from Early Neuroblast Markers") + 
  geom_smooth(method = "lm", level = 0.95, se = TRUE, color = "black", fill = "lightgray") + 
  annotate("text", x = 0.25, y = 1.3, label = r_text, size = 7) +
  theme_classic() +
  theme(
    axis.title.x = element_text(size = 24),
    axis.title.y = element_text(size = 24),
    axis.text.x  = element_text(size = 14, face = "bold"),
    axis.text.y  = element_text(size = 14, face = "bold"),
    legend.title = element_text(size = 18, face = "bold"),
    legend.text  = element_text(size = 16),
    legend.position = c(-0.5, 1.3)
  )



## scatterplot: late neuroblast gsva vs. EXTEND.
source("DataCleaning.R")

late_neuroblast <- list("lateNeuroblast" = c("GAP43", "STMN2", "ISL1", "SYN3", "IL7"))

gsvapar <- gsvaParam(as.matrix(Expression), late_neuroblast, kcdf = "Gaussian")
gsva_result <- gsva(gsvapar)
GSVA_df <- as.data.frame(gsva_result)
GSVA_long_neuroblast <- pivot_longer(GSVA_df, cols = everything(), names_to = "SampleID", 
                                     values_to = "GSVA_Score_neuroblast")


extendScores <- RunEXTEND(as.matrix(Expression))
telomeraseScores <- read_delim("TelomeraseScores.txt")
telomeraseScores <- telomeraseScores[, c("SampleID", "NormEXTENDScores")]
telomeraseScores <- as.data.frame(telomeraseScores)
telomeraseScores <- left_join(telomeraseScores, metadata[, c("SampleID", "TMM", "COG.Risk.Group")], by = "SampleID")

# joining the two results.
GSVA_long_neuroblast <- left_join(GSVA_long_neuroblast, telomeraseScores, by = "SampleID")

# for scatterplot first computing Pearson correlation value.
cor_val <- cor(GSVA_long_neuroblast$GSVA_Score_neuroblast, GSVA_long_neuroblast$NormEXTENDScores, method = "pearson")
r_text <- paste("r =", round(cor_val, 3))

ggplot(GSVA_long_neuroblast, aes(y = GSVA_Score_neuroblast, x = NormEXTENDScores, colour = TMM)) + 
  scale_color_manual(values = c("ALT" = "darkblue", "Telomerase" = "darkred")) + 
  geom_point(size = 4, alpha = 0.9) + 
  labs(x = "EXTEND Scores", y = "GSVA Signature from Late Neuroblast Markers") + 
  geom_smooth(method = "lm", level = 0.95, se = TRUE, color = "black", fill = "lightgray") + 
  annotate("text", x = 0.25, y = 1.3, label = r_text, size = 7) +
  theme_classic() +
  theme(
    axis.title.x = element_text(size = 24),
    axis.title.y = element_text(size = 24),
    axis.text.x  = element_text(size = 14, face = "bold"),
    axis.text.y  = element_text(size = 14, face = "bold"),
    legend.title = element_text(size = 18, face = "bold"),
    legend.text  = element_text(size = 16),
    legend.position = c(-0.5, 1.3)
  )

ggplot(GSVA_long_neuroblast, aes(y = GSVA_Score_neuroblast, x = NormEXTENDScores, colour = COG.Risk.Group)) + 
  scale_color_manual(values = c("Intermediate Risk" = "darkblue", "High Risk" = "darkred", "Low Risk" = "gray")) + 
  geom_point(size = 4, alpha = 0.9) + 
  labs(x = "EXTEND Scores", y = "GSVA Signature from Late Neuroblast Markers") + 
  geom_smooth(method = "lm", level = 0.95, se = TRUE, color = "black", fill = "lightgray") + 
  annotate("text", x = 0.25, y = 1.3, label = r_text, size = 7) +
  theme_classic() +
  theme(
    axis.title.x = element_text(size = 24),
    axis.title.y = element_text(size = 24),
    axis.text.x  = element_text(size = 14, face = "bold"),
    axis.text.y  = element_text(size = 14, face = "bold"),
    legend.title = element_text(size = 18, face = "bold"),
    legend.text  = element_text(size = 16),
    legend.position = c(-0.5, 1.3)
  )



#### scatterplot: extend vs. scp.
source("DataCleaning.R")
gene_set_list_scp <- list("SCP" = c("PLP1", "MPZ", "CDH19", "ERBB3", "ERBB4"))


gsvapar <- gsvaParam(as.matrix(Expression), gene_set_list_scp, kcdf = "Gaussian")
gsva_result <- gsva(gsvapar)
GSVA_df <- as.data.frame(gsva_result)
GSVA_long_scp <- pivot_longer(GSVA_df, cols = everything(), names_to = "SampleID", values_to = "GSVA_Score_scp")


extendScores <- RunEXTEND(as.matrix(Expression))
telomeraseScores <- read_delim("TelomeraseScores.txt")
telomeraseScores <- telomeraseScores[, c("SampleID", "NormEXTENDScores")]
telomeraseScores <- as.data.frame(telomeraseScores)
telomeraseScores <- left_join(telomeraseScores, metadata[, c("SampleID", "TMM", "COG.Risk.Group")], by = "SampleID")

# joining the two results.
GSVA_long_scp <- left_join(GSVA_long_scp, telomeraseScores, by = "SampleID")

# for scatterplot first computing Pearson correlation value.
cor_val <- cor(GSVA_long_scp$GSVA_Score_scp, GSVA_long_scp$NormEXTENDScores, method = "pearson")
r_text <- paste("r =", round(cor_val, 3))

ggplot(GSVA_long_scp, aes(y = GSVA_Score_scp, x = NormEXTENDScores, colour = TMM)) + 
  scale_color_manual(values = c("ALT" = "darkblue", "Telomerase" = "darkred")) + 
  geom_point(size = 4, alpha = 0.9) + 
  labs(x = "EXTEND Scores", y = "GSVA Signature from SCP Markers") + 
  geom_smooth(method = "lm", level = 0.95, se = TRUE, color = "black", fill = "lightgray") + 
  annotate("text", x = 0.25, y = 1.3, label = r_text, size = 7) +
  theme_classic() +
  theme(
    axis.title.x = element_text(size = 24),
    axis.title.y = element_text(size = 24),
    axis.text.x  = element_text(size = 14, face = "bold"),
    axis.text.y  = element_text(size = 14, face = "bold"),
    legend.title = element_text(size = 18, face = "bold"),
    legend.text  = element_text(size = 16),
    legend.position = c(-0.5, 1.3)
  )

ggplot(GSVA_long_scp, aes(y = GSVA_Score_scp, x = NormEXTENDScores, colour = COG.Risk.Group)) + 
  scale_color_manual(values = c("Intermediate Risk" = "darkblue", "High Risk" = "darkred", "Low Risk" = "gray")) + 
  geom_point(size = 4, alpha = 0.9) + 
  labs(x = "EXTEND Scores", y = "GSVA Signature from SCP Markers") + 
  geom_smooth(method = "lm", level = 0.95, se = TRUE, color = "black", fill = "lightgray") + 
  annotate("text", x = 0.25, y = 1.3, label = r_text, size = 7) +
  theme_classic() +
  theme(
    axis.title.x = element_text(size = 24),
    axis.title.y = element_text(size = 24),
    axis.text.x  = element_text(size = 14, face = "bold"),
    axis.text.y  = element_text(size = 14, face = "bold"),
    legend.title = element_text(size = 18, face = "bold"),
    legend.text  = element_text(size = 16),
    legend.position = c(-0.5, 1.3)
  )


#### EXTEND vs. chromaffin.
source("DataCleaning.R")
gene_set_list_chromaffin <- list("Chromaffin" = c("DBH", "TH", "CHGA", "DDC", "PNMT"))


gsvapar <- gsvaParam(as.matrix(Expression), gene_set_list_chromaffin, kcdf = "Gaussian")
gsva_result <- gsva(gsvapar)
GSVA_df <- as.data.frame(gsva_result)
GSVA_long_chromaffin <- pivot_longer(GSVA_df, cols = everything(), names_to = "SampleID", values_to = "GSVA_Score_chromaffin")


extendScores <- RunEXTEND(as.matrix(Expression))
telomeraseScores <- read_delim("TelomeraseScores.txt")
telomeraseScores <- telomeraseScores[, c("SampleID", "NormEXTENDScores")]
telomeraseScores <- as.data.frame(telomeraseScores)
telomeraseScores <- left_join(telomeraseScores, metadata[, c("SampleID", "TMM", "COG.Risk.Group")], by = "SampleID")

# joining the two results.
GSVA_long_chromaffin <- left_join(GSVA_long_chromaffin, telomeraseScores, by = "SampleID")

# for scatterplot first computing Pearson correlation value.
cor_val <- cor(GSVA_long_chromaffin$GSVA_Score_chromaffin, GSVA_long_chromaffin$NormEXTENDScores, method = "pearson")
r_text <- paste("r =", round(cor_val, 3))

ggplot(GSVA_long_chromaffin, aes(y = GSVA_Score_chromaffin, x = NormEXTENDScores, colour = TMM)) + 
  scale_color_manual(values = c("ALT" = "darkblue", "Telomerase" = "darkred")) + 
  geom_point(size = 4, alpha = 0.9) + 
  labs(x = "EXTEND Scores", y = "GSVA Signature from Chromaffin Markers") + 
  geom_smooth(method = "lm", level = 0.95, se = TRUE, color = "black", fill = "lightgray") + 
  annotate("text", x = 0.25, y = 1.3, label = r_text, size = 7) +
  theme_classic() +
  theme(
    axis.title.x = element_text(size = 24),
    axis.title.y = element_text(size = 24),
    axis.text.x  = element_text(size = 14, face = "bold"),
    axis.text.y  = element_text(size = 14, face = "bold"),
    legend.title = element_text(size = 18, face = "bold"),
    legend.text  = element_text(size = 16),
    legend.position = c(-0.5, 1.3)
  )

ggplot(GSVA_long_chromaffin, aes(y = GSVA_Score_chromaffin, x = NormEXTENDScores, colour = COG.Risk.Group)) + 
  scale_color_manual(values = c("Intermediate Risk" = "darkblue", "High Risk" = "darkred", "Low Risk" = "gray")) + 
  geom_point(size = 4, alpha = 0.9) + 
  labs(x = "EXTEND Scores", y = "GSVA Signature from SCP Markers") + 
  geom_smooth(method = "lm", level = 0.95, se = TRUE, color = "black", fill = "lightgray") + 
  annotate("text", x = 0.25, y = 1.3, label = r_text, size = 7) +
  theme_classic() +
  theme(
    axis.title.x = element_text(size = 24),
    axis.title.y = element_text(size = 24),
    axis.text.x  = element_text(size = 14, face = "bold"),
    axis.text.y  = element_text(size = 14, face = "bold"),
    legend.title = element_text(size = 18, face = "bold"),
    legend.text  = element_text(size = 16),
    legend.position = c(-0.5, 1.3)
  )


#### NO_TMM gsva vs. chromaffin.
source("DataCleaning.R")
gene_set_list_chromaffin <- list("Chromaffin" = c("DBH", "TH", "CHGA", "DDC", "PNMT"))

gsvapar <- gsvaParam(as.matrix(Expression), gene_set_list_chromaffin, kcdf = "Gaussian")
gsva_result <- gsva(gsvapar)
GSVA_df <- as.data.frame(gsva_result)
GSVA_long_chromaffin <- pivot_longer(GSVA_df, cols = everything(), names_to = "SampleID", values_to = "GSVA_Score_chromaffin")

gene_list_noTMM <- list("NO_TMM Genes" = c("CPNE8", "PGM2L1", "CNR1", "LIFR", "HOXC9", "SNX16", "HECW2", "ALDH3A2",
                                           "THSD7A", "CPNE3", "IGSF10"))
Expression <- as.matrix(Expression)
noTMM_gsva <- gsvaParam(Expression, gene_list_noTMM, kcdf = "Gaussian")
GSVA_result_noTMM <- gsva(noTMM_gsva)
GSVA_df_noTMM <- as.data.frame(GSVA_result_noTMM)
GSVA_long_noTMM <- pivot_longer(GSVA_df_noTMM, cols = everything(), names_to = "SampleID", values_to = "GSVA_Score")


# joining the two results.
GSVA_long_chromaffin <- left_join(GSVA_long_chromaffin, GSVA_long_noTMM, by = "SampleID")
GSVA_long_chromaffin <- left_join(GSVA_long_chromaffin, metadata[, c("SampleID", "COG.Risk.Group", "TMM")], by = "SampleID")

# for scatterplot first computing Pearson correlation value.
cor_val <- cor(GSVA_long_chromaffin$GSVA_Score_chromaffin, GSVA_long_chromaffin$GSVA_Score, method = "pearson")
r_text <- paste("r =", round(cor_val, 3))

ggplot(GSVA_long_chromaffin, aes(y = GSVA_Score_chromaffin, x = GSVA_Score, colour = TMM)) + 
  scale_color_manual(values = c("ALT" = "darkblue", "Telomerase" = "darkred")) + 
  geom_point(size = 4, alpha = 0.9) + 
  labs(x = "GSVA Score from NO_TMM Markers", y = "GSVA Signature from Chromaffin Markers") + 
  geom_smooth(method = "lm", level = 0.95, se = TRUE, color = "black", fill = "lightgray") + 
  annotate("text", x = 0.25, y = 1.3, label = r_text, size = 7) +
  theme_classic() +
  theme(
    axis.title.x = element_text(size = 24),
    axis.title.y = element_text(size = 24),
    axis.text.x  = element_text(size = 14, face = "bold"),
    
    axis.text.y  = element_text(size = 14, face = "bold"),
    legend.title = element_text(size = 18, face = "bold"),
    legend.text  = element_text(size = 16),
    legend.position = c(-0.5, 1.3)
  )

## EXTEND vs. MES GSVA.
source("DataCleaning.R")


extendScores <- RunEXTEND(as.matrix(Expression))
telomeraseScores <- read_delim("TelomeraseScores.txt")
telomeraseScores <- telomeraseScores[, c("SampleID", "NormEXTENDScores")]
telomeraseScores <- as.data.frame(telomeraseScores)
telomeraseScores <- left_join(telomeraseScores, metadata[, c("SampleID", "TMM", "COG.Risk.Group")], by = "SampleID")

gene_list_mesenchymal <- list("MES" = c("VIM", "FN1", "YAP1", "SNAI2", "PRRX1", "WWTR1"))


gsvapar <- gsvaParam(as.matrix(Expression), gene_list_mesenchymal, kcdf = "Gaussian")
gsva_result <- gsva(gsvapar)
GSVA_df <- as.data.frame(gsva_result)
GSVA_long_mes <- pivot_longer(GSVA_df, cols = everything(), names_to = "SampleID", values_to = "GSVA_Score_Mes")

# joining the two results.
GSVA_long_mes <- left_join(GSVA_long_mes, telomeraseScores, by = "SampleID")

# for scatterplot first computing Pearson correlation value.
cor_val <- cor(GSVA_long_mes$GSVA_Score_Mes, GSVA_long_mes$NormEXTENDScores, method = "pearson")
r_text <- paste("r =", round(cor_val, 3))

ggplot(GSVA_long_mes, aes(y = GSVA_Score_Mes, x = NormEXTENDScores, colour = TMM)) + 
  scale_color_manual(values = c("ALT" = "darkblue", "Telomerase" = "darkred")) + 
  geom_point(size = 4, alpha = 0.9) + 
  labs(y = "GSVA Score from MES Markers", x = "EXTEND Scores") + 
  geom_smooth(method = "lm", level = 0.95, se = TRUE, color = "black", fill = "lightgray") + 
  annotate("text", x = 0.25, y = 1.3, label = r_text, size = 7) +
  theme_classic() +
  theme(
    axis.title.x = element_text(size = 24),
    axis.title.y = element_text(size = 24),
    axis.text.x  = element_text(size = 14, face = "bold"),
    
    axis.text.y  = element_text(size = 14, face = "bold"),
    legend.title = element_text(size = 18, face = "bold"),
    legend.text  = element_text(size = 16),
    legend.position = c(-0.5, 1.3)
  )


## adrenergic vs. NO_TMM GSVA.
source("DataCleaning.R")
gene_list_adrenergic <- list("ADR" = c("DLK1", "DBH", "PHOX2A", "PHOX2B", "GATA2", "GATA3"))

gene_list_noTMM <- list("NO_TMM Genes" = c("CPNE8", "PGM2L1", "CNR1", "LIFR", "HOXC9", "SNX16", "HECW2", "ALDH3A2",
                                           "THSD7A", "CPNE3", "IGSF10"))

Expression <- as.matrix(Expression)
noTMM_gsva <- gsvaParam(Expression, gene_list_noTMM, kcdf = "Gaussian")
GSVA_result_noTMM <- gsva(noTMM_gsva)
GSVA_df_noTMM <- as.data.frame(GSVA_result_noTMM)
GSVA_long_noTMM <- pivot_longer(GSVA_df_noTMM, cols = everything(), names_to = "SampleID", values_to = "GSVA_Score")


gsvapar <- gsvaParam(as.matrix(Expression), gene_list_adrenergic, kcdf = "Gaussian")
gsva_result <- gsva(gsvapar)
GSVA_df <- as.data.frame(gsva_result)
GSVA_long_adr <- pivot_longer(GSVA_df, cols = everything(), names_to = "SampleID", values_to = "GSVA_Score_Adr")

# joining the two results.
GSVA_long_adr <- left_join(GSVA_long_adr, GSVA_long_noTMM, by = "SampleID")
GSVA_long_adr <- left_join(GSVA_long_adr, metadata[, c("SampleID", "COG.Risk.Group", "TMM")], by = "SampleID")

# for scatterplot first computing Pearson correlation value.
cor_val <- cor(GSVA_long_adr$GSVA_Score_Adr, GSVA_long_adr$GSVA_Score, method = "pearson")
r_text <- paste("r =", round(cor_val, 3))

ggplot(GSVA_long_adr, aes(y = GSVA_Score_Adr, x = GSVA_Score, colour = TMM)) + 
  scale_color_manual(values = c("ALT" = "darkblue", "Telomerase" = "darkred")) + 
  geom_point(size = 4, alpha = 0.9) + 
  labs(y = "GSVA Score from ADR Markers", x = "NO_TMM GSVA Scores") + 
  geom_smooth(method = "lm", level = 0.95, se = TRUE, color = "black", fill = "lightgray") + 
  annotate("text", x = 0.25, y = 1.3, label = r_text, size = 7) +
  theme_classic() +
  theme(
    axis.title.x = element_text(size = 24),
    axis.title.y = element_text(size = 24),
    axis.text.x  = element_text(size = 14, face = "bold"),
    
    axis.text.y  = element_text(size = 14, face = "bold"),
    legend.title = element_text(size = 18, face = "bold"),
    legend.text  = element_text(size = 16),
    legend.position = c(-0.5, 1.3)
  )

## adrenergic vs. chromaffin.
source("DataCleaning.R")
gene_list_adrenergic <- list("ADR" = c("DLK1", "DBH", "PHOX2A", "PHOX2B", "GATA2", "GATA3"))
gene_set_list_chromaffin <- list("Chromaffin" = c("DBH", "TH", "CHGA", "DDC", "PNMT"))

gsvapar <- gsvaParam(as.matrix(Expression), gene_list_adrenergic, kcdf = "Gaussian")
gsva_result <- gsva(gsvapar)
GSVA_df <- as.data.frame(gsva_result)
GSVA_long_adr <- pivot_longer(GSVA_df, cols = everything(), names_to = "SampleID", values_to = "GSVA_Score_Adr")

gsvapar <- gsvaParam(as.matrix(Expression), gene_set_list_chromaffin, kcdf = "Gaussian")
gsva_result <- gsva(gsvapar)
GSVA_df <- as.data.frame(gsva_result)
GSVA_long_chromaffin <- pivot_longer(GSVA_df, cols = everything(), names_to = "SampleID", values_to = "GSVA_Score_chromaffin")

# joining the two results.
GSVA_long_adr <- left_join(GSVA_long_adr, GSVA_long_chromaffin, by = "SampleID")
GSVA_long_adr <- left_join(GSVA_long_adr, metadata[, c("SampleID", "COG.Risk.Group", "TMM")], by = "SampleID")

# for scatterplot first computing Pearson correlation value.
cor_val <- cor(GSVA_long_adr$GSVA_Score_Adr, GSVA_long_adr$GSVA_Score_chromaffin, method = "pearson")
r_text <- paste("r =", round(cor_val, 3))

ggplot(GSVA_long_adr, aes(y = GSVA_Score_Adr, x = GSVA_Score_chromaffin, colour = TMM)) + 
  scale_color_manual(values = c("ALT" = "darkblue", "Telomerase" = "darkred")) + 
  geom_point(size = 4, alpha = 0.9) + 
  labs(y = "GSVA Score from ADR Markers", x = "NO_TMM GSVA Scores") + 
  geom_smooth(method = "lm", level = 0.95, se = TRUE, color = "black", fill = "lightgray") + 
  annotate("text", x = 0.25, y = 1.3, label = r_text, size = 7) +
  theme_classic() +
  theme(
    axis.title.x = element_text(size = 24),
    axis.title.y = element_text(size = 24),
    axis.text.x  = element_text(size = 14, face = "bold"),
    
    axis.text.y  = element_text(size = 14, face = "bold"),
    legend.title = element_text(size = 18, face = "bold"),
    legend.text  = element_text(size = 16),
    legend.position = c(-0.5, 1.3))



ggplot(GSVA_long_adr, aes(y = GSVA_Score_Adr, x = GSVA_Score_chromaffin, colour = COG.Risk.Group)) + 
  scale_color_manual(values = c("Intermediate Risk" = "darkblue", "High Risk" = "darkred")) + 
  geom_point(size = 4, alpha = 0.9) + 
  labs(y = "GSVA Score from ADR Markers", x = "Chromaffin GSVA Scores") + 
  geom_smooth(method = "lm", level = 0.95, se = TRUE, color = "black", fill = "lightgray") + 
  annotate("text", x = 0.25, y = 1.3, label = r_text, size = 7) +
  theme_classic() +
  theme(
    axis.title.x = element_text(size = 24),
    axis.title.y = element_text(size = 24),
    axis.text.x  = element_text(size = 14, face = "bold"),
    
    axis.text.y  = element_text(size = 14, face = "bold"),
    legend.title = element_text(size = 18, face = "bold"),
    legend.text  = element_text(size = 16),
    legend.position = c(-0.5, 1.3))


## adrenergic vs. scp.
source("DataCleaning.R")

gene_list_adrenergic <- list("ADR" = c("DLK1", "DBH", "PHOX2A", "PHOX2B", "GATA2", "GATA3"))
gene_set_list_scp <- list("SCP" = c("PLP1", "MPZ", "CDH19", "ERBB3", "ERBB4"))

gsvapar <- gsvaParam(as.matrix(Expression), gene_list_adrenergic, kcdf = "Gaussian")
gsva_result <- gsva(gsvapar)
GSVA_df <- as.data.frame(gsva_result)
GSVA_long_adr <- pivot_longer(GSVA_df, cols = everything(), names_to = "SampleID", values_to = "GSVA_Score_Adr")

gsvapar <- gsvaParam(as.matrix(Expression), gene_set_list_scp, kcdf = "Gaussian")
gsva_result <- gsva(gsvapar)
GSVA_df <- as.data.frame(gsva_result)
GSVA_long_scp <- pivot_longer(GSVA_df, cols = everything(), names_to = "SampleID", values_to = "GSVA_Score_scp")

# joining the two results.
GSVA_long_adr <- left_join(GSVA_long_adr, GSVA_long_scp, by = "SampleID")
GSVA_long_adr <- left_join(GSVA_long_adr, metadata[, c("SampleID", "COG.Risk.Group", "TMM")], by = "SampleID")

# for scatterplot first computing Pearson correlation value.
cor_val <- cor(GSVA_long_adr$GSVA_Score_Adr, GSVA_long_adr$GSVA_Score_scp, method = "pearson")
r_text <- paste("r =", round(cor_val, 3))

ggplot(GSVA_long_adr, aes(y = GSVA_Score_Adr, x = GSVA_Score_scp, colour = TMM)) + 
  scale_color_manual(values = c("ALT" = "darkblue", "Telomerase" = "darkred")) + 
  geom_point(size = 4, alpha = 0.9) + 
  labs(y = "GSVA Score from ADR Markers", x = "SCP GSVA Scores") + 
  geom_smooth(method = "lm", level = 0.95, se = TRUE, color = "black", fill = "lightgray") + 
  annotate("text", x = 0.25, y = 1.3, label = r_text, size = 7) +
  theme_classic() +
  theme(
    axis.title.x = element_text(size = 24),
    axis.title.y = element_text(size = 24),
    axis.text.x  = element_text(size = 14, face = "bold"),
    
    axis.text.y  = element_text(size = 14, face = "bold"),
    legend.title = element_text(size = 18, face = "bold"),
    legend.text  = element_text(size = 16),
    legend.position = c(-0.5, 1.3))


ggplot(GSVA_long_adr, aes(y = GSVA_Score_Adr, x = GSVA_Score_chromaffin, colour = TMM)) + 
  scale_color_manual(values = c("ALT" = "darkblue", "Telomerase" = "darkred")) + 
  geom_point(size = 4, alpha = 0.9) + 
  labs(y = "GSVA Score from ADR Markers", x = "Chromaffin GSVA Scores") + 
  geom_smooth(method = "lm", level = 0.95, se = TRUE, color = "black", fill = "lightgray") + 
  annotate("text", x = 0.25, y = 1.3, label = r_text, size = 7) +
  theme_classic() +
  theme(
    axis.title.x = element_text(size = 24),
    axis.title.y = element_text(size = 24),
    axis.text.x  = element_text(size = 14, face = "bold"),
    
    axis.text.y  = element_text(size = 14, face = "bold"),
    legend.title = element_text(size = 18, face = "bold"),
    legend.text  = element_text(size = 16),
    legend.position = c(-0.5, 1.3))


# adrenergic vs. extend score.

source("DataCleaning.R")


extendScores <- RunEXTEND(as.matrix(Expression))
telomeraseScores <- read_delim("TelomeraseScores.txt")
telomeraseScores <- telomeraseScores[, c("SampleID", "NormEXTENDScores")]
telomeraseScores <- as.data.frame(telomeraseScores)
telomeraseScores <- left_join(telomeraseScores, metadata[, c("SampleID", "TMM", "COG.Risk.Group")], by = "SampleID")

gene_list_adrenergic <- list("MES" = c("DLK1", "DBH", "PHOX2A", "PHOX2B", "GATA2", "GATA3"))


gsvapar <- gsvaParam(as.matrix(Expression), gene_list_adrenergic, kcdf = "Gaussian")
gsva_result <- gsva(gsvapar)
GSVA_df <- as.data.frame(gsva_result)
GSVA_long_adr <- pivot_longer(GSVA_df, cols = everything(), names_to = "SampleID", values_to = "GSVA_Score_Adr")

# joining the two results.
GSVA_long_adr <- left_join(GSVA_long_adr, telomeraseScores, by = "SampleID")

# for scatterplot first computing Pearson correlation value.
cor_val <- cor(GSVA_long_adr$GSVA_Score_Adr, GSVA_long_adr$NormEXTENDScores, method = "pearson")
r_text <- paste("r =", round(cor_val, 3))

ggplot(GSVA_long_adr, aes(y = GSVA_Score_Adr, x = NormEXTENDScores, colour = TMM)) + 
  scale_color_manual(values = c("ALT" = "darkblue", "Telomerase" = "darkred")) + 
  geom_point(size = 4, alpha = 0.9) + 
  labs(y = "GSVA Score from ADR Markers", x = "EXTEND Scores") + 
  geom_smooth(method = "lm", level = 0.95, se = TRUE, color = "black", fill = "lightgray") + 
  annotate("text", x = 0.25, y = 1.3, label = r_text, size = 7) +
  theme_classic() +
  theme(
    axis.title.x = element_text(size = 24),
    axis.title.y = element_text(size = 24),
    axis.text.x  = element_text(size = 14, face = "bold"),
    
    axis.text.y  = element_text(size = 14, face = "bold"),
    legend.title = element_text(size = 18, face = "bold"),
    legend.text  = element_text(size = 16),
    legend.position = c(-0.5, 1.3)
  )


## adrenergic vs. NO_TMM GSVA.
source("DataCleaning.R")
gene_list_adrenergic <- list("ADR" = c("DLK1", "DBH", "PHOX2A", "PHOX2B", "GATA2", "GATA3"))

gene_list_noTMM <- list("NO_TMM Genes" = c("CPNE8", "PGM2L1", "CNR1", "LIFR", "HOXC9", "SNX16", "HECW2", "ALDH3A2",
                                           "THSD7A", "CPNE3", "IGSF10"))

Expression <- as.matrix(Expression)
noTMM_gsva <- gsvaParam(Expression, gene_list_noTMM, kcdf = "Gaussian")
GSVA_result_noTMM <- gsva(noTMM_gsva)
GSVA_df_noTMM <- as.data.frame(GSVA_result_noTMM)
GSVA_long_noTMM <- pivot_longer(GSVA_df_noTMM, cols = everything(), names_to = "SampleID", values_to = "GSVA_Score")


gsvapar <- gsvaParam(as.matrix(Expression), gene_list_adrenergic, kcdf = "Gaussian")
gsva_result <- gsva(gsvapar)
GSVA_df <- as.data.frame(gsva_result)
GSVA_long_adr <- pivot_longer(GSVA_df, cols = everything(), names_to = "SampleID", values_to = "GSVA_Score_Adr")

# joining the two results.
GSVA_long_adr <- left_join(GSVA_long_adr, GSVA_long_noTMM, by = "SampleID")
GSVA_long_adr <- left_join(GSVA_long_adr, metadata[, c("SampleID", "COG.Risk.Group", "TMM")], by = "SampleID")

# for scatterplot first computing Pearson correlation value.
cor_val <- cor(GSVA_long_adr$GSVA_Score_Adr, GSVA_long_adr$GSVA_Score, method = "pearson")
r_text <- paste("r =", round(cor_val, 3))

ggplot(GSVA_long_adr, aes(y = GSVA_Score_Adr, x = GSVA_Score, colour = TMM)) + 
  scale_color_manual(values = c("ALT" = "darkblue", "Telomerase" = "darkred")) + 
  geom_point(size = 4, alpha = 0.9) + 
  labs(y = "GSVA Score from ADR Markers", x = "NO_TMM GSVA Scores") + 
  geom_smooth(method = "lm", level = 0.95, se = TRUE, color = "black", fill = "lightgray") + 
  annotate("text", x = 0.25, y = 1.3, label = r_text, size = 7) +
  theme_classic() +
  theme(
    axis.title.x = element_text(size = 24),
    axis.title.y = element_text(size = 24),
    axis.text.x  = element_text(size = 14, face = "bold"),
    
    axis.text.y  = element_text(size = 14, face = "bold"),
    legend.title = element_text(size = 18, face = "bold"),
    legend.text  = element_text(size = 16),
    legend.position = c(-0.5, 1.3)
  )

ggplot(GSVA_long_adr, aes(y = GSVA_Score_Adr, x = GSVA_Score, colour = COG.Risk.Group)) + 
  scale_color_manual(values = c("Intermediate Risk" = "darkblue", "High Risk" = "darkred")) + 
  geom_point(size = 4, alpha = 0.9) + 
  labs(y = "GSVA Score from ADR Markers", x = "NO_TMM GSVA Scores") + 
  geom_smooth(method = "lm", level = 0.95, se = TRUE, color = "black", fill = "lightgray") + 
  annotate("text", x = 0.25, y = 1.3, label = r_text, size = 7) +
  theme_classic() +
  theme(
    axis.title.x = element_text(size = 24),
    axis.title.y = element_text(size = 24),
    axis.text.x  = element_text(size = 14, face = "bold"),
    
    axis.text.y  = element_text(size = 14, face = "bold"),
    legend.title = element_text(size = 18, face = "bold"),
    legend.text  = element_text(size = 16),
    legend.position = c(-0.5, 1.3)
  )



#######################################################################################



### ADRN score vs. MES score.
source("DataCleaning.R")



gene_sets_phenotype <- list(
  ADR = c("DLK1", "DBH", "PHOX2A", "PHOX2B", "GATA2", "GATA3"),
  MES = c("VIM", "FN1", "YAP1", "SNAI2", "PRRX1", "WWTR1")
)

# gene_sets_phenotype <- list(
#   ADR = c("KLF7", "GATA3", "HAND2", "PHOX2A", "ISL1", "HAND1",
#           "PHOX2B", "TFAP2B", "GATA2", "SATB1", "SIX3", "EYA1",
#           "SOX11", "DACH1", "ASCL1", "HEY1", "KLF13", "PBX3"),
#   MES = c("VIM", "FN1", "MEOX2", "ID1", "EGR3", "AEBP1", "CBFB", "IRF3", "IRF2", "IRF1", 
#           "TBX18", "MAFF", "RUNX2", "ZFP36L1", "NR3C1", "BHLHE41", "GLIS3", "RUNX1", "FOSL1",
#           "FOSL2", "ELK4", "IFI16", "SIX4", "FLI1", "MAML2", "SMAD3", "DCAF6", "WWTR1", "SOX9", 
#           "MEF2D", "ZNF217", "PRRX1", "CREG1", "NOTCH2", "SIX1", "MEOX1"))


gsvapar <- gsvaParam(as.matrix(Expression), gene_sets_phenotype, kcdf = "Gaussian")
gsva_result <- gsva(gsvapar)
GSVA_df <- as.data.frame(gsva_result)
GSVA_df <- as.data.frame(t(GSVA_df))


# for scatterplot first computing Pearson correlation value.
cor_val <- cor(GSVA_df$ADR, GSVA_df$MES, method = "pearson")
r_text <- paste("r =", round(cor_val, 3))

ggplot(GSVA_df, aes(y = MES, x = ADR)) + 
  geom_point(size = 4, alpha = 0.9) + 
  labs(y = "GSVA Score from MES Markers", x = "GSVA Score from ADR Markers") + 
  geom_smooth(method = "lm", level = 0.95, se = TRUE, color = "black", fill = "lightgray") + 
  annotate("text", x = 0.25, y = 1.3, label = r_text, size = 7) +
  theme_classic() +
  theme(
    axis.title.x = element_text(size = 24),
    axis.title.y = element_text(size = 24),
    axis.text.x  = element_text(size = 14, face = "bold"),
    
    axis.text.y  = element_text(size = 14, face = "bold"),
    legend.title = element_text(size = 18, face = "bold"),
    legend.text  = element_text(size = 16),
    legend.position = c(-0.5, 1.3))



