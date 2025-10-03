### dge-limma again.


library(EnhancedVolcano)
library(tidyverse)
library(dplyr)
library(ggpubr)
library(rstatix)
library(GSVA)
library(cmapR)


##### 1) differential expression: NO_TMM vs. ALT.
source("DataCleaning.R")
source("0532Data.R")

# only high risk samples and only genes present in both TARGET and 0532 set.
metadata <- metadata[metadata$COG.Risk.Group == "High Risk", ]
Expression <- Expression[rownames(Expression) %in% rownames(tmm_lcpm), ]

metadata_ALT <- metadata[metadata$TMM == "ALT" | metadata$TMM == "NO_TMM", ]

# creating design matrix for limma.
metadata_ALT$TMM <- as.factor(metadata_ALT$TMM)
metadata_ALT <- metadata_ALT %>%
  arrange(TMM)

design <- model.matrix(~ 0 + TMM, data = metadata_ALT)
colnames(design) <- levels(metadata_ALT$TMM)


# Estimating array weights to account for sample-specific variability in library sizes.
ExpressionALT <- Expression[, colnames(Expression) %in% metadata_ALT$SampleID]
ExpressionALT <- ExpressionALT[, match(metadata_ALT$SampleID, colnames(ExpressionALT))]

weights <- arrayWeights(ExpressionALT, design = design)

# Fitting linear model using limma.
fit <- lmFit(ExpressionALT, design, weights = weights)
fit <- eBayes(fit, trend = TRUE)

# Contrasts: NO_TMM vs ALT && ALT vs NO_TMM.
contrast.matrix <- makeContrasts(
  NO_TMMvsALT = NO_TMM-ALT,
  ALTvsNO_TMM = ALT-NO_TMM,
  levels = design
)

fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2, trend = TRUE)


# top results for ALT.
ALT_results <- topTable(fit2, coef = "ALTvsNO_TMM", number = Inf, adjust = "fdr")


# Filtering for FDR < 0.01.
ALT_candidates <- ALT_results[round(ALT_results$adj.P.Val, 2) <= 0.05 & (ALT_results$logFC >= 0.5 | ALT_results$logFC <= -0.5), ]


##### 2) t-test of the DEGs.
# Initializing an empty data frame for storing t-test results.
ALT_candidates <- rownames_to_column(ALT_candidates, var = "Gene")
dge_gene <- ExpressionALT[rownames(ExpressionALT) %in% ALT_candidates$Gene, ,drop = FALSE]

t_test_results_ALT <- data.frame(
  Gene = character(),
  p_value_t_test = numeric(),
  stringsAsFactors = FALSE
)

# all p-values for FDR adjustment later.
all_p_values <- numeric()


# Looping through each gene in dge_gene.
for (i in 1:nrow(ALT_candidates)) {
  
  gene_id <- ALT_candidates$Gene[i]
  
  
  # Extracting the expression values for this gene, keeping it as a matrix. We are only selecting significant genes (all samples for this gene).
  gene_Expression <- dge_gene[rownames(dge_gene) == gene_id, , drop = FALSE]
  
  # Creating C1 and C2 group. For the gene we are working with, all NO_TMM sample expression data placed in c1_group and all ALT sample expression data placed in c2_group.
  c1_group <- gene_Expression[, metadata_ALT$TMM == "NO_TMM", drop = FALSE]
  c2_group <- gene_Expression[, metadata_ALT$TMM == "ALT", drop = FALSE]
  
  # for a gene, compare c1 and c2 samples.
  t_test <- t.test(c1_group, c2_group)
  p_value_t_test <- t_test$p.value
  all_p_values <- c(all_p_values, p_value_t_test)
  
  
  # Storing the results.
  t_test_results_ALT <- rbind(t_test_results_ALT, data.frame(
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


# Calculating the FDR-adjusted p-values
t_test_results_ALT$fdr_t_test <- p.adjust(all_p_values, method = "BH")


# Storing significant t-test results where FDR is less than 0.05.
t_test_results_sig_ALT <- t_test_results_ALT[round(t_test_results_ALT$fdr_t_test, 2) <= 0.01, ]



##### 3) Linear regression test of candidate genes.
regression_results_ALT <- data.frame(
  Gene = character(),
  estimate = numeric(),
  p_value = numeric(),
  R.squared = numeric(),
  stringsAsFactors = FALSE
)

# for each gene in the candidate list, updating regression_results.
for (gene in t_test_results_sig_ALT$Gene) {
  
  data <- data.frame(gene_Expression = as.numeric(dge_gene[gene, ]), 
                     TMMstatus = metadata_ALT$TMM)
  
  # Fit linear model.
  model <- lm(gene_Expression ~ TMMstatus, data = data)
  summary_model <- summary(model)
  
  
  # Store regression results.
  regression_results_ALT <- rbind(regression_results_ALT, data.frame(
    Gene = gene,
    estimate = summary_model$coefficients[2, 1],
    p_value = summary_model$coefficients[2, 4],
    R.Squared = summary_model$r.squared
  ))
  
  # removing intermediates.
  rm(data)
  rm(summary_model)
  rm(gene)
  rm(estimate)
  rm(p_value)
  rm(R.squared)
}

# adjusted p-values.
regression_results_ALT$Adj.P.Value <- p.adjust(regression_results_ALT$p_value, method = "fdr")

# Filtering for r-squared > 0.3 & fdr < 0.01.
regression_results_sig_ALT <- regression_results_ALT %>%
  filter(round(Adj.P.Value, 2) <= 0.01 & round(R.Squared, 1) >= 0.2)



##########################################################################################



##### 1) differential expression: NO_TMM vs. Telomerase.

source("DataCleaning.R")
source("0532Data.R")

# only high risk samples and only genes present in both TARGET and 0532 set.
metadata <- metadata[metadata$COG.Risk.Group == "High Risk", ]
Expression <- Expression[rownames(Expression) %in% rownames(tmm_lcpm), ]



metadata_Telomerase <- metadata[metadata$TMM == "Telomerase" | metadata$TMM == "NO_TMM", ]

# creating design matrix for limma.
metadata_Telomerase$TMM <- as.factor(metadata_Telomerase$TMM)

metadata_Telomerase <- metadata_Telomerase %>%
  arrange(TMM)

# Estimating array weights to account for sample-specific variability in library sizes.
ExpressionTelomerase <- Expression[, colnames(Expression) %in% metadata_Telomerase$SampleID]
ExpressionTelomerase <- ExpressionTelomerase[, match(metadata_Telomerase$SampleID, colnames(ExpressionTelomerase))]


design <- model.matrix(~ 0 + TMM, data = metadata_Telomerase)
colnames(design) <- levels(metadata_Telomerase$TMM)




weights <- arrayWeights(ExpressionTelomerase, design = design)

# Fitting linear model using limma.
fit <- lmFit(ExpressionTelomerase, design, weights = weights)
fit <- eBayes(fit, trend = TRUE)

# Contrasts: NO_TMM vs ALT && ALT vs NO_TMM.
contrast.matrix <- makeContrasts(
  NO_TMMvsTelomerase = NO_TMM - Telomerase,
  TelomerasevsNO_TMM = Telomerase - NO_TMM,
  levels = design
)

fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2, trend = TRUE)


# top results for Telomerase.
telomerase_results <- topTable(fit2, coef = "TelomerasevsNO_TMM", number = Inf, adjust = "fdr")

# Filtering for p-value < 0.01.
telomerase_candidates <- telomerase_results[round(telomerase_results$adj.P.Val, 2) <= 0.05 & (telomerase_results$logFC >= 0.5 | telomerase_results$logFC <= -0.5), ]



##### 2) t-test of the DEGs.
# Initializing an empty data frame for storing t-test results.
telomerase_candidates <- rownames_to_column(telomerase_candidates, var = "Gene")
dge_gene <- ExpressionTelomerase[rownames(ExpressionTelomerase) %in% telomerase_candidates$Gene, ,drop = FALSE]

t_test_results_telomerase <- data.frame(
  Gene = character(),
  p_value_t_test = numeric(),
  stringsAsFactors = FALSE
)

# all p-values for FDR adjustment later.
all_p_values <- numeric()


# Looping through each gene in dge_gene.
for (i in 1:nrow(telomerase_candidates)) {
  
  gene_id <- telomerase_candidates$Gene[i]
  
  
  # Extracting the expression values for this gene, keeping it as a matrix. We are only selecting significant genes (all samples for this gene).
  gene_Expression <- dge_gene[rownames(dge_gene) == gene_id, , drop = FALSE]
  
  # Creating C1 and C2 group. For the gene we are working with, all NO_TMM sample expression data placed in c1_group and all Telomerase sample expression data placed in c2_group.
  c1_group <- gene_Expression[, metadata_Telomerase$TMM == "NO_TMM", drop = FALSE]
  c2_group <- gene_Expression[, metadata_Telomerase$TMM == "Telomerase", drop = FALSE]
  
  # for a gene, compare c1 and c2 samples.
  t_test <- t.test(c1_group, c2_group)
  p_value_t_test <- t_test$p.value
  all_p_values <- c(all_p_values, p_value_t_test)
  
  
  # Storing the results.
  t_test_results_telomerase <- rbind(t_test_results_telomerase, data.frame(
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

# Calculating the FDR-adjusted p-values
t_test_results_telomerase$fdr_t_test <- p.adjust(all_p_values, method = "BH")


# Storing significant t-test results where p-value is less than 0.01.
t_test_results_sig_telomerase <- t_test_results_telomerase[round(t_test_results_telomerase$fdr_t_test, 2) <= 0.01, ]


##### 3) Linear regression test of candidate genes.
regression_results_telomerase <- data.frame(
  Gene = character(),
  estimate = numeric(),
  p_value = numeric(),
  R.squared = numeric(),
  stringsAsFactors = FALSE
)

# for each gene in the candidate list, updating regression_results.
for (gene in rownames(dge_gene)) {
  
  data <- data.frame(gene_Expression = as.numeric(dge_gene[gene, ]), 
                     TMMstatus = metadata_Telomerase$TMM)
  
  # Fit linear model.
  model <- lm(gene_Expression ~ TMMstatus, data = data)
  summary_model <- summary(model)
  
  
  # Store regression results.
  regression_results_telomerase <- rbind(regression_results_telomerase, data.frame(
    Gene = gene,
    estimate = summary_model$coefficients[2, 1],
    p_value = summary_model$coefficients[2, 4],
    R.Squared = summary_model$r.squared
  ))
  
  # removing intermediates.
  rm(data)
  rm(summary_model)
  rm(gene)
  rm(estimate)
  rm(p_value)
  rm(R.squared)
}

# adjusted p-values.
regression_results_telomerase$Adj.P.Value <- p.adjust(regression_results_telomerase$p_value, method = "fdr")

# Filtering for r-squared >= 0.3 & fdr <= 0.01.
regression_results_sig_telomerase <- regression_results_telomerase %>%
  filter(round(Adj.P.Value, 2) <= 0.01 & round(R.Squared, 1) >= 0.2)



##### 

# combining together.

source("DataCleaning.R")
source("0532Data.R")

# only high risk samples and only genes present in both TARGET and 0532 set.
metadata <- metadata[metadata$COG.Risk.Group == "High Risk", ]
Expression <- Expression[rownames(Expression) %in% rownames(tmm_lcpm), ]

metadata <- metadata %>%
  arrange(TMM_Case, TMM)
Expression <- Expression[, match(metadata$SampleID, colnames(Expression))]

candidates <- rbind(regression_results_sig_telomerase, regression_results_sig_ALT)

candidates_tmm_upregulated_ALT <- candidates[candidates$Gene %in% ALT_candidates[ALT_candidates$logFC > 0, ]$Gene, ]
candidates_tmm_upregulated_ALT <- candidates_tmm_upregulated_ALT %>%
  arrange(Adj.P.Value)

candidates_tmm_upregulated_telomerase <- candidates[candidates$Gene %in% telomerase_candidates[telomerase_candidates$logFC > 0, ]$Gene, ]
candidates_tmm_upregulated_telomerase <- candidates_tmm_upregulated_telomerase %>%
  arrange(Adj.P.Value)



candidates_notmm_upregulated_ALT <- candidates[candidates$Gene %in% ALT_candidates[ALT_candidates$logFC < 0, ]$Gene, ]
candidates_notmm_upregulated_ALT <- candidates_notmm_upregulated_ALT %>%
  arrange(Adj.P.Value)


candidates_notmm_upregulated_telomerase <- candidates[candidates$Gene %in% telomerase_candidates[telomerase_candidates$logFC < 0, ]$Gene, ]
candidates_notmm_upregulated_telomerase <- candidates_notmm_upregulated_telomerase %>%
  arrange(Adj.P.Value)



#####

# finding best signature.

# signature for genes upregulated in ALT over no_tmm.
tmm_target1 <- run_exhaustive_forward_auc(expr_matrix = ExpressionALT, metadata = metadata_ALT, 
                                         candidate_genes = candidates_tmm_upregulated_ALT$Gene, 
                                         phenotype_col = "TMM",
                                         label_one = "ALT", label_two = "NO_TMM",
                                         max_genes = 20,
                                         pivot_gene = "LCN15")

# signature for genes upregulated in Telomerase over no_tmm.
tmm_target2 <- run_exhaustive_forward_auc(expr_matrix = ExpressionTelomerase, metadata = metadata_Telomerase, 
                                         candidate_genes = candidates_tmm_upregulated_telomerase$Gene, 
                                         phenotype_col = "TMM",
                                         label_one = "Telomerase", label_two = "NO_TMM",
                                         max_genes = 20,
                                         pivot_gene = "TERT")


# signature for genes upregulated in no_tmm over ALT.
tmm_target3 <- run_exhaustive_forward_auc(expr_matrix = ExpressionALT, metadata = metadata_ALT, 
                                         candidate_genes = candidates_notmm_upregulated_ALT$Gene, 
                                         phenotype_col = "TMM",
                                         label_one = "NO_TMM", label_two = "ALT",
                                         max_genes = 20,
                                         pivot_gene = "ERICH5")

# signature for genes upregulated in no_tmm over Telomerase.
tmm_target4 <- run_exhaustive_forward_auc(expr_matrix = ExpressionTelomerase, metadata = metadata_Telomerase, 
                                          candidate_genes = candidates_notmm_upregulated_telomerase$Gene, 
                                          phenotype_col = "TMM",
                                          label_one = "NO_TMM", label_two = "Telomerase",
                                          max_genes = 20,
                                          pivot_gene = "ERICH5")

# genes positively upregulated in TMM.
candidate_genes <- list(TMM = c("TERT", "TUT1", "PLCXD1", "PRR22", "ZDHHC11B", "SLC38A5", "ECSIT", "LCN15", "CRIP1", "SLCO1A2"))
gsvapar <- gsvaParam(as.matrix(Expression), candidate_genes, kcdf = "Gaussian")
gsva_result <- gsva(gsvapar)
gsva_result <- as.data.frame(gsva_result)
rownames(gsva_result) <- NULL


gsva_long <- gsva_result %>%
  pivot_longer(cols = everything(),
               names_to = "SampleID",
               values_to = "GSVA_Score"
  )


gsva_long <- left_join(gsva_long, metadata[, c("SampleID", "TMM", "TMM_Case")], by = "SampleID")


ggplot(gsva_long, aes(x = TMM, y = GSVA_Score, fill = TMM, color = TMM)) +
  geom_boxplot(size = 0.2, alpha = 0.5, outlier.shape = NA) +
  geom_point(position = position_jitter(width = 0.2), size = 3) +
  scale_fill_manual(values = c("Telomerase" = "lightpink2", 
                               "NO_TMM" = "lightgreen",
                               "ALT" = "blue")) +
  scale_color_manual(values = c("Telomerase"="darkred", 
                                "NO_TMM" = "darkgreen",
                                "ALT" = "blue")) +
  theme_classic() +
  labs(x = "TMM Group", y = "GSVA Score") +
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
  

# genes positively upregulated in NO_TMM.
candidate_genes2 <- list(TMM = c("ERICH5", "MFAP3L", "FAT4", "MMP16", "TENM4", "SDC1", "LRRN1", "SECISBP2L", "IPP", "TRAK2", "SLC38A2",
                                 "OPHN1", "MT1F"
                                 ))
gsvapar <- gsvaParam(as.matrix(Expression), candidate_genes2, kcdf = "Gaussian")
gsva_result <- gsva(gsvapar)
gsva_result <- as.data.frame(gsva_result)
rownames(gsva_result) <- NULL


gsva_long <- gsva_result %>%
  pivot_longer(cols = everything(),
               names_to = "SampleID",
               values_to = "GSVA_Score"
  )


gsva_long <- left_join(gsva_long, metadata[, c("SampleID", "TMM", "TMM_Case")], by = "SampleID")


ggplot(gsva_long, aes(x = TMM, y = GSVA_Score, fill = TMM, color = TMM)) +
  geom_boxplot(size = 0.2, alpha = 0.5, outlier.shape = NA) +
  geom_point(position = position_jitter(width = 0.2), size = 3) +
  scale_fill_manual(values = c("Telomerase" = "lightpink2", 
                               "NO_TMM" = "lightgreen",
                               "ALT" = "blue")) +
  scale_color_manual(values = c("Telomerase"="darkred", 
                                "NO_TMM" = "darkgreen",
                                "ALT" = "blue")) +
  theme_classic() +
  labs(x = "TMM Group", y = "GSVA Score") +
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



#######################


## doing same thing in Ackerman data.
gct_file <- parse_gctx("Neuroblastoma_208Samples.gct")
ackerman_NB <- gct_file@mat

# Loading metadata.
ackerman_metadata <- read.table(file = 'NBL_Ackerman_CompleteMeta.txt', header = TRUE, sep = '\t')
# ackerman_metadata <- ackerman_metadata %>%
#   filter(Risk == "YES")

ackerman_metadata <- ackerman_metadata %>%
  filter(
    !(TERTRearrangement == "+" & TMM_Category != "Telomerase"),
    
    !( (ATRXMutation != "-" & ATRXMutation != "<NA>" & !is.na(ATRXMutation)) 
       & TMM_Category != "ALT"),
    
    !(APB == "+" & TMM_Category != "ALT")
  )

# only including SampleID in microarray data present in metadata.
ackerman_NB <- ackerman_NB[, colnames(ackerman_NB) %in% ackerman_metadata$SampleID]
ackerman_metadata <- ackerman_metadata[ackerman_metadata$SampleID %in% colnames(ackerman_NB), ]

ackerman_metadata <- ackerman_metadata %>%
  arrange(TMM_Case, TMM_Category)

ackerman_NB <- ackerman_NB[, match(ackerman_metadata$SampleID, colnames(ackerman_NB))]

# genes positively upregulated in TMM.
gsvapar <- gsvaParam(as.matrix(ackerman_NB), candidate_genes, kcdf = "Gaussian")
gsva_result <- gsva(gsvapar)
gsva_result <- as.data.frame(gsva_result)
rownames(gsva_result) <- NULL


gsva_long <- gsva_result %>%
  pivot_longer(cols = everything(),
               names_to = "SampleID",
               values_to = "GSVA_Score"
  )


gsva_long <- left_join(gsva_long, ackerman_metadata[, c("SampleID", "TMM_Category", "TMM_Case")], by = "SampleID")


ggplot(gsva_long, aes(x = TMM_Category, y = GSVA_Score, fill = TMM_Category, color = TMM_Category)) +
  geom_boxplot(size = 0.2, alpha = 0.5, outlier.shape = NA) +
  geom_point(position = position_jitter(width = 0.2), size = 3) +
  scale_fill_manual(values = c("Telomerase" = "lightpink2", 
                               "NO_TMM" = "lightgreen",
                               "ALT" = "blue")) +
  scale_color_manual(values = c("Telomerase"="darkred", 
                                "NO_TMM" = "darkgreen",
                                "ALT" = "blue")) +
  theme_classic() +
  labs(x = "TMM Group", y = "GSVA Score") +
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

# genes positively upregulated in NO_TMM.
gsvapar <- gsvaParam(as.matrix(ackerman_NB), candidate_genes2, kcdf = "Gaussian")
gsva_result <- gsva(gsvapar)
gsva_result <- as.data.frame(gsva_result)
rownames(gsva_result) <- NULL


gsva_long <- gsva_result %>%
  pivot_longer(cols = everything(),
               names_to = "SampleID",
               values_to = "GSVA_Score"
  )


gsva_long <- left_join(gsva_long, ackerman_metadata[, c("SampleID", "TMM_Category", "TMM_Case")], by = "SampleID")


ggplot(gsva_long, aes(x = TMM_Category, y = GSVA_Score, fill = TMM_Category, color = TMM_Category)) +
  geom_boxplot(size = 0.2, alpha = 0.5, outlier.shape = NA) +
  geom_point(position = position_jitter(width = 0.2), size = 3) +
  scale_fill_manual(values = c("Telomerase" = "lightpink2", 
                               "NO_TMM" = "lightgreen",
                               "ALT" = "blue")) +
  scale_color_manual(values = c("Telomerase"="darkred", 
                                "NO_TMM" = "darkgreen",
                                "ALT" = "blue")) +
  theme_classic() +
  labs(x = "TMM Group", y = "GSVA Score") +
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

#######################################


## doing same thing in 0532.

source("0532Data.R")

# genes upregulated in TMM.
gsvapar <- gsvaParam(as.matrix(tmm_lcpm), candidate_genes, kcdf = "Gaussian")

gsva_result <- gsva(gsvapar)
gsva_result <- as.data.frame(gsva_result)
rownames(gsva_result) <- NULL


gsva_long <- gsva_result %>%
  pivot_longer(cols = everything(),
               names_to = "SampleID",
               values_to = "GSVA_Score"
  )


gsva_long <- left_join(gsva_long, metadata_0532[, c("RNAseq_SampleID", "TMM", "TMMCase")], 
                       by = c("SampleID" = "RNAseq_SampleID"))


ggplot(gsva_long, aes(x = TMM, y = GSVA_Score, fill = TMM, color = TMM)) +
  geom_boxplot(size = 0.2, alpha = 0.5, outlier.shape = NA) +
  geom_point(position = position_jitter(width = 0.2), size = 3) +
  scale_fill_manual(values = c("TERT+" = "lightpink2", 
                               "TMM-" = "lightgreen",
                               "ALT+" = "blue")) +
  scale_color_manual(values = c("TERT+"="darkred", 
                                "TMM-" = "darkgreen",
                                "ALT+" = "blue")) +
  theme_classic() +
  labs(x = "TMM Group", y = "GSVA Score") +
  theme(
    axis.text.x = element_text(vjust = 1, hjust = 1),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 12, face = "bold"),
    legend.position = "none"
  ) +
  stat_compare_means(comparisons = list(c("ALT+","TMM-")), method= "t.test",
                     method.args = list(alternative ="two.sided"), size = 6, tip.length = 0.01,
                     label.y = 1.2) + 
  stat_compare_means(comparisons = list(c("TERT+","TMM-")), method= "t.test",
                     method.args = list(alternative ="two.sided"), size = 6, tip.length = 0.01,
                     label.y = 1.4)

svapar <- gsvaParam(as.matrix(tmm_lcpm), candidate_genes, kcdf = "Gaussian")

gsva_result <- gsva(gsvapar)
gsva_result <- as.data.frame(gsva_result)
rownames(gsva_result) <- NULL


gsva_long <- gsva_result %>%
  pivot_longer(cols = everything(),
               names_to = "SampleID",
               values_to = "GSVA_Score"
  )


gsva_long <- left_join(gsva_long, metadata_0532[, c("RNAseq_SampleID", "TMM", "TMMCase")], 
                       by = c("SampleID" = "RNAseq_SampleID"))


ggplot(gsva_long, aes(x = TMM, y = GSVA_Score, fill = TMM, color = TMM)) +
  geom_boxplot(size = 0.2, alpha = 0.5, outlier.shape = NA) +
  geom_point(position = position_jitter(width = 0.2), size = 3) +
  scale_fill_manual(values = c("TERT+" = "lightpink2", 
                               "TMM-" = "lightgreen",
                               "ALT+" = "blue")) +
  scale_color_manual(values = c("TERT+"="darkred", 
                                "TMM-" = "darkgreen",
                                "ALT+" = "blue")) +
  theme_classic() +
  labs(x = "TMM Group", y = "GSVA Score") +
  theme(
    axis.text.x = element_text(vjust = 1, hjust = 1),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 12, face = "bold"),
    legend.position = "none"
  ) +
  stat_compare_means(comparisons = list(c("ALT+","TMM-")), method= "t.test",
                     method.args = list(alternative ="two.sided"), size = 6, tip.length = 0.01,
                     label.y = 1.2) + 
  stat_compare_means(comparisons = list(c("TERT+","TMM-")), method= "t.test",
                     method.args = list(alternative ="two.sided"), size = 6, tip.length = 0.01,
                     label.y = 1.4)



# genes upregulated in NO_TMM.
gsvapar <- gsvaParam(as.matrix(tmm_lcpm), candidate_genes2, kcdf = "Gaussian")

gsva_result <- gsva(gsvapar)
gsva_result <- as.data.frame(gsva_result)
rownames(gsva_result) <- NULL


gsva_long <- gsva_result %>%
  pivot_longer(cols = everything(),
               names_to = "SampleID",
               values_to = "GSVA_Score"
  )


gsva_long <- left_join(gsva_long, metadata_0532[, c("RNAseq_SampleID", "TMM", "TMMCase")], 
                       by = c("SampleID" = "RNAseq_SampleID"))


ggplot(gsva_long, aes(x = TMM, y = GSVA_Score, fill = TMM, color = TMM)) +
  geom_boxplot(size = 0.2, alpha = 0.5, outlier.shape = NA) +
  geom_point(position = position_jitter(width = 0.2), size = 3) +
  scale_fill_manual(values = c("TERT+" = "lightpink2", 
                               "TMM-" = "lightgreen",
                               "ALT+" = "blue")) +
  scale_color_manual(values = c("TERT+"="darkred", 
                                "TMM-" = "darkgreen",
                                "ALT+" = "blue")) +
  theme_classic() +
  labs(x = "TMM Group", y = "GSVA Score") +
  theme(
    axis.text.x = element_text(vjust = 1, hjust = 1),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 12, face = "bold"),
    legend.position = "none"
  ) +
  stat_compare_means(comparisons = list(c("ALT+","TMM-")), method= "t.test",
                     method.args = list(alternative ="two.sided"), size = 6, tip.length = 0.01,
                     label.y = 1.2) + 
  stat_compare_means(comparisons = list(c("TERT+","TMM-")), method= "t.test",
                     method.args = list(alternative ="two.sided"), size = 6, tip.length = 0.01,
                     label.y = 1.4)


#######################################################################




#########################################################################################

# boxplots of only the candidate genes.

# layout for subplots.
num_genes <- length(candidate_genes)

# Saving the plot as a PDF.
pdf("rna-seq-regression_results.pdf", width = 4, height = 4)

regression_test_candidates <- Expression[rownames(Expression) %in% candidate_genes, ]
regression_test_candidates <- regression_test_candidates[, match(metadata$SampleID, colnames(regression_test_candidates))]

for (gene in rownames(regression_test_candidates)) {
  plot_data <- data.frame(gene_Expression = as.numeric(regression_test_candidates[gene, ]), 
                          TMM = metadata$TMM, TMM_Case = metadata$TMM_Case)
  
  # Fit linear model.
  model <- lm(gene_Expression ~ TMM_Case, data = plot_data)
  summary_model <- summary(model)
  r2_label <- paste0("R² = ", round(summary_model$r.squared, 3))
  
  
  # boxplot.
  p <- ggplot(plot_data, aes(x = TMM, y = gene_Expression, fill = TMM, color = TMM)) +
    geom_boxplot(size = 0.2, alpha = 0.5, outlier.shape = NA) +
    geom_point(position = position_jitter(width = 0.2), size = 3) +
    scale_fill_manual(values = c("ALT"="lightblue", 
                                 "NO_TMM" = "lightgreen", 
                                 "Telomerase"="lightpink2")) +
    scale_color_manual(values = c("ALT"="blue", 
                                  "Telomerase"="darkred", 
                                  "NO_TMM" = "darkgreen")) +
    theme_classic() +
    labs(title = gene, x = "Class", y = "Expression") +
    theme(
      axis.text.x = element_text(vjust = 1, hjust = 1),
      axis.title = element_text(size = 6),
      axis.text = element_text(size = 6, face = "bold"),
      legend.position = "none"
    ) +
    annotate("text", label = r2_label, y = max(plot_data$gene_Expression, na.rm = TRUE) * 1.05, x = 2,
             hjust = 0, size = 4) +
    stat_compare_means(comparisons = list(c("Telomerase","NO_TMM")), method= "t.test",
                       method.args = list(alternative ="two.sided"), size = 4, tip.length = 0.01, label.y = max(plot_data$gene_Expression, na.rm = TRUE) * 0.75) +
    stat_compare_means(comparisons = list(c("ALT","NO_TMM")), method= "t.test",
                       method.args = list(alternative ="two.sided"), size = 4, tip.length = 0.01, label.y = max(plot_data$gene_Expression, na.rm = TRUE) * 0.90)
  
  print(p)
  
  rm(summary_model)
  rm(plot_data)
  rm(model)
  rm(r2_label)
  
}

dev.off()

######################################################################################



###############################################################################################################
######################################################################################################################



###############################################################################################









