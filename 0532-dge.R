#############


source("0532DataCleaning.R")

# t-test of the genes.
##### 2) t-test of the DEGs.

# first looking at difference between NO_TMM and Telomerase
noTMM_candidates <- tmm_lcpm[
  rownames(tmm_lcpm) %in% rownames(candidate_genes_Telomerase),
  colnames(tmm_lcpm) %in% metadata_0532$RNAseq_SampleID[metadata_0532$TMM %in% c("TERT+", "TMM-")],
  drop = FALSE
]

# Initializing an empty data frame for storing t-test results.
t_test_results_Telomerase <- data.frame(
  Gene = character(),
  p_value_t_test = numeric(),
  stringsAsFactors = FALSE
)



# all p-values for FDR adjustment later.
all_p_values <- numeric()

# Looping through each gene in dge_gene.
for (i in 1:nrow(noTMM_candidates)) {
  
  gene_id <- rownames(noTMM_candidates)[i]
  
  
  # Extracting the expression values for this gene, keeping it as a matrix. We are only selecting significant genes (all samples for this gene).
  Expression <- noTMM_candidates[rownames(noTMM_candidates) == gene_id, , drop = FALSE]
  
  # Creating C1 and C2 group. For the gene we are working with, all NO_TMM sample expression data placed in c1_group and all TMM sample expression data placed in c2_group.
  c1_group <- Expression[, colnames(Expression) %in% metadata_0532$RNAseq_SampleID[metadata_0532$TMM == "TMM-"], drop = FALSE]
  c2_group <- Expression[, colnames(Expression) %in% metadata_0532$RNAseq_SampleID[metadata_0532$TMM == "TERT+"], drop = FALSE]
  
  
  # for a gene, compare c1 and c2 samples.
  t_test <- t.test(c1_group, c2_group)
  p_value_t_test <- t_test$p.value
  all_p_values <- c(all_p_values, p_value_t_test)
  
  
  # Storing the results.
  t_test_results_Telomerase <- rbind(t_test_results_Telomerase, data.frame(
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
t_test_results_Telomerase$fdr_t_test <- p.adjust(all_p_values, method = "BH")


# Storing significant t-test results where FDR is less than 0.01.
t_test_results_Telomerase_sig <- t_test_results_Telomerase[round(t_test_results_Telomerase$fdr_t_test, 2) <= 0.01, ]

#####

## t-test between NO_TMM and ALT
noTMM_candidates <- tmm_lcpm[
  rownames(tmm_lcpm) %in% rownames(candidate_genes_ALT),
  colnames(tmm_lcpm) %in% metadata_0532$RNAseq_SampleID[metadata_0532$TMM %in% c("ALT+", "TMM-")],
  drop = FALSE
]


# Initializing an empty data frame for storing t-test results.
t_test_results_ALT <- data.frame(
  Gene = character(),
  p_value_t_test = numeric(),
  stringsAsFactors = FALSE
)



# all p-values for FDR adjustment later.
all_p_values <- numeric()

# Looping through each gene in dge_gene.
for (i in 1:nrow(noTMM_candidates)) {
  
  gene_id <- rownames(noTMM_candidates)[i]
  
  
  # Extracting the expression values for this gene, keeping it as a matrix. We are only selecting significant genes (all samples for this gene).
  Expression <- noTMM_candidates[rownames(noTMM_candidates) == gene_id, , drop = FALSE]
  
  # Creating C1 and C2 group. For the gene we are working with, all NO_TMM sample expression data placed in c1_group and all TMM sample expression data placed in c2_group.
  c1_group <- Expression[, colnames(Expression) %in% metadata_0532$RNAseq_SampleID[metadata_0532$TMM == "TMM-"], drop = FALSE]
  c2_group <- Expression[, colnames(Expression) %in% metadata_0532$RNAseq_SampleID[metadata_0532$TMM == "ALT+"], drop = FALSE]
  
  
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


# Storing significant t-test results where FDR is less than 0.01.
t_test_results_ALT_sig <- t_test_results_ALT[round(t_test_results_ALT$fdr_t_test, 2) <= 0.05, ]



######################################################################################################


##### 3) Linear regression test of candidate genes.

# first ALT.

regression_results_ALT <- data.frame(
  Gene = character(),
  estimate = numeric(),
  p_value = numeric(),
  R.squared = numeric(),
  stringsAsFactors = FALSE
)

dge_gene_ALT <-  tmm_lcpm[
  rownames(tmm_lcpm) %in% t_test_results_ALT_sig$Gene,
  colnames(tmm_lcpm) %in% metadata_0532$RNAseq_SampleID[metadata_0532$TMM %in% c("ALT+", "TMM-")],
  drop = FALSE]


meta_sub <- metadata_0532[metadata_0532$RNAseq_SampleID %in% colnames(dge_gene_ALT), ]
meta_sub <- meta_sub[match(colnames(dge_gene_ALT), meta_sub$RNAseq_SampleID), ]


# for each gene in the candidate list, updating regression_results.
for (gene in rownames(dge_gene_ALT)) {
  
  data <- data.frame(gene_Expression = as.numeric(dge_gene_ALT[gene, ]), 
                     TMMstatus = meta_sub$TMM)
  
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

# Filtering for r-squared >= 0.3 & fdr <= 0.01.
regression_results_ALT_sig <- regression_results_ALT %>%
  filter(round(Adj.P.Value, 2) <= 0.05 & round(R.Squared, 1) >= 0.2)

regression_results_ALT_sig <- regression_results_ALT_sig %>%
  arrange(desc(R.Squared))

write.table(regression_results_ALT_sig, file = "0532-ALT.tsv", 
            sep = '\t', quote = FALSE, row.names = TRUE)


### Now, turn for Telomerase.

regression_results_Telomerase <- data.frame(
  Gene = character(),
  estimate = numeric(),
  p_value = numeric(),
  R.squared = numeric(),
  stringsAsFactors = FALSE
)

dge_gene_Telomerase <-  tmm_lcpm[
  rownames(tmm_lcpm) %in% t_test_results_Telomerase_sig$Gene,
  colnames(tmm_lcpm) %in% metadata_0532$RNAseq_SampleID[metadata_0532$TMM %in% c("TERT+", "TMM-")],
  drop = FALSE]


meta_sub <- metadata_0532[metadata_0532$RNAseq_SampleID %in% colnames(dge_gene_Telomerase), ]
meta_sub <- meta_sub[match(colnames(dge_gene_Telomerase), meta_sub$RNAseq_SampleID), ]


# for each gene in the candidate list, updating regression_results.
for (gene in rownames(dge_gene_Telomerase)) {
  
  data <- data.frame(gene_Expression = as.numeric(dge_gene_Telomerase[gene, ]), 
                     TMMstatus = meta_sub$TMM)
  
  # Fit linear model.
  model <- lm(gene_Expression ~ TMMstatus, data = data)
  summary_model <- summary(model)
  
  
  # Store regression results.
  regression_results_Telomerase <- rbind(regression_results_Telomerase, data.frame(
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
regression_results_Telomerase$Adj.P.Value <- p.adjust(regression_results_Telomerase$p_value, method = "fdr")

# Filtering for r-squared >= 0.2 & fdr <= 0.01.
regression_results_Telomerase_sig <- regression_results_Telomerase %>%
  filter(round(Adj.P.Value, 2) <= 0.05 & round(R.Squared, 1) >= 0.2)

regression_results_Telomerase_sig <- regression_results_Telomerase_sig %>%
  arrange(desc(R.Squared))

write.table(regression_results_Telomerase_sig, file = "0532-Telomerase.tsv", 
            sep = '\t', quote = FALSE, row.names = TRUE)



############
candidates_tmm_upregulated_telomerase_0532 <- regression_results_Telomerase_sig %>%
  filter(estimate < 0)

candidates_notmm_upregulated_telomerase_0532 <- regression_results_Telomerase_sig %>%
  filter(estimate > 0)

candidates_tmm_upregulated_ALT_0532 <- regression_results_ALT_sig %>%
  filter(estimate < 0)

candidates_notmm_upregulated_ALT_0532 <- regression_results_ALT_sig %>%
  filter(estimate > 0)


###########

metadata_0532_ALT <- metadata_0532[metadata_0532$TMM %in% c("ALT+", "TMM-"), ]
metadata_0532_ALT <- metadata_0532_ALT %>%
  arrange(TMM)
geneExpressionALT <- tmm_lcpm[, colnames(tmm_lcpm) %in% metadata_0532_ALT$RNAseq_SampleID]
geneExpressionALT <- geneExpressionALT[, match(metadata_0532_ALT$RNAseq_SampleID, colnames(geneExpressionALT))]

metadata_0532_Telomerase <- metadata_0532[metadata_0532$TMM %in% c("TERT+", "TMM-"), ]
metadata_0532_Telomerase <- metadata_0532_Telomerase %>%
  arrange(TMM)
geneExpressionTelomerase <- tmm_lcpm[, colnames(tmm_lcpm) %in% metadata_0532_Telomerase$RNAseq_SampleID]
geneExpressionTelomerase <- geneExpressionTelomerase[, match(metadata_0532_Telomerase$RNAseq_SampleID, colnames(geneExpressionTelomerase))]


# signature for genes upregulated in TMM over no_tmm.

candidates_tmm_upregulated_0532 <- rbind(candidates_tmm_upregulated_telomerase_0532, candidates_tmm_upregulated_ALT_0532)
a <- intersect(candidates_tmm_upregulated_telomerase_0532$Gene, candidates_tmm_upregulated_ALT_0532$Gene)

tmm_target_upregulated <- run_exhaustive_forward_auc(expr_matrix = tmm_lcpm, metadata = metadata_0532, 
                                                     candidate_genes = unique(candidates_tmm_upregulated_0532$Gene), 
                                                     phenotype_col = "TMMCase",
                                                     label_one = "TMM", label_two = "NO_TMM",
                                                     max_genes = 20,
                                                     pivot_gene = "TERT")

# signature for genes upregulated in TMM over no_tmm.

candidates_notmm_upregulated_0532 <- rbind(candidates_notmm_upregulated_telomerase_0532, candidates_notmm_upregulated_ALT_0532)
b <- intersect(candidates_notmm_upregulated_telomerase_0532$Gene, candidates_notmm_upregulated_ALT_0532$Gene)

tmm_target_downregulated <- run_exhaustive_forward_auc(expr_matrix = tmm_lcpm, metadata = metadata_0532, 
                                                     candidate_genes = unique(candidates_notmm_upregulated_0532$Gene), 
                                                     phenotype_col = "TMMCase",
                                                     label_one = "NO_TMM", label_two = "TMM",
                                                     max_genes = 20,
                                                     pivot_gene = "KCTD21")


# tmm_target5 <- run_exhaustive_forward_auc(expr_matrix = geneExpressionALT, metadata = metadata_0532_ALT, 
#                                           candidate_genes = candidates_tmm_upregulated_ALT_0532$Gene, 
#                                           phenotype_col = "TMM",
#                                           label_one = "ALT+", label_two = "TMM-",
#                                           max_genes = 20,
#                                           pivot_gene = "FOXJ1")
# 
# # signature for genes upregulated in Telomerase over no_tmm.
# tmm_target6 <- run_exhaustive_forward_auc(expr_matrix = geneExpressionTelomerase, metadata = metadata_0532_Telomerase, 
#                                           candidate_genes = candidates_tmm_upregulated_telomerase_0532$Gene, 
#                                           phenotype_col = "TMM",
#                                           label_one = "TERT+", label_two = "TMM-",
#                                           max_genes = 20,
#                                           pivot_gene = "TERT")
# 
# # signature of genes upregulated in no_tmm over ALT.
# tmm_target7 <- c("AMOTL1", "SEC63")
# 
# # signature for genes upregulated in no_tmm over Telomerase.
# tmm_target8 <- run_exhaustive_forward_auc(expr_matrix = geneExpressionTelomerase, metadata = metadata_0532_Telomerase, 
#                                           candidate_genes = candidates_notmm_upregulated_telomerase_0532$Gene, 
#                                           phenotype_col = "TMM",
#                                           label_one = "TMM-", label_two = "TERT+",
#                                           max_genes = 20,
#                                           pivot_gene = "ARMC2")

####

# genes positively upregulated in TMM.
candidate_genes <- list(TMM = c("TERT", "MAGEC2", "MFSD11"))
gsvapar <- gsvaParam(as.matrix(tmm_lcpm), candidate_genes, kcdf = "Gaussian")
gsva_result <- gsva(gsvapar)
gsva_result <- as.data.frame(gsva_result)
rownames(gsva_result) <- NULL


gsva_long <- gsva_result %>%
  pivot_longer(cols = everything(),
               names_to = "SampleID",
               values_to = "GSVA_Score"
  )


gsva_long <- left_join(gsva_long, metadata_0532[, c("RNAseq_SampleID", "TMM", "TMMCase")], by = c("SampleID" = "RNAseq_SampleID"))

gsva_long$TMM <- factor(gsva_long$TMM, levels = c("ALT+", "TMM-", "TERT+"))


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

# genes positively upregulated in NO_TMM.
candidate_genes2 <- list(TMM = c("KCTD21", "FAM124A", "SEC63", "PAK3", "PPP6R2", "ARMC2"))
gsvapar <- gsvaParam(as.matrix(tmm_lcpm), candidate_genes2, kcdf = "Gaussian")
gsva_result <- gsva(gsvapar)
gsva_result <- as.data.frame(gsva_result)
rownames(gsva_result) <- NULL


gsva_long <- gsva_result %>%
  pivot_longer(cols = everything(),
               names_to = "SampleID",
               values_to = "GSVA_Score"
  )


gsva_long <- left_join(gsva_long, metadata_0532[, c("RNAseq_SampleID", "TMM", "TMMCase")], by = c("SampleID" = "RNAseq_SampleID"))

gsva_long$TMM <- factor(gsva_long$TMM, levels = c("ALT+", "TMM-", "TERT+"))


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


####### doing the same thing for TARGET data.


# genes positively upregulated in TMM.
gsvapar <- gsvaParam(as.matrix(tmm_lcpm_target), candidate_genes, kcdf = "Gaussian")
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
gsvapar <- gsvaParam(as.matrix(tmm_lcpm_target), candidate_genes2, kcdf = "Gaussian")
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


##### doing the same for the Ackerman data.
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




######################### subtracting gsva from no_tmm upregulated by tmm upregulated gene score.

gsvapar <- gsvaParam(as.matrix(tmm_lcpm_target), candidate_genes, kcdf = "Gaussian")
gsva_result <- gsva(gsvapar)
gsva_result <- as.data.frame(gsva_result)
rownames(gsva_result) <- NULL
gsva_long <- gsva_result %>%
  pivot_longer(cols = everything(),
               names_to = "SampleID",
               values_to = "GSVA_Score_down"
  )

gsvapar2 <- gsvaParam(as.matrix(tmm_lcpm_target), candidate_genes2, kcdf = "Gaussian")
gsva_result2 <- gsva(gsvapar2)
gsva_result2 <- as.data.frame(gsva_result2)
rownames(gsva_result2) <- NULL
gsva_long2 <- gsva_result2 %>%
  pivot_longer(cols = everything(),
               names_to = "SampleID",
               values_to = "GSVA_Score_up"
  )

gsva_long <- left_join(gsva_long, gsva_long2, by = "SampleID")
gsva_long <- left_join(gsva_long, metadata[, c("SampleID", "TMM", "TMM_Case")], by = "SampleID")
gsva_long$GSVA_Score <- gsva_long$GSVA_Score_up - gsva_long$GSVA_Score_down

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

### doing the same for 0532 data.
gsvapar <- gsvaParam(as.matrix(tmm_lcpm), candidate_genes, kcdf = "Gaussian")
gsva_result <- gsva(gsvapar)
gsva_result <- as.data.frame(gsva_result)
rownames(gsva_result) <- NULL
gsva_long <- gsva_result %>%
  pivot_longer(cols = everything(),
               names_to = "SampleID",
               values_to = "GSVA_Score_down"
  )

gsvapar2 <- gsvaParam(as.matrix(tmm_lcpm), candidate_genes2, kcdf = "Gaussian")
gsva_result2 <- gsva(gsvapar2)
gsva_result2 <- as.data.frame(gsva_result2)
rownames(gsva_result2) <- NULL
gsva_long2 <- gsva_result2 %>%
  pivot_longer(cols = everything(),
               names_to = "SampleID",
               values_to = "GSVA_Score_up"
  )

gsva_long <- left_join(gsva_long, gsva_long2, by = "SampleID")
gsva_long <- left_join(gsva_long, metadata_0532[, c("RNAseq_SampleID", "TMM", "TMMCase")], by = c("SampleID" = "RNAseq_SampleID"))
gsva_long$GSVA_Score <- gsva_long$GSVA_Score_up - gsva_long$GSVA_Score_down

gsva_long$TMM <- factor(gsva_long$TMM, levels = c("ALT+", "TMM-", "TERT+"))


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


#### doing the same thing for Ackerman.
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
               values_to = "gsva_Score_down"
  )

# genes positively upregulated in NO_TMM.
gsvapar2 <- gsvaParam(as.matrix(ackerman_NB), candidate_genes2, kcdf = "Gaussian")
gsva_result2 <- gsva(gsvapar2)
gsva_result2 <- as.data.frame(gsva_result2)
rownames(gsva_result2) <- NULL


gsva_long2 <- gsva_result2 %>%
  pivot_longer(cols = everything(),
               names_to = "SampleID",
               values_to = "gsva_Score_up"
  )

gsva_long <- left_join(gsva_long, gsva_long2, by = "SampleID")
gsva_long <- left_join(gsva_long, ackerman_metadata[, c("SampleID", "TMM_Category", "TMM_Case")], by = "SampleID")
gsva_long$GSVA_Score <- gsva_long$gsva_Score_up - gsva_long$gsva_Score_down

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

#######################################################################################################

###########################################################################################################

## trying on 



# 
# ### clustering using these genes.
# 
# # trying with gene rankings.
# ranked_tmm_lcpm <- apply(tmm_lcpm, 2, function(x) rank(x, ties.method = "average"))
# 
# cluster_genes <- ranked_tmm_lcpm[
#   rownames(ranked_tmm_lcpm) %in% c("AMOTL1", "ADBR2", "DLGAP4", "CCBE1", "KIFAP3", "TERT", "MAGEC2", "MFSD11", "LRRC56"), ,
#   drop = FALSE
# ]
# 
# 
# # Scale and transpose first.
# pca_scale <- scale(t(cluster_genes))
# 
# # Run PCA
# pca_result <- prcomp(pca_scale, center = TRUE, scale. = TRUE)
# 
# # plotting the visualize the variance.
# pca_var <- pca_result$sdev^2
# pca_var_explained <- pca_var / sum(pca_var)
# 
# 
# # Make full scree plot to identify plateau.
# plot(pca_var_explained * 100, type = "b", pch = 19,
#      xlab = "Principal Component",
#      ylab = "Variance Explained (%)",
#      main = "Scree Plot: Full PCA")
# 
# # From the elbow plot, graph seems to plateau at PC:1-10.
# pc_scores <- pca_result$x[, 1:3]
# 
# ##### Clustering the whole expression data.
# 
# # k-means clustering.
# set.seed(123)
# umap_res <- umap(pc_scores)
# 
# km_res <- kmeans(umap_res, centers = 2)
# 
# # Plotting with clusters.
# tmm_colors <- c("TERT+" = "red", "TMM-" = "green", "ALT+" = "blue")
# 
# plot(umap_res, col = tmm_colors[metadata_0532$TMM], pch = 19,
#      xlab = "UMAP 1", ylab = "UMAP 2", main = "k-means Clusters on UMAP")
# 
# ########
# 
# ## gsva t-test to see the difference.
# pos_0532 <- list(NO_TMM = c("AMOTL1", "ADRB2", "DLGAP4", "CCBE1", "KIFAP3"))
# neg_0532 <- list(MO_TMM = c("TERT", "MAGEC2", "MFSD11", "LRRC56"))
# 
# 
# # genes positively upregulated in NO_TMM.
# gsvapar <- gsvaParam(as.matrix(tmm_lcpm), pos_0532, kcdf = "Gaussian")
# gsva_result <- gsva(gsvapar)
# gsva_result <- as.data.frame(gsva_result)
# rownames(gsva_result) <- NULL
# 
# 
# gsva_long <- gsva_result %>%
#   pivot_longer(cols = everything(),
#                names_to = "SampleID",
#                values_to = "GSVA_Score"
#   )
# 
# 
# gsva_long <- left_join(gsva_long, metadata_0532[, c("RNAseq_SampleID", "TMM", "TMMCase")], by = c("SampleID" = "RNAseq_SampleID"))
# 
# 
# ggplot(gsva_long, aes(x = TMMCase, y = GSVA_Score, fill = TMMCase, color = TMMCase)) +
#   geom_boxplot(size = 0.2, alpha = 0.5, outlier.shape = NA) +
#   geom_point(position = position_jitter(width = 0.2), size = 3) +
#   scale_fill_manual(values = c("TMM" = "lightpink2", 
#                                "NO_TMM" = "lightgreen")) +
#   scale_color_manual(values = c("TMM"="darkred", 
#                                 "NO_TMM" = "darkgreen")) +
#   theme_classic() +
#   labs(x = "TM Group", y = "GSVA Score") +
#   theme(
#     axis.text.x = element_text(vjust = 1, hjust = 1),
#     axis.title = element_text(size = 18),
#     axis.text = element_text(size = 12, face = "bold"),
#     legend.position = "none"
#   ) +
#   stat_compare_means(comparisons = list(c("TMM","NO_TMM")), method= "t.test",
#                      method.args = list(alternative ="two.sided"), size = 6, tip.length = 0.01,
#                      label.y = 1.2)
# 
# # genes upregulated in TMM.
# gsvapar <- gsvaParam(as.matrix(tmm_lcpm), neg_0532, kcdf = "Gaussian")
# gsva_result <- gsva(gsvapar)
# gsva_result <- as.data.frame(gsva_result)
# rownames(gsva_result) <- NULL
# 
# 
# gsva_long <- gsva_result %>%
#   pivot_longer(cols = everything(),
#                names_to = "SampleID",
#                values_to = "GSVA_Score"
#   )
# 
# 
# gsva_long <- left_join(gsva_long, metadata_0532[, c("RNAseq_SampleID", "TMM", "TMMCase")], by = c("SampleID" = "RNAseq_SampleID"))
# 
# 
# ggplot(gsva_long, aes(x = TMMCase, y = GSVA_Score, fill = TMMCase, color = TMMCase)) +
#   geom_boxplot(size = 0.2, alpha = 0.5, outlier.shape = NA) +
#   geom_point(position = position_jitter(width = 0.2), size = 3) +
#   scale_fill_manual(values = c("TMM" = "lightpink2", 
#                                "NO_TMM" = "lightgreen")) +
#   scale_color_manual(values = c("TMM"="darkred", 
#                                 "NO_TMM" = "darkgreen")) +
#   theme_classic() +
#   labs(x = "TM Group", y = "GSVA Score") +
#   theme(
#     axis.text.x = element_text(vjust = 1, hjust = 1),
#     axis.title = element_text(size = 18),
#     axis.text = element_text(size = 12, face = "bold"),
#     legend.position = "none"
#   ) +
#   stat_compare_means(comparisons = list(c("TMM","NO_TMM")), method= "t.test",
#                      method.args = list(alternative ="two.sided"), size = 6, tip.length = 0.01,
#                      label.y = 1.2)
# 
# 
# #####################################################################################
# ## using signature from log counts data in 0532 data.
# 
# no_tmm_upregulated <- list(NO_TMM = c("CPNE8", "PGM2L1", "LIFR", "CNR1", "HECW2", "CPNE3", "HOXC9", "SNX16", "IGSF10"))
# tmm_upregulated <- list(NO_TMM = c("PRR7", "SAC3D1", "CCDC86", "DDN"))
# 
# # genes positively upregulated in NO_TMM.
# gsvapar <- gsvaParam(as.matrix(tmm_lcpm), tmm_upregulated, kcdf = "Gaussian")
# gsva_result <- gsva(gsvapar)
# gsva_result <- as.data.frame(gsva_result)
# rownames(gsva_result) <- NULL
# 
# 
# gsva_long <- gsva_result %>%
#   pivot_longer(cols = everything(),
#                names_to = "SampleID",
#                values_to = "GSVA_Score"
#   )
# 
# 
# gsva_long <- left_join(gsva_long, metadata_0532[, c("RNAseq_SampleID", "TMM", "TMMCase")], by = c("SampleID" = "RNAseq_SampleID"))
# 
# 
# ggplot(gsva_long, aes(x = TMMCase, y = GSVA_Score, fill = TMMCase, color = TMMCase)) +
#   geom_boxplot(size = 0.2, alpha = 0.5, outlier.shape = NA) +
#   geom_point(position = position_jitter(width = 0.2), size = 3) +
#   scale_fill_manual(values = c("TMM" = "lightpink2", 
#                                "NO_TMM" = "lightgreen")) +
#   scale_color_manual(values = c("TMM"="darkred", 
#                                 "NO_TMM" = "darkgreen")) +
#   theme_classic() +
#   labs(x = "TM Group", y = "GSVA Score") +
#   theme(
#     axis.text.x = element_text(vjust = 1, hjust = 1),
#     axis.title = element_text(size = 18),
#     axis.text = element_text(size = 12, face = "bold"),
#     legend.position = "none"
#   ) +
#   stat_compare_means(comparisons = list(c("TMM","NO_TMM")), method= "t.test",
#                      method.args = list(alternative ="two.sided"), size = 6, tip.length = 0.01,
#                      label.y = 1.2)
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# ######################################################################################################################
# # clustering based on these genes.
# 
# ranked_tmm_lcpm <- apply(tmm_lcpm, 2, function(x) rank(x, ties.method = "average"))
# 
# 
# cluster_genes <- ranked_tmm_lcpm[
#   rownames(tmm_lcpm) %in% c("CPNE8", "PGM2L1", "LIFR", "CNR1", "HECW2", "CPNE3", "HOXC9", "SNX16", "IGSF10",
#                             "PRR7", "SAC3D1", "CCDC86", "DDN"), ,
#   drop = FALSE
# ]
# 
# 
# # Scale and transpose first.
# pca_scale <- scale(t(cluster_genes))
# 
# # Run PCA
# pca_result <- prcomp(pca_scale, center = TRUE, scale. = TRUE)
# 
# # plotting the visualize the variance.
# pca_var <- pca_result$sdev^2
# pca_var_explained <- pca_var / sum(pca_var)
# 
# 
# # Make full scree plot to identify plateau.
# plot(pca_var_explained * 100, type = "b", pch = 19,
#      xlab = "Principal Component",
#      ylab = "Variance Explained (%)",
#      main = "Scree Plot: Full PCA")
# 
# # From the elbow plot, graph seems to plateau at PC:1-10.
# pc_scores <- pca_result$x[, 1:6]
# 
# ##### Clustering the whole expression data.
# 
# # k-means clustering.
# set.seed(123)
# umap_res <- umap(pc_scores)
# 
# km_res <- kmeans(umap_res, centers = 2)
# 
# # Plotting with clusters.
# tmm_colors <- c("TERT+" = "red", "TMM-" = "green", "ALT+" = "blue")
# 
# plot(umap_res, col = tmm_colors[metadata_0532$TMM], pch = 19,
#      xlab = "UMAP 1", ylab = "UMAP 2", main = "k-means Clusters on UMAP")
