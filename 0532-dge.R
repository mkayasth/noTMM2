#############

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
t_test_results_Telomerase_sig <- t_test_results_Telomerase[round(t_test_results_Telomerase$fdr_t_test, 2) <= 0.05, ]

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
  filter(round(Adj.P.Value, 2) <= 0.05 & round(R.Squared, 1) >= 0.3)

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

# Filtering for r-squared >= 0.3 & fdr <= 0.01.
regression_results_Telomerase_sig <- regression_results_Telomerase %>%
  filter(round(Adj.P.Value, 2) <= 0.05 & round(R.Squared, 1) >= 0.3)
write.table(regression_results_Telomerase_sig, file = "0532-Telomerase.tsv", 
            sep = '\t', quote = FALSE, row.names = TRUE)


#######################################################################################################

###########################################################################################################

### clustering using these genes.

# trying with gene rankings.
ranked_tmm_lcpm <- apply(tmm_lcpm, 2, function(x) rank(x, ties.method = "average"))

cluster_genes <- ranked_tmm_lcpm[
  rownames(ranked_tmm_lcpm) %in% c("AMOTL1", "ADBR2", "DLGAP4", "CCBE1", "KIFAP3", "TERT", "MAGEC2", "MFSD11", "LRRC56"), ,
  drop = FALSE
]


# Scale and transpose first.
pca_scale <- scale(t(cluster_genes))

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

# From the elbow plot, graph seems to plateau at PC:1-10.
pc_scores <- pca_result$x[, 1:3]

##### Clustering the whole expression data.

# k-means clustering.
set.seed(123)
umap_res <- umap(pc_scores)

km_res <- kmeans(umap_res, centers = 2)

# Plotting with clusters.
tmm_colors <- c("TERT+" = "red", "TMM-" = "green", "ALT+" = "blue")

plot(umap_res, col = tmm_colors[metadata_0532$TMM], pch = 19,
     xlab = "UMAP 1", ylab = "UMAP 2", main = "k-means Clusters on UMAP")

########

## gsva t-test to see the difference.
pos_0532 <- list(NO_TMM = c("AMOTL1", "ADRB2", "DLGAP4", "CCBE1", "KIFAP3"))
neg_0532 <- list(MO_TMM = c("TERT", "MAGEC2", "MFSD11", "LRRC56"))


# genes positively upregulated in NO_TMM.
gsvapar <- gsvaParam(as.matrix(tmm_lcpm), pos_0532, kcdf = "Gaussian")
gsva_result <- gsva(gsvapar)
gsva_result <- as.data.frame(gsva_result)
rownames(gsva_result) <- NULL


gsva_long <- gsva_result %>%
  pivot_longer(cols = everything(),
               names_to = "SampleID",
               values_to = "GSVA_Score"
  )


gsva_long <- left_join(gsva_long, metadata_0532[, c("RNAseq_SampleID", "TMM", "TMMCase")], by = c("SampleID" = "RNAseq_SampleID"))


ggplot(gsva_long, aes(x = TMMCase, y = GSVA_Score, fill = TMMCase, color = TMMCase)) +
  geom_boxplot(size = 0.2, alpha = 0.5, outlier.shape = NA) +
  geom_point(position = position_jitter(width = 0.2), size = 3) +
  scale_fill_manual(values = c("TMM" = "lightpink2", 
                               "NO_TMM" = "lightgreen")) +
  scale_color_manual(values = c("TMM"="darkred", 
                                "NO_TMM" = "darkgreen")) +
  theme_classic() +
  labs(x = "TM Group", y = "GSVA Score") +
  theme(
    axis.text.x = element_text(vjust = 1, hjust = 1),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 12, face = "bold"),
    legend.position = "none"
  ) +
  stat_compare_means(comparisons = list(c("TMM","NO_TMM")), method= "t.test",
                     method.args = list(alternative ="two.sided"), size = 6, tip.length = 0.01,
                     label.y = 1.2)

# genes upregulated in TMM.
gsvapar <- gsvaParam(as.matrix(tmm_lcpm), neg_0532, kcdf = "Gaussian")
gsva_result <- gsva(gsvapar)
gsva_result <- as.data.frame(gsva_result)
rownames(gsva_result) <- NULL


gsva_long <- gsva_result %>%
  pivot_longer(cols = everything(),
               names_to = "SampleID",
               values_to = "GSVA_Score"
  )


gsva_long <- left_join(gsva_long, metadata_0532[, c("RNAseq_SampleID", "TMM", "TMMCase")], by = c("SampleID" = "RNAseq_SampleID"))


ggplot(gsva_long, aes(x = TMMCase, y = GSVA_Score, fill = TMMCase, color = TMMCase)) +
  geom_boxplot(size = 0.2, alpha = 0.5, outlier.shape = NA) +
  geom_point(position = position_jitter(width = 0.2), size = 3) +
  scale_fill_manual(values = c("TMM" = "lightpink2", 
                               "NO_TMM" = "lightgreen")) +
  scale_color_manual(values = c("TMM"="darkred", 
                                "NO_TMM" = "darkgreen")) +
  theme_classic() +
  labs(x = "TM Group", y = "GSVA Score") +
  theme(
    axis.text.x = element_text(vjust = 1, hjust = 1),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 12, face = "bold"),
    legend.position = "none"
  ) +
  stat_compare_means(comparisons = list(c("TMM","NO_TMM")), method= "t.test",
                     method.args = list(alternative ="two.sided"), size = 6, tip.length = 0.01,
                     label.y = 1.2)


#####################################################################################
## using signature from log counts data in 0532 data.

no_tmm_upregulated <- list(NO_TMM = c("CPNE8", "PGM2L1", "LIFR", "CNR1", "HECW2", "CPNE3", "HOXC9", "SNX16", "IGSF10"))
tmm_upregulated <- list(NO_TMM = c("PRR7", "SAC3D1", "CCDC86", "DDN"))

# genes positively upregulated in NO_TMM.
gsvapar <- gsvaParam(as.matrix(tmm_lcpm), tmm_upregulated, kcdf = "Gaussian")
gsva_result <- gsva(gsvapar)
gsva_result <- as.data.frame(gsva_result)
rownames(gsva_result) <- NULL


gsva_long <- gsva_result %>%
  pivot_longer(cols = everything(),
               names_to = "SampleID",
               values_to = "GSVA_Score"
  )


gsva_long <- left_join(gsva_long, metadata_0532[, c("RNAseq_SampleID", "TMM", "TMMCase")], by = c("SampleID" = "RNAseq_SampleID"))


ggplot(gsva_long, aes(x = TMMCase, y = GSVA_Score, fill = TMMCase, color = TMMCase)) +
  geom_boxplot(size = 0.2, alpha = 0.5, outlier.shape = NA) +
  geom_point(position = position_jitter(width = 0.2), size = 3) +
  scale_fill_manual(values = c("TMM" = "lightpink2", 
                               "NO_TMM" = "lightgreen")) +
  scale_color_manual(values = c("TMM"="darkred", 
                                "NO_TMM" = "darkgreen")) +
  theme_classic() +
  labs(x = "TM Group", y = "GSVA Score") +
  theme(
    axis.text.x = element_text(vjust = 1, hjust = 1),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 12, face = "bold"),
    legend.position = "none"
  ) +
  stat_compare_means(comparisons = list(c("TMM","NO_TMM")), method= "t.test",
                     method.args = list(alternative ="two.sided"), size = 6, tip.length = 0.01,
                     label.y = 1.2)


# clustering based on these genes.

ranked_tmm_lcpm <- apply(tmm_lcpm, 2, function(x) rank(x, ties.method = "average"))


cluster_genes <- ranked_tmm_lcpm[
  rownames(tmm_lcpm) %in% c("CPNE8", "PGM2L1", "LIFR", "CNR1", "HECW2", "CPNE3", "HOXC9", "SNX16", "IGSF10",
                            "PRR7", "SAC3D1", "CCDC86", "DDN"), ,
  drop = FALSE
]


# Scale and transpose first.
pca_scale <- scale(t(cluster_genes))

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

# From the elbow plot, graph seems to plateau at PC:1-10.
pc_scores <- pca_result$x[, 1:6]

##### Clustering the whole expression data.

# k-means clustering.
set.seed(123)
umap_res <- umap(pc_scores)

km_res <- kmeans(umap_res, centers = 2)

# Plotting with clusters.
tmm_colors <- c("TERT+" = "red", "TMM-" = "green", "ALT+" = "blue")

plot(umap_res, col = tmm_colors[metadata_0532$TMM], pch = 19,
     xlab = "UMAP 1", ylab = "UMAP 2", main = "k-means Clusters on UMAP")