### run DataCleaning.R first :)


library(EnhancedVolcano)
library(tidyverse)
library(dplyr)


##### 1) Doing differential gene expression using Limma for our log raw counts data. NO_TMM vs. TMM.

# creating design matrix for limma.
metadata$TMM_Case <- as.factor(metadata$TMM_Case)
metadata <- metadata %>%
  arrange(TMM_Case)

design <- model.matrix(~ 0 + TMM_Case, data = metadata)
colnames(design) <- levels(metadata$TMM_Case)


# Estimating array weights to account for sample-specific variability in library sizes.
Expression <- Expression[, match(metadata$SampleID, colnames(Expression))]

weights <- arrayWeights(Expression, design = design)

# Fitting linear model using limma.
fit <- lmFit(Expression, design, weights = weights)
fit <- eBayes(fit, trend = TRUE)

# Contrasts: NO_TMM vs TMM && TMM Vs. NO_TMM.
contrast.matrix <- makeContrasts(
  NO_TMMvsTMM = NO_TMM-TMM,
  TMMvsNO_TMM = TMM-NO_TMM,
  levels = design
)

fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2, trend = TRUE)

# top results for NO_TMM.
noTMM_results <- topTable(fit2, coef = "NO_TMMvsTMM", number = Inf, adjust = "fdr")

# top results for TMM.
TMM_results <- topTable(fit2, coef = "TMMvsNO_TMM", number = Inf, adjust = "fdr")

# Filtering for FDR < 0.01.
noTMM_candidates <- noTMM_results[round(noTMM_results$adj.P.Val, 2) <= 0.01 & (noTMM_results$logFC >= 0.5 | noTMM_results$logFC <= -0.5), ]
TMM_candidates <- TMM_results[round(TMM_results$adj.P.Val, 2) <= 0.01 & (TMM_results$logFC >= 0.5 | TMM_results$logFC <= -0.5), ]


# in noTMM-TMM, positive genes will have gene status NO_TMM and negative logFC genes will have gene status TMM.
noTMM_results <- noTMM_results %>%
  mutate("Gene Status" = case_when(
   rownames(noTMM_results) %in% rownames(noTMM_candidates) & noTMM_results$logFC > 0 ~ "NO_TMM",
   rownames(noTMM_results) %in% rownames(noTMM_candidates) & noTMM_results$logFC < 0 ~ "TMM",
   TRUE ~ "Not significant"))

noTMM_candidates <- noTMM_candidates %>%
  mutate("Gene Status" = case_when(
    noTMM_candidates$logFC > 0 ~ "NO_TMM",
     noTMM_candidates$logFC < 0 ~ "TMM",
    TRUE ~ "Not significant"))



##### 2) t-test of the DEGs.
# Initializing an empty data frame for storing t-test results.
noTMM_candidates <- rownames_to_column(noTMM_candidates, var = "Gene")
dge_gene <- Expression[rownames(Expression) %in% noTMM_candidates$Gene, ,drop = FALSE]

t_test_results <- data.frame(
  Gene = character(),
  p_value_t_test = numeric(),
  stringsAsFactors = FALSE
)

# Looping through each gene in dge_gene.
for (i in 1:nrow(noTMM_candidates)) {
  
  gene_id <- noTMM_candidates$Gene[i]
  
  
  # Extracting the expression values for this gene, keeping it as a matrix. We are only selecting significant genes (all samples for this gene).
  gene_Expression <- dge_gene[rownames(dge_gene) == gene_id, , drop = FALSE]
  
  # Creating C1 and C2 group. For the gene we are working with, all NO_TMM sample expression data placed in c1_group and all TMM sample expression data placed in c2_group.
  c1_group <- gene_Expression[, metadata$TMM_Case == "NO_TMM", drop = FALSE]
  c2_group <- gene_Expression[, metadata$TMM_Case == "TMM", drop = FALSE]
  
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
  rm(p_value_t_test)
}


t_test_results <- merge(t_test_results, noTMM_candidates[, c("Gene", "Gene Status")],  by = "Gene", all.x = TRUE)


# Storing significant t-test results where p-value is less than 0.01.
t_test_results_sig <- t_test_results[round(t_test_results$p_value_t_test, 2) <= 0.01, ]

##### Making volcano plot for NO_TMM vs. TMM.

top_labels <- t_test_results_sig$Gene

EnhancedVolcano(noTMM_results,
                lab = rownames(noTMM_results),
                x = 'logFC',
                y = 'adj.P.Val',
                #selectLab = top_labels,  # highlighting signature genes.
                xlab = bquote(~Log[2]~ 'fold change'),
                ylab = bquote(~-Log[10]~ 'FDR'),
                title = NULL,
                subtitle = NULL,
                pCutoff = 0.01,
                FCcutoff = 0.5,
                pointSize = 2.0,
                arrowheads = FALSE,
                max.overlaps = 13,
                labFace = 'bold',
                boxedLabels = TRUE,
                labSize = 5.0,
                drawConnectors = TRUE,
                widthConnectors = 0.3,
                colAlpha = 0.8,
                legendLabels = c('Not Significant','Significant logFC','Significant P-value ','Significant P-value & LogFC'),
                col = c('grey80', 'grey50', 'grey25', 'purple'),
                ylim = c(0, 5),
                caption = "Cutoffs: FDR <= 0.01, |Log2FC| >= 0.5"
) + theme_classic() + 
  theme(axis.title = element_text(size = 18),
        axis.text = element_text(size = 15),
        legend.text = element_text(size = 12),
        legend.title = element_blank(),
        plot.caption = element_text(size = 14))


##########################################################################################

##### 1) differential expression: NO_TMM vs. ALT.

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

# top results for NO_TMM.
noTMM_results <- topTable(fit2, coef = "NO_TMMvsALT", number = Inf, adjust = "fdr")

# top results for ALT.
ALT_results <- topTable(fit2, coef = "ALTvsNO_TMM", number = Inf, adjust = "fdr")


# Filtering for FDR < 0.01.
noTMM_candidates <- noTMM_results[round(noTMM_results$adj.P.Val, 2) <= 0.01 & (noTMM_results$logFC >= 0.5 | noTMM_results$logFC <= -0.5), ]
ALT_candidates <- ALT_results[round(ALT_results$adj.P.Val, 2) <= 0.01 & (ALT_results$logFC >= 0.5 | ALT_results$logFC <= -0.5), ]


# in noTMM-TMM, positive genes will have gene status NO_TMM and negative logFC genes will have gene status TMM.
noTMM_results <- noTMM_results %>%
  mutate("Gene Status" = case_when(
    rownames(noTMM_results) %in% rownames(noTMM_candidates) & noTMM_results$logFC > 0 ~ "NO_TMM",
    rownames(noTMM_results) %in% rownames(noTMM_candidates) & noTMM_results$logFC < 0 ~ "ALT",
    TRUE ~ "Not significant"))

noTMM_candidates <- noTMM_candidates %>%
  mutate("Gene Status" = case_when(
    noTMM_candidates$logFC > 0 ~ "NO_TMM",
    noTMM_candidates$logFC < 0 ~ "ALT",
    TRUE ~ "Not significant"))




##### 2) t-test of the DEGs.
# Initializing an empty data frame for storing t-test results.
noTMM_candidates <- rownames_to_column(noTMM_candidates, var = "Gene")
dge_gene <- ExpressionALT[rownames(ExpressionALT) %in% noTMM_candidates$Gene, ,drop = FALSE]

t_test_results2 <- data.frame(
  Gene = character(),
  p_value_t_test = numeric(),
  stringsAsFactors = FALSE
)

# Looping through each gene in dge_gene.
for (i in 1:nrow(noTMM_candidates)) {
  
  gene_id <- noTMM_candidates$Gene[i]
  
  
  # Extracting the expression values for this gene, keeping it as a matrix. We are only selecting significant genes (all samples for this gene).
  gene_Expression <- dge_gene[rownames(dge_gene) == gene_id, , drop = FALSE]
  
  # Creating C1 and C2 group. For the gene we are working with, all NO_TMM sample expression data placed in c1_group and all ALT sample expression data placed in c2_group.
  c1_group <- gene_Expression[, metadata_ALT$TMM == "NO_TMM", drop = FALSE]
  c2_group <- gene_Expression[, metadata_ALT$TMM == "ALT", drop = FALSE]
  
  # for a gene, compare c1 and c2 samples.
  t_test <- t.test(c1_group, c2_group)
  p_value_t_test <- t_test$p.value
  
  
  # Storing the results.
  t_test_results2 <- rbind(t_test_results2, data.frame(
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
  rm(p_value_t_test)
}


t_test_results2 <- merge(t_test_results2, noTMM_candidates[, c("Gene", "Gene Status")],  by = "Gene", all.x = TRUE)


# Storing significant t-test results where p-value is less than 0.01.
t_test_results_sig2 <- t_test_results2[round(t_test_results2$p_value_t_test, 2) <= 0.01, ]

##### Making volcano plot for NO_TMM vs. TMM.

top_labels <- t_test_results_sig2$Gene

EnhancedVolcano(noTMM_results,
                lab = rownames(noTMM_results),
                x = 'logFC',
                y = 'adj.P.Val',
                # selectLab = top_labels,  # highlighting signature genes.
                xlab = bquote(~Log[2]~ 'fold change'),
                ylab = bquote(~-Log[10]~ 'FDR'),
                title = NULL,
                subtitle = NULL,
                pCutoff = 0.01,
                FCcutoff = 0.5,
                pointSize = 2.0,
                arrowheads = FALSE,
                max.overlaps = 17,
                labFace = 'bold',
                boxedLabels = TRUE,
                labSize = 5.0,
                drawConnectors = TRUE,
                widthConnectors = 0.3,
                colAlpha = 0.8,
                legendLabels = c('Not Significant','Significant logFC','Significant P-value ','Significant P-value & LogFC'),
                col = c('grey80', 'grey50', 'grey25', 'purple'),
                ylim = c(0, 5),
                caption = "Cutoffs: FDR <= 0.01, |Log2FC| >= 0.5"
) + theme_classic() + 
  theme(axis.title = element_text(size = 18),
        axis.text = element_text(size = 15),
        legend.text = element_text(size = 12),
        legend.title = element_blank(),
        plot.caption = element_text(size = 14))



##########################################################################################

##### 1) differential expression: NO_TMM vs. Telomerase.

metadata_Telomerase <- metadata[metadata$TMM == "Telomerase" | metadata$TMM == "NO_TMM", ]

# creating design matrix for limma.
metadata_Telomerase$TMM <- as.factor(metadata_Telomerase$TMM)

metadata_Telomerase <- metadata_Telomerase %>%
  arrange(TMM)


design <- model.matrix(~ 0 + TMM, data = metadata_Telomerase)
colnames(design) <- levels(metadata_Telomerase$TMM)


# Estimating array weights to account for sample-specific variability in library sizes.
ExpressionTelomerase <- Expression[, colnames(Expression) %in% metadata_Telomerase$SampleID]
ExpressionTelomerase <- ExpressionTelomerase[, match(metadata_Telomerase$SampleID, colnames(ExpressionTelomerase))]

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

# top results for NO_TMM.
noTMM_results <- topTable(fit2, coef = "NO_TMMvsTelomerase", number = Inf, adjust = "fdr")

# top results for Telomerase.
telomerase_results <- topTable(fit2, coef = "TelomerasevsNO_TMM", number = Inf, adjust = "fdr")

### no significant fdr for either case.

# Filtering for p-value < 0.01.
noTMM_candidates <- noTMM_results[round(noTMM_results$adj.P.Val, 2) <= 0.01 & (noTMM_results$logFC >= 0.5 | noTMM_results$logFC <= -0.5), ]
telomerase_candidates <- telomerase_results[round(telomerase_results$adj.P.Val, 2) <= 0.01 & (telomerase_results$logFC >= 0.5 | telomerase_results$logFC <= -0.5), ]


# in noTMM-TMM, positive genes will have gene status NO_TMM and negative logFC genes will have gene status TMM.
noTMM_results <- noTMM_results %>%
  mutate("Gene Status" = case_when(
    rownames(noTMM_results) %in% rownames(noTMM_candidates) & noTMM_results$logFC > 0 ~ "NO_TMM",
    rownames(noTMM_results) %in% rownames(noTMM_candidates) & noTMM_results$logFC < 0 ~ "Telomerase",
    TRUE ~ "Not significant"))

noTMM_candidates <- noTMM_candidates %>%
  mutate("Gene Status" = case_when(
    noTMM_candidates$logFC > 0 ~ "NO_TMM",
    noTMM_candidates$logFC < 0 ~ "Telomerase",
    TRUE ~ "Not significant"))




##### 2) t-test of the DEGs.
# Initializing an empty data frame for storing t-test results.
noTMM_candidates <- rownames_to_column(noTMM_candidates, var = "Gene")
dge_gene <- ExpressionTelomerase[rownames(ExpressionTelomerase) %in% noTMM_candidates$Gene, ,drop = FALSE]

t_test_results3 <- data.frame(
  Gene = character(),
  p_value_t_test = numeric(),
  stringsAsFactors = FALSE
)

# Looping through each gene in dge_gene.
for (i in 1:nrow(noTMM_candidates)) {
  
  gene_id <- noTMM_candidates$Gene[i]
  
  
  # Extracting the expression values for this gene, keeping it as a matrix. We are only selecting significant genes (all samples for this gene).
  gene_Expression <- dge_gene[rownames(dge_gene) == gene_id, , drop = FALSE]
  
  # Creating C1 and C2 group. For the gene we are working with, all NO_TMM sample expression data placed in c1_group and all Telomerase sample expression data placed in c2_group.
  c1_group <- gene_Expression[, metadata_Telomerase$TMM == "NO_TMM", drop = FALSE]
  c2_group <- gene_Expression[, metadata_Telomerase$TMM == "Telomerase", drop = FALSE]
  
  # for a gene, compare c1 and c2 samples.
  t_test <- t.test(c1_group, c2_group)
  p_value_t_test <- t_test$p.value
  
  
  # Storing the results.
  t_test_results3 <- rbind(t_test_results3, data.frame(
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
  rm(p_value_t_test)
}


t_test_results3 <- merge(t_test_results3, noTMM_candidates[, c("Gene", "Gene Status")],  by = "Gene", all.x = TRUE)


# Storing significant t-test results where p-value is less than 0.01.
t_test_results_sig3 <- t_test_results3[round(t_test_results3$p_value_t_test, 2) <= 0.01, ]

##### Making volcano plot for NO_TMM vs. TMM.

top_labels <- t_test_results_sig3$Gene

EnhancedVolcano(noTMM_results,
                lab = rownames(noTMM_results),
                x = 'logFC',
                y = 'adj.P.Value',
                #selectLab = top_labels,  # highlighting signature genes.
                xlab = bquote(~Log[2]~ 'fold change'),
                ylab = bquote(~-Log[10]~ 'FDR'),
                title = NULL,
                subtitle = NULL,
                pCutoff = 0.01,
                FCcutoff = 0.5,
                pointSize = 2.0,
                arrowheads = FALSE,
                max.overlaps = 17,
                labFace = 'bold',
                boxedLabels = TRUE,
                labSize = 5.0,
                drawConnectors = TRUE,
                widthConnectors = 0.3,
                colAlpha = 0.8,
                legendLabels = c('Not Significant','Significant logFC','Significant P-value ','Significant P-value & LogFC'),
                col = c('grey80', 'grey50', 'grey25', 'purple'),
                ylim = c(0, 5),
                caption = "Cutoffs: FDR <= 0.01, |Log2FC| >= 0.5"
) + theme_classic() + 
  theme(axis.title = element_text(size = 18),
        axis.text = element_text(size = 15),
        legend.text = element_text(size = 12),
        legend.title = element_blank(),
        plot.caption = element_text(size = 14))








