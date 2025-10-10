### run DGE-LimmaLogCounts first so that Expression and geneExpression has the same dimensions. :)

library(EnhancedVolcano)
library(tidyverse)
library(dplyr)
library(ggpubr)
library(biomaRt)
library(edgeR)

source("DataCleaning.R")

# un-log the data.
gene_Expression <- round(2^Expression - 1)


metadata <- metadata %>%
  filter(COG.Risk.Group == "High Risk")

# filtering for protein coding genes.
mart <- useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl")
protein_coding_genes <- getBM(attributes = c("ensembl_gene_id", "gene_biotype", "hgnc_symbol"),
                              filters="biotype",
                              values = "protein_coding",
                              mart = mart)

gene_Expression <- gene_Expression[rownames(gene_Expression) %in% protein_coding_genes$hgnc_symbol, ]
gene_Expression <- gene_Expression[, colnames(gene_Expression) %in% metadata$SampleID]


# setting metadata order.
metadata <- metadata %>%
  arrange(TMM_Case, TMM)


gene_Expression <- gene_Expression[, match(metadata$SampleID, colnames(gene_Expression)), drop = FALSE]


##### Now, running edgeR for DGE //

########################################################

# building model matrix.

# First, determining the factors of TMM.
group1 <- as.factor(metadata$TMM)

# model matrix ~ without an intercept term.
design <- model.matrix(~group1+0)

# creating differential gene expression object.
dge_TMM <- DGEList(counts=gene_Expression,group=group1)

# removing lowly expressed genes with cpm < 1 in 5% of the samples.
keep <- filterByExpr(dge_TMM, design = design)
dge_TMM <- dge_TMM[keep, , keep.lib.sizes = FALSE]

# TMM normalization.
dge_TMM <- calcNormFactors(dge_TMM, method = "TMM")


# Calculating dispersion and fitting the model.
d <- estimateDisp(dge_TMM, design, verbose=TRUE)
fit <- glmQLFit(d, design, robust = TRUE)

# contrast parameter (ALT-NO_TMM).
contrast <- makeContrasts(altVSnotmm = group1ALT - group1NO_TMM,
                          TelomeraseVSnotmm = group1Telomerase - group1NO_TMM,
                          ALTvsTelomerase = group1ALT - group1Telomerase,
                          levels = design
)


# differential expression test.
fitALT <- glmQLFTest(fit, contrast = contrast[, "altVSnotmm"])
fitTelomerase <- glmQLFTest(fit, contrast = contrast[, "TelomeraseVSnotmm"])

# results
top_ALT_target <- topTags(fitALT, n = Inf)
top_Telomerase_target <- topTags(fitTelomerase, n = Inf)



# filtering for candidate genes.
candidate_genes_ALT_target <- subset(top_ALT_target$table, FDR <= 0.05 & abs(logFC) >= 0.5)
candidate_genes_Telomerase_target <- subset(top_Telomerase_target$table, FDR <= 0.05 & abs(logFC) >= 0.5)



tmm_cpm_target  <- cpm(dge_TMM, normalized.lib.sizes = TRUE)            
tmm_lcpm_target <- cpm(dge_TMM, log = TRUE, prior.count = 1)   


# only including genes present in both TARGET and 0532 data.
tmm_lcpm_0532 <- read_tsv("tmm_lcpm_0532.tsv")
tmm_lcpm_0532 <- as.data.frame(tmm_lcpm_0532)
rownames(tmm_lcpm_0532) <- tmm_lcpm_0532$...1 
tmm_lcpm_0532$...1 <- NULL

tmm_lcpm_target <- tmm_lcpm_target[rownames(tmm_lcpm_target) %in% rownames(tmm_lcpm_0532), ,drop = FALSE]
write.table(tmm_lcpm_target, file = "tmm_lcpm_target.tsv", sep = "\t",
            row.names = TRUE, col.names = NA)

candidate_genes_ALT_target <- candidate_genes_ALT_target[rownames(candidate_genes_ALT_target) %in% rownames(tmm_lcpm_target), ]
candidate_genes_Telomerase_target <- candidate_genes_Telomerase_target[rownames(candidate_genes_Telomerase_target) %in% rownames(tmm_lcpm_target), ]

#############

# t-test of the genes.
##### 2) t-test of the DEGs.

# first looking at difference between NO_TMM and Telomerase
noTMM_candidates <- tmm_lcpm_target[
  rownames(tmm_lcpm_target) %in% rownames(candidate_genes_Telomerase_target),
  colnames(tmm_lcpm_target) %in% metadata$SampleID[metadata$TMM %in% c("Telomerase", "NO_TMM")],
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
  c1_group <- Expression[, colnames(Expression) %in% metadata$SampleID[metadata$TMM == "NO_TMM"], drop = FALSE]
  c2_group <- Expression[, colnames(Expression) %in% metadata$SampleID[metadata$TMM == "Telomerase"], drop = FALSE]
  
  
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


# Storing significant t-test results where FDR is less than 0.05.
t_test_results_Telomerase_sig <- t_test_results_Telomerase[round(t_test_results_Telomerase$fdr_t_test, 2) <= 0.01, ]

#####

## t-test between NO_TMM and ALT
noTMM_candidates <- tmm_lcpm_target[
  rownames(tmm_lcpm_target) %in% rownames(candidate_genes_ALT_target),
  colnames(tmm_lcpm_target) %in% metadata$SampleID[metadata$TMM %in% c("ALT", "NO_TMM")],
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
  c1_group <- Expression[, colnames(Expression) %in% metadata$SampleID[metadata$TMM == "NO_TMM"], drop = FALSE]
  c2_group <- Expression[, colnames(Expression) %in% metadata$SampleID[metadata$TMM == "ALT"], drop = FALSE]
  
  
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
t_test_results_ALT_sig <- t_test_results_ALT[round(t_test_results_ALT$fdr_t_test, 2) <= 0.01, ]



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

dge_gene_ALT <-  tmm_lcpm_target[
  rownames(tmm_lcpm_target) %in% t_test_results_ALT_sig$Gene,
  colnames(tmm_lcpm_target) %in% metadata$SampleID[metadata$TMM %in% c("ALT", "NO_TMM")],
  drop = FALSE]


meta_sub <- metadata[metadata$SampleID %in% colnames(dge_gene_ALT), ]
meta_sub <- meta_sub[match(colnames(dge_gene_ALT), meta_sub$SampleID), ]


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
  filter(round(Adj.P.Value, 2) <= 0.01 & round(R.Squared, 1) >= 0.2)
write.table(regression_results_ALT_sig, file = "TARGET-ALT.tsv", sep = '\t', quote = FALSE, row.names = TRUE)




### Now, turn for Telomerase.

regression_results_Telomerase <- data.frame(
  Gene = character(),
  estimate = numeric(),
  p_value = numeric(),
  R.squared = numeric(),
  stringsAsFactors = FALSE
)

dge_gene_Telomerase <-  tmm_lcpm_target[
  rownames(tmm_lcpm_target) %in% t_test_results_Telomerase_sig$Gene,
  colnames(tmm_lcpm_target) %in% metadata$SampleID[metadata$TMM %in% c("Telomerase", "NO_TMM")],
  drop = FALSE]


meta_sub <- metadata[metadata$SampleID %in% colnames(dge_gene_Telomerase), ]
meta_sub <- meta_sub[match(colnames(dge_gene_Telomerase), meta_sub$SampleID), ]


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
  filter(round(Adj.P.Value, 2) <= 0.01 & round(R.Squared, 1) >= 0.2)
write.table(regression_results_Telomerase_sig, file = "TARGET-Telomerase.tsv", sep = '\t', quote = FALSE, row.names = TRUE)


############## 

candidates_tmm_upregulated_telomerase_target <- regression_results_Telomerase_sig %>%
  filter(estimate > 0)

candidates_notmm_upregulated_telomerase_target <- regression_results_Telomerase_sig %>%
  filter(estimate < 0)

candidates_tmm_upregulated_ALT_target <- regression_results_ALT_sig %>%
  filter(estimate < 0)

candidates_notmm_upregulated_ALT_target <- regression_results_ALT_sig %>%
  filter(estimate > 0)



candidates_tmm_upregulated_ALT_target <- candidates_tmm_upregulated_ALT_target %>%
  arrange(desc(R.Squared))
candidates_tmm_upregulated_telomerase_target <- candidates_tmm_upregulated_telomerase_target %>%
  arrange(desc(R.Squared))
candidates_notmm_upregulated_ALT_target <- candidates_notmm_upregulated_ALT_target %>%
  arrange(desc(R.Squared))
candidates_notmm_upregulated_telomerase_target <- candidates_notmm_upregulated_telomerase_target %>%
  arrange(desc(R.Squared))


########################################################################################
# finding best signature.

gene_Expression_ALT <- tmm_lcpm_target[, colnames(tmm_lcpm_target) %in% metadata_ALT$SampleID]
gene_Expression_ALT <- tmm_lcpm_target[, match(metadata_ALT$SampleID, colnames(tmm_lcpm_target))]

gene_Expression_Telomerase <- tmm_lcpm_target[, colnames(tmm_lcpm_target) %in% metadata_Telomerase$SampleID]
gene_Expression_Telomerase <- tmm_lcpm_target[, match(metadata_Telomerase$SampleID, colnames(tmm_lcpm_target))]


# # signature for genes upregulated in ALT over no_tmm.
# tmm_target9 <- run_exhaustive_forward_auc(expr_matrix = gene_Expression_ALT, metadata = metadata_ALT, 
#                                           candidate_genes = candidates_tmm_upregulated_ALT_target$Gene, 
#                                           phenotype_col = "TMM",
#                                           label_one = "ALT", label_two = "NO_TMM",
#                                           max_genes = 20,
#                                           pivot_gene = "LCN15")

# signature for genes upregulated in Telomerase over no_tmm.

# tmm_target10 <- run_exhaustive_forward_auc(expr_matrix = gene_Expression_Telomerase, metadata = metadata_Telomerase, 
#                                           candidate_genes = candidates_tmm_upregulated_telomerase_target$Gene, 
#                                           phenotype_col = "TMM",
#                                           label_one = "Telomerase", label_two = "NO_TMM",
#                                           max_genes = 20,
#                                           pivot_gene = "TERT")
# 
# 
# # signature for genes upregulated in no_tmm over ALT.
# tmm_target11 <- run_exhaustive_forward_auc(expr_matrix = gene_Expression_ALT, metadata = metadata_ALT, 
#                                           candidate_genes = candidates_notmm_upregulated_ALT_target$Gene, 
#                                           phenotype_col = "TMM",
#                                           label_one = "NO_TMM", label_two = "ALT",
#                                           max_genes = 20,
#                                           pivot_gene = "CPNE8")
# 
# # signature for genes upregulated in no_tmm over Telomerase.
# tmm_target12 <- run_exhaustive_forward_auc(expr_matrix = gene_Expression_Telomerase, metadata = metadata_Telomerase, 
#                                           candidate_genes = candidates_notmm_upregulated_telomerase_target$Gene, 
#                                           phenotype_col = "TMM",
#                                           label_one = "NO_TMM", label_two = "Telomerase",
#                                           max_genes = 20,
#                                           pivot_gene = "FAT4")

candidates_tmm_upregulated_target <- rbind(candidates_tmm_upregulated_telomerase_target, candidates_tmm_upregulated_ALT_target)
a <- intersect(candidates_tmm_upregulated_telomerase_target$Gene, candidates_tmm_upregulated_ALT_target$Gene)

tmm_target_upregulated <- run_exhaustive_forward_auc(expr_matrix = tmm_lcpm_target, metadata = metadata, 
                                                                                                candidate_genes = unique(candidates_tmm_upregulated_target$Gene), 
                                                                                                phenotype_col = "TMM_Case",
                                                                                                label_one = "TMM", label_two = "NO_TMM",
                                                                                                max_genes = 20,
                                                                                                pivot_gene = "WDR24")


candidates_notmm_upregulated_target <- rbind(candidates_notmm_upregulated_ALT_target, candidates_notmm_upregulated_telomerase_target)
b <- intersect(candidates_notmm_upregulated_telomerase_target$Gene, candidates_notmm_upregulated_ALT_target$Gene)

tmm_target_downregulated <- run_exhaustive_forward_auc(expr_matrix = tmm_lcpm_target, metadata = metadata, 
                                                       candidate_genes = candidates_notmm_upregulated_target$Gene, 
                                                       phenotype_col = "TMM_Case",
                                                       label_one = "NO_TMM", label_two = "TMM",
                                                       max_genes = 20,
                                                       pivot_gene = "SLC38A2")

########################################################################


########## genes positively upregulated in TMM.

candidate_genes <- list(TMM = c("WDR24", "ARL17A", "STOML2", "MACROD1", "LTB4R",
                                "TSEN54", "TTLL12", "CLUH", "SLC38A5", "C5orf47", "PES1"))
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


# ### trying the same thing in 0532.
source("0532DataCleaning.R")

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

#### trying the same thing in Ackerman.

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


####### genes positively upregulated in NO_TMM.

candidate_genes2 <- list(TMM = c("SLC38A2", "SECISBP2L", "FMO5", "MFAP3L"))
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


# ### trying the same thing in 0532.
source("0532DataCleaning.R")

# genes upregulated in TMM.
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

#### trying the same thing in Ackerman.

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

candidate_genes <- list(TMM = c("WDR24", "ARL17A", "STOML2", "MACROD1", "LTB4R",
                                "TSEN54", "TTLL12", "CLUH", "SLC38A5", "C5orf47", "PES1"))
candidate_genes2 <- list(TMM = c("SECISBP2L", "SLC38A2", "FMO5", "MFAP3L"))

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
candidate_genes <- list(TMM = c("WDR24", "ARL17A", "STOML2", "MACROD1", "LTB4R",
                                "TSEN54", "TTLL12", "CLUH", "SLC38A5", "C5orf47", "PES1"))
candidate_genes2 <- list(TMM = c("SECISBP2L", "SLC38A2", "FMO5", "MFAP3L"))

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
candidate_genes <- list(TMM = c("WDR24", "ARL17A", "STOML2", "MACROD1", "LTB4R",
                                "TSEN54", "TTLL12", "CLUH", "SLC38A5", "C5orf47", "PES1"))
candidate_genes2 <- list(TMM = c("SECISBP2L", "SLC38A2", "FMO5", "MFAP3L"))

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






#####################################################################################
# boxplots of only the candidate genes.

# layout for subplots.
num_genes <- length(candidates_notmm_upregulated_ALT_0532$Gene)

# Saving the plot as a PDF.
pdf("rna-seq-regression_results.pdf", width = 4, height = 4)

regression_test_candidates <- tmm_lcpm[rownames(tmm_lcpm) %in% candidates_notmm_upregulated_ALT_0532$Gene, ]
regression_test_candidates <- regression_test_candidates[, match(metadata_0532$RNAseq_SampleID, colnames(regression_test_candidates))]

metadata_0532$TMM <- factor(metadata_0532$TMM,
                            levels = c("ALT+", "TMM-", "TERT+"))

for (gene in rownames(regression_test_candidates)) {
  plot_data <- data.frame(gene_Expression = as.numeric(regression_test_candidates[gene, ]), 
                          TMM = metadata_0532$TMM, TMM_Case = metadata_0532$TMMCase)
  
  
  # Fit linear model.
  model <- lm(gene_Expression ~ TMM_Case, data = plot_data)
  summary_model <- summary(model)
  r2_label <- paste0("R² = ", round(summary_model$r.squared, 3))
  
  
  # boxplot.
  p <- ggplot(plot_data, aes(x = TMM, y = gene_Expression, fill = TMM, color = TMM)) +
    geom_boxplot(size = 0.2, alpha = 0.5, outlier.shape = NA) +
    geom_point(position = position_jitter(width = 0.2), size = 3) +
    scale_fill_manual(values = c("ALT+"="lightblue", 
                                 "TMM-" = "lightgreen", 
                                 "TERT+"="lightpink2")) +
    scale_color_manual(values = c("ALT+"="blue", 
                                  "TERT+"="darkred", 
                                  "TMM-" = "darkgreen")) +
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
    stat_compare_means(comparisons = list(c("TERT+","TMM-")), method= "t.test",
                       method.args = list(alternative ="two.sided"), size = 4, tip.length = 0.01, label.y = max(plot_data$gene_Expression, na.rm = TRUE) * 0.75) +
    stat_compare_means(comparisons = list(c("ALT+","TMM-")), method= "t.test",
                       method.args = list(alternative ="two.sided"), size = 4, tip.length = 0.01, label.y = max(plot_data$gene_Expression, na.rm = TRUE) * 0.90)
  
  print(p)
  
  rm(summary_model)
  rm(plot_data)
  rm(model)
  rm(r2_label)
  
}

dev.off()



# genes positively upregulated in TMM.
# 
# candidate_genes <- list(TMM = c("WDR24", "PGP", "USH1G", "LRRC74B", "GTDC1", "KCHN4", "ENGASE", "WWP1", "TERT", "ENTR1",
#                                 "CLN6", "MFSD3", "LTB4R", "HOOK2", "EXOSC5", "DNAJB5"))
# gsvapar <- gsvaParam(as.matrix(tmm_lcpm_target), candidate_genes, kcdf = "Gaussian")
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
# gsva_long <- left_join(gsva_long, metadata[, c("SampleID", "TMM", "TMM_Case")], by = "SampleID")
# 
# 
# ggplot(gsva_long, aes(x = TMM, y = GSVA_Score, fill = TMM, color = TMM)) +
#   geom_boxplot(size = 0.2, alpha = 0.5, outlier.shape = NA) +
#   geom_point(position = position_jitter(width = 0.2), size = 3) +
#   scale_fill_manual(values = c("Telomerase" = "lightpink2", 
#                                "NO_TMM" = "lightgreen",
#                                "ALT" = "blue")) +
#   scale_color_manual(values = c("Telomerase"="darkred", 
#                                 "NO_TMM" = "darkgreen",
#                                 "ALT" = "blue")) +
#   theme_classic() +
#   labs(x = "TMM Group", y = "GSVA Score") +
#   theme(
#     axis.text.x = element_text(vjust = 1, hjust = 1),
#     axis.title = element_text(size = 18),
#     axis.text = element_text(size = 12, face = "bold"),
#     legend.position = "none"
#   ) +
#   stat_compare_means(comparisons = list(c("ALT","NO_TMM")), method= "t.test",
#                      method.args = list(alternative ="two.sided"), size = 6, tip.length = 0.01,
#                      label.y = 1.2) +
#   stat_compare_means(comparisons = list(c("Telomerase","NO_TMM")), method= "t.test",
#                      method.args = list(alternative ="two.sided"), size = 6, tip.length = 0.01,
#                      label.y = 1.4) 
# 
# 
# # genes positively upregulated in NO_TMM.
# candidate_genes2 <- list(TMM = c("CPNE8", "DPY19L4", "SDC1", "FMO5", "WWP1", "FLRT2", 
#                                  "FAT4", "LYSMD3", "MINDY2", "FGFR1OP2", "RAB6D", "SLC38A2", "OPHN1", "MYO9A"))
# gsvapar <- gsvaParam(as.matrix(tmm_lcpm_target), candidate_genes2, kcdf = "Gaussian")
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
# gsva_long <- left_join(gsva_long, metadata[, c("SampleID", "TMM", "TMM_Case")], by = "SampleID")
# 
# 
# ggplot(gsva_long, aes(x = TMM, y = GSVA_Score, fill = TMM, color = TMM)) +
#   geom_boxplot(size = 0.2, alpha = 0.5, outlier.shape = NA) +
#   geom_point(position = position_jitter(width = 0.2), size = 3) +
#   scale_fill_manual(values = c("Telomerase" = "lightpink2", 
#                                "NO_TMM" = "lightgreen",
#                                "ALT" = "blue")) +
#   scale_color_manual(values = c("Telomerase"="darkred", 
#                                 "NO_TMM" = "darkgreen",
#                                 "ALT" = "blue")) +
#   theme_classic() +
#   labs(x = "TMM Group", y = "GSVA Score") +
#   theme(
#     axis.text.x = element_text(vjust = 1, hjust = 1),
#     axis.title = element_text(size = 18),
#     axis.text = element_text(size = 12, face = "bold"),
#     legend.position = "none"
#   ) +
#   stat_compare_means(comparisons = list(c("ALT","NO_TMM")), method= "t.test",
#                      method.args = list(alternative ="two.sided"), size = 6, tip.length = 0.01,
#                      label.y = 1.2) + 
#   stat_compare_means(comparisons = list(c("Telomerase","NO_TMM")), method= "t.test",
#                      method.args = list(alternative ="two.sided"), size = 6, tip.length = 0.01,
#                      label.y = 1.4)
# 
# #### trying the same thing in log TARGET.
# 
# source("DataCleaning.R")
# 
# # only high risk samples and only genes present in both TARGET and 0532 set.
# metadata <- metadata[metadata$COG.Risk.Group == "High Risk", ]
# metadata <- metadata %>%
#   arrange(TMM_Case, TMM)
# 
# Expression <- Expression[rownames(Expression) %in% rownames(tmm_lcpm), colnames(Expression) %in% metadata$SampleID]
# Expression <- Expression[, match(metadata$SampleID, colnames(Expression))]
# 
# 
# 
# 
# # genes positively upregulated in TMM.
# gsvapar <- gsvaParam(as.matrix(Expression), candidate_genes, kcdf = "Gaussian")
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
# gsva_long <- left_join(gsva_long, metadata[, c("SampleID", "TMM", "TMM_Case")], by = "SampleID")
# 
# 
# ggplot(gsva_long, aes(x = TMM, y = GSVA_Score, fill = TMM, color = TMM)) +
#   geom_boxplot(size = 0.2, alpha = 0.5, outlier.shape = NA) +
#   geom_point(position = position_jitter(width = 0.2), size = 3) +
#   scale_fill_manual(values = c("Telomerase" = "lightpink2", 
#                                "NO_TMM" = "lightgreen",
#                                "ALT" = "blue")) +
#   scale_color_manual(values = c("Telomerase"="darkred", 
#                                 "NO_TMM" = "darkgreen",
#                                 "ALT" = "blue")) +
#   theme_classic() +
#   labs(x = "TMM Group", y = "GSVA Score") +
#   theme(
#     axis.text.x = element_text(vjust = 1, hjust = 1),
#     axis.title = element_text(size = 18),
#     axis.text = element_text(size = 12, face = "bold"),
#     legend.position = "none"
#   ) +
#   stat_compare_means(comparisons = list(c("ALT","NO_TMM")), method= "t.test",
#                      method.args = list(alternative ="two.sided"), size = 6, tip.length = 0.01,
#                      label.y = 1.2) +
#   stat_compare_means(comparisons = list(c("Telomerase","NO_TMM")), method= "t.test",
#                      method.args = list(alternative ="two.sided"), size = 6, tip.length = 0.01,
#                      label.y = 1.4) 
# 
# 
# # genes positively upregulated in NO_TMM.
# gsvapar <- gsvaParam(as.matrix(Expression), candidate_genes2, kcdf = "Gaussian")
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
# gsva_long <- left_join(gsva_long, metadata[, c("SampleID", "TMM", "TMM_Case")], by = "SampleID")
# 
# 
# ggplot(gsva_long, aes(x = TMM, y = GSVA_Score, fill = TMM, color = TMM)) +
#   geom_boxplot(size = 0.2, alpha = 0.5, outlier.shape = NA) +
#   geom_point(position = position_jitter(width = 0.2), size = 3) +
#   scale_fill_manual(values = c("Telomerase" = "lightpink2", 
#                                "NO_TMM" = "lightgreen",
#                                "ALT" = "blue")) +
#   scale_color_manual(values = c("Telomerase"="darkred", 
#                                 "NO_TMM" = "darkgreen",
#                                 "ALT" = "blue")) +
#   theme_classic() +
#   labs(x = "TMM Group", y = "GSVA Score") +
#   theme(
#     axis.text.x = element_text(vjust = 1, hjust = 1),
#     axis.title = element_text(size = 18),
#     axis.text = element_text(size = 12, face = "bold"),
#     legend.position = "none"
#   ) +
#   stat_compare_means(comparisons = list(c("ALT","NO_TMM")), method= "t.test",
#                      method.args = list(alternative ="two.sided"), size = 6, tip.length = 0.01,
#                      label.y = 1.2) +
#   stat_compare_means(comparisons = list(c("Telomerase","NO_TMM")), method= "t.test",
#                      method.args = list(alternative ="two.sided"), size = 6, tip.length = 0.01,
#                      label.y = 1.4) 
# 
# ### trying the same thing in 0532.
# source("0532DataCleaning.R")
# 
# # genes upregulated in TMM.
# gsvapar <- gsvaParam(as.matrix(tmm_lcpm), candidate_genes, kcdf = "Gaussian")
# 
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
# gsva_long <- left_join(gsva_long, metadata_0532[, c("RNAseq_SampleID", "TMM", "TMMCase")], 
#                        by = c("SampleID" = "RNAseq_SampleID"))
# 
# gsva_long$TMM <- factor(gsva_long$TMM, levels = c("ALT+", "TMM-", "TERT+"))
# 
# 
# ggplot(gsva_long, aes(x = TMM, y = GSVA_Score, fill = TMM, color = TMM)) +
#   geom_boxplot(size = 0.2, alpha = 0.5, outlier.shape = NA) +
#   geom_point(position = position_jitter(width = 0.2), size = 3) +
#   scale_fill_manual(values = c("TERT+" = "lightpink2", 
#                                "TMM-" = "lightgreen",
#                                "ALT+" = "blue")) +
#   scale_color_manual(values = c("TERT+"="darkred", 
#                                 "TMM-" = "darkgreen",
#                                 "ALT+" = "blue")) +
#   theme_classic() +
#   labs(x = "TMM Group", y = "GSVA Score") +
#   theme(
#     axis.text.x = element_text(vjust = 1, hjust = 1),
#     axis.title = element_text(size = 18),
#     axis.text = element_text(size = 12, face = "bold"),
#     legend.position = "none"
#   ) +
#   stat_compare_means(comparisons = list(c("ALT+","TMM-")), method= "t.test",
#                      method.args = list(alternative ="two.sided"), size = 6, tip.length = 0.01,
#                      label.y = 1.2) + 
#   stat_compare_means(comparisons = list(c("TERT+","TMM-")), method= "t.test",
#                      method.args = list(alternative ="two.sided"), size = 6, tip.length = 0.01,
#                      label.y = 1.4)
# 
# # candidate genes upregulated in TMM-.
# gsvapar <- gsvaParam(as.matrix(tmm_lcpm), candidate_genes2, kcdf = "Gaussian")
# 
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
# gsva_long <- left_join(gsva_long, metadata_0532[, c("RNAseq_SampleID", "TMM", "TMMCase")], 
#                        by = c("SampleID" = "RNAseq_SampleID"))
# 
# gsva_long$TMM <- factor(gsva_long$TMM, levels = c("ALT+", "TMM-", "TERT+"))
# 
# 
# ggplot(gsva_long, aes(x = TMM, y = GSVA_Score, fill = TMM, color = TMM)) +
#   geom_boxplot(size = 0.2, alpha = 0.5, outlier.shape = NA) +
#   geom_point(position = position_jitter(width = 0.2), size = 3) +
#   scale_fill_manual(values = c("TERT+" = "lightpink2", 
#                                "TMM-" = "lightgreen",
#                                "ALT+" = "blue")) +
#   scale_color_manual(values = c("TERT+"="darkred", 
#                                 "TMM-" = "darkgreen",
#                                 "ALT+" = "blue")) +
#   theme_classic() +
#   labs(x = "TMM Group", y = "GSVA Score") +
#   theme(
#     axis.text.x = element_text(vjust = 1, hjust = 1),
#     axis.title = element_text(size = 18),
#     axis.text = element_text(size = 12, face = "bold"),
#     legend.position = "none"
#   ) +
#   stat_compare_means(comparisons = list(c("ALT+","TMM-")), method= "t.test",
#                      method.args = list(alternative ="two.sided"), size = 6, tip.length = 0.01,
#                      label.y = 1.2) + 
#   stat_compare_means(comparisons = list(c("TERT+","TMM-")), method= "t.test",
#                      method.args = list(alternative ="two.sided"), size = 6, tip.length = 0.01,
#                      label.y = 1.4)
# 
# 
# 
# 
# #### trying the same thing in Ackerman.
# 
# gct_file <- parse_gctx("Neuroblastoma_208Samples.gct")
# ackerman_NB <- gct_file@mat
# 
# # Loading metadata.
# ackerman_metadata <- read.table(file = 'NBL_Ackerman_CompleteMeta.txt', header = TRUE, sep = '\t')
# # ackerman_metadata <- ackerman_metadata %>%
# #   filter(Risk == "YES")
# 
# ackerman_metadata <- ackerman_metadata %>%
#   filter(
#     !(TERTRearrangement == "+" & TMM_Category != "Telomerase"),
#     
#     !( (ATRXMutation != "-" & ATRXMutation != "<NA>" & !is.na(ATRXMutation)) 
#        & TMM_Category != "ALT"),
#     
#     !(APB == "+" & TMM_Category != "ALT")
#   )
# 
# # only including SampleID in microarray data present in metadata.
# ackerman_NB <- ackerman_NB[, colnames(ackerman_NB) %in% ackerman_metadata$SampleID]
# ackerman_metadata <- ackerman_metadata[ackerman_metadata$SampleID %in% colnames(ackerman_NB), ]
# 
# ackerman_metadata <- ackerman_metadata %>%
#   arrange(TMM_Case, TMM_Category)
# 
# ackerman_NB <- ackerman_NB[, match(ackerman_metadata$SampleID, colnames(ackerman_NB))]
# 
# # genes positively upregulated in TMM.
# gsvapar <- gsvaParam(as.matrix(ackerman_NB), candidate_genes, kcdf = "Gaussian")
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
# gsva_long <- left_join(gsva_long, ackerman_metadata[, c("SampleID", "TMM_Category", "TMM_Case")], by = "SampleID")
# 
# 
# ggplot(gsva_long, aes(x = TMM_Category, y = GSVA_Score, fill = TMM_Category, color = TMM_Category)) +
#   geom_boxplot(size = 0.2, alpha = 0.5, outlier.shape = NA) +
#   geom_point(position = position_jitter(width = 0.2), size = 3) +
#   scale_fill_manual(values = c("Telomerase" = "lightpink2", 
#                                "NO_TMM" = "lightgreen",
#                                "ALT" = "blue")) +
#   scale_color_manual(values = c("Telomerase"="darkred", 
#                                 "NO_TMM" = "darkgreen",
#                                 "ALT" = "blue")) +
#   theme_classic() +
#   labs(x = "TMM Group", y = "GSVA Score") +
#   theme(
#     axis.text.x = element_text(vjust = 1, hjust = 1),
#     axis.title = element_text(size = 18),
#     axis.text = element_text(size = 12, face = "bold"),
#     legend.position = "none"
#   ) +
#   stat_compare_means(comparisons = list(c("ALT","NO_TMM")), method= "t.test",
#                      method.args = list(alternative ="two.sided"), size = 6, tip.length = 0.01,
#                      label.y = 1.2) +
#   stat_compare_means(comparisons = list(c("Telomerase","NO_TMM")), method= "t.test",
#                      method.args = list(alternative ="two.sided"), size = 6, tip.length = 0.01,
#                      label.y = 1.4)
# 
# 
# 
# 
# # genes positively upregulated in NO_TMM.
# gsvapar <- gsvaParam(as.matrix(ackerman_NB), candidate_genes2, kcdf = "Gaussian")
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
# gsva_long <- left_join(gsva_long, ackerman_metadata[, c("SampleID", "TMM_Category", "TMM_Case")], by = "SampleID")
# 
# 
# ggplot(gsva_long, aes(x = TMM_Category, y = GSVA_Score, fill = TMM_Category, color = TMM_Category)) +
#   geom_boxplot(size = 0.2, alpha = 0.5, outlier.shape = NA) +
#   geom_point(position = position_jitter(width = 0.2), size = 3) +
#   scale_fill_manual(values = c("Telomerase" = "lightpink2", 
#                                "NO_TMM" = "lightgreen",
#                                "ALT" = "blue")) +
#   scale_color_manual(values = c("Telomerase"="darkred", 
#                                 "NO_TMM" = "darkgreen",
#                                 "ALT" = "blue")) +
#   theme_classic() +
#   labs(x = "TMM Group", y = "GSVA Score") +
#   theme(
#     axis.text.x = element_text(vjust = 1, hjust = 1),
#     axis.title = element_text(size = 18),
#     axis.text = element_text(size = 12, face = "bold"),
#     legend.position = "none"
#   ) +
#   stat_compare_means(comparisons = list(c("ALT","NO_TMM")), method= "t.test",
#                      method.args = list(alternative ="two.sided"), size = 6, tip.length = 0.01,
#                      label.y = 1.2) +
#   stat_compare_means(comparisons = list(c("Telomerase","NO_TMM")), method= "t.test",
#                      method.args = list(alternative ="two.sided"), size = 6, tip.length = 0.01,
#                      label.y = 1.4)
# 
# 
# 
# #############################
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
# 
# 
# 
# 
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
# km_res <- kmeans(umap_res, centers = 3)
# 
# # Plotting with clusters.
# tmm_colors <- c("Telomerase" = "red", "NO_TMM" = "green", "ALT" = "blue")
# 
# plot(umap_res, col = tmm_colors[metadata$TMM], pch = 19,
#      xlab = "UMAP 1", ylab = "UMAP 2", main = "k-means Clusters on UMAP")

