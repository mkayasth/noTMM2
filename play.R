############### making cls and gct files for running GSEA.

##### first for log counts TARGET data.

# run DataCleaning.R first.

##### making gct file for DGE.
metadata_ALT <- metadata_ALT %>%
  arrange(TMM_Case, TMM)
ExpressionALT <- Expression[, colnames(Expression) %in% metadata_ALT$SampleID]
ExpressionALT <- ExpressionALT[, match(metadata_ALT$SampleID, colnames(ExpressionALT))]

metadata_Telomerase <- metadata_Telomerase %>%
  arrange(TMM_Case, TMM)
ExpressionTelomerase <- Expression[, colnames(Expression) %in% metadata_Telomerase$SampleID]
ExpressionTelomerase <- ExpressionTelomerase[, match(metadata_Telomerase$SampleID, colnames(ExpressionTelomerase))]


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

write_gct(ExpressionTelomerase, "ExpressionTelomeraseTargetLog.gct")
write_gct(ExpressionALT, "ExpressionAltTargetLog.gct")


 ### cls files.

# cls file for alt first.
out_cls <- "classesALTTarget.cls"


common_samples <- intersect(colnames(ExpressionALT), metadata_ALT$SampleID)
ExpressionALT <- ExpressionALT[, common_samples, drop = FALSE]
metadata_ALT <- metadata_ALT[match(common_samples, metadata_ALT$SampleID), ]

# Extract phenotype classes in the exact order of ExpressionTelomerase columns
classes <- as.character(metadata_ALT$TMM)


levels <- unique(classes)
k <- length(levels)
class_index <- match(classes, levels) - 1  # 0-based: first class = 0, second = 1, etc.

# Write CLS file
con <- file(out_cls, open = "wt")
on.exit(close(con), add = TRUE)

writeLines(sprintf("%d %d 1", length(classes), k), con)
writeLines(paste("#", paste(levels, collapse = " ")), con)
writeLines(paste(class_index, collapse = " "), con)

close(con)

# cls file for telomerase.
out_cls <- "classesTelomeraseTarget.cls"

# Ensure ExpressionTelomerase and metadata_Telomerase are aligned
common_samples <- intersect(colnames(ExpressionTelomerase), metadata_Telomerase$SampleID)
ExpressionTelomerase <- ExpressionTelomerase[, common_samples, drop = FALSE]
metadata_Telomerase <- metadata_Telomerase[match(common_samples, metadata_Telomerase$SampleID), ]

# Extract phenotype classes in the exact order of ExpressionTelomerase columns
classes <- as.character(metadata_Telomerase$TMM)


levels <- unique(classes)
k <- length(levels)
class_index <- match(classes, levels) - 1  # 0-based: first class = 0, second = 1, etc.

# Write CLS file
con <- file(out_cls, open = "wt")
on.exit(close(con), add = TRUE)

writeLines(sprintf("%d %d 1", length(classes), k), con)
writeLines(paste("#", paste(levels, collapse = " ")), con)
writeLines(paste(class_index, collapse = " "), con)

close(con)

########################################################################

##### gct and cls file for 0532 dataset.
metadata_Telomerase <- metadata_0532 %>%
  filter(TMM %in% c("TMM-", "TERT+"))

metadata_Telomerase <- metadata_Telomerase %>%
  arrange(TMMCase, TMM)

metadata_ALT <- metadata_0532 %>%
  filter(TMM %in% c("ALT+", "TMM-"))

metadata_ALT <- metadata_ALT %>%
  arrange(TMMCase, TMM)

geneExpressionTelomerase <- tmm_lcpm[,  colnames(tmm_lcpm) %in% metadata_Telomerase$RNAseq_SampleID]
geneExpressionTelomerase <- geneExpressionTelomerase[, match(metadata_Telomerase$RNAseq_SampleID, 
                                                             colnames(geneExpressionTelomerase))]

geneExpressionALT <- tmm_lcpm[, colnames(tmm_lcpm) %in% metadata_ALT$RNAseq_SampleID]
geneExpressionALT <- geneExpressionALT[, match(metadata_ALT$RNAseq_SampleID, 
                                               colnames(geneExpressionALT))]

write_gct(geneExpressionTelomerase, "ExpressionTelomerase0532.gct")
write_gct(geneExpressionALT, "ExpressionAlt0532.gct")

out_cls   <- "classesALT0532.cls"
classes <- as.character(metadata_ALT[["TMM"]])
levels <- unique(classes)
k <- length(levels)
class_index <- match(classes, levels) - 1  # 0-based.
con <- file(out_cls, open = "wt")
on.exit(close(con), add = TRUE)
writeLines(sprintf("%d %d 1", ncol(geneExpressionALT), k), con)
writeLines(paste("#", paste(levels, collapse = " ")), con)
writeLines(paste(class_index, collapse = " "), con)
close(con)


out_cls   <- "classesTelomerase0532.cls"
classes <- as.character(metadata_Telomerase[["TMM"]])
levels <- unique(classes)
k <- length(levels)
class_index <- match(classes, levels) - 1  # 0-based.
con <- file(out_cls, open = "wt")
on.exit(close(con), add = TRUE)
writeLines(sprintf("%d %d 1", ncol(geneExpressionTelomerase), k), con)
writeLines(paste("#", paste(levels, collapse = " ")), con)
writeLines(paste(class_index, collapse = " "), con)
close(con)



#################################################################################
##### making gct file for target-raw counts.
metadata_ALT <- metadata_ALT %>%
  arrange(TMM_Case, TMM)
tmm_lcpm_target_ALT <- tmm_lcpm_target[, colnames(tmm_lcpm_target) %in% metadata_ALT$SampleID]
tmm_lcpm_target_ALT <- tmm_lcpm_target_ALT[, match(metadata_ALT$SampleID, colnames(tmm_lcpm_target_ALT))]

metadata_Telomerase <- metadata_Telomerase %>%
  arrange(TMM_Case, TMM)
tmm_lcpm_target_telomerase <- tmm_lcpm_target[, colnames(tmm_lcpm_target) %in% metadata_Telomerase$SampleID]
tmm_lcpm_target_telomerase <- tmm_lcpm_target_telomerase[, match(metadata_Telomerase$SampleID, colnames(tmm_lcpm_target_telomerase))]

write_gct(tmm_lcpm_target_ALT, "ExpressionALTTargetRaw.gct")
write_gct(tmm_lcpm_target_telomerase, "ExpressionTelomeraseTargetRaw.gct")


#############################################################################################
# TARGET raw counts: candidate_genes_ALT_target, candidate_genes_Telomerase_target; TARGET log counts: ALT_candidates, telomerase_candidates.
# 0532 raw counts: candidate_genes_ALT, candidate_genes_Telomerase.


## finding common genes: DGE genes in TARGET log counts vs. TARGET raw counts.
length(intersect(rownames(candidate_genes_ALT_target), rownames(ALT_candidates)))
length(intersect(rownames(candidate_genes_Telomerase_target), rownames(telomerase_candidates)))


## finding common genes: DGE genes in TARGET log counts vs. 0532.
length(intersect(rownames(candidate_genes_ALT), rownames(ALT_candidates)))
length(intersect(rownames(candidate_genes_Telomerase), rownames(telomerase_candidates)))


## finding common genes: DGE genes in TARGET raw counts vs. 0532.
length(intersect(rownames(candidate_genes_ALT), rownames(candidate_genes_ALT_target)))
length(intersect(rownames(candidate_genes_Telomerase), rownames(candidate_genes_Telomerase_target)))

########################

# tmm_lcpm_target, metadata.
# Expression, metadata.
# tmm_lcpm, metadata_0532.


#### t-test: seeing if TARGET raw counts DGE are significant in TARGET log counts.
tmm_lcpm_target <- tmm_lcpm_target[rownames(tmm_lcpm_target) %in% rownames(Expression), ]
Expression <- Expression[rownames(Expression) %in% rownames(tmm_lcpm_target), ]


### now, telomerase candidates.

t_test_results <- data.frame(
  Gene = character(),
  p_value_t_test = numeric(),
  mean_NO_TMM = numeric(),
  mean_Telomerase = numeric(),
  stringsAsFactors = FALSE
)

for (gene_id in rownames(candidate_genes_Telomerase_target)) {
  
  if (!(gene_id %in% rownames(Expression))) next
  
  # extract numeric vector
  gene_expr <- as.numeric(Expression[gene_id, ])
  names(gene_expr) <- colnames(Expression)
  
  # group by phenotype
  c1 <- gene_expr[metadata_Telomerase$TMM == "NO_TMM"]
  c2 <- gene_expr[metadata_Telomerase$TMM == "Telomerase"]
  
  # skip if groups too small
  if (length(c1) < 2 || length(c2) < 2) next
  
  t_res <- t.test(c1, c2)
  
  t_test_results <- rbind(
    t_test_results,
    data.frame(
      Gene = gene_id,
      p_value_t_test = t_res$p.value,
      mean_NO_TMM = mean(c1, na.rm = TRUE),
      mean_Telomerase = mean(c2, na.rm = TRUE)
    )
  )
}

# Adjust FDR and add log2FC
t_test_results$FDR <- p.adjust(t_test_results$p_value_t_test, method = "fdr")
t_test_results$log2FC <- t_test_results$mean_Telomerase - t_test_results$mean_NO_TMM

t_test_results <- t_test_results[order(t_test_results$FDR), ]
t_test_results_sig <- t_test_results[t_test_results$FDR < 0.05, ]



