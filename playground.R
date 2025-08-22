
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
    axis.title = element_text(size = 6),
    axis.text = element_text(size = 6, face = "bold"),
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
    axis.title = element_text(size = 6),
    axis.text = element_text(size = 6, face = "bold"),
    legend.position = "none"
  ) +
  stat_compare_means(comparisons = list(c("TMM","NO_TMM")), method= "t.test",
                     method.args = list(alternative ="two.sided"), size = 4, tip.length = 0.01,
                     label.y = 1)



###################################################################################################

# heatmap based on adrenergic and mesenchymal elements from single cell.

library(ComplexHeatmap)
library(grid)

source("DataCleaning.R")

metadata <- metadata %>%
  arrange(TMM_Case, TMM)

Expression <- Expression[, match(metadata$SampleID, colnames(Expression))]

heatmap_genes <- c("PHOX2A", "KLF7", "PHOX2B", "TH", "DBH", "TBX2", "ISL1", "GATA3", "HAND2", "GATA2", "ZNF536",
                   "CD44", "FN1", "VIM", "IRF1", "IRF2", "RUNX1", "RUNX2", "MEOX1", "MEOX2", "SIX1", "SIX4", "SOX9", "SMAD3", "WWTR1", "PRRX1")

heatmap_genes <- c("KLF7", "GATA3", "HAND2", "PHOX2A", "ISL1", "HAND1", "PHOX2B", "TFAP2B", "GATA2", "SATB1", "SIX3", "EYA1", "SOX11", "DACH1", "ASCL1", "HEY1", "KLF13", "PBX3",
                   "VIM", "FN1", "MEOX2", "ID1", "EGR3", "AEBP1", "CBFB", "IRF3", "IRF2", "IRF1", "TBX18", "MAFF", "RUNX2", "ZFP36L1", "NR3C1", "BHLHE41", "GLIS3", "RUNX1", "FOSL1", "FOSL2", "ELK4", "IFI16", "SIX4", "FLI1", "MAML2", "SMAD3", "DCAF6", "WWTR1", "SOX9", "MEF2D", "ZNF217", "PRRX1", "CREG1", "NOTCH2", "SIX1", "MEOX1")
  
  
  
expression_matrix <- t(scale(t(Expression[rownames(Expression) %in% heatmap_genes, ,drop = FALSE])))

# building split vector aligned to rownames(expression_matrix).
grp_map <- setNames(c(rep("Adrenergic", 18), rep("Mesenchymal", 36)), heatmap_genes)
row_groups <- factor(grp_map[rownames(expression_matrix)], levels = c("Adrenergic","Mesenchymal"))

ht1 <- Heatmap(
  expression_matrix,
  row_split = row_groups,
  cluster_columns = TRUE,
  cluster_rows = FALSE,
  show_column_names = TRUE,
  height = unit(17, "cm"),
  width = unit(15, "cm"),
  row_names_gp = gpar(fontsize = 8, fontface = "bold"),
  column_names_gp = gpar(fontsize = 2),
  column_names_rot = 45
)

pdf("adr-mesHeatmap.pdf", width = 8, height = 12)
draw(ht1)
dev.off()



