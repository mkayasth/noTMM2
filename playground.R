
### looking at TERT distributions across the phenotype (Expression Levels.)

library(tidyverse)
library(dplyr)
library(ggpubr)
library(survival)
library(survminer)


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

write_gct(Expression, "Expression.gct")

### .cls file from metadata.

out_cls   <- "classes.cls"
 
classes <- as.character(metadata[["TMM_Case"]])

levels <- unique(classes)
k <- length(levels)
class_index <- match(classes, levels) - 1  # 0-based

con <- file(out_cls, open = "wt")
on.exit(close(con), add = TRUE)
writeLines(sprintf("%d %d 1", ncol(Expression), k), con)
writeLines(paste("#", paste(levels, collapse = " ")), con)
writeLines(paste(class_index, collapse = " "), con)
close(con)



