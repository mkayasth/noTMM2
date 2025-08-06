###

# We have FPKM data and log TPM, raw counts data from TARGET NB xenahub.
# fpkm not very appropriate for DGE because of normalization of sequence depth and gene length within sample, but factor not comparable between samples.
# proceeding with log(raw counts + 1) data and using Limma for dge. Limma, however, wont recover count-based dispersion info and log transformed data may not be perfectly linear.

###

library(limma)
library(tidyverse)
library(biomaRt)
library(dplyr)

# Loading bulk-seq log counts RNA and metadata file.
Expression <- read_delim("TARGET-NBL.star_counts.tsv")
metadata <- read_delim("Metadata_TARGETFinal_08012024.txt")


