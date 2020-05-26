#remotes::install_github("rnabioco/clustifyr")
<<<<<<< HEAD
library(dplyr)
library(Seurat)
library(patchwork)
=======
>>>>>>> 508c844ae9ac5118d3e0c2ad1246c73cdb9672dd
library(clustifyr)
library(tidyverse)
library(usethis)

mat_humanMelanomaDC <- read_tsv("ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE137nnn/GSE137710/suppl/GSE137710_human_melanoma_counts_normalized_9315x19445.tsv.gz") %>%
  select(-X1) %>% t()
colnames(mat_humanMelanomaDC) <- seq_len(ncol(mat_humanMelanomaDC))
meta_humanMelanomaDC <- read_tsv("ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE137nnn/GSE137710/suppl/GSE137710_human_melanoma_cell_metadata_9315x14.tsv.gz") %>%
  mutate("rn" = row_number())
ref_humanMelanomaDC <- average_clusters(
  mat_humanMelanomaDC,
  meta_humanMelanomaDC,
  "cell_type"
)

usethis::use_data(ref_humanMelanomaDC, compress = "xz", overwrite = TRUE)
<<<<<<< HEAD

melanoma <- CreateSeuratObject(counts = mat_humanMelanomaDC, project = "Melanoma", min.cells = 3, min.features = 200) 
=======
>>>>>>> 508c844ae9ac5118d3e0c2ad1246c73cdb9672dd
