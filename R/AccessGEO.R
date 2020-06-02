if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("GEOquery")
GEOquery::getGEOSuppFiles("GSM1514028", fetch_files = F)