if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("GEOquery")
GEOquery::getGEOSuppFiles("GSM1514028", fetch_files = F)
GEOquery::getGEOSuppFiles("GSE137710", fetch_files = F)
GEOquery::getGEOSuppFiles("GSE139829", fetch_files = F)
