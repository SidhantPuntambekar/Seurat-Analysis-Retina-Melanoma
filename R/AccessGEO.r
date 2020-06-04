if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("GEOquery")
GEOquery::getGEOSuppFiles("GSM1514028", fetch_files = F)
GEOquery::getGEOSuppFiles("GSE137710", fetch_files = F)
GEOquery::getGEOSuppFiles("GSE139829", fetch_files = F)
browseVignettes("GEOquery")

library(GEOquery)
gds <- getGEO("GDS507")
gsm <- getGEO("GSM11805")
head(Meta(gsm))

Table(gsm)[1:5,]
Columns(gsm)
Columns(gds)[,1:3]

gse <- getGEO("GSE781", GSEMatrix = FALSE)
head(Meta(gse))
names(GSMList(gse))
GSMList(gse)[[1]]
names(GPLList(gse))
