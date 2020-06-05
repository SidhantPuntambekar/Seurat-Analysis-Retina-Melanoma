#Query for GEO ref matrix
#"expression profiling by high throughput sequencing"[DataSet Type] AND ("mtx"[Supplementary Files] OR "tsv"[Supplementary Files] OR "csv"[Supplementary Files] OR "txt"[Supplementary Files]) AND "single cell"[All Fields] AND "metadata"[All Fields]
records <- read.csv(file = "~/scRNA-seq-Practice/data_Reference_Matrix/records.csv")
head(records)
tail(records)
for (i in 1:36)
{
  print(GEOquery::getGEOSuppFiles(records[i,1,1], fetch_files = F))
  print(i)
}
