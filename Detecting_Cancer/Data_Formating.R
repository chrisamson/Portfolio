library(GEOquery)
library(limma)
library(umap)
library(tidyr)
library(stringr)
library(RSQLite)

#Data is from Asakura K et al. studies of circulating microRNA in cancer
#cases and controls

#Get initial data
DATA <- getGEO("GSE137140", GSEMatrix = TRUE)

#Extract samples
samples <- DATA$GSE137140_series_matrix.txt.gz@phenoData@data$geo_accession

#Get data for each sample
TABLE <- data.frame()
for (x in 1:length(samples)){
  print(samples[x])
  sample_data <- getGEO(samples[x], GSEMatrix = TRUE)
  info <- data.frame(X = sample_data@header$characteristics_ch1) %>% separate(X, into = c('A','B'), sep = ': ')
  tmp <- as.data.frame(t(rbind(data.frame(row.names = c('sample_id'), values = c(samples[x])),
                               data.frame(row.names = str_replace_all(gsub(r"{\s*\([^\)]+\)}","",info$A), ' ','_'), values = info$B),
                               data.frame(row.names = sample_data@dataTable@table$ID_REF, values = sample_data@dataTable@table$VALUE))))
  rm(info)
  TABLE[x,] = tmp
  if (x == 1){
    TABLE <- tmp
  }else{
    TABLE[x,] = tmp
  }
  rm(tmp)
}
row.names(TABLE) <- NULL

#Connect to SQLite database
#Note: due to the size of this table in an ideal situation this data would be
#stored in a SQL database other than SQLite, however SQLite was used due to the
#ease of reproducibility of this project
con <- dbConnect(RSQLite::SQLite(), "../Databases/Lung_Cancer_miRNA_DB.sqlite")
#Deposit as three tables due to the number of columns involved

#Table Part A
dbWriteTable(con, "miRNA_Table_A" ,TABLE[,1:2000])

#Table Part B
dbWriteTable(con, "miRNA_Table_B" ,cbind(sample_id = TABLE$sample_id,TABLE[,2001:2570]))

#Disconect
dbDisconnect(con)
