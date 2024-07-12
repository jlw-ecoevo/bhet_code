# Helper script for running gRodon on a large batch of genomes

# JLW 2020

# Load Packages ----------------------------------------------------------------

library(Biostrings)
library(gRodon)
library(dplyr)
library(parallel)
library(data.table)

# Functions --------------------------------------------------------------------

merge.easy <- function(df1,df2,key){
  df1 <- data.table(df1,key=key)
  df2 <- data.table(df2,key=key)
  return(unique(merge(df1,df2,all.x=TRUE,by=.EACHI,allow.cartesian=TRUE)))
}

predictGenome <- function(i,metadata){
  path <- metadata$Path[i]
  genome <- metadata$Accession[i]
  #print(genome)
  #print(path)
  
  genes <- readDNAStringSet(path)
  #ribo <- try(readLines(paste0(path,".ribo")))
  #highly_expressed <- gsub(" .*","",names(genes)) %in% ribo
  highly_expressed <- grepl("ribosomal protein",names(genes),ignore.case=T)
  
  #print(genome)
  print(path)
  print(sum(highly_expressed))
  
  growth99 <- try(predictGrowth(genes,
                            highly_expressed,
                            temperature = metadata$OGT99[i],
                            mode="full",
                            training_set="madin"))
  growth95 <- try(predictGrowth(genes,
                            highly_expressed,
                            temperature = metadata$OGT95[i],
                            mode="full",
                            training_set="madin"))
  growth20C <- try(predictGrowth(genes,
                                highly_expressed,
                                temperature = 20,
                                mode="full",
                                training_set="madin"))
  
  if((!inherits(growth99,"try-error")) & (!inherits(growth95,"try-error")) & (!inherits(growth20C,"try-error"))){
    return(data.frame(Genome = genome,
                      ASV = metadata$ASV[i],
                      OGT99 = growth99$OGT,
                      OGT95 = growth95$OGT,
                      CUBHE = growth99$CUBHE,
                      CUB = growth99$CUB,
                      dCUB = (growth99$CUB-growth99$CUBHE)/growth99$CUB,
                      d99 = growth99$d,
                      d95 = growth95$d,
                      d20C = growth20C$d,
                      nHE = growth99$nHE,
                      nGenes = length(genes),
                      stringsAsFactors = F))
  } else {
    return(data.frame(Genome = genome,
                      ASV = metadata$ASV[i],
                      OGT99 = NA,
                      OGT95 = NA,
                      CUBHE = NA,
                      CUB = NA,
                      dCUB = NA,
                      d99 = NA,
                      d95 = NA,
                      d20C =NA,
                      nHE = NA,
                      nGenes = length(genes),
                      stringsAsFactors = F))
  }
}


# Load Data --------------------------------------------------------------------
  
setwd("/project/jakeweis_896/gtdb220/gtdb_genomes_reps_r220/")
load("asv_OGT99.rda") # OGT metadata, precalculated in temp99.R
paths <- readLines("ffn.files") # Should hold the filepaths of the genomes you want to run prediction on

path_df <- data.frame(Accession = basename(paths) %>% substr(1,15),
                      Path = paths,
                      stringsAsFactors = F)

print(head(paths))
print(head(path_df))
asv_df <- merge.easy(asv_df,path_df,key="Accession") # Assuming the metadata file matches your paths file

# Run --------------------------------------------------------------------------

print(head(asv_df))
pred_list <- mclapply(1:nrow(asv_df),
                    predictGenome,
                    metadata = asv_df,
                    mc.cores = 63)

pred_df <- do.call("rbind",pred_list)

save(pred_df,file="gtdb220_genome_predictions.rda")
