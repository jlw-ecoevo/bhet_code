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
  
  growth <- try(predictGrowth(genes,
                            highly_expressed,
                            temperature = 20,
                            mode="full",
                            training_set="madin"))
  
  if((!inherits(growth,"try-error"))){
    return(data.frame(Genome = genome,
                      OGT = growth$OGT,
                      CUBHE = growth$CUBHE,
                      CUB = growth$CUB,
                      dCUB = (growth$CUB-growth$CUBHE)/growth$CUB,
                      CPB = growth$CPB,
                      ConsistencyHE = growth$ConsistencyHE,
                      LowerCI = growth$LowerCI,
                      UpperCI = growth$UpperCI,
                      d = growth$d,
                      nHE = growth$nHE,
                      nGenes = length(genes),
                      stringsAsFactors = F))
  } else {
    return(data.frame(Genome = genome,
                      OGT = NA,
                      CUBHE = NA,
                      CUB = NA,
                      dCUB = NA,
                      CPB = NA,
                      ConsistencyHE = NA,
                      LowerCI = NA,
                      UpperCI = NA,
                      d = NA,
                      nHE = NA,
                      nGenes = length(genes),
                      stringsAsFactors = F))
  }
}


# Load Data --------------------------------------------------------------------
  
setwd("/project/jakeweis_896/gtdb220/gtdb_genomes_reps_r220/")

paths <- readLines("ffn.files") #filepaths of the genomes you want to run prediction on

path_df <- data.frame(Accession = basename(paths) %>% substr(1,15),
                      Path = paths,
                      stringsAsFactors = F)

# Run --------------------------------------------------------------------------


pred_list <- mclapply(1:nrow(path_df),
                    predictGenome,
                    metadata = path_df,
                    mc.cores = 63)

pred_df <- do.call("rbind",pred_list)

save(pred_df,file="gtdb220_genome_predictions_20C.rda")
