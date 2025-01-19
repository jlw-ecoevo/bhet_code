#SUmmarize growth rate predictions by ASV

# Load packages and data -------------------------------------------------------

library(dplyr)
library(data.table)

merge.easy <- function(df1,df2,key){
  df1 <- data.table(df1,key=key)
  df2 <- data.table(df2,key=key)
  return(unique(merge(df1,df2,all.x=TRUE,by=.EACHI,allow.cartesian=TRUE)))
}

boxcoxTransform <- function(x, lambda, back_transform = F) {
  if (back_transform == TRUE) {
    (x*lambda +1)^(1/lambda)  %>% return()
  } else {
    (((x^lambda) - 1) / lambda) %>% return()
  }
}

setwd("C:/Users/jlwei/Documents/bhet_code/Data")
load("gtdb220_genome_predictions_20C.rda")
names(pred_df)[1] <- "Accession"
pred_df <- pred_df %>%
  mutate(d20C=as.numeric(d)) %>%
  subset(d20C<1e3) %>%
  subset(select=c(Accession,
                  CUBHE,
                  CUB,
                  dCUB,
                  CPB,
                  ConsistencyHE,
                  nHE,
                  nGenes,
                  d20C))
load("gRodon_sysdata.rda")
pred_mod <- gRodon_model_meta_temp_madin

tax_bins <- read.delim("ASVs_aggname_toJL_Feb2023.csv",sep=",") %>%
  subset(select=-c(X,taxonomy))
names(tax_bins)[1] <- "ASV"

load("asv_OGT99.rda")
P16_df <- merge.easy(asv_df,pred_df,key="Accession") %>%
  subset(!is.na(dCUB)) %>%
  mutate(OGT=OGT99) %>%
  subset(ASV %in% tax_bins$ASV)

length(unique(P16_df$Accession))

tax <- read.delim("bac120_taxonomy_r220.tsv",header=F)
tax$Accession <- substr(tax$V1,4,18)

P16_spp <- tax$V2[tax$Accession %in% unique(P16_df$Accession)] %>% unique()
P16_genomes <- setdiff(tax$Accession[tax$V2 %in% P16_spp],pred_df$Accession)

writeLines(P16_genomes,"P16_nonrep_genomes.txt")

