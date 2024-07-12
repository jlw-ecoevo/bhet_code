# JLW - 2023
#Estimate optimal growth temperature for an ASV from distributional data

# Load Packages ----------------------------------------------------------------

library(dplyr)
library(data.table)
library(ggplot2)
library(ggpubr)

# Functions --------------------------------------------------------------------

merge.easy <- function(df1,df2,key){
  df1 <- data.table(df1,key=key)
  df2 <- data.table(df2,key=key)
  return(unique(merge(df1,df2,all.x=TRUE,by=.EACHI,allow.cartesian=TRUE)))
}

# Load Data --------------------------------------------------------------------

setwd("C:/Users/jlwei/Documents/bhet_code/Data")
asv <- read.delim("221216-1006_I8-I9_4.95-fold-18S-correction_merged_16S_18S_proportions.QCd.prok-nonphoautototrophic.tsv")
meta <- read.delim("grump.metadata.I8S_I9N.tsv") %>%
  mutate(SID = SampleID)
asv_genome <- 
  read.delim("gtdb220/240705_I8-I9.prok-nonphoautototrophic.BLAST-95pcID-vs-GTDB-r220-allssu.tsv",
             header=F) %>%
  group_by(V1) %>%
  subset(select = c(V1,V2)) %>%
  mutate(V2 = substr(V2,4,18))
names(asv_genome) <- c("ASV","Accession")

names(asv) <- gsub("[.]","-",names(asv))
asv <- asv %>% subset(OTU_ID %in% asv_genome$ASV)


# Estimate OGT w/ q99 Method ---------------------------------------------------

temp_df <- data.frame(ASV=character(),
                      OGT99=numeric())
for(i in 1:nrow(asv)){
  q99 <- quantile(unlist(asv[i,-c(1:2)]),0.99)[1]
  s99 <- setdiff(names(asv[,asv[i,]>q99]),c("OTU_ID","taxonomy"))
  t99 <- meta$Temperature[meta$SID %in% s99]
  
  q95 <- quantile(unlist(asv[i,-c(1:2)]),0.95)[1]
  s95 <- setdiff(names(asv[,asv[i,]>q95]),c("OTU_ID","taxonomy"))
  t95 <- meta$Temperature[meta$SID %in% s95]
  
  temp_df <- rbind(temp_df,
                   data.frame(ASV = asv[i,1],
                              OGT99 = mean(t99),
                              nt99 = length(t99),
                              OGT95 = mean(t95),
                              nt95 = length(t95)))
}

plot(temp_df$OGT99,temp_df$OGT95)
# Merge w/ Genome Info ---------------------------------------------------------

asv_df <- merge.easy(asv_genome,temp_df,key="ASV") %>%
  subset(!is.na(OGT99))
save(asv_df,file="asv_OGT99_I8-I9.rda")

