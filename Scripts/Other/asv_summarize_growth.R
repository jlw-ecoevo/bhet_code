#SUmmarize growth rate predictions by ASV

# Load packages and data -------------------------------------------------------

library(dplyr)
library(data.table)

merge.easy <- function(df1,df2,key){
  df1 <- data.table(df1,key=key)
  df2 <- data.table(df2,key=key)
  return(unique(merge(df1,df2,all.x=TRUE,by=.EACHI,allow.cartesian=TRUE)))
}

setwd("~/bhet_code/Data/")
load("gtdb_genome_predictions.rda")
pred_df <- pred_df %>% subset(d99<1e3)
asv <- read.delim("221118-1030_P16N-S_3.84-fold-18S-correction_merged_16S_18S_proportions.QCd.prok-nonphoautototrophic.tsv") %>%
  subset(select=c(OTU_ID,taxonomy)) %>%
  subset(grepl("g__",taxonomy)) %>%
  mutate(Genus = taxonomy %>%
           gsub(pattern = "s__.*", replace = ""))
names(asv)[1] <- c("ASV")
pred_df <- merge.easy(pred_df,asv,key="ASV")

# Growth Rates (ASVs) ----------------------------------------------------------

pred_df_asv <- pred_df %>%
  subset(nHE>=10) %>%
  group_by(ASV) %>%
  summarise(OGT99=mean(OGT99),
            CUBHE=mean(CUBHE),
            CUB=mean(CUB),
            dCUB=mean(dCUB),
            d99=mean(d99),
            nHE=mean(nHE),
            nGenes=mean(nGenes))

# Save ASV Data ----------------------------------------------------------------

setwd("~/bhet_code/Data/")
write.csv(pred_df_asv,file="asv_growth_data.csv")

#------------------------------------------------------------------------------
# 20C -------------------------------------------------------------------------
#------------------------------------------------------------------------------



setwd("~/bhet_code/Data/")
load("gtdb_genome_predictions_20C.rda")
pred_df <- pred_df %>% subset(d99<1e3)
asv <- read.delim("221118-1030_P16N-S_3.84-fold-18S-correction_merged_16S_18S_proportions.QCd.prok-nonphoautototrophic.tsv") %>%
  subset(select=c(OTU_ID,taxonomy)) %>%
  subset(grepl("g__",taxonomy)) %>%
  mutate(Genus = taxonomy %>%
           gsub(pattern = "s__.*", replace = ""))
names(asv)[1] <- c("ASV")
pred_df <- merge.easy(pred_df,asv,key="ASV")


pred_df_asv <- pred_df %>%
  subset(nHE>=10) %>%
  group_by(ASV) %>%
  summarise(OGT99=mean(OGT99),
            CUBHE=mean(CUBHE),
            CUB=mean(CUB),
            dCUB=mean(dCUB),
            d99=mean(d99),
            nHE=mean(nHE),
            nGenes=mean(nGenes))

setwd("~/bhet_code/Data/")
write.csv(pred_df_asv,file="asv_growth_data_20C.csv")
