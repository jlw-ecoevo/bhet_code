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

load("asv_OGT99_I8-I9.rda")
I8I9_df <- merge.easy(asv_df,pred_df,key="Accession") %>%
  subset(!is.na(dCUB)) %>%
  mutate(OGT=OGT99)
load("asv_OGT99_GA02.rda")
GA02_df <- merge.easy(asv_df,pred_df,key="Accession") %>%
  subset(!is.na(dCUB)) %>%
  mutate(OGT=OGT99)
load("asv_OGT99.rda")
P16_df <- merge.easy(asv_df,pred_df,key="Accession") %>%
  subset(!is.na(dCUB)) %>%
  mutate(OGT=OGT99)

# temperature corrected growth -------------------------------------------------

I8I9_df <- I8I9_df %>%
  mutate(CUBHE=as.numeric(CUBHE),
         dCUB=as.numeric(dCUB),
         CPB=as.numeric(CPB),
         ConsistencyHE=as.numeric(ConsistencyHE),
         nHE=as.numeric(nHE),
         nGenes=as.numeric(nGenes),
         d20C=as.numeric(d20C),
         OGT=as.numeric(OGT))
GA02_df <- GA02_df %>%
  mutate(CUBHE=as.numeric(CUBHE),
         dCUB=as.numeric(dCUB),
         CPB=as.numeric(CPB),
         ConsistencyHE=as.numeric(ConsistencyHE),
         nHE=as.numeric(nHE),
         nGenes=as.numeric(nGenes),
         d20C=as.numeric(d20C),
         OGT=as.numeric(OGT))
P16_df <- P16_df %>%
  mutate(CUBHE=as.numeric(CUBHE),
         dCUB=as.numeric(dCUB),
         CPB=as.numeric(CPB),
         ConsistencyHE=as.numeric(ConsistencyHE),
         nHE=as.numeric(nHE),
         nGenes=as.numeric(nGenes),
         d20C=as.numeric(d20C),
         OGT=as.numeric(OGT))


# #gut-check
# x <- boxcoxTransform(x = predict.lm(object = pred_mod, newdata = pred_df),
#                 lambda = lambda_milc_madin,
#                 back_transform = T)
# plot(x,pred_df$d,log="xy")      

I8I9_df$d99 <- boxcoxTransform(x = predict.lm(object = pred_mod, 
                                              newdata = I8I9_df),
                             lambda = lambda_milc_madin,
                             back_transform = T)
GA02_df$d99 <- boxcoxTransform(x = predict.lm(object = pred_mod, 
                                              newdata = GA02_df),
                               lambda = lambda_milc_madin,
                               back_transform = T)
P16_df$d99 <- boxcoxTransform(x = predict.lm(object = pred_mod, 
                                             newdata = P16_df),
                               lambda = lambda_milc_madin,
                               back_transform = T)


# Growth Rates (ASVs) ----------------------------------------------------------

I8I9_df_asv <- I8I9_df %>%
  subset(nHE>=10) %>%
  group_by(ASV) %>%
  summarise(OGT99=mean(OGT99),
            CUBHE=mean(CUBHE),
            CUB=mean(CUB),
            dCUB=mean(dCUB),
            d99=mean(d99),
            d20C=mean(d20C),
            nHE=mean(nHE),
            nGenes=mean(nGenes))

GA02_df_asv <- GA02_df %>%
  subset(nHE>=10) %>%
  group_by(ASV) %>%
  summarise(OGT99=mean(OGT99),
            CUBHE=mean(CUBHE),
            CUB=mean(CUB),
            dCUB=mean(dCUB),
            d99=mean(d99),
            d20C=mean(d20C),
            nHE=mean(nHE),
            nGenes=mean(nGenes))

P16_df_asv <- P16_df %>%
  subset(nHE>=10) %>%
  group_by(ASV) %>%
  summarise(OGT99=mean(OGT99),
            CUBHE=mean(CUBHE),
            CUB=mean(CUB),
            dCUB=mean(dCUB),
            d99=mean(d99),
            d20C=mean(d20C),
            nHE=mean(nHE),
            nGenes=mean(nGenes))

# Save ASV Data ----------------------------------------------------------------

write.csv(I8I9_df_asv,file="asv_growth_data_I8-I9.csv")
write.csv(GA02_df_asv,file="asv_growth_data_GA02.csv")
write.csv(P16_df_asv,file="asv_growth_data_P16.csv")

I8I9_df$transect <- "I8I9"
GA02_df$transect <- "GA02"
P16_df$transect <- "P16"

pred_df <- rbind(I8I9_df,GA02_df,P16_df)
save(pred_df,file="gtdb220_genome_predictions.rda")

# Hits per genome --------------------------------------------------------------

P16_count <- P16_df %>%
  group_by(ASV) %>%
  summarize(nGenomes=length(Accession))

ggplot(P16_count,aes(x=nGenomes)) +
  geom_histogram() +
  scale_x_log10()

sum(P16_count$nGenomes>1)/nrow(P16_count)
