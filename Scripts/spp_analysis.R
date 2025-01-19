
library(dplyr)
library(ggplot2)
library(data.table)
library(ggpubr)

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
P16_genomes <- tax$Accession[tax$V2 %in% P16_spp]


setwd("C:/Users/jlwei/Documents/bhet_code/Data/")
load("gtdb220_genome_predictions_20C.rda")
pred_rep <- pred_df
load("gtdb220_spp_predictions.rda")
pred_df$Genome <- substr(pred_df$Genome,3,17)
pred_rep <- pred_rep[,names(pred_df)]
pred_all <- rbind(pred_rep,pred_df) %>%
  subset(!grepl("Error",Genome)) %>%
  subset(nHE>=10)
pred_P16 <- pred_all %>%
  subset(Genome %in% P16_genomes)

names(tax) <- c("acc","Taxonomy","Genome")
pred_P16 <- merge.easy(pred_P16,tax,key="Genome")

dbcan_rep <- read.delim("dbcan_counts2_gtdb220.tsv",header = F,sep="\t")[,c(2,1)]
names(dbcan_rep) <- c("nCAZy","Genome")
dbcan_rep$Genome <- dbcan_rep$Genome %>%
  substr(4,18)
dbcan_spp <- read.delim("dbcan_counts_spp.tsv",header = F,sep=" ")
names(dbcan_spp) <- c("nCAZy","Genome")
dbcan_spp$Genome <- dbcan_spp$Genome %>%
  substr(1,15)
dbcan <- rbind(dbcan_rep,dbcan_spp)
pred_P16 <- merge.easy(pred_P16,dbcan,key="Genome")


setwd("C:/Users/jlwei/Documents/bhet_code/Data")
load("PCA_model_P16.rda")
# Run PCA, take PC1
pred_P16 <- pred_P16 %>%
  mutate(dCUB=as.numeric(dCUB),
         d=as.numeric(d),
         nHE=as.numeric(nHE),
         nGenes=as.numeric(nGenes))
pred_P16$AllN <- pred_P16$nCAZy/pred_P16$nGenes
x <- predict(x,pred_P16[,c("AllN","dCUB","nGenes")])
pred_P16$PC1 <- x[,1]
#Make it positive for plotting
pred_P16$IoC <- pred_P16$PC1-min(pred_P16$PC1)+1

pred_P16$Genus <- pred_P16$Taxonomy %>%
  gsub(pattern=";s__.*",replace="")

tax_keep <- names(table(pred_P16$Taxonomy))[table(pred_P16$Taxonomy)>=10]
tax_keep2 <- names(table(pred_P16$Taxonomy))[table(pred_P16$Taxonomy)>=100]
genus_keep <- names(table(pred_P16$Genus))[table(pred_P16$Genus)>=10]

ggplot(pred_P16 %>% subset(Taxonomy %in% tax_keep),
       aes(x=reorder(Taxonomy,dCUB),y=dCUB)) +
  geom_boxplot() +
  theme_pubclean() +
  geom_hline(yintercept=-0.08,lty=2)

ggplot(pred_P16 %>% subset(Taxonomy %in% tax_keep),
       aes(x=reorder(Taxonomy,nCAZy),y=nCAZy)) +
  geom_boxplot() +
  theme_pubclean() +
  geom_hline(yintercept=-0.08,lty=2)

p3 <- ggplot(pred_P16 %>% subset(Taxonomy %in% tax_keep),
       aes(x=reorder(Taxonomy,IoC),y=IoC)) +
  geom_boxplot() +
  theme_pubclean() +
  geom_hline(yintercept=c(2.12,3.03),lty=2) +
  ylab("Index of Copiotrophy") +
  xlab("Species") +
  theme(axis.text.x=element_blank(),axis.ticks.x = element_blank()) 

p3b <- ggplot(pred_P16 %>% subset(Genus %in% genus_keep),
             aes(x=reorder(Genus,IoC),y=IoC)) +
  geom_boxplot() +
  theme_pubclean() +
  geom_hline(yintercept=c(2.12,3.03),lty=2) +
  ylab("Index of Copiotrophy") +
  xlab("Species") +
  theme(axis.text.x=element_blank(),axis.ticks.x = element_blank()) 
p3b

p1 <- ggplot(pred_P16 %>% subset(Taxonomy %in% tax_keep2),
       aes(x=nCAZy,y=dCUB,
           color=Taxonomy %>% gsub(pattern=".*;s__",replace=""))) +
  geom_point(pch=19,size=3,alpha=0.25) +
  theme_pubclean() +
  stat_ellipse(type="norm",level=0.5,lwd=1.1) +
  theme(legend.position = "none") +
  xlab("Number of CAZymes") +
  ylab("Codon Usage Bias (dCUB)") +
  # scale_color_brewer(palette = "Set2") +
  geom_hline(yintercept=-0.08,lty=2)

p2 <- ggplot(pred_P16 %>% subset(Taxonomy %in% tax_keep2),
             aes(x=reorder(Taxonomy %>% gsub(pattern=".*;s__",replace=""),IoC),
                 color=Taxonomy %>% gsub(pattern=".*;s__",replace=""),
                 y=IoC)) +
  geom_boxplot(alpha=1,fill="lightgray") +
  geom_jitter(alpha=0.25,width=0.1,pch=19) +
  theme_pubclean() +
  theme(axis.text.x = element_text(angle = 80, hjust=1)) +
  xlab("") +
  ylab("Index of Copiotrophy")  +
  theme(legend.position = "none")+
  # scale_color_brewer(palette = "Set3") +
  geom_hline(yintercept=c(2.12,3.03),lty=2)

setwd("C:/Users/jlwei/Documents/bhet_code/Figures")
pdf("spp_dCUB_nCAZy.pdf",width=10,height=10)
ggarrange(p1,p2,ncol=2,labels=c("(a)","(b)"))
dev.off()

setwd("C:/Users/jlwei/Documents/bhet_code/Figures")
pdf("spp_IOC.pdf",width=12,height=10)
p3
dev.off()

setwd("C:/Users/jlwei/Documents/bhet_code/Figures")
pdf("genus_IOC.pdf",width=12,height=10)
p3b
dev.off()

