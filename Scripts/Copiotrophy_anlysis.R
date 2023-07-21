# JLW - 2023
# Generate figures for heterotroph paper

# Load Packages & Helper Functions ---------------------------------------------

library(ggplot2)
library(ggpubr)
library(pheatmap)
library(dplyr)
library(data.table)
library(mclust)
library(matrixStats)
library(ggridges)
library(ggfortify)
library(corrplot)
library(ggcorrplot)
library(GGally)
library(patchwork)

merge.easy <- function(df1,df2,key){
  df1 <- data.table(df1,key=key)
  df2 <- data.table(df2,key=key)
  return(unique(merge(df1,df2,all.x=TRUE,by=.EACHI,allow.cartesian=TRUE)))
}

# Load Datasets ----------------------------------------------------------------

#ASV Transect abunmdance data
setwd("~/bhet_code/Data")
asv <- read.delim("221118-1030_P16N-S_3.84-fold-18S-correction_merged_16S_18S_proportions.QCd.prok-nonphoautototrophic.tsv")
#Make relative abundance
asv[,3:197] <- t(t(asv[,-c(1:2)])/colSums(asv[,-c(1:2)]))

# Match ASV to genome based on blast best hit
asv_genome <- 
  read.delim("221201_P16N-S.prok-nonphoautototrophic.BLAST-95pcID-vs-GTDB-r207-allproks.tsv",
             header=F) %>%
  group_by(V1) %>%
  slice_max(V3, with_ties=F) %>%
  subset(select = c(V1,V2)) %>%
  mutate(V2 = substr(V2,4,18))
names(asv_genome) <- c("OTU_ID","Accession")
asvg <- asv_genome
asv <- asv %>% subset(OTU_ID %in% asv_genome$OTU_ID)
asv$taxonomy <- NULL
genome_ab <- merge.easy(asv_genome,asv,key="OTU_ID") %>%
  subset(select=-c(OTU_ID)) %>%
  group_by(Accession) %>%
  summarise_all(sum)
genome_meanab <- data.frame(Genome=genome_ab$Accession,
                            MeanRelAb=rowMeans(genome_ab[,-1]))
genome_maxab <- data.frame(Genome=genome_ab$Accession,
                            MaxRelAb=rowMaxs(as.matrix(genome_ab[,-1])))

#Taxonomic clusters based on expert annotation
setwd("~/bhet_code/Data")
tax_bins <- read.delim("ASVs_aggname_toJL_Feb2023.csv",sep=",") %>%
  subset(select=-c(X,taxonomy))
names(tax_bins)[1] <- "ASV"

#METABOLIC output for genomes
setwd("~/bhet_code/Data/METABOLIC_result_each_spreadsheet/")
x1 <- read.table("METABOLIC_result_worksheet5_batch1.tsv",sep="\t",head=T)
x2 <- read.table("METABOLIC_result_worksheet5_batch2.tsv",sep="\t",head=T)
x <- merge.easy(x1,x2,key="CAZyme.ID")
xm <- x %>% subset(select=names(x)[grep("Hit.numbers",names(x))])
rownames(xm) <- x$CAZyme.ID

#Save CAZy data organized
setwd("~/bhet_code/Data")
write.csv(xm,file="CAZy_METABOLIC_BestHitGenomes.csv")

#CUB/growth data
setwd("~/bhet_code/Data")
load("gtdb_genome_predictions.rda")

asv_genome <- pred_df %>%
  subset(select=c(Genome,ASV)) %>%
  merge.easy(.,tax_bins,key="ASV") %>%
  as.data.frame() %>%
  unique() %>%
  group_by(Genome,Jesse_aggname) %>%
  tally() %>%
  subset(!is.na(Jesse_aggname)) %>% 
  group_by(Genome) %>%
  filter(n == max(n)) %>%
  subset(select=-c(n))
xml <- xm>0
gdf <- merge.easy(pred_df %>%
                    group_by(Genome) %>%
                    summarise_all(mean),
                  asv_genome,key="Genome")
gdf_cazy <- data.frame(Genome = names(xm) %>% 
                         substr(4,18),
                       All = colSums(xm),
                       GH = xm %>% 
                         subset(grepl("GH",rownames(xm))) %>% 
                         colSums,
                       PL = xm %>% 
                         subset(grepl("PL",rownames(xm))) %>% 
                         colSums)

#Put it together for plotting
setwd("~/bhet_code/Data")
asv_cazy <- merge.easy(pred_df,gdf_cazy,key="Genome") %>%
  # subset(nHE>=10) %>%
  group_by(ASV) %>%
  summarise(CAZy=mean(All,na.rm=T),
            CAZy_GH=mean(GH,na.rm=T),
            CAZy_PL=mean(PL,na.rm=T))
write.csv(asv_cazy,file="asv_cazy.csv")
gdf_cazy <- merge.easy(gdf_cazy,gdf,key="Genome")
plot_df <- gdf_cazy %>%
  subset(d99<1e3 &
           nHE>=10 &
           !is.na(Jesse_aggname)) %>%
  mutate(All1=All+1,
         dCUBp=-dCUB + 2*max(dCUB))
plot_df <- merge.easy(plot_df,genome_meanab,key="Genome")
plot_df <- merge.easy(plot_df,genome_maxab,key="Genome")
plot_df$Jesse_aggname <- gsub(".__","",plot_df$Jesse_aggname) %>%
  gsub(pattern="Marinimicrobia_[(]SAR406_clade[)]",replace="SAR406_clade") %>%
  gsub(pattern="Pseudomonadales_multiplegrep",replace="Pseudomonadales, excl. SUP05 and SAR86")

plot_df_means <- plot_df %>% 
  group_by(Jesse_aggname) %>%
  summarise(All1=weightedMean(All1,w=MeanRelAb),
            dCUBp=weightedMean(dCUBp,w=MeanRelAb),
            All_2sd=2*sd(All1),
            dCUBp_2sd=2*sd(dCUBp))

#Filter non-het ASV
plot_df <- plot_df %>% subset(!Jesse_aggname %in%c("Nitroso","Nitrosp"))


# Plot & Analyze ---------------------------------------------------------------

### Metabolism Summary per guild -----------------------------------------------

setwd("~/bhet_code/Figures")
pdf("guild_metabolism.pdf",width=14,height=12)
ggplot(plot_df,
       aes(x=All1,y=dCUBp,size=MaxRelAb)) +
  geom_point(pch=21,alpha=0.5,fill="black") +
  geom_point(data=plot_df_means,
             aes(x=All1,y=dCUBp),
             color="red",alpha=0.5,
             size=5) +
  # geom_density_2d(na.rm=T,n=20,color="red",lwd=0.1) +
  scale_fill_viridis_c() +
  scale_x_log10() +
  scale_y_log10() +
  theme_bw() +
  geom_hline(yintercept=0.08+2*max(plot_df$dCUB),lty=2) +
  geom_vline(xintercept = 8,lty=2) +
  xlab("Number of Carbohydrate Active Enzymes") +
  theme(legend.position = "right") +
  ylab("Codon Usage Bias") +
  labs(size="Maximum Relative\nAbundance Among\nHeterotrophs") +
  facet_wrap(~Jesse_aggname) 
dev.off()

setwd("~/bhet_code/Figures")
pdf("guild_metabolism_all.pdf",width=14,height=12)
ggplot(plot_df,
       aes(x=All1,y=dCUBp,size=MaxRelAb,fill=Jesse_aggname)) +
  geom_point(pch=21,alpha=0.75) +
  geom_point(data=plot_df_means,
             aes(x=All1,y=dCUBp),
             color="red",
             fill="red",
             alpha=0.5,
             size=5) +
  scale_x_log10() +
  scale_y_log10() +
  theme_bw() +
  geom_hline(yintercept=0.08+2*max(plot_df$dCUB),lty=2) +
  scale_color_brewer(palette="Set2") +
  geom_vline(xintercept = 8,lty=2) +
  xlab("Number of Carbohydrate Active Enzymes") +
  theme(legend.position = "right") +
  ylab("Codon Usage Bias") +
  labs(size="Maximum Relative\nAbundance Among\nHeterotrophs")
dev.off()

### Generate Index of Copiotropy -----------------------------------------------

# Run PCA, take PC1
plot_df$AllN <- plot_df$All/plot_df$nGenes
x <- prcomp(plot_df[,c("AllN","dCUB","nGenes")],scale=T)
autoplot(x,loadings=T,loadings.label=T,loadings.label.vjust = 1.3)
plot_df$PC1 <- x$x[,1]
#Make it positive for plotting
plot_df$PC1t <- plot_df$PC1-min(plot_df$PC1)+1

#Cut out anything with missing nHE or predictions
plot_df <- plot_df %>% subset(d99<1e3 &
                            nHE>=10 & 
                            !is.na(dCUB))

### Index of Copiotrophy -------------------------------------------------------

#Clustering
mC <- Mclust(log10(plot_df$PC1t),
             prior = priorControl(functionName="defaultPrior", shrinkage=10),
             verbose=T)
mC$parameters$mean
mC$parameters$variance$sigmasq
p <- (colSums(mC$z)/sum(mC$z))
x_seq <- seq(-.5,1,0.01)
cl_df <- data.frame(x=10^x_seq,
                     cl1=dnorm(x_seq,
                               mean=mC$parameters$mean[1],
                               sd=sqrt(mC$p$variance$sigmasq[1]))*p[1],
                     cl2=dnorm(x_seq,
                               mean=mC$parameters$mean[2],
                               sd=sqrt(mC$p$variance$sigmasq[2]))*p[2])


#Confidence limits
plot(10^mC$data,mC$uncertainty)
abline(v=1.96)
abline(v=2.95)
abline(h=0.05)
c_l <- max(10^mC$data[mC$uncertainty<0.05 & 10^mC$data<2])
c_h <- min(10^mC$data[mC$uncertainty<0.05 & 10^mC$data>2])

#Plot: Cluster Plot for IoC
pC <- ggplot() + 
  scale_x_log10() +
  geom_density(data=plot_df,aes(x=PC1t),lwd=2) + 
  geom_polygon(data=cl_df,aes(x=x, y=cl1, fill="Cluster 1"),
               alpha=0.75, color="black") +
  geom_polygon(data=cl_df,aes(x=x, y=cl2, fill="Cluster 2"), 
               alpha=0.75, color="black") +
  theme_pubclean() + scale_color_brewer(palette = "Set1")  +
  xlab("Index of Copiotrophy")  + 
  ylab("Density") +
  theme(legend.position = "right") + 
  labs(fill="") +
  # geom_vline(xintercept=1.75,lty=2) +
  geom_vline(xintercept=c_l,lty=2) +
  geom_vline(xintercept=c_h,lty=2) +
  scale_fill_brewer(palette="Set2")
pC

#Plot: Guild IoC values
pCg <- ggplot(plot_df,aes(y=PC1t,x=reorder(Jesse_aggname,PC1,FUN=median))) +
  geom_violin(width=1.25,color="lightblue",fill="lightblue") +
  # geom_boxplot(width=0.1) +
  geom_point(aes(size=MaxRelAb),alpha=0.25) +
  theme_pubclean() +
  scale_y_log10() +
  # geom_hline(yintercept=1.75,lty=2) +
  geom_hline(yintercept=c_h,lty=2) +
  geom_hline(yintercept=c_l,lty=2) +
  ylab("Index of Copiotrophy")  + 
  xlab("")  + 
  labs(size="Maximum Relative\nAbundance Among\nHeterotrophs") +
  theme(legend.position = "right") +
  coord_flip()
pCg

#Plot: Cluster confidence limits
pCu <- ggplot(data.frame(x=10^mC$data,
                         y=mC$uncertainty),
              aes(x=x,y=y)) +
  geom_line() +
  scale_x_log10() +
  ylab("Uncertainty of Classification") +
  xlab("Index of Copiotrophy") +
  theme_pubclean() +
  geom_hline(yintercept=0.05,lty=2) +
  geom_vline(xintercept=c_h,lty=2) +
  geom_vline(xintercept=c_l,lty=2)
pCu

#Plot: Arrange IoC plots
setwd("~/bhet_code/Figures")
pdf("Copiotrophy_classification_uncertainty.pdf",width=10,height=15)
ggarrange(pC,pCu,pCg,nrow=3,labels=c("(a)","(b)","(c)"))
dev.off()


#Plot: Ridgeline plot with quartiles for IoC of guilds
pCr <- ggplot(plot_df,
              aes(x=PC1t, y=reorder(Jesse_aggname,PC1,FUN=median), 
                  fill = factor(stat(quantile)))) + 
  stat_density_ridges(
    geom = "density_ridges_gradient", calc_ecdf = TRUE,
    quantiles = 4, quantile_lines = TRUE,
    jittered_points = TRUE,
    position = position_points_jitter(width = 0.05, height = 0),
    scale=0.9,
    point_shape = '|', point_size = 3, point_alpha = 1, alpha = 0.5
  ) +
  scale_fill_manual(values=c("white","lightgray","darkgray","black")) +
  labs(fill="Quartile") + 
  theme_pubclean() + 
  scale_x_log10() +
  xlab("Index of Copiotrophy") + 
  ylab("") +
  # geom_vline(xintercept=1.75,lty=2)
  geom_vline(xintercept=c_h,lty=2) +
  geom_vline(xintercept=c_l,lty=2)
pCr

setwd("~/bhet_code/Figures")
pdf("Copiotrophy_ridges.pdf",width=10,height=10)
pCr
dev.off()

#Plot: Arrange clustering with ridges
setwd("~/bhet_code/Figures")
pdf("Copiotrophy_ridges_plusdist.pdf",width=10,height=12)
ggarrange(pC + theme(legend.position=c(0.2,0.8)),
          pCr + theme(legend.position = c(-0.3,0.9)),
          nrow=2,heights=c(1,2),labels=c("(a)","(b)"))
dev.off()


#Write output ASV and IoC data, along with IoC inputs
CI <- plot_df %>% subset(select=c(Genome,PC1t,dCUB,nGenes,All))
names(CI)[5] <- "nCAZy"
names(asvg) <- c("ASV","Genome")
asv_CI <- merge.easy(asvg,CI,key="Genome") %>%
  subset(!is.na(PC1t)) %>%
  subset(select=c(ASV,PC1t,dCUB,nGenes,nCAZy))
names(asv_CI)[2] <- "CopiotrophyIndex"
setwd("~/bhet_code/Data")
write.table(asv_CI,file="Copiotrophy_index.tsv")
write.csv(asv_CI,file="Copiotrophy_index.csv")

### IoC Components -------------------------------------------------------------


cor_df <- plot_df %>% subset(select=c(dCUBp,All1,nGenes)) %>%
                               mutate_all(log10)
names(cor_df) <- c("Codon Usage Bias (transformed)",
                   "Number of CAZymes (transformed)",
                   "Number of Genes  (transformed)")

ggpairs(cor_df) + theme_pubclean()

plot_ls <- ggpairs(cor_df) + theme_pubclean()


setwd("~/bhet_code/Figures")
pdf("Copiotrophy_dist.pdf",width=20,height=10,onefile=FALSE)
wrap_elements(ggmatrix_gtable(plot_ls)) + 
  wrap_elements( ggarrange(pC,pCg,nrow=2)) + 
  plot_annotation(tag_levels = 'a')
dev.off()

setwd("~/bhet_code/Figures")
pdf("Copiotrophy_dist_a.pdf",width=10,height=10,onefile=FALSE)
wrap_elements(ggmatrix_gtable(plot_ls))
dev.off()

setwd("~/bhet_code/Figures")
pdf("Copiotrophy_dist_b.pdf",width=10,height=5,onefile=FALSE)
pC
dev.off()

setwd("~/bhet_code/Figures")
pdf("Copiotrophy_dist_c.pdf",width=10,height=5,onefile=FALSE)
pCg
dev.off()
