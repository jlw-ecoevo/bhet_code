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
setwd("C:/Users/jlwei/Documents/bhet_code/Data")
asv <- read.delim("221118-1030_P16N-S_3.84-fold-18S-correction_merged_16S_18S_proportions.QCd.prok-nonphoautototrophic.tsv")
#Make relative abundance
asv[,3:197] <- t(t(asv[,-c(1:2)])/colSums(asv[,-c(1:2)]))

# Match ASV to genome based on blast best hit
asv_genome <- 
  read.delim("gtdb220/240705_P16N-S.prok-nonphoautototrophic.BLAST-95pcID-vs-GTDB-r220-allssu.tsv",
             header=F) %>%
  group_by(V1) %>%
  subset(V3>=99) %>%
  # slice_max(V3, with_ties=F) %>%
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
genome_meanab <- data.frame(Accession=genome_ab$Accession,
                            MeanRelAb=rowMeans(genome_ab[,-1]))
genome_maxab <- data.frame(Accession=genome_ab$Accession,
                            MaxRelAb=rowMaxs(as.matrix(genome_ab[,-1])))

#Taxonomic clusters based on expert annotation
setwd("C:/Users/jlwei/Documents/bhet_code/Data")
tax_bins <- read.delim("ASVs_aggname_toJL_Feb2023.csv",sep=",") %>%
  subset(select=-c(X,taxonomy))
names(tax_bins)[1] <- "ASV"

load("gtdb220_genome_predictions_99PID.rda")
pred_df <- pred_df %>% subset(transect=="P16")
dbcan <- read.delim("dbcan_counts2_gtdb220.tsv",header=F,sep="\t") %>%
  mutate(nCAZy=V2,
         Accession=substr(V1,4,18)) %>%
  subset(select=c(Accession,nCAZy))

asv_genome <- pred_df %>%
  subset(select=c(Accession,ASV)) %>%
  merge.easy(.,tax_bins,key="ASV") %>%
  as.data.frame() %>%
  unique() %>%
  group_by(Accession,Jesse_aggname) %>%
  tally() %>%
  subset(!is.na(Jesse_aggname)) %>% 
  group_by(Accession) %>%
  filter(n == max(n)) %>%
  subset(select=-c(n))
gdf <- merge.easy(pred_df %>%
                    group_by(Accession) %>%
                    summarise_all(mean),
                  asv_genome,key="Accession")
gdf_cazy <- merge.easy(gdf,dbcan,key="Accession")

#Put it together for plotting
asv_cazy <- merge.easy(pred_df,dbcan,key="Accession") %>%
  # subset(nHE>=10) %>%
  group_by(ASV) %>%
  summarise(CAZy=mean(nCAZy,na.rm=T))
write.csv(asv_cazy,file="asv_cazy.csv")
plot_df <- gdf_cazy %>%
  subset(d99<1e3 &
           nHE>=10 &
           !is.na(Jesse_aggname)) %>%
  mutate(All1=nCAZy+1,
         dCUBp=-dCUB + 2*max(dCUB))
plot_df <- merge.easy(plot_df,genome_meanab,key="Accession")
plot_df <- merge.easy(plot_df,genome_maxab,key="Accession")
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

#calculate percent mags
genome_acc <- unique(plot_df$Accession)
mag_acc <- c(readLines("ar53_mags_marine.txt"),
             readLines("bac120_mags_marine.txt")) %>%
  substr(4,18)
table(genome_acc %in% mag_acc)/length(genome_acc)
table(genome_acc %in% mag_acc)

# Plot & Analyze ---------------------------------------------------------------

### Generate Index of Copiotropy -----------------------------------------------

setwd("C:/Users/jlwei/Documents/bhet_code/Data")
load("PCA_model_P16.rda")

# Run PCA, take PC1
plot_df$AllN <- plot_df$All/plot_df$nGenes
x <- predict(x,plot_df[,c("AllN","dCUB","nGenes")])
autoplot(x,loadings=T,loadings.label=T,loadings.label.vjust = 1.3)
plot_df$PC1 <- x[,1]
#Make it positive for plotting
plot_df$PC1t <- plot_df$PC1-min(plot_df$PC1)+1

#Cut out anything with missing nHE or predictions
plot_df <- plot_df %>% subset(d99<1e3 &
                                nHE>=10 & 
                                !is.na(dCUB))


### Index of Copiotrophy -------------------------------------------------------

#Clustering
mC <- Mclust(log10(plot_df$PC1t),
             prior = priorControl(functionName="defaultPrior", shrinkage=100),
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
                               sd=sqrt(mC$p$variance$sigmasq[1]))*p[2])


#Confidence limits
plot(10^mC$data,mC$uncertainty)
points(10^mC$data,mC$classification/10,col="red")
abline(v=2.12)
abline(v=3.03)
abline(h=0.1)
c_l <- max(10^mC$data[mC$uncertainty<0.1 & 10^mC$data<2.5])
c_h <- min(10^mC$data[mC$uncertainty<0.1 & 10^mC$data>2.5])

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
  geom_hline(yintercept=0.1,lty=2) +
  geom_vline(xintercept=c_h,lty=2) +
  geom_vline(xintercept=c_l,lty=2)
pCu

#Plot: Arrange IoC plots

ggarrange(pC,pCu,pCg,nrow=3,labels=c("(a)","(b)","(c)"))



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

setwd("C:/Users/jlwei/Documents/bhet_code/Figures")
pdf("Copiotrophy_ridges_99PID.pdf",width=10,height=10)
pCr
dev.off()
