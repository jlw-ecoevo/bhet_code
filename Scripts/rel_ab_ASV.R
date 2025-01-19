library(ggplot2)
library(ggpubr)

setwd("C:/Users/jlwei/Documents/bhet_code/Data")

asv <- read.delim("221118-1030_P16N-S_3.84-fold-18S-correction_merged_16S_18S_proportions.QCd.prok-nonphoautototrophic.tsv")
rownames(asv) <- asv$OTU_ID

asv_genome <- 
  read.delim("gtdb220/240705_P16N-S.prok-nonphoautototrophic.BLAST-95pcID-vs-GTDB-r220-allssu.tsv",
             header=F) %>%
  group_by(V1) %>%
  slice_max(V3, with_ties=F) %>%
  subset(select = c(V1,V2)) %>%
  mutate(V2 = substr(V2,4,18))
names(asv_genome) <- c("OTU_ID","Accession")

x <- data.frame(x=colSums(asv[unique(asv_genome$OTU_ID),-c(1:2)]))
median(x$x)


P16 <- ggplot(x,aes(x=x)) +
  geom_density() +
  theme_pubclean() +
  xlab("Proportion of Community") +
  geom_vline(xintercept=median(x$x),lty=2) +
  scale_x_continuous(breaks = round(seq(0, 1, by = 0.1),1)) +
  ggtitle("P16")




asv <- read.delim("221223-1135_GA02_7.37-fold-18S-correction_merged_16S_18S_proportions.QCd.prok-nonphoautototrophic.tsv")
rownames(asv) <- asv$OTU_ID

asv_genome <- 
  read.delim("gtdb220/240705_GA02.prok-nonphoautototrophic.BLAST-95pcID-vs-GTDB-r220-allssu.tsv",
             header=F) %>%
  group_by(V1) %>%
  subset(select = c(V1,V2)) %>%
  mutate(V2 = substr(V2,4,18))
names(asv_genome) <- c("ASV","Accession")

x <- data.frame(x=colSums(asv[unique(asv_genome$ASV),-c(1:2)]))
median(x$x)


GA02 <- ggplot(x,aes(x=x)) +
  geom_density() +
  theme_pubclean() +
  xlab("Proportion of Community") +
  geom_vline(xintercept=median(x$x),lty=2) +
  scale_x_continuous(breaks = round(seq(0, 1, by = 0.1),1))+
  ggtitle("GA02")


asv <- read.delim("221216-1006_I8-I9_4.95-fold-18S-correction_merged_16S_18S_proportions.QCd.prok-nonphoautototrophic.tsv")
rownames(asv) <- asv$OTU_ID

asv_genome <- 
  read.delim("gtdb220/240705_I8-I9.prok-nonphoautototrophic.BLAST-95pcID-vs-GTDB-r220-allssu.tsv",
             header=F) %>%
  group_by(V1) %>%
  subset(select = c(V1,V2)) %>%
  mutate(V2 = substr(V2,4,18))
names(asv_genome) <- c("ASV","Accession")

x <- data.frame(x=colSums(asv[unique(asv_genome$ASV),-c(1:2)]))
median(x$x)


I8I9 <- ggplot(x,aes(x=x)) +
  geom_density() +
  theme_pubclean() +
  xlab("Proportion of Community") +
  geom_vline(xintercept=median(x$x),lty=2) +
  scale_x_continuous(breaks = round(seq(0, 1, by = 0.1),1))+
  ggtitle("I8-I9")

setwd("C:/Users/jlwei/Documents/bhet_code/Figures")
pdf("Community_Covered_ASV_Match.pdf",width=10,height=5)
ggarrange(P16,GA02,I8I9,ncol=3,labels=c("(a)","(b)","(c)"))
dev.off()
