#!/usr/bin/Rscript

#library(reshape)
library(RColorBrewer)
library("seqinr")
library(ggtree)
library(magrittr)
library(BiocManager)
library(lmerTest)
library(lme4)
library(cowplot)
library(microshades)
library(ggpmisc)
library(MuMIn)

source("R/1_Lab_filter.R")

# Do ASV match beween MA and SA?
rownames(Eim2@tax_table) %in% rownames(Eim@tax_table)
colnames(Eim@otu_table)[1] == colnames(Eim2@otu_table)[1]
colnames(Eim@otu_table)[2] == colnames(Eim2@otu_table)[2]

#plotting individual ASV's
MA.e <- psmelt(Eim)
SA.e <- psmelt(Eim2)

Eim2.g <- tax_glom(Eim2, taxrank="Genus")
SA.e.g <- psmelt(Eim2.g)
Eim.g <- tax_glom(Eim, taxrank="Genus")
MA.e.g <- psmelt(Eim.g)

SA.e$ASV <- "ASV"
SA.e$ASV[which(SA.e$OTU==colnames(Eim2@otu_table)[4])] <- "ASV4"
SA.e$ASV[which(SA.e$OTU==colnames(Eim2@otu_table)[3])] <- "ASV3"
SA.e$ASV[which(SA.e$OTU==colnames(Eim2@otu_table)[2])] <- "ASV2"
SA.e$ASV[which(SA.e$OTU==colnames(Eim2@otu_table)[1])] <- "ASV1"

MA.e$ASV <- "ASV"
MA.e$ASV[which(MA.e$OTU==colnames(Eim@otu_table)[2])] <- "ASV2"
MA.e$ASV[which(MA.e$OTU==colnames(Eim@otu_table)[1])] <- "ASV1"

SA.e2 <- SA.e[which(SA.e$OTU==colnames(Eim2@otu_table)[2]),]
SA.e1 <- SA.e[which(SA.e$OTU==colnames(Eim2@otu_table)[1]),]

## dataset with asv1 and 2
MA.e2 <- MA.e[which(MA.e$OTU==colnames(Eim@otu_table)[2]),]
MA.e1 <- MA.e[which(MA.e$OTU==colnames(Eim@otu_table)[1]),]

# No GC and ASV abundance zeros
SA.e.0 <- SA.e[SA.e$Genome_copies_ngDNA>0,]
SA.e.0 <- SA.e.0[SA.e.0$Abundance>0,]
SA.e.0 <-SA.e.0[!is.na(SA.e.0$Abundance),]

MA.e.0 <- MA.e[MA.e$Genome_copies_ngDNA>0,]
MA.e.0 <- MA.e.0[MA.e.0$Abundance>0,]
MA.e.0 <-MA.e.0[!is.na(MA.e.0$Abundance),]

## Let's do the correlations first
SA.ASV <- SA.e[SA.e$Genome_copies_ngDNA>0,]
SA.e.g0 <- SA.e.g[SA.e.g$Genome_copies_ngDNA>0,]
SA.e.g0 <- SA.e.g0[SA.e.g0$Abundance>0,]
SA.e.g0 <- SA.e.g0[!is.na(SA.e.g0$Genome_copies_ngDNA),]

MA.e.g0 <- MA.e.g[MA.e.g$Genome_copies_ngDNA>0,]
MA.e.g0 <- MA.e.g0[MA.e.g0$Abundance>0,]
MA.e.g0 <- MA.e.g0[!is.na(MA.e.g0$Genome_copies_ngDNA),]

# sanity test here
cor.test(log(SA.e.g0$Abundance), log(SA.e.g0$Genome_copies_ngDNA))
cor.test(log(MA.e.g0$Abundance), log(MA.e.g0$Genome_copies_ngDNA))

SA.e10 <- SA.e1[SA.e1$Genome_copies_ngDNA>0,]
SA.e10 <- SA.e10[!is.na(SA.e10$Genome_copies_ngDNA),]
SA.e10 <- SA.e10[SA.e10$Abundance>0,]

SA.e20 <- SA.e2[SA.e2$Genome_copies_ngDNA>0,]
SA.e20 <- SA.e20[!is.na(SA.e10$Genome_copies_ngDNA),]
SA.e20 <- SA.e20[SA.e20$Abundance>0,]

SA.e30 <- SA.e3[SA.e3$Genome_copies_ngDNA>0,]
SA.e30 <- SA.e30[!is.na(SA.e30$Genome_copies_ngDNA),]
SA.e30 <- SA.e30[SA.e30$Abundance>0,]

SA.e40 <- SA.e4[SA.e4$Genome_copies_ngDNA>0,]
SA.e40 <- SA.e40[!is.na(SA.e40$Genome_copies_ngDNA),]
SA.e40 <- SA.e40[SA.e40$Abundance>0,]

### now same for MA
MA.e10 <- MA.e1[MA.e1$Genome_copies_ngDNA>0,]
MA.e10 <- MA.e10[!is.na(MA.e10$Genome_copies_ngDNA),]
MA.e10 <- MA.e10[MA.e10$Abundance>0,]

MA.e20 <- MA.e2[MA.e2$Genome_copies_ngDNA>0,]
MA.e20 <- MA.e20[!is.na(MA.e10$Genome_copies_ngDNA),]
MA.e20 <- MA.e20[MA.e20$Abundance>0,]

####
#### Checking specificity
summary(SA.e.g$Abundance[SA.e.g$Genome_copies_ngDNA==0]>0) # 13 false positives
summary(MA.e.g$Abundance[MA.e.g$Genome_copies_ngDNA==0]>0) # 3 false positives

# relative abundances of true positives and false positives
SA.TP <- SA.e.g$Abundance[SA.e.g$Genome_copies_ngDNA>0]
SA.TP <- SA.TP[SA.TP>0]
mean(SA.TP, na.rm=TRUE)
sd(SA.TP, na.rm=TRUE)

SA.FP <- SA.e.g$Abundance[SA.e.g$Genome_copies_ngDNA==0]
SA.FP <- SA.FP[SA.FP>0]
mean(SA.FP, na.rm=TRUE)
sd(SA.FP, na.rm=TRUE)

MA.TP <- MA.e.g$Abundance[MA.e.g$Genome_copies_ngDNA>0]
MA.TP <- MA.TP[MA.TP>0]
mean(MA.TP, na.rm=TRUE)
sd(MA.TP, na.rm=TRUE)

MA.FP <- MA.e.g$Abundance[MA.e.g$Genome_copies_ngDNA==0]
MA.FP <- MA.FP[MA.FP>0]
mean(MA.FP, na.rm=TRUE)
sd(MA.FP, na.rm=TRUE)

################ Plotting
## defining colour for dpi
coul <- c(microshades_palette("micro_purple", 4, lightest=FALSE), microshades_palette("micro_orange", 3, lightest=FALSE), microshades_palette("micro_brown", 3, lightest=FALSE))
coul2 <- coul[-1]

# plotting with abundance and prevalence filter correlation
cor.sa <- cor.test(log(SA.e.g0$Abundance), log(SA.e.g0$Genome_copies_ng), method="pearson")
plot_SA_all <- ggplot(SA.e.g0, aes(x=log(Abundance), y=log(Genome_copies_ngDNA)))+
    geom_point(shape=21, size=2.5, aes(fill= dpi), color= "white", alpha=0.8)+
    scale_fill_manual(values=coul)+
    ylab("Eimeria genome copies (log)")+
    xlab("Eimeria ASV abundance (log)")+
    ggtitle("Standard PCR: Eimeria quantification")+
    guides(fill=guide_legend(ncol=2))+
    annotate(geom="text", x=min(log(SA.e.g0$Abundance)), y=max(log(SA.e.g0$Genome_copies_ngDNA)), hjust=0.05, label=paste("Pearson's rho=", round(cor.sa$estimate, 2), ", p<0.001, ", "df= ", cor.sa$parameter, sep=""))+
        theme_bw(base_size=10)+
    theme(axis.title.x = element_text(vjust = 0, size = 12),
          axis.title.y = element_text(vjust = 2, size = 12),
          plot.title=element_text(face="bold",
                                  margin=margin(10,0,10,0),
                                  size=12))

ASV.col <-  c("dodgerblue4",
              "darkolivegreen4",
              "darkorchid3",
               "goldenrod1")
ASV_ab <- ggplot(SA.e, aes(x=Abundance, y=ASV))+
    geom_jitter( position=position_jitter(height=0.2), shape=21, size=2.5, aes(fill= ASV), color= "white", alpha=0.7)+
    scale_fill_manual(values=ASV.col)+
    stat_summary(fun=mean, geom="point", colour="black", shape=23, size=5, aes(fill= ASV))+
    ylab("")+
    xlab("Eimeria ASV abundance")+
    ggtitle("Standard PCR: Abundance of individual ASVs")+
    coord_cartesian(xlim=c(0,1))+
    theme_bw(base_size=10)+
    theme(axis.title.x = element_text(vjust = 0, size = 12),
          axis.title.y = element_text(vjust = 2, size = 12),
          plot.title=element_text(face="bold",
                                  margin=margin(10,0,10,0),
                                  size=12),
          legend.position="none")


Figure2 <- cowplot::plot_grid(plot_SA_all, ASV_ab, nrow=2, labels="auto")
ggplot2::ggsave("fig/Figure2.pdf", Figure2, width = 6, height = 8, dpi = 300)
ggplot2::ggsave("fig/Figure2.png", Figure2, width = 6, height = 8, dpi = 300)

############### Do all ASVs explain Eimeria genome copies?
################ preparing dataset for regressions
# removeing zeros from qPCR
# can't remove zeros from indidual ASV's since then I end up with a handful of samples
Eim2.0 = phyloseq::prune_samples(Eim2@sam_data$Genome_copies_ngDNA>0, Eim2)
Eim.0 = phyloseq::prune_samples(Eim@sam_data$Genome_copies_ngDNA>0, Eim)

## OK, I will now do a linear model with only postive PCR samples and positive sequencing reads (total Eimeria sums) and include the sum of all Eimeria ASV's.

SA.df <- data.frame(Eim2.0@sam_data$Genome_copies_ngDNA, Eim2.0@otu_table[,1], Eim2.0@otu_table[,2],Eim2.0@otu_table[,3],Eim2.0@otu_table[,4],Eim2.0@sam_data$dpi, Eim2.0@sam_data$EH_ID)
names(SA.df) <- c("Genome_copies", "ASV1", "ASV2", "ASV3", "ASV4", "dpi", "EH_ID")
SA.df$TotalE <- SA.df$ASV1+SA.df$ASV2+SA.df$ASV3+SA.df$ASV4
SA.df <- SA.df[SA.df$TotalE>0,]
nrow(SA.df)
lmm.sa <- lmer(data=SA.df, log(Genome_copies)~log(1+ASV1)+log(1+ASV2)+log(1+ASV3)+log(1+ASV4)+ (1|dpi))
lmm.sa0 <- lmer(data=SA.df, log(Genome_copies)~1+ (1|dpi))

ranova(lmm.sa)
summary(lmm.sa)

# calculating r-squared following Nakagawa and Schielzeth 2012
# Extraction of fitted value for the alternative model
Fixed <- fixef(lmm.sa)[2]*getME(lmm.sa, "X")[,2]+fixef(lmm.sa)[3]*getME(lmm.sa, "X")[,3]+fixef(lmm.sa)[4]*getME(lmm.sa, "X")[,4]++fixef(lmm.sa)[5]*getME(lmm.sa, "X")[,5]

# Calculation of the variance in fitted values
VarF <- var(Fixed)
# R2GLMM(m) - marginal R2GLMM
# VarCorr() extracts variance components
# attr(VarCorr(lmer.model),'sc')^2 extracts the residual variance
VarF/(VarF + VarCorr(lmm.sa)$dpi[1] + attr(VarCorr(lmm.sa), "sc")^2)

# R2GLMM(c) - conditional R2GLMM for full model
(VarF + VarCorr(lmm.sa)$dpi[1])/(VarF + VarCorr(lmm.sa)$dpi[1] + (attr(VarCorr(lmm.sa), "sc")^2))

### or we do it with a package...
r.squaredGLMM(lmm.sa) ## well, that was easy

#plot(lmm.sa)
anova(lmm.sa, test="LRT")
anova(lmm.sa, lmm.sa0)

## plotting estimated random effects for each dpi and its interval estimate
#library(merTools)
#plotREsim(REsim(lmm.sa))
#library(sjPlot)
#sjPlot::plot_model(lmm.sa)

#### for multiamplicon
MA.df <- data.frame(Eim.0@sam_data$Genome_copies_ngDNA, Eim.0@otu_table[,1], Eim.0@otu_table[,2],Eim.0@sam_data$dpi)
names(MA.df) <- c("Genome_copies", "ASV1", "ASV2", "dpi")
MA.df$TotalE <- MA.df$ASV1+MA.df$ASV2
MA.df <- MA.df[MA.df$TotalE>0,]
nrow(MA.df)

lmm.ma <- lmer(data=MA.df, log(Genome_copies)~log(1+ASV1)+log(1+ASV2)+(1|dpi))
summary(lmm.ma)
lmm.ma0 <- lmer(data=MA.df, log(Genome_copies)~1+(1|dpi))

# LR tests for significance
anova(lmm.ma, test="LRT")
anova(lmm.ma, lmm.ma0)
ranova(lmm.ma)

## marginal and conditional r2
r.squaredGLMM(lmm.ma)

#### Plot by EH_ID, only positive animals
col <- c(microshades_palette("micro_purple", 4, lightest=FALSE), microshades_palette("micro_orange", 4, lightest=FALSE), microshades_palette("micro_brown", 3, lightest=FALSE))

SA2 <- ggplot(SA.e2, aes(x=dpi, y=log(1+Abundance)))+
    geom_point(shape=21, size=2.5, aes(fill= dpi), color= "white", alpha=0.7)+
    scale_fill_manual(values=col)+
    stat_summary(fun=mean, geom="point", colour="black", shape=23, size=5, aes(fill= dpi))+
    geom_line(aes(group=EH_ID), alpha=0.05)+
    ylab("Eimeria ASV2 abundance, log(+1)")+
    xlab("Day post infection")+
    ggtitle("Standard PCR ASV2")+
    coord_cartesian(ylim=c(0,0.1))+
    theme_bw(base_size=10)+
    theme(axis.title.x = element_text(vjust = 0, size = 12),
          axis.title.y = element_text(vjust = 2, size = 12),
          plot.title=element_text(face="bold",
                                  margin=margin(10,0,10,0),
                                  size=12),
          legend.position="none")


SA1 <- ggplot(SA.e1, aes(x=dpi, y=log(1+Abundance)))+
    geom_point(shape=21, size=2.5, aes(fill= dpi), color= "white", alpha=0.7)+
    scale_fill_manual(values=col)+
    stat_summary(fun=mean, geom="point", colour="black", shape=23, size=5, aes(fill= dpi))+
    geom_line(aes(group=EH_ID), alpha=0.05)+
    ylab("Eimeria ASV1 abundance, log(+1)")+
    xlab("Day post infection")+
    ggtitle("Standard PCR ASV1")+
    coord_cartesian(ylim=c(0,1))+
    theme_bw(base_size=10)+
    theme(axis.title.x = element_text(vjust = 0, size = 12),
          axis.title.y = element_text(vjust = 2, size = 12),
          plot.title=element_text(face="bold",
                                  margin=margin(10,0,10,0),
                                  size=12),
          legend.position="none")

MA2 <- ggplot(MA.e2, aes(x=dpi, y=log(1+Abundance)))+
    geom_point(shape=21, size=2.5, aes(fill= dpi), color= "white", alpha=0.7)+
    scale_fill_manual(values=col)+
    stat_summary(fun=mean, geom="point", colour="black", shape=23, size=5, aes(fill= dpi))+
    geom_line(aes(group=EH_ID), alpha=0.05)+
    ylab("Eimeria ASV2 abundance, log(+1)")+
    xlab("Day post infection")+
    ggtitle("Microfluidics PCR ASV2")+
    coord_cartesian(ylim=c(0,0.1))+
    theme_bw(base_size=10)+
    theme(axis.title.x = element_text(vjust = 0, size = 12),
          axis.title.y = element_text(vjust = 2, size = 12),
          plot.title=element_text(face="bold",
                                  margin=margin(10,0,10,0),
                                  size=12),
          legend.position="none")

MA1 <- ggplot(MA.e1, aes(x=dpi, y=log(1+Abundance)))+
    geom_point(shape=21, size=2.5, aes(fill= dpi), color= "white", alpha=0.7)+
    scale_fill_manual(values=col)+
    stat_summary(fun=mean, geom="point", colour="black", shape=23, size=5, aes(fill= dpi))+
    geom_line(aes(group=EH_ID), alpha=0.05)+
    ylab("Eimeria ASV1 abundance, log(+1)")+
    xlab("Day post infection")+
    ggtitle("Microfluidics PCR ASV1")+
        coord_cartesian(ylim=c(0,1))+
    theme_bw(base_size=10)+
    theme(axis.title.x = element_text(vjust = 0, size = 12),
          axis.title.y = element_text(vjust = 2, size = 12),
          plot.title=element_text(face="bold",
                                  margin=margin(10,0,10,0),
                                  size=12),
          legend.position="none")

FigureS3 <-cowplot::plot_grid(SA1, SA2, MA1, MA2,
              align="vh",
              labels="auto",
              nrow=2)

ggplot2::ggsave("fig/FigureS3.pdf", FigureS3, width = 10, height = 6, dpi = 300)
ggplot2::ggsave("fig/FigureS3.png", FigureS3, width = 10, height = 6, dpi = 300)


############### comparison between standard and microfluidics ASVs
ma.sa_1 <- SA.e[, c("ASV", "Abundance", "EH_ID", "dpi", "labels")]
ma.sa_2 <- MA.e[, c("ASV", "Abundance", "labels")]
ma.sa <- merge(ma.sa_1, ma.sa_2, by=c("ASV", "labels"))

cor.test(ma.sa[ma.sa$ASV=="ASV1","Abundance.x"],ma.sa[ma.sa$ASV=="ASV1","Abundance.y"], method="pearson")

cor.test(ma.sa[ma.sa$ASV=="ASV2","Abundance.x"],ma.sa[ma.sa$ASV=="ASV2","Abundance.y"], method="pearson")

ASV1.c <-ggplot(ma.sa[ma.sa$ASV=="ASV1",], aes(x=Abundance.x, y=Abundance.y))+
    geom_point(colour="gray40", alpha=0.7)+
    xlab("Standard PCR ASV1 abundance")+
    ylab("Microfluidics PCR ASV1 abundance")+
    stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")),
                                          parse = TRUE) +  
    stat_poly_line(method="lm", color="black", fill="firebrick")+
    coord_cartesian(ylim=c(0,1), xlim=c(0,1))+
    theme_bw(base_size=10)+
    theme(axis.title.x = element_text(vjust = 0, size = 12),
          axis.title.y = element_text(vjust = 2, size = 12))

ASV2.c <-ggplot(ma.sa[ma.sa$ASV=="ASV2",], aes(x=Abundance.x, y=Abundance.y))+
    geom_point(colour="gray40", alpha=0.7)+
#    scale_fill_manual(values=c("#009E73", "mediumvioletred"), name="")+
    xlab("Standard PCR ASV2 abundance")+
    ylab("Microfluidics PCR ASV2 abundance")+
    stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")),
                                          parse = TRUE) +  
    stat_poly_line(method="lm", color="black", fill="firebrick")+
    coord_cartesian(ylim=c(0,0.1), xlim=c(0,0.1))+
    theme_bw(base_size=10)+
    theme(axis.title.x = element_text(vjust = 0, size = 12),
          axis.title.y = element_text(vjust = 2, size = 12))

Fig1 <- plot_grid(ASV1.c, ASV2.c, labels="auto", nrow=1, align="hv")
Fig1

ggplot2::ggsave(file="fig/Figure1.pdf", Fig1, width = 8, height = 4, dpi = 300)
ggplot2::ggsave(file="fig/Figure1.png", Fig1, width = 8, height = 4, dpi = 300)


