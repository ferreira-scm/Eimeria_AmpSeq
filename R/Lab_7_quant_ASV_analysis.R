#!/usr/bin/Rscript

#library(reshape)
library(RColorBrewer)
library("seqinr")
library(ggtree)
library(magrittr)
library(BiocManager)
library(lmerTest)
library(lme4)

source("R/Lab_5_filtering.R")

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

SA.e4 <- SA.e[which(SA.e$OTU==colnames(Eim2@otu_table)[4]),]
SA.e3 <- SA.e[which(SA.e$OTU==colnames(Eim2@otu_table)[3]),]
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

### they are correlate pretty well, ASV1 being the best
cor.test(log(SA.e10$Genome_copies_ngDNA), log(SA.e10$Abundance))
cor.test(log(SA.e20$Genome_copies_ngDNA), log(SA.e20$Abundance))
cor.test(log(SA.e30$Genome_copies_ngDNA), log(SA.e30$Abundance))
cor.test(log(SA.e40$Genome_copies_ngDNA), log(SA.e40$Abundance))

### now same for MA
MA.e10 <- MA.e1[MA.e1$Genome_copies_ngDNA>0,]
MA.e10 <- MA.e10[!is.na(MA.e10$Genome_copies_ngDNA),]
MA.e10 <- MA.e10[MA.e10$Abundance>0,]

MA.e20 <- MA.e2[MA.e2$Genome_copies_ngDNA>0,]
MA.e20 <- MA.e20[!is.na(MA.e10$Genome_copies_ngDNA),]
MA.e20 <- MA.e20[MA.e20$Abundance>0,]

## ASV1 is also best for MA
cor.test(log(MA.e10$Genome_copies_ngDNA), log(MA.e10$Abundance))
cor.test(log(MA.e20$Genome_copies_ngDNA), log(MA.e20$Abundance))

####
## Let's do the correlations first with OPG
SA.opg <- SA.e.g[SA.e.g$OPG>0,]
SA.opg <- SA.opg[SA.opg$Abundance>0,]

### OPG correlation
cor.test(log(SA.opg$OPG), log(SA.opg$Abundance))
col7 <- c("#F1B6DA", "#C51B7D", "#DE77AE","#01665E", "#35978F", "#80CDC1","#C7EAE5")
OPG_Abundance <- ggplot(SA.opg, aes(y=log(OPG), x=log(Abundance), fill=dpi))+
    geom_point(shape=21, size=4, alpha=0.8)+
    scale_fill_manual(values=col7)+
    ylab("Oocysts/g faeces (log)")+
#    guides(fill=guide_legend(ncol=3, byrow=TRUE))+
    xlab("Eimeria ASV abundance (log)")+
    annotate(geom="text", x=min(log(SA.opg$Abundance)), y=max(log(SA.opg$OPG)), hjust=0.05, label="Pearson rho=0.67, p<0.001", size=3)+
    theme_bw(base_size=10)+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.position = "none",
          axis.line = element_line(colour = "black"))

## now same but for MA
MA.opg <- MA.e.g[MA.e.g$OPG>0,]
MA.opg <- MA.opg[MA.opg$Abundance>0,]
### OPG correlation
cor.test(log(MA.opg$OPG), log(MA.opg$Abundance))

OPG_Abundance_MA <- ggplot(MA.opg, aes(y=log(OPG), x=log(Abundance), fill=dpi))+
    geom_point(shape=21, size=4, alpha=0.8)+
    scale_fill_manual(values=col7)+
    ylab("Oocysts/g faeces (log)")+
#    guides(fill=guide_legend(ncol=3, byrow=TRUE))+
    xlab("Eimeria ASV abundance (log)")+
    annotate(geom="text", x=min(log(MA.opg$Abundance)), y=max(log(MA.opg$OPG)), hjust=0.05, label="Pearson rho=0.59, p<0.001", size=3)+
    theme_bw(base_size=10)+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
           legend.position = "none",
          axis.line = element_line(colour = "black"))

OPG_Abundance_MA

legend <- get_legend(OPG_Abundance_MA+
          guides(fill=guide_legend(nrow=1, byrow=TRUE))+
          theme(legend.position="top"))

OPG_ab <- plot_grid(OPG_Abundance, OPG_Abundance_MA, labels="auto")

OPG_ab_panel <- plot_grid(legend, OPG_ab,  nrow=2, rel_heights=c(0.1, 0.8))


# saving object so that we can plot them together with wild
#saveRDS(OPG_Abundance, "tmp/OPG_Abundance.R")
# saving object so that we can plot them together with wild
#saveRDS(OPG_Abundance_MA, "tmp/OPG_Abundance_MA.R")
#saveRDS(OPG_ab_panel, "tmp/OPG_Abundance_MA_panel.R")


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
plot_SA_all <- ggplot(SA.e.g0, aes(x=log(Abundance), y=log(Genome_copies_ngDNA)))+
    geom_point(shape=21, size=2.5, aes(fill= dpi), color= "white", alpha=0.8)+
    scale_fill_manual(values=coul)+
    ylab("Eimeria genome copies (log)")+
    xlab("Eimeria ASV abundance (log)")+
    ggtitle("Microfluidics PCR Eimeria quantification")+
#    annotate(geom="text", x=min(log(SA.e.g0$Abundance)), y=max(log(SA.e.g0$Genome_copies_ngDNA)), hjust=0.05, label="Conditional R2=0.89")+
    theme_bw(base_size=10)+
    theme(axis.title.x = element_text(vjust = 0, size = 12),
          axis.title.y = element_text(vjust = 2, size = 12),
          plot.title=element_text(face="bold",
                                  margin=margin(10,0,10,0),
                                  size=12),
          legend.position="none")

#          legend.position= "none")
plot_MA_all <-ggplot(MA.e.g0, aes(x=log(Abundance), y=log(Genome_copies_ngDNA)))+
        geom_point(shape=21, size=2.5, aes(fill= dpi), color= "white", alpha=0.8)+
    scale_fill_manual(values=coul2)+
    ylab("Eimeria genome copies (log)")+
    xlab("Eimeria ASV abundance (log)")+
    ggtitle("Microfluidics PCR Eimeria quantification")+
#    annotate(geom="text", x=min(log(MA.e.g0$Abundance)), y=max(log(MA.e.g0$Genome_copies_ngDNA)), hjust=0.05, label="Conditional R2=0.79")+
    theme_bw(base_size=10)+
    theme(axis.title.x = element_text(vjust = 0, size = 12),
          axis.title.y = element_text(vjust = 2, size = 12),
          plot.title=element_text(face="bold",
                                  margin=margin(10,0,10,0),
                                  size=12),
    legend.position= "none")

plot_MA_all

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
library(MuMIn)
r.squaredGLMM(lmm.sa) ## well, that was easy

#plot(lmm.sa)
anova(lmm.sa, test="LRT")
anova(lmm.sa, lmm.sa0)

# saving model results
sink("fig/Eimeira_asv_SA.txt")
summary(lmm.sa)
sink()

## plotting estimated random effects for each dpi and its interval estimate
library(merTools)
plotREsim(REsim(lmm.sa))
confint(lmm.sa)

library(sjPlot)
sjPlot::plot_model(lmm.sa)

## is using negative binomial better than transforming?
gau.m.0.nb <- glmer.nb(data=SA.df, Genome_copies~ASV1+ASV2+ASV3+ASV4+(1|dpi))
summary(gau.m.0.nb)


#### for multiamplicon
MA.df <- data.frame(Eim.0@sam_data$Genome_copies_ngDNA, Eim.0@otu_table[,1], Eim.0@otu_table[,2],Eim.0@sam_data$dpi)
names(MA.df) <- c("Genome_copies", "ASV1", "ASV2", "dpi")
MA.df$TotalE <- MA.df$ASV1+MA.df$ASV2
MA.df <- MA.df[MA.df$TotalE>0,]
nrow(MA.df)

lmm.ma <- lmer(data=MA.df, log(Genome_copies)~log(1+ASV1)+log(1+ASV2)+(1|dpi))
summary(lmm.ma)
lmm.ma0 <- lmer(data=MA.df, log(Genome_copies)~1+(1|dpi))

#sink("fig/Eimeira_asv_MA.txt")
#print(summary(lmm.ma))
#sink()

# LR tests for significance
anova(lmm.ma, test="LRT")

anova(lmm.ma, lmm.ma0)

ranova(lmm.ma)

## marginal and conditional r2
library(MuMIn)
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

MA2

FigureS3 <-cowplot::plot_grid(SA1, SA2, MA1, MA2,
              align="vh",
              labels="auto",
              nrow=2)

ggplot2::ggsave("fig/FigureS3.pdf", FigureS3, width = 10, height = 8, dpi = 300)
#ggplot2::ggsave("fig/FigureS3.png", fCor, width = 8, height = 12, dpi = 300)


MA.all <- ggplot(MA.e.g, aes(x=dpi, y=log(1+Abundance), fill=EH_ID))+
    geom_point(shape=21, position=position_jitter(0.1), size=2, alpha=0.7)+
    geom_line(aes(group=EH_ID), alpha=0.2)+
#    scale_fill_manual(values=coul)+
    ggtitle("Multi-amplicon")+
    labs(y="Eimeria ASV abundance( log(1+)", x="Days post infection")+
    theme_bw(base_size=10)+
    theme(plot.title=element_text(hjust=0.5, face="bold"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          text=element_text(size=10),
          legend.position = "none",
          axis.line = element_line(colour = "black"))

SA.all <- ggplot(SA.e.g, aes(x=dpi, y=log(1+Abundance), fill=EH_ID))+
    geom_point(shape=21, position=position_jitter(0.1), size=2, alpha=0.7)+
    geom_line(aes(group=EH_ID), alpha=0.2)+
#    scale_fill_manual(values=coul)+
    ggtitle("Single amplicon")+
    labs(y="Eimeria ASV abundance log(1+)", x="Days post infection")+
    theme_bw(base_size=10)+
    theme(plot.title=element_text(hjust=0.5, face="bold"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          text=element_text(size=10),
          legend.position = "none",
          axis.line = element_line(colour = "black"))

MA.all

SA.all


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

library(ggpmisc)

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

ASV2.c

ASV1.c

Fig1 <- plot_grid(ASV1.c, ASV2.c, labels="auto", nrow=1, align="hv")

Fig1

ggplot2::ggsave(file="fig/Figure1.pdf", Fig1, width = 8, height = 4, dpi = 300)
ggplot2::ggsave(file="fig/Figure1.png", Fig1, width = 8, height = 4, dpi = 300)

sa <- SA.e[, c("ASV", "Abundance", "EH_ID", "dpi", "labels")]
sa$amp <- "sa"
ma <- MA.e[, c("ASV", "Abundance", "EH_ID", "dpi", "labels")]
ma$amp <- "ma"
sama <- rbind(sa, ma)

ASV_sama <- ggplot(sama, aes(x=(Abundance), y=ASV, fill=amp))+
    geom_point(size=4, shape=21, alpha=0.7, position=position_jitterdodge(dodge.width=0.75, jitter.width=0.1))+
    geom_boxplot(alpha=0.3, colour="black", outlier.shape = NA)+
    scale_fill_manual(values=c("#CC6677", "#DDCC77"), name="", labels=c("multi-amplicon", "single-amplicon"))+
    xlab("Eimeria ASV abundance")+
    ylab("")+
    guides(fill = guide_legend(override.aes = list(linetype = 0)),
           color = guide_legend(override.aes = list(linetype = 0)))+
    theme_bw(base_size=10)+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          text=element_text(size=10),
          legend.position = "top",
          axis.line = element_line(colour = "black"))

ASV_sama



ggplot2::ggsave(file="fig/Figure1-ASV_concordance.pdf", Fig1, width=8, height=4, dpi=300)
ggplot2::ggsave(file="fig/Figure1-ASV_concordance.png", Fig1, width=8, height=4, dpi=300)


sama_1 <- ma.sa[ma.sa$ASV=="ASV1",]

sama_2 <- ma.sa[ma.sa$ASV=="ASV2",]

t.test(log(1+sama_1$Abundance.x), log(1+sama_1$Abundance.y))
wilcox.test(sama_1$Abundance.x, sama_1$Abundance.y)

t.test(log(1+sama_2$Abundance.x), log(1+sama_2$Abundance.y))
wilcox.test(sama_2$Abundance.x, sama_2$Abundance.y)

ma.sa[ma.sa$ASV=="ASV1",]

SA_Eimeria.ASVs <- ggplot(SA.e.0, aes(y=log(Genome_copies_ngDNA), x=log(Abundance), fill=ASV))+
    geom_point(shape=21, size=4, alpha=0.7)+
    scale_fill_manual(values=c("#009E73", "#F0E442", "#0072B2",
                               "#D55E00", "#E63C9A"))+
    geom_smooth(method=lm, colour="black", aes(colour="ASV"))+
    ylab("Genome copies/ng DNA (log)")+
#    guides(fill=guide_legend(ncol=3, byrow=TRUE))+
    xlab("ASV abundance (log)")+
#    annotate(geom="text", x=12, y=19, label="Pearson rho=0.92, p<0.001", size=3)+
    theme_bw()+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          text=element_text(size=12),
          legend.position = "top",
          axis.line = element_line(colour = "black"))

MA_Eimeria.ASVs <- ggplot(MA.e.0, aes(x=log(Genome_copies_ngDNA), y=log(Abundance), fill=ASV))+
    geom_point(shape=21, size=4, alpha=0.7)+
    scale_fill_manual(values=c("#009E73", "#F0E442", "#0072B2"))+
    geom_smooth(method=lm, colour="black", aes(colour="ASV"))+
    ylab("Genome copies/ng DNA (log)")+
    xlab("ASV abundance (log)")+
#    annotate(geom="text", x=12, y=19, label="Pearson rho=0.92, p<0.001", size=3)+
    theme_bw()+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          text=element_text(size=11),
          legend.position = "top",
          axis.line = element_line(colour = "black"))

ASV.SA.MA <- plot_grid(SA_Eimeria.ASVs, MA_Eimeria.ASVs, nrow=1, labels="auto")

plot_grid(SA1,SA2,SA3,SA4) -> SA.asv
plot_grid(MA1,MA2, nrow=2) -> MA.asv

row1 <- plot_grid(SA1,MA1, nrow=1, labels=c("b", "c"))
row2 <- plot_grid(SA2,MA2, nrow=1, labels=c("d", "e"))

SAMA.asv <- plot_grid(qpcr, row1, row2, labels=c("a", "", ""), ncol=1)

SAMA.all <- plot_grid(qpcr, SA.all, MA.all, labels="auto", nrow=1)

SAMA.asv

SAMA.all

#ggplot2::ggsave(file="fig/Eimeria_ASVs_dpi.pdf", SAMA.asv, width = 10, height = 15, dpi = 300)
#ggplot2::ggsave(file="fig/Eimeria_ASVs_dpi.png", SAMA.asv, width = 10, height = 15, dpi = 300)

#ggplot2::ggsave(file="fig/Sup1_Eimeria_SA_MA_all_dpi.pdf", SAMA.all, width = 10, height = 4, dpi = 300)

#ggplot2::ggsave(file="fig/Sup1_Eimeria_SA_MA_all_dpi.png", SAMA.all, width = 10, height = 4, dpi = 300)

#ggplot2::ggsave(file="fig/MA_SA_Eimeria_ASVs.pdf", ASV.SA.MA, width = 13, height = 5, dpi = 300)
#ggplot2::ggsave(file="fig/MA_SA_Eimeria_ASVs.png", ASV.SA.MA, width = 13, height = 5, dpi = 300)


#ggplot2::ggsave(file="fig/ASV_abundace_SAMA.pdf", ASV_sama, width=4, height=4, dpi=300)
#ggplot2::ggsave(file="fig/ASV_abundace_SAMA.png", ASV_sama, width=4, height=4, dpi=300)



#################### plotting individuals by ASV
library(cowplot) # to plot a list of plots

cl <- colorRampPalette(brewer.pal(8, "Accent"))(6)

length(levels(as.factor(SA.e$EH_ID)))

SA.e$EH_ID

p.ID <- function(i){
    SA.e%>%
    dplyr::filter(EH_ID%in%SA.e$EH_ID[i])%>%
    dplyr::select(EH_ID, dpi, ASV, Genome_copies_ngDNA, Abundance)%>%
    ggplot(aes(x=dpi, y=log(Abundance+1), fill=ASV))+
    geom_point(shape=21, position=position_jitter(0.2), size=4, alpha=0.8, aes(fill=ASV), color="black")+
    geom_line(aes(group=ASV), color="gray", alpha=0.5)+
    scale_fill_manual(values=c("#009E73", "mediumvioletred","#F0E442","#0072B2"), name="")+
    ylab("Single-amplicon ASV abundance /ng DNA (log(1+)")+
    theme_bw()+
    theme(text = element_text(size=16), axis.title.x = element_blank(), legend.position = "top") -> nm
    nm
}

library(dplyr)
p.EH <- list()
for (i in 1:22){
    p <- p.ID(i)
    p.EH[[i]] <- p
}

#pdf("fig/Eimeria_ASVs_ID.pdf")
#for (i in 1:22) {
#    print(p.EH[[i]])
#}
#dev.off()

p.ID.MA <- function(i){
    MA.e%>%
    dplyr::filter(EH_ID%in%MA.e$EH_ID[i])%>%
    dplyr::select(EH_ID, dpi, ASV, Genome_copies_ngDNA, Abundance)%>%
    ggplot(aes(x=dpi, Abundance+1, fill=ASV))+
    geom_point(shape=21, position=position_jitter(0.2), size=4, alpha=0.7, aes(fill=ASV), color="black")+
    ylab("Multi-amplicon ASV abundance /ng DNA (log(1+)")+
    geom_line(aes(group=ASV), color="gray", alpha=0.5)+
    scale_fill_manual(values=c("#009E73", "mediumvioletred"), name="")+
    theme_bw()+
    theme(text = element_text(size=16), axis.title.x = element_blank(), legend.position = "top") -> nm
    nm
}

p.EH.M <- list()
for (i in 1:22){
    p <- p.ID.MA(i)
    p.EH.M[[i]] <- p
}

p.EH.M[1]

#pdf("fig/MA/Eimeria_ASVs_ID.pdf")
#for (i in 1:22) {
#    print(p.EH.M[[i]])
#}
#dev.off()

################################
ggplot(SA.e, aes(x=Abundance, y=ASV))+
    geom_point(aes(fill=ASV), size=4, shape=21, apha=0.5, position=position_jitter(height=0.2))+
#      geom_density(aes(color = ASV), size = 1) +
    scale_fill_manual(values=c("#009E73", "#F0E442", "#0072B2", "#D55E00", "#E63C9A"))+
#        scale_colour_manual(values=c("#009E73", "#F0E442", "#0072B2", "#D55E00", "#E63C9A"))+
#    ylab("SA - ASV abundance (log1+)")+
#    xlab("MA - ASV abundance (log1+)")+
    theme_bw()+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          text=element_text(size=12),
           legend.position = "top",
          axis.line = element_line(colour = "black"))


ASV1 <- ma.sa[ma.sa$ASV=="ASV1",]
ASV2 <- ma.sa[ma.sa$ASV=="ASV2",]
ASV3 <- ma.sa[ma.sa$ASV=="ASV3",]

cor.test(ASV1$Abundance.x, ASV1$Abundance.y)
cor.test(ASV2$Abundance.x, ASV2$Abundance.y)
cor.test(log(1+ASV1$Abundance.x), log(1+ASV1$Abundance.y))
cor.test(log(1+ASV2$Abundance.x), log(1+ASV2$Abundance.y))



