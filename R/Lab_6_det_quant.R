library(lme4)
library(MASS)
library(DescTools)
library(cocor)
library(reshape2)
library(microbiome)
library(phyloseq)

source("R/Lab_5_filtering.R")

############################################
##############################################

#### sensitivity, specificity, positive predicted value and negative predictive value
#### for MA
sensit(Eim)

MA.asv1 <- prune_taxa(rownames(tax_table(Eim2))[1], Eim)
MA.asv2 <- prune_taxa(rownames(tax_table(Eim2))[2], Eim)
sensit(MA.asv1)
sensit(MA.asv2)

#### for SA
sensit(Eim2)

SA.asv1 <- prune_taxa(rownames(tax_table(Eim2))[1], Eim2)
SA.asv2 <- prune_taxa(rownames(tax_table(Eim2))[2], Eim2)
sensit(SA.asv1)
sensit(SA.asv2)

# now plotting and doing correlation analysis and comparisons
### for single amplicon
Plotting_cor_MA.l(ps=sin.PS18S.slv, f.sin18.slv, "FigureS1_single-amplicon", dir="fig/")
#adjusting p-value
# tss,rle,clr,acs, rare
p <- c(0.3243, 0.00001, 0.00001, 0.9665, 0.3759)
p.adjust(p, method="BH")

### for MA
Plotting_cor_MA.l(ps=all.PS.l.slv[[37]], ps.f=f.all.l.slv[[37]], "FigureS2_multi-amplicon", dir="fig/")
# tss,rle,clr,acs, rare
p <- c(0.0046, 0.00001, 0.7266, 0.0002, 0.0009)
round(p.adjust(p, method="BH"),3)

###################################### Plotting SA normalisation

library(microshades)
library("ggpmisc")
eimf <- subset_taxa(f.all.l.slv[[37]], Genus%in%"g__Eimeria")

#create total sums and Eimeria sums data frame
df <- data.frame(eimf@sam_data$labels)
df$eim <-(sample_sums(otu_table(eimf)))
names(df) <- c("labels", "FilEimeriaSums")

sam <- data.frame(sample_data(eimf))
all(sam$labels==df$labels)
df$Genome_copies_ngDNA <- sam$Genome_copies_ngDNA
df$dpi <- sam$dpi
#correlation tests
df$logFilEimeriaSums <- log(df$FilEimeriaSums)
df$logGC <- log(df$Genome_copies_ngDNA)

df <- df[df$FilEimeriaSums>0,] # only positive samples
df <- df[df$Genome_copies_ngDNA>0,]
a.cor <- cor.test(df$logGC, df$logFilEimeriaSums, method="pearson")

print(a.cor)
coul <- c(microshades_palette("micro_purple", 3, lightest=FALSE), microshades_palette("micro_orange", 3, lightest=FALSE), microshades_palette("micro_brown", 3, lightest=FALSE))

# plotting with abundance and prevalence filter correlation
a <- ggplot(df, aes(y=logGC, x=logFilEimeriaSums))+
    geom_jitter(shape=21, position=position_jitter(0.2), size=2.5, aes(fill= dpi), color= "white", alpha=0.8)+
    scale_fill_manual(values=coul)+
    ylab("Eimeria genome copies (log)")+
    xlab("Eimeria ASV abundance (log)")+
    ggtitle("No normalization: ASV read counts")+
    annotate(geom="text", x=min(df$logFilEimeriaSums), y=max(df$logGC), hjust=0.05, label=paste("Spearman rho=", round(a.cor$estimate, 2), ", p<0.001, ", "df= ", a.cor$parameter, sep=""))+
    theme_bw(base_size=10)+
    theme(axis.title.x = element_text(vjust = 0, size = 12),
          axis.title.y = element_text(vjust = 2, size = 12),
          plot.title=element_text(face="bold",
                                  margin=margin(10,0,10,0),
                                  size=12),
          legend.position= "none")

################################################################
#### using relative abundance
#create total sums and Eimeria sums data frame
df.b <-data.frame(sample_sums(otu_table(Eim)))
df.b$labels <- rownames(df.b)
names(df.b) <- c("TSS_Eim", "labels")
df <-merge(df, df.b, by="labels")
df$logTSS_Eim <-log(df$TSS_Eim)

b.cor <- cor.test(df$logGC, df$logTSS_Eim, method="pearson")

# plotting TSS correlation
b <- ggplot(df, aes(y=logGC, x=logTSS_Eim))+
    geom_jitter(shape=21, position=position_jitter(0.2), size=2.5, aes(fill= dpi), color= "white", alpha=0.8)+
    scale_fill_manual(values=coul)+
    ylab("Eimeria genome copies (log)")+
    xlab("Eimeria ASV abundance (log)")+
    ggtitle("Total sum scaling")+
    annotate(geom="text", x=min(df$logTSS_Eim), y=max(df$logGC), hjust=0.05, label=paste("Spearman rho=", round(b.cor$estimate, 2), ", p<0.001, ", "df= ", b.cor$parameter, sep=""))+
    theme_bw(base_size=10)+
    theme(axis.title.x = element_text(vjust = 0, size = 12),
          axis.title.y = element_text(vjust = 2, size = 12),
          plot.title=element_text(face="bold",
                                  margin=margin(10,0,10,0),
                                  size=12),
          legend.position= "none")
b

###############################################################
#### using TMM: Trimmed Mean by M-Values
library(microbial)

TMM.eim <- normalize(f.all.l.slv[[37]], method="TMM")
TMM.eim <- subset_taxa(TMM.eim, Genus%in%"g__Eimeria")

df.c <-data.frame(sample_sums(otu_table(TMM.eim)))
df.c$labels <- rownames(df.c)
names(df.c) <- c("TMM_Eim", "labels")

#merge
df <-merge(df, df.c, by="labels")

df$logTMM <- log(df$TMM_Eim)

#correlation tests
c.cor <- cor.test(df$logGC, df$logTMM, method="pearson")
c.cor

# plotting CLR correlation
c <- ggplot(df, aes(y=logGC, x=logTMM))+
    geom_jitter(shape=21, position=position_jitter(0.2), size=2.5, aes(fill= dpi), color= "white", alpha=0.8)+
    scale_fill_manual(values=coul)+
    ylab("Eimeria genome copies (log)")+
    xlab("Eimeria ASV abundance (log)")+
    ggtitle("Trimmed mean by M-values")+
    annotate(geom="text", x=min(df$logTMM), y=max(df$logGC), hjust=0.05, label=paste("Spearman rho=", round(c.cor$estimate, 2), ", p<0.001, ", "df= ", c.cor$parameter, sep=""))+
    theme_bw(base_size=10)+
    theme(axis.title.x = element_text(vjust = 0, size = 12),
          axis.title.y = element_text(vjust = 2, size = 12),
          plot.title=element_text(face="bold",
                                  margin=margin(10,0,10,0),
                                  size=12),
          legend.position= "none")
c

########### centered log ratio
Eimclr =microbiome::transform(f.all.l.slv[[37]], transform="clr")
Eimclr <- subset_taxa(Eimclr, Genus%in%"g__Eimeria")
df.d <-data.frame(sample_sums(otu_table(Eimclr)))
df.d$labels <- rownames(df.d)
names(df.d) <- c("clr_Eim", "labels")

#merge
df <-merge(df, df.d, by="labels")

#correlation tests
d.cor <- cor.test(df$logGC, df$clr_Eim, method="pearson")
d.cor

# plotting CLR correlation
d <- ggplot(df, aes(y=logGC, x=clr_Eim))+
    geom_jitter(shape=21, position=position_jitter(0.2), size=2.5, aes(fill= dpi), color= "white", alpha=0.8)+
    scale_fill_manual(values=coul)+
    ylab("Eimeria genome copies (log)")+
    xlab("Eimeria ASV abundance (log)")+
    ggtitle("Centered log ratio")+
    annotate(geom="text", x=min(df$clr_Eim), y=max(df$logGC), hjust=0.05, label=paste("Spearman rho=", round(d.cor$estimate, 2), ", p<0.001, ", "df= ", d.cor$parameter, sep=""))+
    theme_bw(base_size=10)+
    theme(axis.title.x = element_text(vjust = 0, size = 12),
          axis.title.y = element_text(vjust = 2, size = 12),
          plot.title=element_text(face="bold",
                                  margin=margin(10,0,10,0),
                                  size=12),
          legend.position= "none")
d

############### rarefaction
rare <- rarefy_even_depth(f.all.l.slv[[37]], rngseed=1234)
rare <- subset_taxa(rare, Genus%in%"g__Eimeria")

#create total sums and Eimeria sums data frame
df.e <- data.frame(sample_sums(otu_table(rare)))
df.e$labels <- rownames(df.e)
names(df.e) <- c("Eim_rare", "labels")

#merge
df <- merge(df,df.e, by="labels")
df.e <- df[df$Eim_rare>0,]
df.e$logEim_rare <- log(df.e$Eim_rare)

e.cor <- cor.test(df.e$logGC, df.e$logEim_rare, method="pearson")
e.cor
# plotting with rarefied abundance
e <- ggplot(df.e, aes(y=logGC, x=logEim_rare))+
    geom_jitter(shape=21, position=position_jitter(0.2), size=2.5, aes(fill= dpi), color= "white", alpha=0.8)+
    scale_fill_manual(values=coul)+
    ylab("Eimeria genome copies (log)")+
    xlab("Eimeria ASV abundance (log)")+
    ggtitle("Rarefied read counts")+
    annotate(geom="text", x=min(df.e$logEim_rare), y=max(df.e$logGC), hjust=0.05, label=paste("Spearman rho=", round(e.cor$estimate, 2), ", p<0.001, ", "df= ", e.cor$parameter, sep=""))+
    theme_bw(base_size=10)+
    theme(axis.title.x = element_text(vjust = 0, size = 12),
          axis.title.y = element_text(vjust = 2, size = 12),
          plot.title=element_text(face="bold",
                                  margin=margin(10,0,10,0),
                                  size=12),
          legend.position= "none")

e
    
# save plots of what we have so far
legend <- get_legend(a+
          guides(fill=guide_legend(nrow=3, byrow=TRUE))+
          theme(legend.position="top"))

fCor <-plot_grid(a, b, c, d, e,legend,
              align="vh",
              labels=c("a", "b", "c", "d", "e"),
              nrow=3)
fCor

ggplot2::ggsave("fig/FigureS2.pdf", fCor, width = 8, height = 12, dpi = 300)
ggplot2::ggsave("fig/FigureS2.png", fCor, width = 8, height = 12, dpi = 300)


# testing differences between correlations
cocor(~logGC + logFilEimeriaSums | logGC + logTSS_Eim, data = df, test = c("hittner2003", "zou2007"))

cocor(~logGC + logFilEimeriaSums | logGC + logTMM, data = df, test = c("hittner2003", "zou2007"))

cocor(~logGC + logFilEimeriaSums | logGC + clr_Eim, data = df, test = c("hittner2003", "zou2007"))

cocor(~logGC + logFilEimeriaSums | logGC + logEim_rare, data = df.e, test = c("hittner2003", "zou2007"))

rm(df)
rm(df.e)
rm(eimf)
#########################################################
##################################################
# plotting MA normalisations
eimf <- subset_taxa(f.sin18.slv, Genus%in%"g__Eimeria")

#create total sums and Eimeria sums data frame
df <- data.frame(eimf@sam_data$labels)
df$eim <-(sample_sums(otu_table(eimf)))
names(df) <- c("labels", "FilEimeriaSums")

sam <- data.frame(sample_data(eimf))
all(sam$labels==df$labels)
df$Genome_copies_ngDNA <- sam$Genome_copies_ngDNA
df$dpi <- sam$dpi
#correlation tests
df$logFilEimeriaSums <- log(df$FilEimeriaSums)
df$logGC <- log(df$Genome_copies_ngDNA)

df <- df[df$FilEimeriaSums>0,] # only positive samples
df <- df[df$Genome_copies_ngDNA>0,]
a.cor <- cor.test(df$logGC, df$logFilEimeriaSums, method="pearson")

a.cor
coul <- c(microshades_palette("micro_purple", 4, lightest=FALSE), microshades_palette("micro_orange", 3, lightest=FALSE), microshades_palette("micro_brown", 3, lightest=FALSE))

# plotting with abundance and prevalence filter correlation
a <- ggplot(df, aes(y=logGC, x=logFilEimeriaSums))+
    geom_jitter(shape=21, position=position_jitter(0.2), size=2.5, aes(fill= dpi), color= "white", alpha=0.8)+
    scale_fill_manual(values=coul)+
    ylab("Eimeria genome copies (log)")+
    xlab("Eimeria ASV abundance (log)")+
    ggtitle("No normalization: ASV read counts")+
    annotate(geom="text", x=min(df$logFilEimeriaSums), y=max(df$logGC), hjust=0.05, label=paste("Spearman rho=", round(a.cor$estimate, 2), ", p<0.001, ", "df= ", a.cor$parameter, sep=""))+
    theme_bw(base_size=10)+
    theme(axis.title.x = element_text(vjust = 0, size = 12),
          axis.title.y = element_text(vjust = 2, size = 12),
          plot.title=element_text(face="bold",
                                  margin=margin(10,0,10,0),
                                  size=12),
          legend.position= "none")

################################################################
#### using relative abundance
#create total sums and Eimeria sums data frame
df.b <-data.frame(sample_sums(otu_table(Eim2)))
df.b$labels <- rownames(df.b)
names(df.b) <- c("TSS_Eim", "labels")
df <-merge(df, df.b, by="labels")
df$logTSS_Eim <-log(df$TSS_Eim)
b.cor <- cor.test(df$logGC, df$logTSS_Eim, method="pearson")
b.cor

# plotting TSS correlation
b <- ggplot(df, aes(y=logGC, x=logTSS_Eim))+
    geom_jitter(shape=21, position=position_jitter(0.2), size=2.5, aes(fill= dpi), color= "white", alpha=0.8)+
    scale_fill_manual(values=coul)+
    ylab("Eimeria genome copies (log)")+
    xlab("Eimeria ASV abundance (log)")+
    ggtitle("Total sum scaling")+
    annotate(geom="text", x=min(df$logTSS_Eim), y=max(df$logGC), hjust=0.05, label=paste("Spearman rho=", round(b.cor$estimate, 2), ", p<0.001, ", "df= ", b.cor$parameter, sep=""))+
    theme_bw(base_size=10)+
    theme(axis.title.x = element_text(vjust = 0, size = 12),
          axis.title.y = element_text(vjust = 2, size = 12),
          plot.title=element_text(face="bold",
                                  margin=margin(10,0,10,0),
                                  size=12),
          legend.position= "none")
###############################################################
#### using TMM: Trimmed Mean by M-Values
TMM.eim <- normalize(f.sin18.slv, method="TMM")
TMM.eim <- subset_taxa(TMM.eim, Genus%in%"g__Eimeria")

df.c <-data.frame(sample_sums(otu_table(TMM.eim)))
df.c$labels <- rownames(df.c)
names(df.c) <- c("TMM_Eim", "labels")

#merge
df <-merge(df, df.c, by="labels")

df$logTMM <- log(df$TMM_Eim)

#correlation tests
c.cor <- cor.test(df$logGC, df$logTMM, method="pearson")
c.cor

# plotting CLR correlation
c <- ggplot(df, aes(y=logGC, x=logTMM))+
    geom_jitter(shape=21, position=position_jitter(0.2), size=2.5, aes(fill= dpi), color= "white", alpha=0.8)+
    scale_fill_manual(values=coul)+
    ylab("Eimeria genome copies (log)")+
    xlab("Eimeria ASV abundance (log)")+
    ggtitle("Trimmed mean by M-values")+
    annotate(geom="text", x=min(df$logTMM), y=max(df$logGC), hjust=0.05, label=paste("Spearman rho=", round(c.cor$estimate, 2), ", p<0.001, ", "df= ", c.cor$parameter, sep=""))+
    theme_bw(base_size=10)+
    theme(axis.title.x = element_text(vjust = 0, size = 12),
          axis.title.y = element_text(vjust = 2, size = 12),
          plot.title=element_text(face="bold",
                                  margin=margin(10,0,10,0),
                                  size=12),
          legend.position= "none")

########### centered log ratio
Eimclr =microbiome::transform(f.sin18.slv, transform="clr")
Eimclr <- subset_taxa(Eimclr, Genus%in%"g__Eimeria")
df.d <-data.frame(sample_sums(otu_table(Eimclr)))
df.d$labels <- rownames(df.d)
names(df.d) <- c("clr_Eim", "labels")

#merge
df <-merge(df, df.d, by="labels")

#correlation tests
d.cor <- cor.test(df$logGC, df$clr_Eim, method="pearson")
d.cor

# plotting CLR correlation
d <- ggplot(df, aes(y=logGC, x=clr_Eim))+
    geom_jitter(shape=21, position=position_jitter(0.2), size=2.5, aes(fill= dpi), color= "white", alpha=0.8)+
    scale_fill_manual(values=coul)+
    ylab("Eimeria genome copies (log)")+
    xlab("Eimeria ASV abundance (log)")+
    ggtitle("Centered log ratio")+
    annotate(geom="text", x=min(df$clr_Eim), y=max(df$logGC), hjust=0.05, label=paste("Spearman rho=", round(d.cor$estimate, 2), ", p<0.001, ", "df= ", d.cor$parameter, sep=""))+
    theme_bw(base_size=10)+
    theme(axis.title.x = element_text(vjust = 0, size = 12),
          axis.title.y = element_text(vjust = 2, size = 12),
          plot.title=element_text(face="bold",
                                  margin=margin(10,0,10,0),
                                  size=12),
          legend.position= "none")

############### rarefaction
rare <- rarefy_even_depth(f.sin18.slv, rngseed=1234)
rare <- subset_taxa(rare, Genus%in%"g__Eimeria")

#create total sums and Eimeria sums data frame
df.e <- data.frame(sample_sums(otu_table(rare)))
df.e$labels <- rownames(df.e)
names(df.e) <- c("Eim_rare", "labels")

#merge
df <- merge(df,df.e, by="labels")

df.e <- df[df$Eim_rare>0,]
df.e$logEim_rare <-log(df.e$Eim_rare)

names(df.e)

e.cor <- cor.test(df.e$logGC, df.e$logEim_rare, method="pearson")
e.cor

# plotting with rarefied abundance
e <- ggplot(df.e, aes(y=logGC, x=logEim_rare))+
    geom_jitter(shape=21, position=position_jitter(0.2), size=2.5, aes(fill= dpi), color= "white", alpha=0.8)+
    scale_fill_manual(values=coul)+
    ylab("Eimeria genome copies (log)")+
    xlab("Eimeria ASV abundance (log)")+
    ggtitle("Rarefied read counts")+
    annotate(geom="text", x=min(df.e$logEim_rare), y=max(df.e$logGC), hjust=0.05, label=paste("Spearman rho=", round(e.cor$estimate, 2), ", p<0.001, ", "df= ", e.cor$parameter, sep=""))+
    theme_bw(base_size=10)+
    theme(axis.title.x = element_text(vjust = 0, size = 12),
          axis.title.y = element_text(vjust = 2, size = 12),
          plot.title=element_text(face="bold",
                                  margin=margin(10,0,10,0),
                                  size=12),
          legend.position= "none")
e
    
# save plots of what we have so far
legend <- get_legend(a+
          guides(fill=guide_legend(nrow=3, byrow=TRUE))+
          theme(legend.position="top"))

fCor <-plot_grid(a, b, c, d, e,legend,
              align="vh",
              labels=c("a", "b", "c", "d", "e"),
              nrow=3)
fCor

ggplot2::ggsave("fig/FigureS1.pdf", fCor, width = 8, height = 12, dpi = 300)
ggplot2::ggsave("fig/FigureS1.png", fCor, width = 8, height = 12, dpi = 300)


# testing differences between correlations
cocor(~logGC + logFilEimeriaSums | logGC + logTSS_Eim, data = df, test = c("hittner2003", "zou2007"))

cocor(~logGC + logFilEimeriaSums | logGC + logTMM, data = df, test = c("hittner2003", "zou2007"))

cocor(~logGC + logFilEimeriaSums | logGC + clr_Eim, data = df, test = c("hittner2003", "zou2007"))

cocor(~logGC + logFilEimeriaSums | logGC + logEim_rare, data = df.e, test = c("hittner2003", "zou2007"))

head(df.e)

nrow(df.e)


e

