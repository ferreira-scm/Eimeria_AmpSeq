library(lme4)
library(MASS)
library(DescTools)
library(cocor)
library(reshape2)
library(microbiome)
library(phyloseq)
library(cowplot)
library(microshades)
library("ggpmisc")
library(microbial)

source("R/1_Lab_filter.R")

############################################
##############################################

#### sensitivity, specificity, positive predicted value and negative predictive value
#### for MA
sensit <- function(Eim_nf){
# GC- Eim -
    tneg <- summary(sample_sums(Eim_nf@otu_table)==0&Eim_nf@sam_data$Genome_copies_gFaeces==0)[[3]]

# GC+ Eim -
    fneg <- summary(sample_sums(Eim_nf@otu_table)==0&Eim_nf@sam_data$Genome_copies_gFaeces>0)[[3]]

# GC+ Eim+
    tpos <- summary(sample_sums(Eim_nf@otu_table)>0&Eim_nf@sam_data$Genome_copies_gFaeces>0)[[3]]

# GC- Eim+
    fpos <- try(summary(sample_sums(Eim_nf@otu_table)>0&Eim_nf@sam_data$Genome_copies_gFaeces==0)[[3]])

    if (class(fpos)=="try-error"){
        fpos <- "0"
    }
print(paste("false positives:", fpos,sep=" "))
print(paste("false negatives:", fneg,sep=" "))


# GC+ GC-
gcpos <- summary(Eim_nf@sam_data$Genome_copies_gFaeces==0)[[2]]
gcneg <- summary(Eim_nf@sam_data$Genome_copies_gFaeces==0)[[3]]

# Eim + Eim -
epos <- summary(sample_sums(Eim_nf@otu_table)==0)[[2]]
eneg <- summary(sample_sums(Eim_nf@otu_table)==0)[[3]]

e.s <- matrix(as.numeric(c(tneg, fpos, fneg, tpos)), ncol=2, byrow=TRUE)
colnames(e.s) <- c("Eim-","Eim+")
rownames(e.s) <- c("GC-", "GC+")
margin1 <- margin.table(e.s, margin=1)
margin2 <- margin.table(e.s, margin=2)

#sensitivity
print("Sensitivity:")
print(e.s[2,2]/margin1[2]*100)
#specificity
print("Specificity:")
print(e.s[1,1]/margin1[1]*100)
#ppv
print("positive predictive value")
print(e.s[2,2]/margin2[2]*100)
#npv
print("negative predictive value")
print(e.s[1,1]/margin2[1]*100)
}


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
###################################### Plotting SA normalisation

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
coul <- c(microshades_palette("micro_purple", 3, lightest=FALSE), microshades_palette("micro_orange", 4, lightest=FALSE), microshades_palette("micro_brown", 3, lightest=FALSE))

coul2 <- coul[-1]

# plotting with abundance and prevalence filter correlation
a <- ggplot(df, aes(y=logGC, x=logFilEimeriaSums))+
    geom_jitter(shape=21, position=position_jitter(0.2), size=2.5, aes(fill= dpi), color= "white", alpha=0.8)+
    scale_fill_manual(values=coul2)+
    ylab("Eimeria genome copies (log)")+
    xlab("Eimeria ASV abundance (log)")+
    ggtitle("No normalization: ASV read counts")+
        labs(fill="DPI")+
    annotate(geom="text", x=min(df$logFilEimeriaSums), y=max(df$logGC), hjust=0.05, label=paste("Pearson's rho=", round(a.cor$estimate, 2), ", p<0.001, ", "df= ", a.cor$parameter, sep=""))+
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
    scale_fill_manual(values=coul2)+
    ylab("Eimeria genome copies (log)")+
    xlab("Eimeria ASV abundance (log)")+
        labs(fill="DPI")+
    ggtitle("Total sum scaling")+
    annotate(geom="text", x=min(df$logTSS_Eim), y=max(df$logGC), hjust=0.05, label=paste("Pearson's rho=", round(b.cor$estimate, 2), ", p<0.001, ", "df= ", b.cor$parameter, sep=""))+
    theme_bw(base_size=10)+
    theme(axis.title.x = element_text(vjust = 0, size = 12),
          axis.title.y = element_text(vjust = 2, size = 12),
          plot.title=element_text(face="bold",
                                  margin=margin(10,0,10,0),
                                  size=12),
          legend.position= "none")

###############################################################
#### using TMM: Trimmed Mean by M-Values

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
    scale_fill_manual(values=coul2)+
    ylab("Eimeria genome copies (log)")+
    xlab("Eimeria ASV abundance (log)")+
        labs(fill="DPI")+
    ggtitle("Trimmed mean by M-values")+
    annotate(geom="text", x=min(df$logTMM), y=max(df$logGC), hjust=0.05, label=paste("Pearson's rho=", round(c.cor$estimate, 2), ", p<0.001, ", "df= ", c.cor$parameter, sep=""))+
    theme_bw(base_size=10)+
    theme(axis.title.x = element_text(vjust = 0, size = 12),
          axis.title.y = element_text(vjust = 2, size = 12),
          plot.title=element_text(face="bold",
                                  margin=margin(10,0,10,0),
                                  size=12),
          legend.position= "none")

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
    scale_fill_manual(values=coul2)+
    ylab("Eimeria genome copies (log)")+
    xlab("Eimeria ASV abundance (log)")+
        labs(fill="DPI")+
    ggtitle("Centered log ratio")+
    annotate(geom="text", x=min(df$clr_Eim), y=max(df$logGC), hjust=0.05, label=paste("Pearson's rho=", round(d.cor$estimate, 2), ", p<0.001, ", "df= ", d.cor$parameter, sep=""))+
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
    scale_fill_manual(values=coul2)+
    ylab("Eimeria genome copies (log)")+
    xlab("Eimeria ASV abundance (log)")+
        labs(fill="DPI")+
    ggtitle("Rarefied read counts")+
    annotate(geom="text", x=min(df.e$logEim_rare), y=max(df.e$logGC), hjust=0.05, label=paste("Pearson's rho=", round(e.cor$estimate, 2), ", p<0.001, ", "df= ", e.cor$parameter, sep=""))+
    theme_bw(base_size=10)+
    theme(axis.title.x = element_text(vjust = 0, size = 12),
          axis.title.y = element_text(vjust = 2, size = 12),
          plot.title=element_text(face="bold",
                                  margin=margin(10,0,10,0),
                                  size=12),
          legend.position= "none")

# save plots of what we have so far
legend <- get_legend(a+
          guides(fill=guide_legend(nrow=3, byrow=TRUE))+
          theme(legend.position="top"))

fCor <-cowplot::plot_grid(a, b, c, d, e,legend,
              align="vh",
              labels=c("a", "b", "c", "d", "e"),
              nrow=3)

fCor

ggplot2::ggsave("fig/FigureS2.pdf", fCor, width = 170, height = 280, units="mm", dpi = 300)


# testing differences between correlations
cocor(~logGC + logFilEimeriaSums | logGC + logTSS_Eim, data = df, test = c("hittner2003", "zou2007"))

cocor(~logGC + logFilEimeriaSums | logGC + logTMM, data = df, test = c("hittner2003", "zou2007"))

cocor(~logGC + logFilEimeriaSums | logGC + clr_Eim, data = df, test = c("hittner2003", "zou2007"))

cocor(~logGC + logFilEimeriaSums | logGC + logEim_rare, data = df.e, test = c("hittner2003", "zou2007"))

p <- c(0.3243, 0.104, 0.00001, 0.3759)
p.adjust(p, method="BH")


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

# plotting with abundance and prevalence filter correlation
a <- ggplot(df, aes(y=logGC, x=logFilEimeriaSums))+
    geom_jitter(shape=21, position=position_jitter(0.2), size=2.5, aes(fill= dpi), color= "white", alpha=0.8)+
    scale_fill_manual(values=coul)+
    ylab("Eimeria genome copies (log)")+
    xlab("Eimeria ASV abundance (log)")+
    labs(fill="DPI")+
    ggtitle("No normalization: ASV read counts")+
    annotate(geom="text", x=min(df$logFilEimeriaSums), y=max(df$logGC), hjust=0.05, label=paste("Pearson's rho=", round(a.cor$estimate, 2), ", p<0.001, ", "df= ", a.cor$parameter, sep=""))+
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
    labs(fill="DPI")+
    ggtitle("Total sum scaling")+
    annotate(geom="text", x=min(df$logTSS_Eim), y=max(df$logGC), hjust=0.05, label=paste("Pearson's rho=", round(b.cor$estimate, 2), ", p<0.001, ", "df= ", b.cor$parameter, sep=""))+
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
    labs(fill="DPI")+
    ggtitle("Trimmed mean by M-values")+
    annotate(geom="text", x=min(df$logTMM), y=max(df$logGC), hjust=0.05, label=paste("Pearson's rho=", round(c.cor$estimate, 2), ", p<0.001, ", "df= ", c.cor$parameter, sep=""))+
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
    labs(fill="DPI")+
    ggtitle("Centered log ratio")+
    annotate(geom="text", x=min(df$clr_Eim), y=max(df$logGC), hjust=0.05, label=paste("Pearson's rho=", round(d.cor$estimate, 2), ", p<0.001, ", "df= ", d.cor$parameter, sep=""))+
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
    labs(fill="DPI")+
    ggtitle("Rarefied read counts")+
    annotate(geom="text", x=min(df.e$logEim_rare), y=max(df.e$logGC), hjust=0.05, label=paste("Pearson's rho=", round(e.cor$estimate, 2), ", p<0.001, ", "df= ", e.cor$parameter, sep=""))+
    theme_bw(base_size=10)+
    theme(axis.title.x = element_text(vjust = 0, size = 12),
          axis.title.y = element_text(vjust = 2, size = 12),
          plot.title=element_text(face="bold",
                                  margin=margin(10,0,10,0),
                                  size=12),
          legend.position= "none")

# save plots of what we have so far
legend <- get_legend(a+
          guides(fill=guide_legend(nrow=3, byrow=TRUE))+
          theme(legend.position="top"))

fCor <-cowplot::plot_grid(a, b, c, d, e,legend,
              align="vh",
              labels=c("a", "b", "c", "d", "e"),
              nrow=3)
fCor

ggplot2::ggsave("fig/FigureS1.pdf", fCor, width =170, height = 280, units="mm", dpi = 300)


# testing differences between correlations
cocor(~logGC + logFilEimeriaSums | logGC + logTSS_Eim, data = df, test = c("hittner2003", "zou2007"))

cocor(~logGC + logFilEimeriaSums | logGC + logTMM, data = df, test = c("hittner2003", "zou2007"))

cocor(~logGC + logFilEimeriaSums | logGC + clr_Eim, data = df, test = c("hittner2003", "zou2007"))

cocor(~logGC + logFilEimeriaSums | logGC + logEim_rare, data = df.e, test = c("hittner2003", "zou2007"))

p <- c(0.0046, 0.9445, 0.7266, 0.0009)
p.adjust(p, method="BH")

