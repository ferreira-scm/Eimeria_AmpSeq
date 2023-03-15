######

fil <- function(ps){
    x = phyloseq::taxa_sums(ps)
    # abundance filtering at 0.005%
    keepTaxa = (x / sum(x) > 0.00005)
#    keepTaxa = (x / sum(x) > 0.0005)
    summary(keepTaxa)
    ps = phyloseq::prune_taxa(keepTaxa, ps)
# plus prevalnce filter at 1%
    KeepTaxap <- microbiome::prevalence(ps)>0.01
    ps <- phyloseq::prune_taxa(KeepTaxap, ps)
# subset samples based on total read count (500 reads)
#ps <- phyloseq::subset_samples(ps, phyloseq::sample_sums(ps) > 500)
    ps <- phyloseq::prune_samples(sample_sums(ps)>100, ps)
    ps
}

############# to erase
subPS <- function(ps) {
    ps <- transform_sample_counts(ps, function(x) x / sum(x))
    otu_table(ps) <- otu_table(ps)*sample_data(ps)$Total_DNA
    ps <- subset_taxa(ps, Genus%in%"g__Eimeria")
    ps <-aggregate_taxa(ps, level="Genus")
    ps <- psmelt(ps)
    ps <- ps[ps$Abundance>0,]
    ps <- ps[ps$Genome_copies_ngDNA>0,]
    ps$logA <- log(ps$Abundance)
    ps$logGC <- log(ps$Genome_copies_ngDNA)
    return(ps)
}

# to erase
p_tss <- function(df, lb, name){
ggplot(df, aes(x=logA, y=logGC))+
    geom_jitter(shape=21, position=position_jitter(0.002), size=4, aes(fill= dpi), color= "black", alpha=0.7)+
    scale_fill_brewer(palette="Spectral")+
    geom_smooth(method = "lm", se=FALSE, na.rm=TRUE) +
    stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")),
                                          parse = TRUE) +  
    ylab("Genome copies/ng DNA (log)")+
    xlab(paste(name, "Eimeria (log)", sep=" "))+
    ggtitle(name)+
    labs(tag=lb)+
    theme_bw()+
    theme(text = element_text(size=12),
#          axis.title.x = element_blank(),
          legend.position= "bottom")+
    guides(fill=guide_legend(nrow=2, byrow=TRUE))
}

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
    annotate(geom="text", x=min(df$logFilEimeriaSums), y=max(df$logGC), hjust=0.05, label=paste("Spearman rho=", round(a.cor$estimate, 2), ", p<0.001, ", "df= ", e.cor$parameter, sep=""))+
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
    annotate(geom="text", x=min(df$logTSS_Eim), y=max(df$logGC), hjust=0.05, label=paste("Spearman rho=", round(b.cor$estimate, 2), ", p<0.001, ", "df= ", e.cor$parameter, sep=""))+
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
    annotate(geom="text", x=min(df$logTMM), y=max(df$logGC), hjust=0.05, label=paste("Spearman rho=", round(c.cor$estimate, 2), ", p<0.001, ", "df= ", e.cor$parameter, sep=""))+
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
    annotate(geom="text", x=min(df$clr_Eim), y=max(df$logGC), hjust=0.05, label=paste("Spearman rho=", round(d.cor$estimate, 2), ", p<0.001, ", "df= ", e.cor$parameter, sep=""))+
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
e.cor <- cor.test(df.e$logGC, df.e$logEim_rare], method="pearson")
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

ggplot2::ggsave("fig/FigureS1.pdf", fCor, width = 6, height = 10, dpi = 300)

#ggplot2::ggsave(file=paste(dir, name, "TSS-CS.png", sep=""), f, width = 5, height = 5, dpi = 600)

### terrible coding here now...
names(sdt)

cor.df <- sdt[,c("logGC", "logEimeriaSums", "logFilEimeriaSums", "logTSS_Eim", "logACS_Eim", "clr_Eim", "REL_Eim", "logEim_rare")]

cor.df$REL_Eim <- cor.df$REL_Eim*-1

cor.df1 <- na.omit(cor.df[!is.infinite(cor.df$logEimeriaSums),])

#cor.df1 <- na.omit(cor.df1[cor.df1$logFilEimeriaSums>0,])
cor.df2 <- na.omit(cor.df[!is.infinite(cor.df$logTSS_Eim),])
cor.df2 <- cor.df2[!is.infinite(cor.df2$logGC),]
print(cocor(~logGC + logFilEimeriaSums | logGC + logTSS_Eim, data = cor.df2, test = c("hittner2003", "zou2007")))

cor.df3 <- na.omit(cor.df[!is.infinite(cor.df$logFilEimeriaSums),])
cor.df3 <- cor.df3[!is.infinite(cor.df3$logGC),]
print(cocor(~logGC + logFilEimeriaSums | logGC + REL_Eim, data = cor.df3,test = c("hittner2003", "zou2007")))

cor.df4 <- na.omit(cor.df[!is.infinite(cor.df$logTSS_Eim),])
cor.df4 <- cor.df4[!is.infinite(cor.df4$logGC),]
print(cocor(~logGC + logFilEimeriaSums | logGC + clr_Eim, data = cor.df4,
            test = c("hittner2003", "zou2007")))

cor.df5 <- na.omit(cor.df[!is.infinite(cor.df$logACS_Eim),])
cor.df5 <- cor.df5[!is.infinite(cor.df5$logGC),]
print(cocor(~logGC + logFilEimeriaSums | logGC + logACS_Eim, data = cor.df5,
            test = c("hittner2003", "zou2007")))

cor.df6 <- na.omit(cor.df[!is.infinite(cor.df$logEim_rare),])
cor.df6 <- cor.df6[!is.infinite(cor.df6$logGC),]
print(cocor(~logGC + logFilEimeriaSums | logGC + logEim_rare, data = cor.df6,
            test = c("hittner2003", "zou2007")))

#png(filename=paste(dir, name, "COR.png", sep=""),
#    width =14, height = 14, units = "in", res= 300)
#fCor
#dev.off()

}


