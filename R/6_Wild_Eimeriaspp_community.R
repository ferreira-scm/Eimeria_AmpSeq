y# Eimeria species
source("R/5_Wild_Eimeria_spp_assignment.R")

library(cowplot)
library(ggeffects)

######### Preparing phyloseq objects for plotting and so on
Eim.TSSw@tax_table[,6] <- "Eimeria"
Eim_sp <- tax_glom(Eim.TSSw, "Species")

amp_names <- gsub("_ASV.*", "", names18S)
amp_names <- gsub("_R", "", amp_names)

Eim.amp <- Eim.TSSw
Eim.amp@tax_table[,6] <- amp_names

eim.m0 <- psmelt(Eim_sp)
eim.m <- psmelt(Eim.TSSw)
eim.mp <- psmelt(Eim.amp)

# relevel
dist_bc <- (vegdist(Eim.TSSw@otu_table, method="bray"))
res <- pcoa(dist_bc)

EH_sort <- names(sort(res$vectors[,1]))

eim.m$Sample <- factor(eim.m$Sample, levels= EH_sort)
eim.m0$Sample <- factor(eim.m0$Sample, levels= EH_sort) # I will still sort with the distances of all amplicons, so that I can align plots
eim.mp$Sample <- factor(eim.mp$Sample, levels= EH_sort)

# sanity check
EH_sort==levels(eim.m$Sample)

eim.mp$Genus <- as.factor(eim.mp$Genus)


nb.cols <- length(levels(eim.mp$Genus))+1
mycolors <- colorRampPalette(brewer.pal(8, "Dark2"))(nb.cols)

################# now plotting both 18S and 28S genes
Com.m.all <- ggplot(eim.mp, aes(x=Sample, y=Abundance, fill=Species))+
    geom_bar(position="stack", stat="identity")+
#    scale_fill_manual(values=c("forestgreen", "pink", "dodgerblue4",  "darkgray", "darkred"))+
    scale_fill_manual(values=c("darkgoldenrod2","darkolivegreen","cornflowerblue"))+
    labs(fill="Eimeria", x="Sample", y="Eimeria ASV abundance", tag="a")+
    theme_bw(base_size=12)+
    theme(axis.text.x=element_blank(),
      axis.ticks.x=element_blank(),
      legend.key = element_blank(),
      panel.grid.major = element_blank(),
      legend.position="none")
#    coord_flip()

Com.m.all_amp <- ggplot(eim.mp, aes(x=Sample, y=Abundance, fill=Genus))+
    geom_bar(position="stack", stat="identity")+
    scale_fill_manual(values=mycolors)+
    labs(fill="Amplicon", x="Sample", y="Eimeria ASV abundance", tag="b")+
    theme_bw(base_size=12)+
    theme(axis.text.x=element_blank(),
      axis.ticks.x=element_blank(),
      legend.key = element_blank(),
      panel.grid.major = element_blank(),
      legend.position="none")

legend <- get_legend(Com.m.all_amp+
                     guides(fill=guide_legend(override.aes=list(size=6),nrow=3))+
                     theme(legend.text = element_text(size = 9),
      legend.position="bottom"))

legend2 <- get_legend(Com.m.all+
                     guides(fill=guide_legend(override.aes=list(size=4), nrow=1))+
                     theme(legend.text = element_text(face = 'italic'),
                           legend.title=element_text(face="italic"),
                           legend.position="top"))

cowplot::plot_grid(Com.m.all, Com.m.all_amp, nrow=2)

library(cowplot)

Comp_amplicon2 <- cowplot::plot_grid(Com.m.all_amp, legend, Com.m.all, legend2, ncol=1, rel_heights=c(0.8, 0.25,0.8,0.04))

Comp_amplicon2


ggplot2::ggsave(file="fig/FigureS4.pdf", Comp_amplicon2, height=8, width=10, dpi=400)

library(viridis)
library(wesanderson)
library(scales)

pal <- wes_palette("Zissou1", 500, type = "continuous")
Eim_heat_all <- ggplot(eim.m0, aes(Sample, Species, fill=Abundance))+
    geom_tile()+
    labs(y="Eimeria", x="Sample", fill="Eimeria ASV abundance")+
    scale_fill_gradientn(colours = pal)+
    theme_bw(base_size=12)+
    theme(axis.text.y = element_text(face = 'italic'),
      axis.title.y=element_text(face="italic"),
      axis.text.x=element_blank(),
      axis.ticks.x=element_blank(),
      legend.key = element_blank(),
      strip.background = element_rect(colour="black", fill="white"),
      legend.text=element_text(size=10),
#      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position="bottom")

Eim_heat_all

## something seems weird in this plot
Amplicon_ASV <- ggplot(eim.m, aes(Sample, Genus, fill=Abundance))+
    geom_tile()+
      labs(y="Amplicon", x="Sample", fill="Eimeria ASV abundance")+
    scale_fill_gradientn(colours = pal, values=rescale(c(0,0.0005,1)),)+
    theme_bw(base_size=12)+
      theme(axis.title.x=element_blank(),
      axis.text.x=element_blank(),
      axis.ticks.x=element_blank(),
      legend.key = element_blank(),
      strip.background = element_rect(colour="black", fill="white"),
      legend.text=element_text(size=10),
#      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position="bottom")

Sp.m <-ggplot(eim.m[eim.m$Abundance>0,], aes(x=Genus, y=Abundance, fill=Species))+
#    geom_bar(position="stack", stat="identity")+
    geom_point(size=4, colour="gray", shape=21, position=position_jitterdodge(dodge.width=0.8, jitter.width=0.15), alpha=0.4)+
    geom_boxplot(alpha=0.3, colour="black", outlier.shape = NA)+
    scale_fill_manual(values=c("forestgreen", "dodgerblue4", "darkred"))+
    labs(fill="Eimeria", x="", y="Eimeria abundance")+
    theme_bw(base_size=12)+
    guides(fill=guide_legend(ncol=4))+
    theme(axis.text.y = element_text(colour = 'black', size = 10),
       legend.key = element_blank(),
      panel.grid.major = element_blank(),
#      panel.grid.minor = element_blank(),
      legend.text = element_text(face = 'italic'),
      legend.title=element_text(face="italic"),
      legend.position="top")+
    coord_flip()

Sp.m

# sanity check
levels(Eim_heat_all$data$Sample)==levels(Com.m.all_amp$data$Sample)



ggsave("fig/FigureS4_Eimeria_amplicon_sp.pdf", Sp.m, height=4, width=10, dpi=400)
ggsave("fig/FigureS4_Eimeria_amplicon_sp.png", Sp.m, height=4, width=10, dpi=400)

#####################################################################
# Ok, let's try and figure it out what is happening with these co-infections
# removing empty samples
Fer <- subset_taxa(Eim_sp, Species %in% "ferrisi")
sample_data(Eim_sp)$Ferrisi <- sample_sums(Fer)

Fal <- subset_taxa(Eim_sp, Species %in% "falciformis")
sample_data(Eim_sp)$Falciformis <- sample_sums(Fal)

Ver <- subset_taxa(Eim_sp, Species %in% "vermiformis")
sample_data(Eim_sp)$Vermiformis <- sample_sums(Ver)

#Eimra0 <- subset_samples(Eim.ra0, !BMI=="NA")

Eimdf <- sample_data(Eim_sp)
Eimdf$Locality <- as.factor(Eimdf$Locality)
Eimdf <- as.data.frame(Eimdf)
class(Eimdf) <- "data.frame"


dis <- phyloseq::distance(Eim_sp, method="jaccard", type="samples")

dis2 <- phyloseq::distance(prune_samples(rownames(Eimdf[!is.na(Eimdf$BMI),]), Eim_sp), method="jaccard", type="samples")

Eimdf1 <- Eimdf[!is.na(Eimdf$BMI),]


permaPS=adonis2(dis2~
#            Eimdf1$HI+
            Eimdf1$Locality+
            Eimdf1$Year+
            Eimdf1$Sex+
           Eimdf1$BMI,
           permutations = 1000, method = "bray",
           by="margin")

permaPS

library(merTools)
library(MuMIn)

library(lmerTest)

Eimdf1$logFer <- log(1+Eimdf1$Ferrisi)
Eimdf1$logVer <- log(1+Eimdf1$Vermiformis)
Eimdf1$logFal <- log(1+Eimdf1$Falciformis)


summary(Eimdf1$BMI)

BMIm <- lmerTest::lmer(BMI~Ferrisi + Falciformis + Vermiformis +(1|Locality), data=Eimdf1)
BMIm0 <- lmerTest::lmer(BMI~1 + (1|Locality), data=Eimdf1)

summary(BMIm)

#plot(BMIm)
#qqnorm(resid(BMIm))
#qqline(resid(BMIm))
#plotREsim(REsim(BMIm))

r.squaredGLMM(BMIm)

anova(BMIm, BMIm0)
anova(BMIm, test="LRT")

library(lmerTest)
ranova(BMIm)

library(ggeffects)

fe <- ggpredict(BMIm, terms="Ferrisi", type="random")

fa <- ggpredict(BMIm, terms="Falciformis", type="random")

fe$conf.low

fe.eff <- ggplot(Eimdf1, aes(x=Ferrisi, y=BMI))+
    geom_point(colour="gray40", alpha=0.7)+
    geom_line(data=fe, aes(x=x, y=predicted), colour="darkgoldenrod2", size=1.5)+
    geom_ribbon(inherit.aes = FALSE, data=fe, aes(x=x, ymin=conf.low, ymax=conf.high), alpha=0.1, fill="darkgoldenrod2")+
    ylab("Body mass index")+
    xlab("E. ferrisi abundance")+
    theme_bw(base_size=10)+
    theme(axis.title.x = element_text(vjust = 0, size = 12),
          axis.title.y = element_text(vjust = 2, size = 12))


fa.eff <- ggplot(Eimdf1, aes(x=Falciformis, y=BMI))+
    geom_point(colour="gray40", alpha=0.7)+
    geom_line(data=fa, aes(x=x, y=predicted), colour="darkolivegreen", size=1.5)+
    geom_ribbon(inherit.aes = FALSE, data=fa, aes(x=x, ymin=conf.low, ymax=conf.high), alpha=0.1, fill="darkolivegreen")+
    ylab("Body mass index")+
    xlab("E. falciformis abundance")+
    theme_bw(base_size=10)+
    theme(axis.title.x = element_text(vjust = 0, size = 12),
          axis.title.y = element_text(vjust = 2, size = 12))

Figure5 <-cowplot::plot_grid(fe.eff, fa.eff, labels="auto", nrow=1, align="hv")
Figure5
ggplot2::ggsave(file="fig/Figure5.pdf", Figure5, width = 8, height = 4, dpi = 300)

################# now we are plotting Eimeria genotyping from PCR with Eimeria ASV abundance, it is currently broken... Trying to find a way to connect dots from the saame samples with lines
eim_sp$Concatenated <- as.factor(eim_sp$Concatenated)
eim_c <- eim_sp[eim_sp$Abundance>0,]
eim_c <- eim_c[!is.na(eim_c$Concatenated),]


length(unique(as.factor(eim_c$Mouse_ID)))


eim_c[eim_c$Concatenated=="E_vermiformis", c("Abundance", "Mouse_ID", "Group_18S", "Group_COI_1", "Group_COI_2", "Group_ORF", "Concatenated", "Species")]

test=transform(eim_c, dAbundance= ifelse(Concatenated=="ferrisi",
                                          as.numeric(Abundance)-0.25,
                                          as.numeric(Abundance)+0.25))

eim_c$dConcatenated <- 200

eim_c$dConcatenated[eim_c$Species=="falciformis"] <- as.numeric(eim_c$Concatenated)-0.25

eim_c$dConcatenated[eim_c$Species=="ferrisi"] <- as.numeric(eim_c$Concatenated)

eim_c$dConcatenated[eim_c$Species=="vermidormis"] <- as.numeric(eim_c$Concatenated)+0.25

eim_c$dConcatenated[eim_c$Species=="falciformis"]

eim_c$dConcatenated

as.numeric(eim_c$Concatenated)-0.25

eim_c$Concatenated

test$dAbundance

eim_c$Species

                                        #https://stackoverflow.com/questions/44656299/lines-connecting-jittered-points-dodging-by-multiple-groups

TISSUE_MA <- ggplot(eim_c, aes(x=Concatenated, y=Abundance, fill=Species))+
    scale_fill_manual(values=c("forestgreen", "dodgerblue4", "darkred"))+
    geom_point(size=4, shape=21, position=position_jitterdodge(dodge.width=0.75, jitter.width=0.1), alpha=0.7)+
    geom_boxplot(alpha=0.3, colour="black", outlier.shape = NA)+
    geom_line(aes(x=dConcatenated, group=Sample), position=position_jitterdodge(dodge.width=0.75, jitter.width=0.1))+
    labs(x="genotyping PCR from tissue DNA", y="Proportion within all ASVs/ng DNA (log)")+#    geom_bar(position="dodge", stat="identity")+
    guides(fill=guide_legend(ncol=4))+
    theme_classic()+
    theme(axis.text.y = element_text(colour = 'black', size = 14, face = 'italic'),
      axis.title.x=element_blank(),
#      axis.text.x=element_blank(),
      axis.ticks.x=element_blank(),
      legend.key = element_blank(),
#      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      strip.background = element_rect(colour="black", fill="white"),
      legend.text = element_text(colour = 'black', size = 10, face = 'italic'),
      legend.position="none"
      )

#TISSUE_MA

#ggsave("fig/Eimeria_qPCR_MA.pdf", TISSUE_MA, height=4, width=5, dpi=400)
#ggsave("fig/Eimeria_qPCR_MA.png", TISSUE_MA, height=4, width=5, dpi=400)

