# Eimeria species
source("R/5_Wild_Eimeria_spp_assignment.R")

library(cowplot)
library(ggeffects)
library(viridis)
library(wesanderson)
library(scales)
library(merTools)
library(MuMIn)
library(lmerTest)

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
                     guides(fill=guide_legend(override.aes=list(size=6),nrow=4))+
                     theme(legend.text = element_text(size = 9),
      legend.position="bottom"))

legend2 <- get_legend(Com.m.all+
                     guides(fill=guide_legend(override.aes=list(size=4), nrow=1))+
                     theme(legend.text = element_text(face = 'italic'),
                           legend.title=element_text(face="italic"),
                           legend.position="top"))

Comp_amplicon2 <- cowplot::plot_grid(legend, Com.m.all_amp, Com.m.all, legend2, ncol=1, rel_heights=c(0.25,0.8,0.8,0.04))
Comp_amplicon2

ggplot2::ggsave(file="fig/FigureS4.pdf", Comp_amplicon2, height=8, width=10, dpi=400)

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
