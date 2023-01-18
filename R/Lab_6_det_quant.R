library(relaimpo)
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
SA.asv3 <- prune_taxa(rownames(tax_table(Eim2))[3], Eim2)
SA.asv4 <- prune_taxa(rownames(tax_table(Eim2))[4], Eim2)

sensit(SA.asv1)

sensit(SA.asv2)

sensit(SA.asv3)

sensit(SA.asv4)

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

