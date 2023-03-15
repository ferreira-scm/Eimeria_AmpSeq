#!/usr/bin/Rscript

library(ggplot2)
#library(dada2)
#library(MultiAmplicon, lib.loc="/usr/local/lib/R/site-library/")
#library(reshape)
#remotes::install_github("vmikk/metagMisc")
library(metagMisc)
require(grid)
library(Hmisc)
library(phyloseq)
#library(grid_extra)
require(cowplot)
library(RColorBrewer)
library(dplyr)
## using the devel
#devtools::load_all("/SAN/Susanas_den/MultiAmplicon/")
source("R/PlottingCor.R")

## Microfluidics PCR multiple amplicons
all.PS.l.slv <- readRDS("tmp/Lab/PhyloSeqList_All_Tax_New.Rds")

## Standard PCR Single amplicon 18S
sin.PS18S.slv <- readRDS("tmp/Lab/PhyloSeqList18S_SILVA.Rds")
sin.PS18S.slv <- sin.PS18S.slv[[2]]

# let's filter
f.sin18.slv <- fil(sin.PS18S.slv)

## filtering MA by amplicon
f.all.l.slv <- list()
for (i in 1:48) {
    try(f.all.l.slv[[i]] <- fil(all.PS.l.slv[[i]]), silent=TRUE)
}
# And now merge
f.all.lp.slv <- f.all.l.slv[[1]]
for (i in 2:47){
    f.all.lp.slv <- try(merge_phyloseq(f.all.lp.slv,f.all.l.slv[[i]]))
}

# and transform TSS
sin18TSS.slv  <- transform_sample_counts(f.sin18.slv, function(x) x / sum(x)) # single amplicon
TSS.wang <- transform_sample_counts(f.all.l.slv[[37]], function(x) x / sum(x)) #"wang" amplicon from multi amplicon

## OK, now we want all a dataset with only Eimeria ASVs
Eim2 <- subset_taxa(sin18TSS.slv, Genus%in%"g__Eimeria") # for single amplicon
Eim <- subset_taxa(TSS.wang, Genus%in%"g__Eimeria") # for multi amplicon but only one amplicon
