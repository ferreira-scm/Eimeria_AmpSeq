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


## eimeria sampling depth
summary(sample_sums(subset_taxa(f.sin18.slv, Genus%in%"g__Eimeria")))

summary(sample_sums(f.sin18.slv))

summary(sample_sums(subset_taxa(f.all.l.slv[[37]], Genus%in%"g__Eimeria")))



Eim <- subset_taxa(TSS.wang, Genus%in%"g__Eimeria") # for multi amplicon but only one amplicon
