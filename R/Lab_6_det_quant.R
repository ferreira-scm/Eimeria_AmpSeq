library(lme4)
library(MASS)
library(DescTools)
library(cocor)
library(reshape2)
library(microbiome)
library(phyloseq)

#library(dada2)
#library(DECIPHER)
#library(phangorn)
#library(ShortRead)
#library(treeio)
#library(ape)


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

################## now the alignments
############################
#######################################################################################
seqs <- DNAStringSet(getSequences(colnames(Eim@otu_table)))
seqs2 <- DNAStringSet(getSequences(colnames(Eim2@otu_table)))

## save ASV and reference 18S sequences for aligments
## let's get all the sequences used in Jarquín-Díaz et al. 2019
access.l <- c("JQ993669", "JQ993670","JQ993671","JQ993665", "JQ993659", "KU174470", "KU174461", "KU174464", "KU174480", "KU174483", "KU174462", "KU174463", "KU174468", "KU174479", "KU174449", "KU174450", "KU174465", "KU174467", "KU174474", "KU174451", "KU174459", "KU174485", "KU174487", "KU174472", "JQ993666","JQ993649", "JQ993650", "AF080614", "MH751998", "KT360995", "JF304148", "U40263","JQ993651", "AF311643","AF311644", "JQ993654", "JQ993655", "JQ993656", "JQ993657", "JQ993658", "JQ993660", "KU174475", "JQ993661", "JQ993662", "JQ993663", "JQ993664", "JQ993652", "JQ993667", "KU174454", "KU174469", "KU174481", "KU174456","KU174484", "KU174478", "KU174473", "KU174471",  "KU174455", "KU174457", "KU174476","KU174466", "KU174486", "KU174453","KU174458", "KU174460", "KU174452","KU174482", "AF246717", "KT184355","JQ993653", "AF307880", "AF339489","AF307878", "AF307876", "AF324214","AF339490", "AF339491", "AF307879", "AF339492", "AF307877", "AF311642", "KU192965", "KU192958", "KU192936", "KU192961", "KU192956", "KU192931", "KU192938", "KU192916")
eim.db <- read.GenBank(access.l)
eim.db.ID <- paste(names(eim.db), attr(eim.db, "species"), sep="_")
# extra 18S falciformis
Eimf <- read.GenBank("KT184339.1")
Eimf.ID <- paste(names(Eimf), attr(Eimf, "species"), sep="_")

# convert to DNAStringset
eim.DB <- eim.db %>% as.character %>% lapply(.,paste0,collapse="") %>% unlist %>% DNAStringSet
names(eim.DB) <- eim.db.ID

Eim.f <- Eimf %>% as.character %>% lapply(.,paste0,collapse="") %>% unlist %>% DNAStringSet
names(Eim.f) <- Eimf.ID

names(eim.DB)

refEim <- readDNAStringSet("/SAN/Susanas_den/AmpMarkers/wildEimeria18S/Eim_ref.fa")
names(refEim) <- gsub("(\\s)", "_", names(refEim))

#gotta correct here
which(names(refEim)=="MH751946.1_Eimeria_vermiformis")
names(refEim)[22] <- "MH751946.1_Eimeria_ferrisi"

names(seqs2) <- c("Lab_single-multi-amplicon_ASV1", "Lab_single-multi-amplicon_ASV2", "Lab_single-amplicon_ASV3", "Lab_single-amplicon_ASV4")

allSeqs <- c(eim.DB, Eim.f, refEim, seqs2)

writeFasta(allSeqs, "tmp/Eimeria_seqs.fa")

