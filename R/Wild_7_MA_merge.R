#!/usr/bin/Rscript
library(phyloseq)

####Merge phyloseq objects
PS1 <-  readRDS(file="tmp/Wild/PhyloSeqCombi_HMHZ_1_1.Rds")
PS2 <-  readRDS(file="tmp/Wild/PhyloSeqCombi_HMHZ_1_2.Rds")    
PSa <- merge_phyloseq(PS1, PS2)
PS3 <-  readRDS(file="tmp/Wild/PhyloSeqCombi_HMHZ_2_1.Rds")
PS4 <-  readRDS(file="tmp/Wild/PhyloSeqCombi_HMHZ_2_2.Rds")
PSb <- merge_phyloseq(PS3, PS4)
PS <- merge_phyloseq(PSa, PSb)
saveRDS(PS, file="tmp/Wild/PhyloSeqCombi_HMHZ_All.Rds") ###Results from full + test run pool 1
#sample_names(PS)

rm(PS1)
rm(PS2)
rm(MA)
rm(PS.l)
rm(PS3)
rm(PS4)
rm(PSa)
rm(PSb)

##Merge PS.list
PS.l.1.1<- readRDS(file="tmp/Wild/PhyloSeqList_HMHZ_1_1.Rds")
PS.l.1.2<- readRDS(file="tmp/Wild/PhyloSeqList_HMHZ_1_2.Rds")
PS.l.2.1<- readRDS(file="tmp/Wild/PhyloSeqList_HMHZ_2_1.Rds")
PS.l.2.2<- readRDS(file="tmp/Wild/PhyloSeqList_HMHZ_2_2.Rds")

###Eliminate the empty primers first and then merge the phyloseq lists... empty list make the next function bug
PS.l.1.1[sapply(PS.l.1.1, is.null)]<- NULL
PS.l.1.2[sapply(PS.l.1.2, is.null)]<- NULL
PS.l.2.1[sapply(PS.l.2.1, is.null)]<- NULL
PS.l.2.2[sapply(PS.l.2.2, is.null)]<- NULL

##It is necessary to eliminate those that are empty in both
along<- names(PS.l.1.1)
PS.l.1 <- lapply(along, function(i) merge_phyloseq(PS.l.1.1[[i]], PS.l.1.2[[i]])) ##Merge all the information from both experiments
names(PS.l.1) <- names(PS.l.1.1)

along<- names(PS.l.2.1)
PS.l.2 <- lapply(along, function(i) merge_phyloseq(PS.l.2.1[[i]], PS.l.2.2[[i]])) ##Merge all the information from both experiments
names(PS.l.2) <- names(PS.l.2.1)

along<- names(PS.l.2)
PS.l <- lapply(along, function(i) merge_phyloseq(PS.l.1[[i]], PS.l.2[[i]])) ##Merge all the information from both experiments
names(PS.l) <- names(PS.l.2)

#### this is no longer necessary, but better check
#x<- names(PS.l)
#x

#x[6]<-"BGf_132_F.BGr_132_R"
#x[22]<-"LSU_Fwd_2_3Mod_55_F.LSU_Rev_4_54_R"
#x[28]<-"NLF184cw_74_F.NL818cw_74_R"
#names(PS.l)<- x

saveRDS(PS.l, file="tmp/Wild/PhyloSeqList_HMHZ_All.Rds")

rm(PS.l.1)
rm(PS.l.1.1)
rm(PS.l.1.2)
rm(PS.l.2)
rm(PS.l.2.1)
rm(PS.l.2.2)
rm(along)

