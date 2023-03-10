sota <- read.csv("https://raw.githubusercontent.com/derele/Mouse_Eimeria_Field/master/data_products/SOTA_Data_Product.csv") # commit c463775

# removing columns I won't use
torm <- c("mtBamH", "YNPAR", "X332", "X347", "X65", "Tsx", "Btk", "Syap1", "Es1", "Gpd1","Idh1", "Mpi",
          "Np", "Sod1", "Es1C","Gpd1C", "Idh1C", "MpiC", "NpC", "Sod1C", "HI_NLoci", "Zfy2", "Y", "Left_Embryo", "Right_Embryo",
          "Address", "Date_count", "N_oocysts_sq1", "N_oocysts_sq2", "N_oocysts_sq3", "N_oocysts_sq4",
          "N_oocysts_sq5", "N_oocysts_sq6", "N_oocysts_sq7", "N_oocysts_sq8", "mean_neubauer",
          "PBS_dil_in_mL", "Ncells", "Feces_Weight", "Fleas", "Liver", "Right_Ovarium_Weight",
          "Left_Ovarium_Weight", "Seminal_Vesicles_Weight", "Left_Testis", "Right_Testis", "Trap_Date",
          "Ticks", "Host", "Sperm")
sota <- sota[,which(!names(sota)%in%torm)]

##Tiny adjustments to metadata
require(dplyr)
sota$Year<- as.factor(sota$Year)
sota$HI <- as.numeric(sota$HI)
sota$HI<-round(sota$HI, 2)
sota$Sex <- as.factor(sota$Sex)
sota%>%
    mutate(BMI=Body_Weight/((Body_Length)^2)) -> sota
sota$Longitude<-round(sota$Longitude, 4)
sota$Latitude<-round(sota$Latitude, 4)
sota$Locality <- paste(sota$Latitude, sota$Longitude)
sota$Transect[sota$Latitude<50] <- "HZ_BAV"
sota$Transect[sota$Latitude>50] <- "HZ_BR"
sota %>%
    mutate(Pinworms = sota$Aspiculuris_sp+sota$Syphacia_sp) -> sota
sota <- sota[!is.na(sota$Mouse_ID),]

## we load Victor data with sequencing information and some qPCR data until raw data is incorporated into sota
ma.meta <- read.csv("/home/victor/AA_HMHZ/Sample_selection_Metabarcoding_Complete.csv", dec=",", stringsAsFactors=FALSE)
ma.meta$Concentration <- as.numeric(ma.meta$Concentration)
ma.meta$Chip_number <- as.factor(ma.meta$Chip_number)
ma.meta%>%
    mutate(Seq_Run = case_when(Chip_number%in% c(1:8) ~ "Pool_1",
                                     Chip_number%in% c(9:15) ~ "Pool_2")) -> ma.meta
ma.meta$Seq_Run <- as.factor(ma.meta$Seq_Run)

tokeep <- c("Mouse_ID", "Group_18S","Group_COI_1", "Group_COI_2", "Group_ORF", "Concatenated", "Concentration", "Chip_number", "Chip_position", "Library_pool", "Seq_Run")
ma.meta <- ma.meta[,tokeep]

sota <- merge(sota, ma.meta, by="Mouse_ID", all=TRUE)
rownames(sota) <- make.unique(sota$Mouse_ID) ##Works when MA contains single run data

rownames(sota)
