# Eimeria_AmpSeq

Code, figures and intermediate files for "Eimeria quatification using
amplicom sequencing" manuscript

Start running the code from Lab_6* and/or Wild_8*
Better, yet, we renamed the script so you start at 1_*

Here is how the figures and tables in the manuscript are created in
the code.


| File reference | Created in                     | File in ropository                             |
|:---------------|--------------------------------|:-----------------------------------------------|
| Table 1        | R/Lab_6_det_quant.R            | text only                                      |
| Table 2        | R/Lab_7_quant_ASV_analysis.R   | fig/Eimeira_asv_MA.txt                         |
| Table 3        | R/Wild_9_Eimeria_spp_com_str.R | text only                                      |
| Figure 1       | R/Lab_7_quant_ASV_analysis.R   | fig/Figure1-ASV_concordance.pdf                |
| Figure 2       | R/Lab_7_quant_ASV_analysis.R   | fig/Figure2.png                                |
| Figure 3a **   | R/Wild_8_AssignEim_tree_cor.R  | text only                                      |
| Figure 3b      | R/Wild_8_AssignEim_tree_cor.R  | tmp/28S_tree/Eimeira_ASV28S_align.fasta.contree|
| Figure 3c      | R/Wild_8_AssignEim_tree_cor.R  | fig/Figure4_Eimeria_ASVs_Network.pdf           |
| Figure 6       | R/Wild_9_Eimeria_spp_com_str.R | fig/Figure6_BMI_Eimeria.pdf                    |
| Table S1       | R/Lab_6_det_quant.R            | text only                                      |
| Table S2       | R/Wild_9_Eimeria_spp_com_str.R | text only                                      |
| Figure S1      | R/Lab_6_det_quant.R            | fig/FigureS1_single-amplicon                   |
| Figure S2      | R/Lab_6_det_quant.R            | fig/FigureS2_multi-amplicon                    |
| Figure S3      | R/Wild_9_Eimeria_spp_com_str.R | FigureS3_distribution.pdf                      |

** Figure 3a reference tree (not collapsed) is found in tmp/reference_tree/Eimeria18S_ref.fasta.contree, and amplicon sequence variant (ASV) trees are found in Eimeria_AmpSeq/tmp/amplicon_alignments/amplicon*.contree

And a short summary of what each script does

## R/Lab_1_Data_preparation.R and R/Lab_2_qPCR_data_preparation.R prepare the metadata 
R/Lab_1_Data_preparation.R produces sample.data
R/Lab_2_qPCR_data_preparation.R calls sample.data and produces sdt (final metadata file)

## R/Lab_3_multiamplicon.R preprocesses sequencing reads from the multi-amplicon with MA (dada2)
calls for R/toPhyloseq.R
saves the phyloseq objects for each amplicon as a list: tmp/Lab/PhyloSeqList_All_Tax_New.Rds
saves ONE phyloseq object with the ASVs from ALL amplicons in it: tmp/Lab/PhyloSeqData_All_Tax_New.Rds

##R/Lab_4_singleamplicon.R preprocesses sequencing reads from single-amplicon with MA (dada2)
calls for R/toPhyloseq.R
saves the phyloseq objects for each amplicon as a list: tmp/Lab/PhyloSeqList18S_SILVA.Rds
saves ONE phyloseq object with the ASVs from ALL amplicons in it: tmp/Lab/PhyloSeqData18S_SILVA.Rds

##R/Lab_5_filtering.R does quality filtering and transforms phyloseq objects 
calls for R/PlottingCor.R for function "fil"
Does some basic ASV per amplicon exploration
Transforms to total sum scaling (relative abundance)
Rescales to sample DNA concentration
Subsets phyloseq objects for Eimeria only reads for downstream analysis

## R/Lab_6_det_quant.R detection analysis and analysis of effect of normalization on the quantification and plotting
calls for R/Lab_5_filtering.R to use its the phyloseq objects
uses function "sensit" from R/PlottingCor.R

## R/Lab_7_quant_ASV_analysis.R quantification analysis
calls for R/Lab_5_filtering.R to use its the phyloseq objects
linear regression analyses of the association of qPCR estimate with the different eimeria ASVs
correlation with oocyst counts (remove)
lots of plots including figure 1, figure 2, table 2

## R/PlottingCor.R defines functions for filtering and plotting correlations

## R/toPhyloseq.R defines function to convert MA object to phyloseq, while this function is broken in MA package

##### to be continued ####



