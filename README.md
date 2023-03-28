# Eimeria_AmpSeq

Code, figures and intermediate files for "Amplicon sequencing allows differential quantification of closely related parasite species: an example from rodent coccidia (Eimeria)" manuscript

Start running the code from 1_Lab_filter.R

Here is how the figures and tables in the manuscript are created in the code.


| File reference | Created in                         | File in repository                             |
|:---------------|------------------------------------|:-----------------------------------------------|
| Table 1        | R/3_Lab_Fig1_Fig2.R                | text only                                      |
| Table 2        | R/2_Lab_S1_S2_detec.R              | text only                                      |
| Table 3        | R/6_Wild_Eimeriaspp_community.R    | text only                                      |
| Figure 1       | R/3_Lab_Fig1_Fig2.R                | fig/Figure1.pdf                                |
| Figure 2       | R/3_Lab_Fig1_Fig2.R                | fig/Figure2.pdf                                |
| Figure 3 **   | R/5_Wild_Eimeria_spp_assignment.R  | composed image **                              |
| Figure 4      | R/5_Wild_Eimeria_spp_assignment.R  | tmp/28S_tree/Eimeira_ASV28S_align.fasta.contree|
| Figure 5       | R/5_Wild_Eimeria_spp_assignment.R  | fig/Figure5.pdf                                |
| Figure 6       | R/6_Wild_Eimeriaspp_community.R    | fig/Figure6.pdf                                |
| Table S1       | R/2_Lab_S1_S2_detec.R              | text only                                      |
| Table S2       | R/6_Wild_Eimeriaspp_community.R    | text only                                      |
| Figure S1      | R/2_Lab_S1_S2_detec.R              | fig/FigureS1.pfd                               |
| Figure S2      | R/2_Lab_S1_S2_detec.R              | fig/FigureS2.pdf                               |
| Figure S3      | R/3_Lab_Fig1_Fig2.R                | FigureS3.pdf                                   |
| Figure S4      | R/6_Wild_Eimeriaspp_community.R    | FigureS4.pdf                                   |


** Figure 3 reference tree (not collapsed) is found in tmp/reference_tree/Eimeria_reference.fasta.contree, and amplicon sequence variant (ASV) trees are found in Eimeria_AmpSeq/tmp/amplicon_alignments/amplicon*.contree

....
Preprocessing steps are done in scripts starting with "Lab_prepro" and "Wild_prepro". To reproduce these, it is necessary to download it from SRA (insert project number here)
Enumerated scripts can be reproduced as all required files are within the repo.
