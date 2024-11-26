# Analysis of scRNAseq data from patients with T1D and healthy donors
 This repository contains an analysis of the single-cell RNA sequencing of PBMC samples from healthy donors and donors with Type 1 Diabetes Mellitus (T1D). The results will be published in:

An imbalance of naïve and effector T-cell phenotypes in early type 1 diabetes across conventional and regulatory subsets 

Veronika Niederlova1,2, Ales Neuwirth1, Vit Neuman3, Juraj Michalik1, Bela Charvatova1, Martin Modrak4, Zdenek Sumnik3, Ondrej Stepanek1

1 Laboratory of Adaptive Immunity, Institute of Molecular Genetics of the Czech Academy of Sciences, Prague, Czechia
2 Department of Cell Biology, Faculty of Science, Charles University in Prague, Czechia
3 Department of Pediatrics, 2nd Faculty of Medicine, Charles University in Prague & Motol University Hospital, Prague, Czechia
4 Department of Bioinformatics, Second Faculty of Medicine, Charles University, Prague, Czechia

The functions and libraries needed for running the analysis are specified in the following file:
* diabetes_analysis_v07.R
  
This analysis itself contains 23 parts:
* 01_Cohort_characterization.ipynb
* 02_CD8_v06_part1_annotation.ipynb
* 03_CD4_v06_part1_annotation.ipynb
* 04_CD8_part2_tcr_analysis.ipynb
* 05_CD4_part2_tcr_analysis.ipynb
* 06_CD4_v05_part3_DIAvsCTRL.ipynb
* 07_CD8_v05_part3_DIAvsCTRL.ipynb
* 08_Subset_analysis_Treg.ipynb
* 09_Subset_analysis_Tgd.ipynb
* 10_DE_genes_networks.ipynb
* 11_Population_dynamics.ipynb
* 12_FACS_data.ipynb
* 13_HLAfreq.ipynb
* 14_HLA_and_TCR.ipynb
* 15_TCR_HPAP.ipynb
* 16_Kallionpaa.ipynb
* 17_Honardoost.ipynb
* 18_NatComm_Diabetes_Data.ipynb
* 19_ParseBio_Dia1M.ipynb
* 20_HPAP_5p_v01.ipynb
* 21_Other_data_Validation.ipynb
* 22_BTN3A2_expression.ipynb
* 23_Treg_Validaton.ipynb

Processed data required to generate figures for the manuscript are available on Zenodo: 

