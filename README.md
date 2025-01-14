# Analysis of scRNAseq data from patients with T1D and healthy donors
![Summary of the paper](./pics/dia_workflow.jpg) 

This repository contains an analysis of the single-cell RNA sequencing of PBMC samples from healthy donors and donors with Type 1 Diabetes Mellitus (T1D). The results will be published in:

><b><i>An imbalance of naïve and effector T-cell phenotypes in early type 1 diabetes across conventional and regulatory subsets </b></i>
><p>Veronika Niederlova<sup>1,2</sup>, Ales Neuwirth<sup>1</sup>, Vit Neuman<sup>3</sup>, Juraj Michalik<sup>1</sup>, Bela Charvatova<sup>1</sup>, Martin Modrak<sup>4</sup>, Zdenek Sumnik<sup>3</sup>, Ondrej Stepanek<sup>1</sup>
>
><sup>1</sup> Laboratory of Adaptive Immunity, Institute of Molecular Genetics of the Czech Academy of Sciences, Prague, Czechia  
><sup>2</sup> Department of Cell Biology, Faculty of Science, Charles University in Prague, Czechia  
><sup>3</sup> Department of Pediatrics, 2nd Faculty of Medicine, Charles University in Prague & Motol University Hospital, Prague, Czechia  
><sup>4</sup> Department of Bioinformatics, Second Faculty of Medicine, Charles University, Prague, Czechia  
></sup>
<br/>


## Requirements
Running the code requires R version 4.2.1 and higher. Specific package versions might be necessary. For more instructions, please take a look at the SessionInfo included in each analysis file. 

The functions and libraries needed for running the analysis are specified in the following file:
* diabetes_analysis_v07.R

## Instructions
For running the code, please clone the repository and download the preprocessed data deposited on Zenodo: [DOI: 10.5281/zenodo.14222418](https://zenodo.org/records/14222418). After the download completes, data files should be placed in the folder `data` of the cloned repository following the structure specified on Zenodo, i.e. initial files should be placed in the `init` and processed files should be placed in the `processed` folder.

The analysis contains 23 parts, the scripts for which are contained within `code02_analysis_of_init` folder:
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

Each file guides you through the analysis as it was performed. Please note that only some of the results presented in these notebooks were included in the manuscript. Additionally, the folder `code01_raw_to_init` contains codes to generate initial data sets (in `init` subfolder of `data` folder) from raw data. Since `init` data are included on Zenodo site you do not have to (and it is impossible without raw data) run these scripts; they are provided for purely informative purpose on general pre-processing strategy.

## Data availability
Processed data required to generate figures for the manuscript are available on Zenodo: 
[DOI: 10.5281/zenodo.14222418](https://zenodo.org/records/14222418) 
Raw data were not deposited to protect the identity of study participants.
