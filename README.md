# Analysis scripts for Chamberlin et al 2022/23
This repo includes the analysis used in the manuscript "Differences in molecular sampling and data processing explain variation among single-cell and single-nucleus RNA-seq experiments"

In this analysis we re-processed raw data from Yao et al 2021 and Thrupp et al 2020 to obtain UMI counts with and without introns. 

Raw data were acquired from:
Yao et al nuclei: https://data.nemoarchive.org/biccn/grant/u19_zeng/zeng/transcriptome/sncell/10x_v3/mouse/raw/MOp/

Yao et al cells: https://data.nemoarchive.org/biccn/grant/u19_zeng/zeng/transcriptome/scell/10x_v3/mouse/raw/MOp/

Thrupp et al: GEO accessions GSE153807 and GSE137444

# final_bias_analysis.Rmd:
Includes main statistical analysis and plotting for the figures 

# preprocessing.Rmd:
Shows how STARsolo outputs were imported and collated

# prep_QC.Rmd:
Confirms that STARsolo output matches the original CellRanger output from Yao et al.
