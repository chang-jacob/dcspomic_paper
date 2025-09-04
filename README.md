# Differential Colocalization in Spatial Omics (DC-SPOMIC)
[![DOI](https://zenodo.org/badge/1015558473.svg)](https://doi.org/10.5281/zenodo.17058880)

This repository contains the code necessary to reproduce results presented in the DC-SPOMIC manuscript. 

## analysis 
This directory contains all .R scripts. It is split into two subdirectories: (1) simulation and (2) crc_tma. 

### simulation
This directory contains two subdirectories for the two simulation studies performed to assess the validity of the DC-SPOMIC framework on simulated data. The scripts contain code to generate and analyze simulated data. 

### crc_tma
This directory contains the scripts used to preproces and analyze the colorectal cancer (CRC) tissue microarray (TMA) dataset. The preprocessing of the CRC TMA can be found in the CODEX_ColonCancer_Schurch.Rmd markdown file. This markdown was copied from the preprocessing performed by the authors of the spicyR [paper](https://github.com/nickcee/spicyRPaper). 
