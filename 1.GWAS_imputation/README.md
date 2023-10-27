# Imputation of GWAS summary statistics (TECAC) 

Harmonization and imputation of TGCT-GWAS summary data to improve overlap with prediction models based on 1000G LD reference panel.

- Description: imputation of Z-score TECAC sumstats to cover finned-mapped SNPs.
- Software: GWAS tools from MetaXcan suite [summary-gwas-imputation](https://github.com/hakyimlab/summary-gwas-imputation)
- Guidelines: see [MetaXcan wiki](https://github.com/hakyimlab/MetaXcan/wiki/Best-practices-for-integrating-GWAS-and-GTEX-v8-transcriptome-prediction-models)

## Analysis steps:
1. GWAS sumstats harmonization
2. GWAS imputation
3. Imputation post-processing
 
## PREREQUISITES
The basic requirements for running GWAS summary-imputation are python>=3.5.

NOTE: installations for TGCT-TWAS project were done under a virtual environment.
Se scripts in: `0.0.1.Install_required_libraries_v1.sh`

To activate virtual environment:
`source Programs/summary-gwas-imputation/bin/activate`

