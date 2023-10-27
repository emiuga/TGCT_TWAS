#!/bin/bash

# Input GWAS summary statistics (imputed)
# Date: July 26, 2021
# Author: Emilio Ugalde
# Reference: http://gusevlab.org/projects/fusion/
# Description: prepare input sumstats for TWAS analysis with FUSION.

#The primary input is genome-wide summary statistics in LD-score format. At minimum, this is a flat file with a header row containing the following fields:
#    SNP – SNP identifier (rsID)
#    A1 – first allele (effect allele)
#    A2 – second allele (other allele)
#    Z – Z-scores, sign with respect to A1. # To convert p-values to z-scores: qnorm(P-value/2, lower.tail=FALSE)*sign(Effect)

# and subsequent data rows for each SNP (all white-space separated). Additional columns are allowed and will be ignored. 
# We recommend using the LDSC munge_stats.py utility for converting GWAS summary data into this format, which detects and reports many common pitfalls. !!!

cd /nfs/GENETEC/TWAS/

# TECAC susmstats filtered by AF>0.1, heter>0.001, and INFO>0.8
#sum=/nfs/GENETEC/TECAC_GWAS/Data/DerivedData/TECAC_META_8site_validation.INFO_08.txt
sum=/nfs/GENETEC/TECAC_GWAS_sumstats/Data/OriginalData/TECAC_META_8site_validation.INFO_08_AF01.txt

# TECAC imputed sumstats
GWAS_prefix="TECAC_META_8site_validation"
dir_sum=/nfs/GENETEC/TECAC_GWAS_sumstats/Data/DerivedData/imputation/processed_summary_imputation
#sum=$dir_sum/imputed_$GWAS_prefix.chr$i.txt.gz 

#---------------
# Format sumstats
#---------------

mkdir -p Data/tmp/sumstats_imputed

# Extract required columns
for i in `seq 22`; do 
	zcat $dir_sum/imputed_$GWAS_prefix.chr$i.txt.gz | awk '{print $1, $5, $6, $10}'  > Data/tmp/sumstats_imputed/$GWAS_prefix.chr$i.txt	
	# Re-name header
	sed -i 's/variant_id/SNP/g' Data/tmp/sumstats_imputed/$GWAS_prefix.chr$i.txt
	sed -i 's/effect_allele/A1/g' Data/tmp/sumstats_imputed/$GWAS_prefix.chr$i.txt
	sed -i 's/non_A1/A2/g' Data/tmp/sumstats_imputed/$GWAS_prefix.chr$i.txt
	sed -i 's/zscore/Z/g' Data/tmp/sumstats_imputed/$GWAS_prefix.chr$i.txt
done




