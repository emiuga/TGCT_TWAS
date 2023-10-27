#!/bin/bash

# GWAS harmonization (autosomal chr)
# Emilio Ugalde
# Mar 7, 2022
# ssh emiuga@vector.meb.ki.se

# Description: parse GWAS sumstas to the format required by MetaXcan.
# Software: MetaXcan's summary-gwas-imputation (GWAS tools) https://github.com/hakyimlab/summary-gwas-imputation/wiki/GWAS-Harmonization-And-Imputation
# Tutorial: https://github.com/hakyimlab/MetaXcan/wiki/Tutorial:-GTEx-v8-MASH-models-integration-with-a-Coronary-Artery-Disease-GWAS

# Version control:
# 	- v1.1: uses liftover to convert hg19 to hg38 genome built.

# Output: 
#	- Harmonized (parsed) TECAC sumstats.

# To submit job in slurm (not use for this scripts):
# cd /nfs/GENETEC/TECAC_GWAS_sumstats/jobs
# sbatch -n 2 -t 120 --mem=8000 ../Programs/imputation/0.1.2.GWAS_sumstats_harmonization_autosom_v1.1.sh &> GWAStools.sum_harmoni.20220307
# To monitor: scontrol show job 1356361 | less
# To see output when finished: less slurm-1356361.out

cd /nfs/GENETEC/TECAC_GWAS_sumstats

# OUPUT directory
OUTPUT="Data/harmonized_gwas/hg38"
mkdir -p $OUTPUT

# INPUT data
gwas_parsed="Data/DerivedData/imputation/tmp/parsing/parsed_tecac_"
panel_hrc="Data/reference_panel/HRC.r1-1.GRCh37.wgs.mac5.sites.metadata.gz"

# SetUp variables
GWAS_TOOLS="/nfs/GENETEC/TWAS/Methods/summary-gwas-imputation/src"
REPO="/nfs/GENETEC/TWAS/Data/imputation_sumstats/sample_data/data"
gwas_prefix="TECAC_META_8site_validation"


#--------------------------------
# GWAS parsing/harmonization
#--------------------------------
# Call virtual environment with installed libraries
source /nfs/GENETEC/TWAS/Methods/summary-gwas-imputation/bin/activate

#------
# Autosomal
#------
# -snp_reference_metadata $panel_hrc METADATA \
time for i in {1..20}; do
python $GWAS_TOOLS/gwas_parsing.py \
-gwas_file ${gwas_parsed}chr$i.txt.gz \
-snp_reference_metadata $REPO/reference_panel_1000G/variant_metadata.txt.gz METADATA \
-liftover $REPO/liftover/hg19ToHg38.over.chain.gz \
-output_column_map name variant_id \
-output_column_map Allele1 non_effect_allele \
-output_column_map Allele2 effect_allele \
-output_column_map Freq1 frequency \
-output_column_map Effect effect_size \
-output_column_map chr chromosome \
-output_column_map position position \
-output_column_map P-value pvalue \
--chromosome_format \
--insert_value sample_size 189839 --insert_value n_cases 10156 \
-output_order variant_id panel_variant_id chromosome position effect_allele non_effect_allele frequency pvalue zscore effect_size standard_error sample_size n_cases \
-output $OUTPUT/$gwas_prefix.chr$i.parsed_hg38.gz &> $OUTPUT/$gwas_prefix.chr$i.parsed_hg38.log
done







