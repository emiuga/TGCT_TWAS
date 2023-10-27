#!/bin/bash

# GWAS imputation: autosomal chromosomes (1000Genomes panel)
# Emilio Ugalde
# Mar 7, 2022
# ssh emiuga@vector.meb.ki.se

# Description: impute TECAC sumstats to match 1000G reference panel.
# Software: MetaXcan's summary-gwas-imputation (GWAS tools) https://github.com/hakyimlab/summary-gwas-imputation/wiki/GWAS-Harmonization-And-Imputation
# Tutorial: https://github.com/hakyimlab/MetaXcan/wiki/Tutorial:-GTEx-v8-MASH-models-integration-with-a-Coronary-Artery-Disease-GWAS

#https://github.com/hakyimlab/summary-gwas-imputation/wiki/GWAS-Harmonization-And-Imputation
#The script gwas_summary_imputation.py imputes missing summary stats (actually, their zscore) using the covariance of the reference data set and the existing summary stats, within a specification of a region (such as Berissa & Pickrell's approximately independent LD blocks, appropriately lifted over). The underlying algorithm is the one from DIST and IMPG.

cd /nfs/GENETEC/TECAC_GWAS_sumstats

# To submit job (chr1-4):
# cd /nfs/GENETEC/TECAC_GWAS_sumstats/jobs
# sbatch -n 8 -t 600 --mem=32000 ../Programs/imputation/0.1.3.GWAS_sumstats_imputation_v1.1.sh &> GWAStools.sum_imputation.20220307.batch_v1.1
# To monitor: scontrol show job 1356362 | less
# To see output when finished: less slurm-1356362.out

# Re-submit job (chr5-20):
# cd /nfs/GENETEC/TECAC_GWAS_sumstats/jobs
# sbatch -n 8 -t 1080 --mem=32000 ../Programs/imputation/0.1.3.GWAS_sumstats_imputation_v1.1.sh &> GWAStools.sum_imputation.20220307.batch_v1.1b
# To monitor: scontrol show job 1356739 | less
# To see output when finished: less slurm-1356739.out

# INPUT data (Output from GWAS harmonization)
INPUT="Data/harmonized_gwas/hg38"
GWAS_prefix="TECAC_META_8site_validation"

# OUPUT directory
OUTPUT="Data/DerivedData/imputation"
DATE="20220307"
mkdir -p $OUTPUT/results
# Note: log files moved to logs/imputation/$DATE

# SETUP
# Call virtual environment with installed libraries
source /nfs/GENETEC/TWAS/Methods/summary-gwas-imputation/bin/activate

# Paths
GWAS_TOOLS="/nfs/GENETEC/TWAS/Methods/summary-gwas-imputation/src"
REPO="/nfs/GENETEC/TWAS/Data/imputation_sumstats/sample_data/data"

# ANALYSIS

# test on chr 21 and 22
#for i in {1..21}; do
#time python $GWAS_TOOLS/gwas_summary_imputation.py \
#-by_region_file $REPO/eur_ld.bed.gz \
#-gwas_file $INPUT/$GWAS_prefix.chr$i.parsed_hg38.gz \
#-parquet_genotype $REPO/reference_panel_1000G/chr$i.variants.parquet \
#-parquet_genotype_metadata $REPO/reference_panel_1000G/variant_metadata.parquet \
#-window 100000 \
#-parsimony 7 \
#-chromosome $i \
#-regularization 0.1 \
#-frequency_filter 0.01 \
#-sub_batches 10 \
#--standardise_dosages \
#--cache_variants \
#-output $OUTPUT/results/${GWAS_prefix}_chr${i}_reg0.1_ff0.01_by_region_$DATE.txt.gz &> $OUTPUT/results/${GWAS_prefix}_chr${i}_reg0.1_ff0.01_by_region_$DATE.txt.log
#done

# complete analysis
#for i in {1..4}; do
for i in {5..20}; do
time python $GWAS_TOOLS/gwas_summary_imputation.py \
-by_region_file $REPO/eur_ld.bed.gz \
-gwas_file $INPUT/$GWAS_prefix.chr$i.parsed_hg38.gz \
-parquet_genotype $REPO/reference_panel_1000G/chr$i.variants.parquet \
-parquet_genotype_metadata $REPO/reference_panel_1000G/variant_metadata.parquet \
-window 100000 \
-parsimony 7 \
-chromosome $i \
-regularization 0.1 \
-frequency_filter 0.01 \
-sub_batches 10 \
--standardise_dosages \
--cache_variants \
-output $OUTPUT/results/${GWAS_prefix}_chr${i}_reg0.1_ff0.01_by_region_$DATE.txt.gz &> $OUTPUT/results/${GWAS_prefix}_chr${i}_reg0.1_ff0.01_by_region_$DATE.txt.log
done

# re-locate log files to dir: logs/imputation/20220307
mkdir -p logs/imputation/$DATE

for i in {1..20}; do
mv $OUTPUT/results/${GWAS_prefix}_chr${i}_*log logs/imputation/$DATE
done



