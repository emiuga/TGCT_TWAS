#!/bin/bash

# GWAS imputation: post-processing
# Emilio Ugalde
# March 8, 2022
# ssh emiuga@vector.meb.ki.se

# Description: post-processing of imputed sumstats and compiling.
# Software: MetaXcan's summary-gwas-imputation (GWAS tools) https://github.com/hakyimlab/summary-gwas-imputation/wiki/GWAS-Harmonization-And-Imputation
# Tutorial: https://github.com/hakyimlab/MetaXcan/wiki/Tutorial:-GTEx-v8-MASH-models-integration-with-a-Coronary-Artery-Disease-GWAS

cd /nfs/GENETEC/TECAC_GWAS_sumstats

# INPUT data (Output from GWAS harmonization)
DIR="Data/DerivedData/imputation"
dir_gwas="Data/harmonized_gwas/hg38"
GWAS_prefix="TECAC_META_8site_validation"
#ls $DIR/results | less

# OUPUT directory
DATE="20220307"
mkdir -p $DIR/processed_summary_imputation
mkdir -p logs/imputation/$DATE


# SETUP
# Call virtual environment with installed libraries
source /nfs/GENETEC/TWAS/Methods/summary-gwas-imputation/bin/activate
# Paths
GWAS_TOOLS="/nfs/GENETEC/TWAS/Methods/summary-gwas-imputation/src"

# ANALYSIS
time for i in {1..22}; do
python $GWAS_TOOLS/gwas_summary_imputation_postprocess.py \
-gwas_file $dir_gwas/$GWAS_prefix.chr$i.*.gz \
-folder $DIR/results \
-pattern ${GWAS_prefix}_chr${i}_reg0* \
-parsimony 7 \
-output $DIR/processed_summary_imputation/imputed_$GWAS_prefix.chr$i.txt.gz &> $DIR/processed_summary_imputation/imputed_$GWAS_prefix.chr$i.txt.log
done 
#real	5m41.103s


## Summarize number of imputed SNPs

# Total variants (from log files):
awk 'NR==6' $DIR/processed_summary_imputation/imputed_$GWAS_prefix.chr$i.txt.log  
for i in `seq 22`; do awk 'NR==7 {print $3}' $DIR/processed_summary_imputation/imputed_$GWAS_prefix.chr$i.txt.log ; done | awk '{sum+=$1} END {print sum}'	# 8638344
# Total imputed variants:
grep 'imputed variants' $DIR/processed_summary_imputation/imputed_$GWAS_prefix.chr*.txt.log | awk '{sum+=$4} END {print sum}'					# 2953506
# # Total KEPT "original" variants:
grep Kept $DIR/processed_summary_imputation/imputed_$GWAS_prefix.chr*.txt.log | awk '{sum+=$4} END {print sum}'							# 5684838
# Sum up: 5684838+2953506=8638344

# # Total READ "original" variants (correspond to the harmonized/parsed sumstats) were ambiguous variants were removed:
grep Read $DIR/processed_summary_imputation/imputed_$GWAS_prefix.chr*.txt.log | awk '{sum+=$4} END {print sum}'							# 6623671

############
## concatenate files
#mkdir -p $OUTPUT/merged
#gzip -d -c $OUTPUT/imputed_${GWAS_prefix}.chr1.txt.gz | head -1 > $OUTPUT/merged/${GWAS_prefix}.txt
#for i in `seq 22`; do gzip -d -c $OUTPUT/imputed_${GWAS_prefix}.chr$i.txt.gz | awk 'NR>1' >> $OUTPUT/merged/${GWAS_prefix}.txt; done
#wc -l $OUTPUT/merged/TECAC_META_8site_validation.INFO_08_AF01.txt 	# 
#gzip $OUTPUT/merged/TECAC_META_8site_validation.INFO_08_AF01.txt

