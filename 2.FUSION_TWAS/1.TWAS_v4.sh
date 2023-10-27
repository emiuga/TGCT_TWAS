#!/bin/bash

# TWAS analysis 
# Date: October 20, 2022
# Author: Emilio Ugalde
# Software: http://gusevlab.org/projects/fusion/
# Description: TWAS (FUSION analysis) using filtered and formatted TECAC sumstats
# Version control: compared with v1.0, v2.0:
# 	- Latest prediction models based on GTEx_v8 (for Testis tissue)
# 	- Add permutation test to compute "an empirical association statistics conditional on the GWAS effect at the locus."
# 	- Use imputed summary GWAS (imputation of Z-score based on summary data).  

# *** Performing the expression imputation (TWAS) ***

#Under the hood, the analysis steps are: 
#	(1) unify the GWAS and reference SNPs and remove/flip alleles as appropriate; 
#	(2) impute GWAS Z-scores for any reference SNPs that were missing using the IMPG algorithm; 
#	(3) estimate the functional-GWAS association statistic; 
#	(4) report results over all features tested.

# Modification:
#	1) 100K iterations for gene-permutation test (using fixed random seed)

# To run:
# cd /nfs/GENETEC/TWAS/
# sh Programs/1.TWAS_v4.sh &> Programs/1.TWAS_v4.log
	
## Setup	
cd /nfs/GENETEC/TWAS/
fusion='Methods/FUSION/fusion_twas-master' # path to pre-computed weights and LDREF panel.
script='Programs/fusion_twas-master/FUSION.assoc_test_ModifiedPerm_v0.1.R'

## INPUT data
# TECAC susmtats (AF>0.01, INFO>0.8 filtered; snpIDs-parsed, munged, and non-duplicated)
GWAS_prefix="TECAC_META_8site_validation.INFO_08_AF01.parsed.munged"
sum="/nfs/GENETEC/TECAC_GWAS_sumstats/Data/DerivedData/$GWAS_prefix.sumstats.gz"

# using imputed sumstas in hg38 version
INPUT="Data/tmp/sumstats_imputed"
GWAS_prefix="TECAC_META_8site_validation"

# models
models='Methods/FUSION/fusion_twas-master/WEIGHTS/'

# OUTPUT
DATE="20221020"
DIR="Output/imputed_TECAC/perm_100K"
mkdir -p $DIR/compiled

# No. of permutations (Significance conditional on high GWAS effects)
# Test is perfomed on genes below a "Minimum p-value" (e.g "for which to initiate permutation test"), which defatult is set to 0.05
perm=100000
seed=12345

# To Look at results
# NOTE: column 20 == TWAS.P; column 21 == permuted P-value; column 21 = no. of permutations; column 23 = PERM.ANL_PV

#-----------------------------------------
# Normal TESTIS tissue: GTEx_v8 
#-----------------------------------------
OUT=Testis_GTEx_v8
LOG=$DIR/$OUT/logs
mkdir -p $DIR/$OUT/logs

time parallel -j4 \
Rscript $script  \
--sumstats $INPUT/$GWAS_prefix.chr{}.txt \
--weights $fusion/WEIGHTS/GTEx_v8/Testis/GTExv8.ALL.Testis.pos \
--weights_dir $fusion/WEIGHTS/GTEx_v8/Testis \
--ref_ld_chr $fusion/LDREF/1000G.EUR. \
--chr {} \
--GWASN 189839 \
--perm $perm \
--rngseed $seed \
--out $DIR/$OUT/TGCT.Testis_GTEx_v8.chr{}.dat &>> $LOG/TGCT.Testis_GTEx_v8.$DATE.log ::: {1..22}


## COMPILE results
awk 'FNR==1{print}' $DIR/$OUT/TGCT.$OUT.chr1.dat > $DIR/compiled/TWAS_TECAC.$OUT.dat
cat $DIR/$OUT/TGCT.$OUT* | grep -v PANEL | awk '$20!="NA" ' >> $DIR/compiled/TWAS_TECAC.$OUT.dat


#-----------------------------------------
# cannonical correlation cross-tissue weights (sCCA)
#-----------------------------------------
OUT=sCCA
LOG=$DIR/$OUT/logs
mkdir -p $DIR/$OUT/logs

time for i in {1..3}; do
parallel -j4 \
Rscript $script  \
--sumstats $INPUT/$GWAS_prefix.chr{}.txt \
--weights $fusion/WEIGHTS/sCCA_weights_v8/sCCA$i.pos  \
--weights_dir $fusion/WEIGHTS/sCCA_weights_v8 \
--ref_ld_chr $fusion/LDREF/1000G.EUR. \
--chr {} \
--GWASN 189839 \
--perm $perm \
--rngseed $seed \
--out $DIR/$OUT/TGCT.sCCA$i.chr{}.dat &>> $LOG/TGCT.sCCA.$DATE.log ::: {1..22}
done


## COMPILE results (per sCCA feature)
for i in {1..3}; do awk 'FNR==1{print}' $DIR/$OUT/TGCT.$OUT$i.chr1.dat > $DIR/compiled/TWAS_TECAC.$OUT$i.dat; done
for i in {1..3}; do cat $DIR/$OUT/TGCT.$OUT$i* | grep -v PANEL | awk '$20!="NA" ' >> $DIR/compiled/TWAS_TECAC.$OUT$i.dat; done
for i in {1..3}; do sed -i "s/^NA/sCCA$i/g" $DIR/compiled/TWAS_TECAC.$OUT$i.dat; done


#-----------------------------------------
### TGCA_Testis tumors (TGCT)  
#-----------------------------------------
OUT=TCGA_TGCT
LOG=$DIR/$OUT/logs
mkdir -p $DIR/$OUT/logs

time parallel -j4 \
Rscript $script  \
--sumstats $INPUT/$GWAS_prefix.chr{}.txt \
--weights $fusion/WEIGHTS/TCGA/TCGA-TGCT.TUMOR.pos  \
--weights_dir $fusion/WEIGHTS/TCGA \
--ref_ld_chr $fusion/LDREF/1000G.EUR. \
--chr {} \
--GWASN 189839 \
--perm $perm \
--rngseed $seed \
--out $DIR/$OUT/TGCT.TCGA_tgct.chr{}.dat &>> $LOG/TGCT.TCGA_tgct.$DATE.log ::: {1..22}


## COMPILE results
awk 'FNR==1{print}' $DIR/$OUT/TGCT.TCGA_tgct.chr1.dat > $DIR/compiled/TWAS_TECAC.$OUT.dat
cat $DIR/$OUT/TGCT.TCGA_tgct* | grep -v PANEL | awk '$20!="NA" ' >> $DIR/compiled/TWAS_TECAC.$OUT.dat









