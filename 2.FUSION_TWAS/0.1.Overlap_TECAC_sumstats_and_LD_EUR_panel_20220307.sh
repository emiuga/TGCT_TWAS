#!/bin/bash

# Overlap between TECAC sumstats and LD reference panel
# Author: Emilio Ugalde <emilio.ugalde.morales@ki.se>
# Date: March 7, 2022
# ssh emiuga@vector.meb.ki.se

# Description:
# Overlap was of 92.95% (INFO>0.8) and 91.16% (INFO>=0.9)
# Overlap using INFO>0.8-IMPUTED sumstats: 97.59%

cd /nfs/GENETEC/TWAS

## Set output file
mkdir -p Data/tmp/LD
out=overlap_sumstats_1000G.LD.EUR.chr

## INPUT
# LD reference panel
fusion=Methods/FUSION/fusion_twas-master
# ls $fusion/LDREF/1000G.EUR.{1..22}.bim
# wc -l $fusion/LDREF/1000G.EUR.{1..22}.bim
# 1190321 total

#----------------
# TECAC susmstats 
#----------------
# Diltered by AF>0.1, heter>0.001, and INFO>0.8

# deprecated: sum=/nfs/GENETEC/TECAC_GWAS/Data/DerivedData/TECAC_META_8site_validation.INFO_08.txt
GWAS_prefix="TECAC_META_8site_validation.INFO_08_AF01.parsed.munged"
## munged:
sum="/nfs/GENETEC/TECAC_GWAS_sumstats/Data/DerivedData/$GWAS_prefix.sumstats.gz"
## before munge:
sum_p="/nfs/GENETEC/TECAC_GWAS_sumstats/Data/DerivedData/imputation/tmp/parsing/parsed_tecac_ALL.nondup.gz"


## Overlap

# munged:
comm -12 <( zcat $sum | awk 'NR>1{print $1}' | sort) <( awk '{print $2}' $fusion/LDREF/1000G.EUR.*.bim | sort) | wc -l	# 1105156
# before munge:
comm -12 <( zcat $sum_p | awk 'NR>1{print $1}' | sort) <( awk '{print $2}' $fusion/LDREF/1000G.EUR.*.bim | sort) | wc -l	# 1105186


#----------------
# TECAC susmstats: harmonized (GWAS_tools/PrediXcan)
#----------------
prefix="/nfs/GENETEC/TECAC_GWAS_sumstats/Data/harmonized_gwas/hg38/TECAC_META_8site_validation"
ls $prefix.chr{1..22}.parsed_hg38.gz | wc -l # 22
zcat $prefix.chr{1..22}.parsed_hg38.gz | grep -v variant_id  | wc -l		# 6623671

## Overlap
comm -12 <( for i in `seq 22`; do zcat $prefix.chr$i.parsed_hg38.gz | awk 'NR>1{print $1}' ; done | sort) <( awk '{print $2}' $fusion/LDREF/1000G.EUR.*.bim | sort) | wc -l # 1091505


#----------------
# TECAC susmstats: imputed (GWAS_tools/PrediXcan)
#----------------
prefix="Data/tmp/sumstats_imputed/TECAC_META_8site_validation"
#ls $prefix.chr{1..22}.parsed_hg38.gz | wc -l # 22
i=22
grep -v SNP $prefix.chr$i.txt  | wc -l					# 119072
awk '{print $2}' $fusion/LDREF/1000G.EUR.$i.bim | sort | uniq | wc -l 	# 17489
## Overlap
comm -12 <(awk 'NR>1{print $1}' $prefix.chr$i.txt | sort | uniq ) <( awk '{print $2}' $fusion/LDREF/1000G.EUR.$i.bim | sort ) | wc -l 	# 16943


prefix="/nfs/GENETEC/TECAC_GWAS_sumstats/Data/DerivedData/imputation/processed_summary_imputation"
files="$prefix/imputed_TECAC_META_8site_validation"
ls $files.chr{1..22}.txt.gz | wc -l # 22
zcat $files.chr{1..22}.txt.gz | grep -v variant_id  | wc -l			# 8638344

## Overlap
comm -12 <( for i in `seq 22`; do zcat $files.chr$i.txt.gz  | awk 'NR>1{print $1}' ; done | sort) <( awk '{print $2}' $fusion/LDREF/1000G.EUR.*.bim | sort) | wc -l 	# 1161657








############################## DO NOT RUN scripts below!
#----------------

#for i in `seq 22`; do 
#	comm -12 <( zcat | awk 'FNR>1{print $1}' $sum | sort) <( awk '{print $2}' $fusion/LDREF/1000G.EUR.*.bim | sort) > Data/tmp/LD/$out.$i ; 
#done

# wc -l Data/tmp/LD/*
# 1106411 total

#----------------
# Overlap with TGCT loci
#----------------
gwas=/nfs/GENETEC/SWENOTECA/PRS/Documents/SNPs/list_78_gwas_SNPs_testis_IDs.txt

grep -wFf $gwas $sum_p | wc -l
#grep -wFf $gwas $sum | wc -l
# 78

#---------------------
# INFO>=0.9
#---------------------
# Set output file
out=info09_overlap_sumstats_1000G.LD.EUR.chr

# TECAC susmstats filtered by AF>0.1, heter>0.001, and INFO>0.8
sum=/nfs/GENETEC/TECAC_GWAS/Data/DerivedData/TECAC_META_8site_validation.INFO_09.txt

for i in `seq 22`; do 
	comm -12 <( awk 'FNR>1{print $1}' $sum | sort -g | grep ^$i: | cut -d: -f2 | sort) <( awk '{print $4}' $fusion/LDREF/1000G.EUR.$i.bim | sort) > Data/tmp/LD/$out.$i ; 
done

wc -l Data/tmp/LD/info09_overlap_*
#1008550 total

#-----
# Overlap with TGCT loci
#-----
grep -wFf $gwas $sum | wc -l
# 75

grep -wFf $gwas $sum | awk '{print $1}' > tmp1
grep -vwFf tmp1 $gwas > tmp2
16:15521872
17:76691564
20:50711198

## look into INFO scores
gzip -d -c /nfs/GENETEC/TECAC_GWAS/Data/OriginalData/TECAC_info.txt.gz | grep -wFf tmp2
16:15521872 0.81458
17:76691564 0.89655
20:50711198 0.84413















