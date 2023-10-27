#!/bin/bash

# GWAS harmonization: data parsing
# Emilio Ugalde
# March 2-4, 2022
# ssh emiuga@vector.meb.ki.se

# Description: parse GWAS sumstas to the format required by MetaXcan.
# Software: MetaXcan's summary-gwas-imputation (GWAS tools)
# Tutorial: https://github.com/hakyimlab/MetaXcan/wiki/Tutorial:-GTEx-v8-MASH-models-integration-with-a-Coronary-Artery-Disease-GWAS

# Current strategy: extract rsIDs from UCSC matching by position
# Altern. strategy: create var_ID with format: chr22_position_A1_A2_b37; with options --keep_non_rsid --model_db_snp_key varID \
# see: https://lab-notes.hakyimlab.org/post/2021/01/07/spredixcan-harmonization-errors/

cd /nfs/GENETEC/TECAC_GWAS_sumstats

# OUPUT directory
DIR="Data/DerivedData/imputation"
#LOG="logs/gwas_harmonization"

# Temporary file
mkdir -p $DIR/tmp/parsing



## INPUT data

# dbSNP 150hg19 dataset
ucsc="/nfs/GENETEC/TWAS/Data/imputation_sumstats/sample_data/data"
# Files subset on bi-allelic snps 
wc -l $ucsc/ucsc/snp150_hg19/chr*.biallelic | grep total					# 216,786,958
# Use file overlaping TECAC_INFO>0.8
wc -l Data/DerivedData/tmp/rsIDs/overlap_UCSC/chr*.biallelic_TECACinfo08 | grep total		# 11,435,838
  
# Summary statistics
sum="Data/OriginalData/TECAC_META_8site_validation.INFO_08_AF01.txt"
sum_rs="Data/DerivedData/tmp/rsIDs/overlap_UCSC/rsIDs_TECAC_META_8site_validation.info08"

gzip -c -d $sum.gz | head -3
#MarkerName	Allele1	Allele2	Freq1	FreqSE	Effect		StdErr	P-value	Direction	HetISq	HetChiSq	HetDf	HetPVal
#5:85928892	t		c		0.9361	0.0065	0.0018		0.0411	0.9643		---++-+?+	0.0	2.136		7	0.952
#2:170966953	t		c		0.0227	0.0094	-0.3787	0.1947	0.05177	?-??????-	27.3	1.375		1	0.2409
head -3 $sum_rs.chr1
#rs61145527	t	c	0.0112	0.0002	-0.1328	0.2120	0.531	???+?-???	44.0	1.784	1	0.1816
#rs75574240	t	c	0.0175	0.0039	-0.1819	0.2118	0.3906	???--?+??	65.1	5.727	2	0.05708



## DATA parsing

# Check P-value hererogeneity  > 0.001
gzip -c -d $sum.gz | awk '$NF <= 0.001' | wc -l 	# 0
cat $sum_rs.chr* | awk '$NF <= 0.001' | wc -l		# 0

## Extract sumstats per chromosome (for chr:pos SNPs)

## Select following variables: markername (MarkerName: $1); non_effect_allele (Allele1: $2); effect_allele (Allele2: $3); effect_size (Effect: $6); standard_error (StdErr: $7) frequency (Freq1: $4); pvalue (P-value: $8) 
# Autosomal
for i in `seq 22`; do cat <(echo 'chr position Allele1 Allele2 Freq1 Effect StdErr P-value' ) <(gzip -c -d $sum.gz | grep "^$i:" | awk '{print $1, $2, $3, $4, $6, $7, $8}' | sed 's/:/ /g' ) > $DIR/tmp/parsing/tmp_chr$i; done
# X chromosome
cat <(echo 'chr position Allele1 Allele2 Freq1 Effect StdErr P-value' ) <(zcat $sum.gz | grep "^X:" | awk '{print $1, $2, $3, $4, $6, $7, $8}' | sed 's/:/ /g' ) > $DIR/tmp/parsing/tmp_chrX

# Write R script
	# Obtain rsID, chr, and position from UCSC dataset
	# Generate file with rsIDs as SNP "name"
echo '
chr <- commandArgs(trailingOnly = TRUE)
library(data.table)
sum.stat <- paste0("Data/DerivedData/imputation/tmp/parsing/tmp_chr", chr)
sum.rs <- paste0("Data/DerivedData/tmp/rsIDs/overlap_UCSC/rsIDs_TECAC_META_8site_validation.info08.chr", chr)
log <- paste0("Data/DerivedData/imputation/tmp/parsing/parsed_tecac_chr", chr,".log")
ss <- fread(sum.stat)
ss$MarkerName <- paste(ss$chr, ss$position, sep=":")
ss_rs <- fread(sum.rs, select=c(1,2,3,4,6,7,8), col.names=c("name", "Allele1",  "Allele2",  "Freq1" , "Effect", "StdErr", "P-value"))
ss_rs$MarkerName <- ss_rs$name

ucsc <- fread(paste0("Data/DerivedData/tmp/rsIDs/overlap_UCSC/chr", chr,".biallelic_TECACinfo08"), col.names=c("chr", "position", "name") )

# rsID
tmp1 <- merge(ss, ucsc[,2:3], by="position")
tmp2 <- merge(ss_rs, ucsc, by="name")
out <- rbind(tmp2, tmp1)

# change to uppercase
out$Allele2 <- toupper(as.character(out$Allele2))
out$Allele1 <- toupper(as.character(out$Allele1))

# Logs
write.table(paste0("No. susmstats variants: ", nrow(ss)+nrow(ss_rs) ), log, col.names=F, row.names=F, quote=F, sep="\t", append=T)
write.table(paste0("No. UCSC variants: ", nrow(ucsc) ), log, col.names=F, row.names=F, quote=F, sep="\t", append=T)
write.table(paste0("rsIDs extracted: ", nrow(out) ) , log, col.names=F, row.names=F, quote=F, sep="\t", append=T)
write.table(paste0("rsIDs extracted (unique): ", length(unique(out$name)) ), log, col.names=F, row.names=F, quote=F, sep="\t", append=T)
write.table(paste0("Variants with more than one rsID: ", nrow(out[duplicated(out$name),]) ), log, col.names=F, row.names=F, quote=F, sep="\t", append=T)
write.table(paste0("Variants with more than one chr_position: ", nrow(out[duplicated(out$position),]) ), log, col.names=F, row.names=F, quote=F, sep="\t", append=T)
# save
write.table(out, paste0("Data/DerivedData/imputation/tmp/parsing/parsed_tecac_chr", chr, ".txt"), col.names=T, row.names=F, quote=F, sep="\t")
' > Programs/imputation/parse_susmstats_20220302.R

# Parse files
for i in `echo {1..22} X `; do Rscript Programs/imputation/parse_susmstats_20220302.R $i ; done

# Total no. of SNPs
grep -v name Data/DerivedData/imputation/tmp/parsing/parsed_tecac_chr*.txt | wc -l	# 7,174,888

# compress files
for i in `echo {1..22} X `; do gzip Data/DerivedData/imputation/tmp/parsing/parsed_tecac_chr$i.txt; done

# remove intermediate files
rm Data/DerivedData/imputation/tmp/parsing/tmp*
#rm -r $DATA/ucsc/snp150_hg19



## Check number of SNPs extracted
grep 'rsIDs extracted:' Data/DerivedData/imputation/tmp/parsing/parsed_tecac_chr*.log | awk '{sum +=$NF} END {print sum}'					# 7174888
grep 'rsIDs extracted (unique):' Data/DerivedData/imputation/tmp/parsing/parsed_tecac_chr*.log | awk '{sum +=$NF} END {print sum}'				# 7072751

## Check number of duplicated ids
grep 'Variants with more than one rsID' Data/DerivedData/imputation/tmp/parsing/parsed_tecac_chr*.log | awk '{sum +=$NF} END {print sum}'			# 102137
grep 'Variants with more than one chr_position' Data/DerivedData/imputation/tmp/parsing/parsed_tecac_chr*.log | awk '{sum +=$NF} END {print sum}'		# 143772



# Concatenate GWAS 'parsed' files
gwas_parsed="Data/DerivedData/imputation/tmp/parsing/parsed_tecac_"

zcat ${gwas_parsed}chr1.txt.gz | head -1 > ${gwas_parsed}_ALL
for i in `echo {1..22} X`; do zcat ${gwas_parsed}chr$i.txt.gz | grep -v name >> ${gwas_parsed}_ALL; done
gzip ${gwas_parsed}ALL

zcat ${gwas_parsed}ALL.gz | wc -l 		# 7174889

## Remove duplicated SNPs
zcat ${gwas_parsed}ALL.gz | awk '{seen[$1]++; if(seen[$1]==1){ print}}' > ${gwas_parsed}ALL.nondup
gzip ${gwas_parsed}ALL.nondup
zcat ${gwas_parsed}ALL.nondup.gz | wc -l 	# 7072752
# 7174889-7072752 = 102137 removed












