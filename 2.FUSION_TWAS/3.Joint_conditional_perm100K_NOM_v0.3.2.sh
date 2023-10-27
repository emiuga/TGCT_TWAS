#!/bin/bash

# Joint/conditional: Testis+TGCT+sCCA_bestR2 TWAS (gene-permutation significant)
# Date: July 14, 2023
# Author: Emilio Ugalde
# Software: http://gusevlab.org/projects/fusion/

# Description: 
# Note: apply lambda regularization for TWAS loci that failed joint-conditional analysis (e.g. produced NaNs prop. dues to co-linearity issues)
# TWAS scritps: Programs/1.TWAS_v4.sh (including '100K' permutation test)
# Compile scripts: Programs/2.Summarize_TWAS_v4_output_v2.1.R

# Version control: 
# In version 0.3.2: fixed legend overlap with gene boxes.
# In version 0.3.1: improved colors of gene boxex.
# Compared with v0.2, in version 0.3:
# 	- (only to update conditional plots!): new legend features added.
# Compared with v0.1, in version 0.2:
# 	- Plots were updated: horizontal line for soft sig. threshold (i.e. 1/nSNPs) was removed from conditional plots. 


cd /nfs/GENETEC/TWAS/

## Scripts
fusion='Programs/fusion_twas-master/FUSION.post_process_AP-EUM_v1.6.2.R'

## Output directory
dir="Output/imputed_TECAC/Joint_conditional_perm100K_NOM"
mkdir -p $dir/compiled
mkdir -p $dir/cond_plots_v3.1

## Backup prev. analysis
#mv $dir/out $dir/out_backup_v01 

#--------------------------------
## Subset top-genes FDR < 0.01
#--------------------------------
Rscript - <<"EOF"
fname='Output/imputed_TECAC/summary_perm100K/TWAS_testis_tgct_scca_19839uniqGENES'
all <- read.table(paste0(fname, ".dat"), header=T, sep="\t")

# Filter by FDR and get top-genes (165 genes)
all$TWAS.P.fdr <- p.adjust(all$TWAS.P, method="fdr")
out <- all[all$TWAS.P.fdr<0.01, ]

# Filter by PERMUTATION test
out$PERM.PV.fdr <- p.adjust(out$PERM.PV, method="fdr")
out <- out[out$PERM.PV<0.05, ]

write.table(out, paste0(fname, "_",nrow(out),"_permutation_NOM.dat"), col.names=T, row.names=F, quote=F, sep="\t" )
EOF

## INPUT
dir_in="Output/imputed_TECAC/summary_perm100K"
# top genes
file="TWAS_testis_tgct_scca_19839uniqGENES_91_permutation_NOM" 
# full twas
twas='TWAS_testis_tgct_scca_19839uniqGENES.dat'
# LDREF path
LDREF='Methods/FUSION/fusion_twas-master/LDREF/1000G.EUR.'

# using imputed sumstas in hg38 version
INPUT="Data/tmp/sumstats_imputed"
GWAS_prefix="TECAC_META_8site_validation"

## Settings
locus_win=100000
min_r2=0.05
max_r2=0.9

lambda=0.15 



#--------------------------------
# Run analysis 
#--------------------------------

dirOut=$dir/out_v3.1
mkdir -p $dirOut

time for i in `seq 22`; do
Rscript $fusion \
--sumstats $INPUT/$GWAS_prefix.chr$i.txt \
--input $dir_in/$file.dat  \
--full_twas $dir_in/$twas  \
--out $dirOut/chr$i.$file \
--ref_ld_chr $LDREF \
--lambda $lambda \
--max_r2 $max_r2 \
--min_r2 $min_r2 \
--chr $i \
--report --save_loci \
--plot --locus_win $locus_win --gene_legend --gwas_line \
--verbose 2 &> $dirOut/chr$i.$file.log 
done
	# real	5m9.503s

#------------
### assessment
#------------
Rscript Programs/Assesment_conditional_v0.1.R --dir $dirOut --dirOut $dir/compiled --name assessment_conditional.dat

### summarize
print=$dir/compiled/assessment_conditional.dat

# number of loci
Nloc=`grep -v LOCUS $print | wc -l ` 	# 51

# Inflated SNPs: those with conditional p-value lower than the original.
# average number of inflated SNPs (overall)
awk -v nloc=$Nloc '{sum += $8} END {print sum/nloc}' $print	# 
# average number of inflated SNPs (above average Pval)
awk -v nloc=$Nloc '{sum += $9} END {print sum/nloc}' $print 	# 
# average number of inflated SNPs (above SNP-adj. threhold)
awk -v nloc=$Nloc '{sum += $10} END {print sum/nloc}' $print 	# 

#------------
# Compile: 
#------------
# copy plots
cp $dirOut/*loc*.pdf $dir/cond_plots_v3.1

# COMPILE REPORT (loci)
nam_report="Report_${Nloc}_loci"
list=`find $dirOut -name *.report | sort -V`
cat $dirOut/*.report  | awk 'NR==1 {$1=""; print $0}' > $dir/compiled/$nam_report.dat
for i in $list; do cat $i | awk '!/FILE/ {$1=""; print $0}' >> $dir/compiled/$nam_report.dat; done


## GENE-level results (included genes): TWAS output merged with Joint-status and TOP.SNP.COR
nam_file="Joint_$file.dat"
list=`find $dirOut -name *genes | sort -V`
cat $dirOut/*genes | head -1 > $dir/compiled/$nam_file
for i in $list; do cat $i | awk '!/PANEL/'  >> $dir/compiled/$nam_file; done
grep -v PANEL $dir/compiled/$nam_file | wc -l								# 91

# Tab. by join-status (TRUE==prioritized by join/cond. analysis)
awk 'NR>1 {print $30}' $dir/compiled/$nam_file | sort | uniq -c
	#    34 FALSE
	#    57 TRUE
     
# Tab. by join-status (TRUE==prioritized by join/cond. analysis
awk 'NR>1 {print $31}' $dir/compiled/$nam_file | sort | uniq -c
	#     28 conditional
	#      6 correlated
	#     26 joint_best
	#      6 joint_other
	#     25 single


## 2023-01-31
## Compile joint-conditional estimates (at per-locus level)

## For genes found as conditionally-independent.
joint_cond_out=$dir/compiled/LOCI_Joint_conditional_independent_genes.dat
# header
cat `find $dirOut -name *LOCI_joint_included_*.dat | head -1` | awk 'NR==1 {print $2,$5,$6,$7,$8}' > $joint_cond_out
# Populate file
cat `find $dirOut -name *.LOCI_joint_included_*.dat | sort -V` | grep -v 'FILE' | awk '{print $2,$5,$6,$7,$8}'  >> $joint_cond_out

## For conditional genes.
joint_drop_out=$dir/compiled/LOCI_Joint_conditional_dropped_genes.dat
# header
cat `find $dirOut -name *LOCI_joint_dropped_*.dat | head -1` | awk 'NR==1 {print $2,$5,$6,$7,$8}' > $joint_drop_out
# Populate file
cat `find $dirOut -name *.LOCI_joint_dropped_*.dat | sort -V` | grep -v 'FILE' | awk '{print $2,$5,$6,$7,$8}'  >> $joint_drop_out

## Retrieve the Z-score threhold used
## NOTE: zthresh = qnorm( 0.05 / nrow( wgtlist ) / 2,lower.tail=F)
## It was set as nominally significance adjusted for the number of gene in the locus

cat $dirOut/*log | grep "a weight must have Z^2" | head -2
	#1 weights considered, a weight must have Z^2 > 3.84 to be retained in the model >> qnorm( 0.05 / 1 / 2,lower.tail=F) == 1.96
	#6 weights considered, a weight must have Z^2 > 6.96 to be retained in the model >> qnorm( 0.05 / 6 / 2,lower.tail=F) == 2.64








