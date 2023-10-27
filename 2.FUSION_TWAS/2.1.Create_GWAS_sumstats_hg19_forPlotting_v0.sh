#!/bin/bash

# Generate full imputed GWAS sumstats with hg19 coordinates for plotting
# Date: June 29, 2022
# Author: Emilio Ugalde

# Strategy: 
#	- Get hg19 coordinates from LDREF panel. 
#	- Retrieve alleles from LDREF panel, because FUSION aligns sumstats to the panel.

# Background: imputing sumstats was based on hg38 coordinates.
# Output file: "Data/tmp/sumstats_imputed/TECAC_META_8site_validation_LDREF_hg19_aligned.txt"

cd /nfs/GENETEC/TWAS

## INPUT
# Sumstats
dir="Data/tmp/sumstats_imputed"
ss="$dir/TECAC_META_8site_validation.chr"

# LDREF
ref="Methods/FUSION/fusion_twas-master/LDREF/1000G.EUR."


## Parse data

# Compile SS
for i in `seq 22`; do awk 'NR>1 {print $1,$4}' $ss$i.txt | grep -wv NA >> $dir/TECAC_META_8site_validation.txt; done
# Compile REF
for i in `seq 22`; do awk '{print $2,$1,$4,$5,$6}' $ref$i.bim >> $dir/LDREF.txt; done

# Merge to get hg19 positions
Rscript - <<'EOF'
file_out="TECAC_META_8site_validation_LDREF_hg19_aligned.txt"
dir="Data/tmp/sumstats_imputed/"
ss <- read.table(paste0(dir, 'TECAC_META_8site_validation.txt'), header=F)	# 8419022
ref <- read.table(paste0(dir, 'LDREF.txt'), header=F)			# 1190321

names(ss) <- c("SNP", "Z")
names(ref) <- c("SNP", "CHR", "POS", "A1", "A2")

# compute GWAS P-value
ss$P <- 2*pnorm(q=abs(ss$Z), lower.tail=FALSE)

# Merge
out <- merge(ss[,c(1,3)], ref, by="SNP")		# 1161660
# Overlap
nrow(out)/nrow(ref)*100				# 97.59%

# Save
write.table(out, file=paste0(dir, file_out), row.names=F, col.names=T, quote=F, sep="\t" )

EOF


