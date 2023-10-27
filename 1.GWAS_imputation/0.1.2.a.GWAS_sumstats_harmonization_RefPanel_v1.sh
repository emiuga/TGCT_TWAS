#!/bin/bash

# Reference panel needed for GWAS harmonization (parsing)
# Emilio Ugalde
# Mar 02, 2022
# ssh emiuga@vector.meb.ki.se

# Description: Parse the snp_reference_metadata file based on HRC GRch37
# GWAS harmonization: MetaXcan's summary-gwas-imputation (GWAS tools) https://github.com/hakyimlab/summary-gwas-imputation/wiki/GWAS-Harmonization-And-Imputation

cd /nfs/GENETEC/TECAC_GWAS_sumstats

mkdir -p Data/reference_panel
hrc=/nfs/GENETEC/GENETEC_GWAS_HRC/Data/OriginalData/HRC.r1-1.GRCh37.wgs.mac5.sites.tab

Rscript - << "EOF"
library(data.table)
hrc <- fread('/nfs/GENETEC/GENETEC_GWAS_HRC/Data/OriginalData/HRC.r1-1.GRCh37.wgs.mac5.sites.tab')
pars_colmns <- c("chromosome", "position", "id", "allele_0", "allele_1", "allele_1_frequency", "rsid") 

out <- hrc[, c(1,2,4,5,8,3)]
out$id <- paste(out$'#CHROM', out$POS, sep=":")
out <- out[, c(1,2,7,3,4,5,6)]
colnames(out) <- pars_colmns
write.table(out, 'Data/reference_panel/HRC.r1-1.GRCh37.wgs.mac5.sites.metadata', col.names=T, row.names=F, quote=F, sep="\t")
EOF
# compress file
gzip Data/reference_panel/HRC.r1-1.GRCh37.wgs.mac5.sites.metadata


