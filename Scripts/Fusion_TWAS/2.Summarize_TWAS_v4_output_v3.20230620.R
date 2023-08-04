
# Summaryize TWAS output (R script)
# Date: Jun 19, 2023
# Author: Emilio Ugalde
# TWAS software: http://gusevlab.org/projects/fusion/
# Description: Summarize TWAS_v4 results (using imputed sumstats) including 100K permutation test.
# Strategy: start with Testis TWAS, the concatenate subsecuent "new" genes added by TGCT and sCCA with best predictive performance, in that order. 

# NOTE: here hg19_GTExV8 gene annotations are used to extract chr_position (NOTE: see Data/GTEx/README_gtexGene.txt regarding gene annotation data (variable's codebook).)

# VERSION control: 
# - In version (v4), genome coordinates and cytoband for hg19 and hg38 were added. 
# - In version (v3), including permutation test using 100K iterations. 
# - In version (v2), models are filtered by cross-validation R2>0.01 and P-value<0.05.

## SETUP
# cd /nfs/GENETEC/TWAS # working directory
# mv Output/imputed_TECAC/summary_perm100K Output/imputed_TECAC/summary_perm100K_deprecated_v3 # deprecate previous oputput
# mkdir -p Output/imputed_TECAC/summary_perm100K

# setwd('/nfs/GENETEC/TWAS/')
suppressMessages(library('dplyr'))
dir <- "Output/imputed_TECAC/perm_100K/compiled/"

# Set OUPUT
dir_out <- "Output/imputed_TECAC/summary_perm100K/"
file_out <- "TWAS_testis_tgct_scca"

# INPUT

## TWAS results
testis <- read.table(paste0(dir, 'TWAS_TECAC.Testis_GTEx_v8.dat'), header=T, sep="\t") 		# 12683
tgct <- read.table(paste0(dir, 'TWAS_TECAC.TCGA_TGCT.dat'), header=T, sep="\t")			# 1291
scca1 <- read.table(paste0(dir, 'TWAS_TECAC.sCCA1.dat'), header=T, sep="\t")			# 13249
scca2 <- read.table(paste0(dir, 'TWAS_TECAC.sCCA2.dat'), header=T, sep="\t")			# 12522
scca3 <- read.table(paste0(dir, 'TWAS_TECAC.sCCA3.dat'), header=T, sep="\t")			# 12032

## FIlter significant models
testis <- testis[testis$MODELCV.R2>0.01 & testis$MODELCV.PV<0.05, ]	# 12206
tgct <- tgct[tgct$MODELCV.R2>0.01 & tgct$MODELCV.PV<0.05, ]		# 1254
scca1 <- scca1[scca1$MODELCV.R2>0.01 & scca1$MODELCV.PV<0.05, ]		# 13009
scca2 <- scca2[scca2$MODELCV.R2>0.01 & scca2$MODELCV.PV<0.05, ]		# 11981
scca3 <- scca3[scca3$MODELCV.R2>0.01 & scca3$MODELCV.PV<0.05, ]		# 11320


## Add Gene Symbols
#--------------
# Format TGCT ids (see scripts: Programs/2.0.Harmonize*sh )
#--------------
fix <- read.table('Data/GTEx/twasTGCTgenes_NOTin_hg19.gtexGeneV8.fixIDs', header=T, sep="\t")
x <- merge(tgct, fix,by="ID")
x$ID <- x$symbol
tgct <- rbind(tgct[!tgct$ID%in%fix$ID,], x[, names(x)%in%names(tgct)])	# 1257 rows (note: 3 genes had more than one symbol)
rm(x, fix)

#--------------
## Harmonize Gene IDs and coordinates (To convert testis and sCCA, which are built in GTEx_v8 hg38, to gene names and  hg19  positions): 
#--------------
ids <- read.table('Data/GTEx/ids.txt', header=T)				# 56200
colnames(ids) <- c("ID", "GENE")

#-------
# Annotations  (hg38)
#-------
# https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/gtexGeneV8.txt.gz

hg38 <- read.table('Data/GTEx/gtexGeneV8.hg38_annot.txt', header=F)
names(hg38) <- c("chr", "P0_hg38", "P1_hg38", "GENE", "ID", "geneType") 

annot <- merge(ids, hg38[, c("ID", "chr", "P0_hg38", "P1_hg38", "geneType")], by="ID" )
annot$ENSEMBL <- annot$ID

#-------
# Annotations  (hg19)
#-------
# Merge using ENSEMBL id
hg19 <- read.table('Data/GTEx/GTExV8_hg19_IDensembl')				# 55550 (all in ids file)
names(hg19)[c(1,4,5)] <- c("ID", "P0_hg19", "P1_hg19")

# merge(ids, ens[, c(1,4,5)], by.x="ID", by.y="V1", all.x=T ) 

annot <- merge(annot, hg19[, c("ID", "P0_hg19", "P1_hg19")], by="ID", all.x=T )


#-------
# Harmonize symbols to match to TGCT
#-------
genes <- read.table('Data/GTEx/GTExV8_hg19_geneNames')				# 54481
dup_genes <- unique(genes$V4[duplicated(genes$V4)])		# 178
dup_geneIDs <- unique(ids$GENE[duplicated(ids$GENE)]) 	# 199
genes_nodup <- genes[!genes$V4%in%dup_genes, ]		# 53021
cor_genes <- merge(genes_nodup[, c(4,1,2,3)], by.x="V4", ids[!ids$GENE%in%dup_geneIDs,], by.y="GENE") 	# 52907
colnames(cor_genes) <- c("ID", "chr", "start", "end", "ENSEMBL")
# Add anontations (merge using ENSEMBL id)
#cor_genes <- merge(cor_genes, annot[, c(4,8)], by.x="ID", by.y="name" )

annot_2 <- annot[annot$ID%in%cor_genes$ENSEMBL, ]	# [1] 52907

#--------------------
# Format TWAS files (add genome coordinates and annotation)
#-------------------- 
# NOTE: use hg19 are referene coordinates!!!
# TGCT panel: use the coordinates from the models, which are in hg19)
# Testis and sCCA (which are built in GTEx_v8 hg38 start site)

## NOTE, deprecated: replace orignal coordinate for the hg19 version, if not missing in "new" file. 
## Example:
##testis$P0 <- ifelse(is.na(testis$start)==T, testis$P0, testis$start)	
##testis$P1 <- ifelse(is.na(testis$end)==T, testis$P1, testis$end)

# SET columns to be retrieved
# col_nam <- append(colnames(testis)[1:6], append(c( "P0_hg19", "P1_hg19", "P0_hg38", "P1_hg38", "geneType"), colnames(testis)[7:ncol(testis)]) )

# TGCT
tgct <- merge(tgct, annot_2[, c("GENE", "P0_hg38", "P1_hg38", "geneType", "ENSEMBL")], by.x="ID", by.y="GENE", all.x=T) 	
as.numeric(table(is.na(tgct$P0_hg38) )[2])	# Missing: 69
tgct$P0_hg19 <- tgct$P0
tgct$P1_hg19 <- tgct$P1
tgct$GENE <- tgct$ID

# TGCT (add ENSEMBL ids and geneType)
#tgct <- merge(tgct, cor_genes, by="ID", all.x=T)
#as.numeric(table(is.na(tgct$start) )[2])	# Missing: 69
#tgct$GENE <- tgct$ID

# Testis
testis <- merge(testis, annot[, c("ID", "P0_hg19", "P1_hg19", "P0_hg38", "P1_hg38", "GENE","geneType")], by="ID", all.x=T)
as.numeric(table(is.na(testis$P0_hg19) )[2])	# Missing: 145
testis$P0 <- testis$P0_hg19
testis$P1 <- testis$P1_hg19

# scca1 
scca1 <- merge(scca1, annot[, c("ID", "P0_hg19", "P1_hg19", "P0_hg38", "P1_hg38", "GENE","geneType")], by="ID", all.x=T)
as.numeric(table(is.na(scca1$P0_hg19) )[2])	# Missing: 168
scca1$P0 <- scca1$P0_hg19
scca1$P1 <- scca1$P1_hg19

# scca2
scca2 <- merge(scca2, annot[, c("ID", "P0_hg19", "P1_hg19", "P0_hg38", "P1_hg38", "GENE","geneType")], by="ID", all.x=T)
as.numeric(table(is.na(scca2$P0_hg19) )[2])	# Missing: 153
scca2$P0 <- scca2$P0_hg19
scca2$P1 <- scca2$P1_hg19

# scca1 
scca3 <- merge(scca3, annot[, c("ID", "P0_hg19", "P1_hg19", "P0_hg38", "P1_hg38", "GENE","geneType")], by="ID", all.x=T)
as.numeric(table(is.na(scca3$P0_hg19) )[2])	# Missing: 142
scca3$P0 <- scca3$P0_hg19
scca3$P1 <- scca3$P1_hg19




#--------------------
# Select sCCA feature with best performance. NOTE: use R2!!
#--------------------

#---------
## model CV performance (R2)
#---------
cols <- c("ID","MODELCV.R2")
tmp_r2 <- scca1[, names(scca1)%in%c("GENE", "ID","MODELCV.R2")]
tmp_r2 <- merge(tmp_r2, scca2[, names(scca2)%in%cols], by="ID", all=T)
tmp_r2 <- merge(tmp_r2, scca3[, names(scca3)%in%cols], by="ID", all=T)
tmp_r2 <- tmp_r2[, c(1,3,2,4,5)]
names(tmp_r2) <- c("ID","GENE","scca1", "scca2", "scca3")
# object with NAs == 0
tmp_r2_na <- tmp_r2
tmp_r2_na[, c(3:5)][is.na(tmp_r2_na[, c(3:5)])] <- 0

# Identify "best" panel
tmp_r2$best <- unlist(apply(tmp_r2[,3:5], 1, function(x){names(which.max(x))} ))
#table(tmp_r2$best)
#scca1 scca2 scca3 
#11107  1361   793

scca_best_r2 <- rbind(scca1[scca1$ID%in%tmp_r2[tmp_r2$best=="scca1",]$ID, ], scca2[scca2$ID%in%tmp_r2[tmp_r2$best=="scca2",]$ID, ], scca3[scca3$ID%in%tmp_r2[tmp_r2$best=="scca3",]$ID, ])

#---------
## model CV performance (P-value)
#---------
cols <- c("ID","MODELCV.PV")
tmp_pv <- scca1[, names(scca1)%in%cols]
tmp_pv <- merge(tmp_pv, scca2[, names(scca2)%in%cols], by="ID", all=T)
tmp_pv <- merge(tmp_pv, scca3[, names(scca3)%in%cols], by="ID", all=T)
names(tmp_pv) <- c("ID","scca1", "scca2", "scca3")
# object with NAs == 1
tmp_pv_na <- tmp_pv
tmp_pv_na[, c(2:4)][is.na(tmp_pv_na[, c(2:4)])] <- 1

# Identify "best" panel
tmp_pv$best <- unlist(apply(tmp_pv[,2:4], 1, function(x){names(which.min(x))} ))
#table(tmp_pv$best)
#scca1 scca2 scca3 
#10930  1435   896

scca_best_pv <- rbind(scca1[scca1$ID%in%tmp_pv[tmp_pv$best=="scca1",]$ID, ], scca2[scca2$ID%in%tmp_pv[tmp_pv$best=="scca2",]$ID, ], scca3[scca3$ID%in%tmp_pv[tmp_pv$best=="scca3",]$ID, ])


#tmp_md <- merge(tmp_r2, tmp_pv, by="ID" )


#--------------------
# Prepare files to be Concatenated
#--------------------
col_nam <- append(names(testis), c("ENSEMBL"))

testis_out <- testis
testis_out$ENSEMBL <- testis_out$ID
testis_out$ID <- testis_out$GENE

scca_best <- scca_best_r2
scca_best$ENSEMBL <- scca_best$ID
scca_best$ID <- scca_best$GENE


#--------------------
# File with unique genes (testis + tgct + scca)
#--------------------

## Testis + tgct
length(unique(tgct$ID[!tgct$ID%in%testis_out$ID]) )											# 414
comb <- rbind(testis_out[, names(testis_out)%in%col_nam], tgct[!tgct$ID%in%testis_out$ID, names(tgct)%in%col_nam] )	# 12206+414=12620

## + best sCCA (ensembl) 
length(scca_best$ENSEMBL[!scca_best$ENSEMBL%in%comb$ENSEMBL]) 			# 7238 using ensembl id
length(scca_best$GENE[!scca_best$GENE%in%comb$GENE])				# 7230 using gene symbols
comb2 <- rbind(comb, scca_best[!scca_best$ENSEMBL%in%comb$ENSEMBL, names(scca_best)%in%col_nam ])			# 12620+7238=19858

## + best sCCA (gene symbol)
length(scca_best$GENE[!scca_best$GENE%in%comb2$GENE]) 			# 0

table(comb2$PANEL)
#GTExv8.ALL.Testis             sCCA1             sCCA2             sCCA3 	TCGA-TGCT.TUMOR
# 	     12206              5740               945               553   		    412



#--------------------
# Duplicated IDs
#--------------------

dup_genes <- unique(comb2$ID[duplicated(comb2$ID)]) # 16
# chr [1:16] "Y_RNA" "LINC01422" "DNAJC9-AS1" "LINC01481" "LINC01238" ...

#----
## Remove duplicated genes
all <- comb2
## Of duplicated genes, exclude transcript with higher P-value
dup <- unique(all$ID[duplicated(all$ID)])
# Get rownanes
dup_exclude <- NULL
for(i in c(1:length(dup)) ){
row_min_p <- rownames(all[all$ID==dup[i], ][which.min(all[all$ID==dup[i], "TWAS.P"]), ])
dup_exclude <- append(dup_exclude, rownames(all[all$ID==dup[i], ])[!rownames(all[all$ID==dup[i], ])%in%row_min_p ] ) 
}
# Remove rows
comb2 <- all[!rownames(all)%in%dup_exclude, ]
#----

table(comb2$PANEL)
#GTExv8.ALL.Testis             sCCA1             sCCA2             sCCA3 	TCGA-TGCT.TUMOR 
#            12192              5738               944               553		    412


#--------------------
# Save per-panel files
#--------------------
# Duplicated genes to remove
dup_exclude_file <- all[rownames(all)%in%dup_exclude, ]$FILE

testis_out <- testis_out[!testis_out$FILE%in%dup_exclude_file ,]
tgct_out <- tgct[!tgct$FILE%in%dup_exclude_file ,]
scca_best_out <- scca_best[!scca_best$FILE%in%dup_exclude_file ,]

write.table(testis_out, paste0(dir_out, 'TWAS_TECAC.Testis_GTEx_v8_formatted.dat'), col.names=T, row.names=F, quote=F, sep="\t") 
write.table(tgct_out, paste0(dir_out, 'TWAS_TECAC.TCGA_TGCT_formatted.dat'), col.names=T, row.names=F, quote=F, sep="\t")
write.table(scca_best_out, paste0(dir_out, 'TWAS_TECAC.best_sCCA_formatted.dat'), col.names=T, row.names=F, quote=F, sep="\t")	

#---------
# Get cytobands (hg19)
#---------
bands_hg19 <- read.delim("Data/Human_karyotype/cytoBand_hg19.txt", header=F)
names(bands_hg19)[2:3] <- c("start", "end")
bands_hg19$CHR <- sub("chr", "", bands_hg19$V1)
bands_hg19$Cytoband_hg19 <- paste0(bands_hg19$CHR, bands_hg19$V4)
## Overlap, as a fraction of gene lenght
md <- merge(comb2[ , c("ID", "CHR","P0", "P1")], bands_hg19, by="CHR" )
# Gene length
md$gene_l <- md$P1-md$P0
# Consider 3 scenarios: 1)full overlap, 2)overlap at-the-end of band, 3) overlap at-the-start of band.
md$overlap <- ifelse(md$P0>=md$start & md$P1<=md$end, md$gene_l/md$gene_l,
                     ifelse(md$P1>md$end & md$P0<md$end, (md$end-md$P0)/md$gene_l, 
                            ifelse(md$P1>md$start & md$P0<md$start, (md$P1-md$start)/md$gene_l, NA ) ) )
## Handle multiple gene-band overlaps
md_ov <- md[is.na(md$overlap)==F, ]
dup <- md_ov$ID[duplicated(md_ov$ID)]
# Select best overlap for genes with multiple instances
# Get rownanes
dup_exclude <- NULL
# keep the one with max overlap
for(i in c(1:length(dup)) ){
row_max <- rownames(md_ov[md_ov$ID==dup[i], ][which.max(md_ov[md_ov$ID==dup[i], "overlap"]), ])
dup_exclude <- append(dup_exclude, rownames(md_ov[md_ov$ID==dup[i], ])[!rownames(md_ov[md_ov$ID==dup[i], ])%in%row_max ] ) 
}
# Remove duplicated rows with poor overlap
md_ov_clean <- md_ov[!rownames(md_ov)%in%dup_exclude, ]
str(factor(md_ov_clean$ID))	# Factor w/ 19601 levels "A1BG","A1BG-AS1",..:
# Add cytoband to file
comb2 <- merge(comb2, md_ov_clean[, c("ID", "Cytoband_hg19")], by="ID", all.x=T)


#---------
# Get cytobands (hg38)
#---------
bands_hg38 <- read.delim("Data/Human_karyotype/cytoBand_hg38.txt", header=F)
names(bands_hg38)[2:3] <- c("start", "end")
bands_hg38$CHR <- sub("chr", "", bands_hg38$V1)
bands_hg38$Cytoband_hg38 <- paste0(bands_hg38$CHR, bands_hg38$V4)
## Overlap, as a fraction of gene lenght
md <- merge(comb2[ , c("ID", "CHR","P0_hg38", "P1_hg38")], bands_hg38, by="CHR" )
# Gene length
md$gene_l <- md$P1_hg38-md$P0_hg38
# Consider 3 scenarios: 1)full overlap, 2)overlap at-the-end of band, 3) overlap at-the-start of band.
md$overlap <- ifelse(md$P0_hg38>=md$start & md$P1_hg38<=md$end, md$gene_l/md$gene_l,
                     ifelse(md$P1_hg38>md$end & md$P0_hg38<md$end, (md$end-md$P0_hg38)/md$gene_l, 
                            ifelse(md$P1_hg38>md$start & md$P0_hg38<md$start, (md$P1_hg38-md$start)/md$gene_l, NA ) ) )
## Handle multiple gene-band overlaps
md_ov <- md[is.na(md$overlap)==F, ]
dup <- md_ov$ID[duplicated(md_ov$ID)]
# Select best overlap for genes with multiple instances
# Get rownanes
dup_exclude <- NULL
# keep the one with max overlap
for(i in c(1:length(dup)) ){
row_max <- rownames(md_ov[md_ov$ID==dup[i], ][which.max(md_ov[md_ov$ID==dup[i], "overlap"]), ])
dup_exclude <- append(dup_exclude, rownames(md_ov[md_ov$ID==dup[i], ])[!rownames(md_ov[md_ov$ID==dup[i], ])%in%row_max ] ) 
}
# Remove duplicated rows with poor overlap
md_ov_clean <- md_ov[!rownames(md_ov)%in%dup_exclude, ]
str(factor(md_ov_clean$ID))	# Factor w/ 19801 levels "A1BG","A1BG-AS1",..:
# Add cytoband to file
comb2 <- merge(comb2, md_ov_clean[, c("ID", "Cytoband_hg38")], by="ID", all.x=T)


#---------
# Consistent direction of association
#---------
table(sign(comb2$TWAS.Z)==sign(comb2$EQTL.Z)*sign(comb2$EQTL.GWAS.Z) )
#FALSE  TRUE 
# 2260 17579
 
comb2$EQTL_aligned <- ifelse(sign(comb2$TWAS.Z)==sign(comb2$EQTL.Z)*sign(comb2$EQTL.GWAS.Z), "yes", "no")

#---------
# False discovery rate (FDR)
#---------

comb2$TWAS.P.fdr <- p.adjust(comb2$TWAS.P, method="fdr")
    
              
#--------------------
# Save file (unfiltered)
#--------------------
write.table(comb2[, append(c(2,3,1), seq(4,ncol(comb2),1))], paste0(dir_out, file_out, "_", length(unique(comb2$ID)),"uniqGENES.dat"), col.names=T, row.names=F, quote=F, sep="\t")



#--------------------
# TOP associations
#--------------------

# Bonferroni
comb_bonf <- comb2[comb2$TWAS.P<0.05/nrow(comb2), ]	# 96
length(unique(comb_bonf$ID))				# 96

# FDR<0.01
comb_fdr <- comb2[comb2$TWAS.P.fdr<0.01, ]		# 165

#---------
## TWAS aligned to EQTL?
#---------
table(comb_bonf$EQTL_aligned=="yes")
#TRUE 
#  96

table(comb_fdr$EQTL_aligned=="yes")
#FALSE  TRUE 
#    1   164

#---------
# NOVEL loci
#---------
# TWAS.Z larger than BEST.GWAS.Z?
table(ifelse(abs(comb_fdr$TWAS.Z) > abs(comb_fdr$BEST.GWAS.Z), "yes", "no" ) )
# no yes 
# 156   9

# Sig. TWAS where best GWAS.P < sig. threhold (i.e BEST.GWAS.P>5e-8)
comb_fdr$BEST.GWAS.P <- 2*pnorm(q=abs(comb_fdr$BEST.GWAS.Z), lower.tail=FALSE)
comb_fdr$novel <- ifelse(comb_fdr$BEST.GWAS.P<5e-8, "no", "yes" )
table(comb_fdr$novel)
# no yes 
# 115   50

#---------
# Save
#---------
write.table(comb_fdr[, append(c(2,3,1), seq(4,ncol(comb2),1))], paste0(dir_out, file_out, "_", nrow(comb2),"uniqGENES_", nrow(comb_fdr),"_FDR01.dat"), col.names=T, row.names=F, quote=F, sep="\t")







