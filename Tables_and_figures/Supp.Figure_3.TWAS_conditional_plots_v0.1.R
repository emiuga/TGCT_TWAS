#' Scripts to compile conditional plots (TCGT TWAS)
#' Emilio Ugalde 
#' July 10, 2023

#' AIM: to generate Suppl. doc file.
#' Figure title: Supplementary Figure 3. TGCT GWAS signals explained by predicted expression at the corresponding TWAS risk locus.
#' Figure captions: Locus index 1. Locus index 2., ..., Locus index 51. 
#' NOTE: Plots for GWAS conditional on predicted  UCK2 and GINM1 expression are shown in Figure 3.  

#' PREPARE FILES
#' Directory
wd="/run/user/1000/gvfs/sftp:host=vector.meb.ki.se/nfs/GENETEC/TWAS/"
dir=paste0(wd, "Output/imputed_TECAC/Joint_conditional_perm100K_NOM/cond_plots_v3.1")

#' Output Directory
out_dir=paste0(wd, "Manuscript/Figures/")
FILE=paste0(out_dir, "Suppl.Figure_3")

#' Compiled top-TWAS results
twas <- read.delim(file=paste0(wd, "Manuscript/Tables/annot_100KpermNOM_only/TWAS_testis_tgct_scca_19839uniqGENES_91permNOM_51LOCI.master_ANNOT.dat"))
# Get loci
top <- twas[is.na(twas$LOCUS)==F, ]
top <- top[!duplicated(top$LOCUS), ]
top <- top[order(top$LOCUS.index), ]

# INDEPENDENT Genes at locus
top_loc <- twas[twas$JOINT & twas$PERM.nom=="yes", ]
top$JOINT.GENES.IDs <- rep(NA, nrow(top))
for(i in 1:nrow(top)) {
  LOCUS=top[i,]$LOCUS
  for(g in 1:nrow(top_loc)) {
    top$JOINT.GENES.IDs[i] <- ifelse(is.na(top$JOINT.GENES.IDs[i])==T, ifelse(top_loc[g,]$LOCUS==LOCUS & top_loc[g,]$JOINT==T, top_loc[g,]$ID, top$JOINT.GENES.IDs[i]),
                                     ifelse(top_loc[g,]$LOCUS==LOCUS & top_loc[g,]$JOINT==T, paste(top$JOINT.GENES.IDs[i], top_loc[g,]$ID, sep=", "), top$JOINT.GENES.IDs[i] ) )
  }
}


#' Get list of files
p <- list.files(path = dir)
#' Locus variable name
p_loc <- unlist(lapply(p, function(x){paste(strsplit(x, "[.]")[[1]][1], strsplit(x, "[.]")[[1]][3], sep="_") } ))
#' Match to TWAS data set
locus <- top[top$LOCUS%in%p_loc, ]$LOCUS
m <- match(locus, p_loc)
p <- p[m]

#' Files to import
#' NOTE: exclude GINM1 and UCK2 alreday in Main Figure_3
#exc <- which(locus%in%twas[twas$ID%in%c("GINM1", "UCK2"), ]$LOCUS )
#p_in <- p[-exc]

#' LOCUS INDEX
ind <- top[top$LOCUS%in%locus, ]$LOCUS.index
#ind_in <- ind[-exc]

#' Objet with plot's file and locus name
N=length(ind)
out <- as.data.frame(matrix(rep(NA, N*4), ncol=4))
names(out) <- c("file", "locus", "locus_N", "genes")
for(i in 1:length(ind)){
  out$file[i] <- p[i] 
  out$locus[i] <- paste0("Locus index ", ind[i])
  out$genes[i] <- top$JOINT.GENES.IDs[i]
  out$locus_N[i] <- ind[i]
}

text <- paste0(out$locus_N[1], ", ",out$genes[1], ".")
for(i in 2:nrow(out) ){ 
  info <- paste0(out$locus_N[i], ", ",out$genes[i], ".")
  text <- paste(text, info, sep=" " )
  }

text <- paste0(out$locus_N[1], " (",out$genes[1], ").")
for(i in 2:nrow(out) ){ 
  info <- paste0(out$locus_N[i], " (",out$genes[i], ").")
  text <- paste(text, info, sep=" " )
}

legend <- paste0("**", out$locus_N[1], "** (*",out$genes[1], "*).")
for(i in 2:nrow(out) ){ 
  info <- paste0("**", out$locus_N[i], "** (*",out$genes[i], "*).")
  legend <- paste(legend, info, sep=" " )
}

description <- "Locus index 1-51. Each figure consists of an upper and lower sub-panel. The upper sub-panel shows genes located at the risk locus according to genomic position (hg19). Based on the joint modeling, genes are labeled as leading in dark blue (i.e. joint-conditionally independent), conditional on leading genes(s) in light blue, or correlated with leading genes in green (i.e. predicted expression R2>0.9); marginally associated genes (i.e., FDRTWAS≥0.01 or Ppermutation≥0.05) and not tested genes (e.g. lacking a significant prediction model) are indicated with empty boxes with dark and light outlines, respectively. The lower sub-panel shows results from the SNP level GWAS-conditional analysis. Each point corresponds to the association between SNP and TGCT status. Gray points indicate the marginal association of a SNP with TGCT status (i.e., GWAS association). Green points indicate the association of the same SNPs with TGCT after conditioning of predicted expression of the leading gene(s) at each locus. The horizontal dashed line indicates the genome-wide significance threshold (i.e., P = 5 ×10-8). Locus index: "

#' Markdown/pandoc
#' ===============

library(rmarkdown)

PLOTDIR=dir

AUTHOR <- "Ugalde-Morales et al."
TITLE  <- "Supplementary Figure 3. TGCT GWAS signals conditional on predicted expression at 51 TWAS risk locus."
DATE   <- "2023-07-14"

fn <- paste0(FILE, ".md")
header <- c(
  "---",
  paste0("title: \"", TITLE, "\""),
  paste0("author: \"", AUTHOR, "\""),
  paste0("date: \"", DATE, "\""),
  "output: word_document",
  "---", ""
)

con <- file(fn, open = "w")
writeLines(header, con)
writeLines(description, con)
writeLines(legend, con)
lapply(out$file, function(x) {
  figtxt <- c(
    paste0("*", out[out$file%in%x, ]$locus, ": ", out[out$file%in%x, ]$genes ,"*"),
    paste0("![Gene ", out[out$file%in%x, ]$genes, "](", file.path(PLOTDIR, x), "){width=100%}"), ""
  )
  writeLines(figtxt, con)
})
close(con)

rmarkdown::pandoc_convert(fn, to = "docx", output = paste0(FILE, ".docx"))
#file.show(paste0(out_dir, "Suppl.Figure_3.docx")



  
  