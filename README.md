# Code availability for research paper: 
Identification of genes associated with testicular germ cell tumor susceptibility through a transcription wide association study.
Ugalde-Morales E. et al. <emilio.ugalde.morales@ki.se>

## DESCRIPTION

The present repository contains computer code used to perform Trancriptome Wide Association Study (TWAS) on TGCT risk, following FUSION methodology (Gusev et al. 2016)[^1], as well as complementary analysis on external gene expression data sets related to normal testis embryonal development (Human gonadal development)[^2], pre-malignant tissue (i.e GCNIS, germ cell neoplasia in situ)[^3], and TGCT tumor tissue[^4]. 

Code regarding data preparation `1.GWAS_imputation`, data analysis on TWAS `2.FUSION_TWAS` and expression data sets `3.Gene_expression`, and generation of tables and figures for publication `Tables_and_figures`, is provided.

Modified version of FUSION software (R scripts) are located under `2.FUSION_TWAS/Software` folder. 

Detailed methodology is described in Methods section of the manuscript. 


## DATA SOURCES 

**GWAS summary statistics**

TGCT GWAS meta-analysis summary statistics (Pluta et al. 2021) at availabe at https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/study.cgi?study_id=phs001349.v2.p1.
Summary-GWAS harmonization and imputation was done using the MetaXcan software (https://github.com/hakyimlab/summary-gwas-imputation), please se README and installation/requirements files within ./Scripts/GWAS_imputation folder.

**TWAS analysis**

Visit FUSION website for software installation and documentation [FUSION](http://gusevlab.org/projects/fusion/).
Pre-computed prediction models on testis tissue (GTExv8.ALL), TGCT (TCGA), and cross-tissue (GTExv8) were download as follow:

```bash
mkdir -p ./Data/WEIGHTS
cd ./Data/WEIGHTS
wget https://s3.us-west-1.amazonaws.com/gtex.v8.fusion/ALL/GTExv8.ALL.Testis.tar.gz
wget http://gusevlab.org/projects/fusion/weights/TCGA-TGCT.TUMOR.tar.bz2
wget http://gusevlab.org/projects/fusion/weights/sCCA_weights_v8_2.zip
tar xzf GTExv8.ALL.Testis.tar.gz
tar xjf TCGA-TGCT.TUMOR.tar.bz2
unzip sCCA_weights_v8_2.zip
```

**Gene expression analysis**

Data sets can be obtained from:

- Human gonadal development [The Reproductive Cell Atlas](https://www.reproductivecellatlas.org/gonads.html).
- Pre-malignan tissue [GCNIS data set](https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-TABM-488).
- TGCT tumor samples [TGCA data set](http://firebrowse.org/?cohort=TGCT&download_dialog=true).

To download single-cell gonadal data on developmental germ cells `human_germcells.h5ad` and male germ cells `human_germcells.h5ad`:

```bash
mkdir - Data/Gonads
cd Data/Gonads
wget https://cellgeni.cog.sanger.ac.uk/vento/reproductivecellatlas/gonads/human_main_male.h5ad
wget https://cellgeni.cog.sanger.ac.uk/vento/reproductivecellatlas/gonads/human_germcells.h5ad
```

## References
[^1]:Gusev A., et al. Nat Genet. 2016 Mar;48(3):245-52. doi: 10.1038/ng.3506. Epub 2016 Feb 8. 
[^2]:Garcia-Alonso L., et al. Nature. 2022 Jul;607(7919):540-547. doi: 10.1038/s41586-022-04918-4.
[^3]:Sonne S. B. et al. Cancer Res. 2009 Jun 2. doi: 10.1158/0008-5472.CAN-08-4554.
[^4]:Shen H., et al. Cell Reports. 2018 Jun 12;23(11):3392-3406. doi: 10.1016/j.celrep.2018.05.039. 

