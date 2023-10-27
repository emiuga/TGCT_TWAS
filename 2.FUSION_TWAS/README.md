
## Software

Modified FUSION scripts:
- `FUSION.assoc_test_ModifiedPerm_v0.1.R` (TWAS analysis, including permutation test using effective number of iterations)
- `FUSION.post_process_AP-EUM_v1.6.R` (Joint-conditional analysis and customized plots on GWAS summary statistics conditional on predicted expression)

## Analysis flow

- `0.1.Overlap_TECAC_sumstats_and_LD_EUR_panel_20220307.sh` (evaluate imput summary data)
- `0.2.Prepare_Input_sumstats_imputed_v1.sh (formatting of` TGCT GWAS summary statistics for FUSION-TWAS analysis)
- `1.TWAS_v4.sh` (bash scripts to run TWAS analysis on the different expression panels, and compile output)
- `2.Summarize_TWAS_v4_output_v3.20230620.R` (scripts to select single gene-trait assocition prioritizing normal testis, followed by additional gene association from TGCT, and cross-tissue prediction models)
- `2.1.Create_GWAS_sumstats_hg19_forPlotting_v0.sh`
- `3.Joint_conditional_perm100K_NOM_v0.3.2.sh` (bash scripts to run joint-conditional analyis on TWAS results)
