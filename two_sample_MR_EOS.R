# =============================================================================
# Two-sample Mendelian randomization
#   Exposure : brain / CSF / plasma protein QTL (pQTL) instruments
#   Outcome  : Early-onset ischemic stroke (EOS)
#
# Pipeline:
#   1. Load pQTL instruments and select the target protein (here: CSF APOE)
#   2. Format exposure and outcome data
#   3. Harmonise
#   4. MR-PRESSO global test + iterative outlier removal
#   5. MR (8 estimators) + OR
#   6. Heterogeneity, pleiotropy (MR-Egger intercept), I^2
#   7. Leave-one-out
#   8. Visualisation (scatter / forest / funnel / leave-one-out)
#   9. Export all results
# =============================================================================

# ---- Libraries --------------------------------------------------------------
library(TwoSampleMR)
library(MRPRESSO)
library(tidyverse)
library(data.table)
library(ggplot2)
library(MendelianRandomization)
library(ieugwasr)
library(plinkbinr)
library(metafor)
library(openxlsx)

# ---- User settings ----------------------------------------------------------
# >>> CHECK THESE PATHS AND COLUMN NAMES AGAINST YOUR OWN FILES <<<

# --- input files ---
path_pqtl   <- "data/brain-csf-plasma-pqtl.xlsx"   # pQTL instruments
path_eos    <- "data/EOS_gwas.txt"                 # EOS outcome summary statistics
out_root    <- "output/EPHX2_result"               # results root directory

# --- target instrument ---
target_tissue  <- "CSF"      # one of: Brain / CSF / Plasma
target_protein <- "APOE"     # protein to analyse
set.seed(5201314)

# --- pQTL (exposure) column names in path_pqtl ---
pqtl_tissue_col   <- "Tissue"
pqtl_protein_col  <- "Protein"
pqtl_snp_col      <- "SNP"
pqtl_beta_col     <- "beta"
pqtl_se_col       <- "se"
pqtl_ea_col       <- "effect_allele"
pqtl_oa_col       <- "other_allele"
pqtl_pval_col     <- "p"

# --- EOS (outcome) column names in path_eos ---
# Defaults follow a typical GWAS summary-statistics layout; edit to match your
# EOS file (e.g. EOSC / Jaworek 2022).
eos_phenotype     <- "EOS"            # label written into the outcome data
eos_snp_col       <- "SNP"
eos_beta_col      <- "beta"
eos_se_col        <- "se"
eos_ea_col        <- "effect_allele"
eos_oa_col        <- "other_allele"
eos_pval_col      <- "pval"
eos_eaf_col       <- "eaf"            # set to NULL below if not available
eos_chr_col       <- "chr"
eos_pos_col       <- "pos"
eos_samplesize    <- 16927            # EOS total sample size (cases + controls); edit as appropriate

# ---- Load pQTL instruments and select target --------------------------------
pqtl <- read.xlsx(path_pqtl)
exposure_tissue <- pqtl[pqtl[[pqtl_tissue_col]] == target_tissue, ]
exp_dat_clump   <- exposure_tissue[exposure_tissue[[pqtl_protein_col]] == target_protein, ]

t1 <- Sys.time()

# ---- Outcome (EOS) ----------------------------------------------------------
outcome <- as.data.frame(fread(path_eos))
outcome$phenotype  <- eos_phenotype
outcome$samplesize <- eos_samplesize

out_dat <- format_data(
  dat               = outcome,
  type              = "outcome",
  snps              = exp_dat_clump[[pqtl_snp_col]],
  phenotype_col     = "phenotype",
  snp_col           = eos_snp_col,
  beta_col          = eos_beta_col,
  se_col            = eos_se_col,
  effect_allele_col = eos_ea_col,
  other_allele_col  = eos_oa_col,
  pval_col          = eos_pval_col,
  eaf_col           = eos_eaf_col,
  chr_col           = eos_chr_col,
  pos_col           = eos_pos_col,
  samplesize_col    = "samplesize"
)
out_dat <- out_dat %>% subset(., !duplicated(SNP))

# ---- Exposure ---------------------------------------------------------------
exp_data_clump <- format_data(
  dat               = exp_dat_clump,
  type              = "exposure",
  phenotype_col     = pqtl_protein_col,
  snp_col           = pqtl_snp_col,
  beta_col          = pqtl_beta_col,
  se_col            = pqtl_se_col,
  effect_allele_col = pqtl_ea_col,
  other_allele_col  = pqtl_oa_col,
  pval_col          = pqtl_pval_col
)

# ---- Harmonise --------------------------------------------------------------
mydata <- harmonise_data(exposure_dat = exp_data_clump, outcome_dat = out_dat, action = 2)
mydata <- mydata[which(mydata$mr_keep == TRUE), ]

mydata_clump         <- mydata
mydata_filteroutcome <- mydata_clump
mydata_outcome_n     <- dim(mydata_filteroutcome)[1]
Nbd                  <- 10000   # NbDistribution for MR-PRESSO

# ---- MR-PRESSO: global test and iterative outlier removal -------------------
if (mydata_outcome_n <= 3) {
  presso_pval   <- "NA"
  mydata_presso <- mydata_filteroutcome
} else {
  presso <- mr_presso(
    BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure",
    SdOutcome = "se.outcome", SdExposure = "se.exposure",
    OUTLIERtest = TRUE, DISTORTIONtest = TRUE,
    data = mydata_filteroutcome,
    NbDistribution = Nbd, SignifThreshold = 0.05
  )
  presso_pval <- presso$`MR-PRESSO results`$`Global Test`$Pvalue
  presso_snp  <- presso$`MR-PRESSO results`$`Outlier Test`

  if (presso_pval >= 0.05) {
    # No significant global pleiotropy: keep all variants
    mydata_presso <- mydata_filteroutcome
  } else {
    # Significant pleiotropy: remove outliers one by one (in order of outlier
    # p-value) until the global test is no longer significant
    out_order <- order(presso_snp$Pvalue)
    for (ii in seq_along(out_order)) {
      snp_ii        <- out_order[1:ii]
      mydata_presso <- mydata_filteroutcome[-snp_ii, ]
      mydata_pres_n <- dim(mydata_presso)[1]

      if (mydata_pres_n <= 3) {
        presso_pval <- "NA"
      } else {
        presso2 <- mr_presso(
          BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure",
          SdOutcome = "se.outcome", SdExposure = "se.exposure",
          OUTLIERtest = TRUE, DISTORTIONtest = TRUE,
          data = mydata_presso,
          NbDistribution = Nbd, SignifThreshold = 0.05
        )
        presso_pval <- presso2$`MR-PRESSO results`$`Global Test`$Pvalue
      }
      message("Outlier removal iteration: ", ii)
      if (presso_pval >= 0.05) break
    }
  }
}

# ---- Instrument-count bookkeeping through the pipeline ----------------------
exp_n            <- dim(exp_dat_clump)[1]
out_n            <- dim(out_dat)[1]
mydata_har_n     <- dim(mydata)[1]
mydata_clu_n     <- dim(mydata_clump)[1]
mydata_outcome_n <- dim(mydata_filteroutcome)[1]
mydata_pres_n    <- dim(mydata_presso)[1]
filter_data <- data.frame(exp_n, out_n, mydata_clu_n, mydata_har_n,
                          mydata_outcome_n, mydata_pres_n, presso_pval)

# ---- MR analysis (8 estimators) ---------------------------------------------
res <- mr(mydata_presso, method_list = c(
  "mr_egger_regression", "mr_ivw", "mr_ivw_mre", "mr_ivw_fe",
  "mr_wald_ratio", "mr_weighted_median", "mr_two_sample_ml",
  "mr_simple_mode", "mr_weighted_mode"
))
# Reduced estimator set used for the scatter plot
res_scatter <- mr(mydata_presso, method_list = c(
  "mr_egger_regression", "mr_ivw", "mr_weighted_median",
  "mr_weighted_mode", "mr_simple_mode"
))
res_or <- generate_odds_ratios(res)

# ---- Heterogeneity, pleiotropy ----------------------------------------------
het <- mr_heterogeneity(mydata_presso,
                        method_list = c("mr_egger_regression", "mr_ivw"))
plt <- mr_pleiotropy_test(mydata_presso)

# ---- I^2 (via single-SNP estimates + fixed-effect meta-analysis) ------------
res_single1 <- mr_singlesnp(mydata_presso, all_method = "mr_ivw")
res_single2 <- res_single1[grep("^rs", res_single1$SNP), ]
res_meta <- metafor::rma(yi = res_single2$b, sei = res_single2$se,
                         weights = 1 / mydata_presso$se.outcome^2,
                         data = res_single2, method = "FE")
I2 <- as.data.frame(cbind(res_meta$I2, res_meta$H2))
colnames(I2) <- c("I2", "H2")   # I2 = percentage; H2 = Cochran's Q-based statistic

# ---- Leave-one-out ----------------------------------------------------------
res_loo <- mr_leaveoneout(mydata_presso)

# ---- Visualisation ----------------------------------------------------------
p1   <- mr_scatter_plot(res, mydata_presso)
p1_2 <- mr_scatter_plot(res_scatter, mydata_presso)
res_single <- mr_singlesnp(mydata_presso)
p2 <- mr_forest_plot(res_single)
p3 <- mr_funnel_plot(res_single)
p4 <- mr_leaveoneout_plot(res_loo)

# ---- Export -----------------------------------------------------------------
out_dir <- file.path(out_root, paste0(target_tissue, "_", target_protein, "_", eos_phenotype))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

ggsave(p1[[1]],   file = file.path(out_dir, "mr_scatter_plot.pdf"),     height = 7, width = 7)
ggsave(p1_2[[1]], file = file.path(out_dir, "mr_scatter_plot2.pdf"),    height = 7, width = 7)
ggsave(p2[[1]],   file = file.path(out_dir, "mr_forest_plot.pdf"),      height = 7, width = 7)
ggsave(p3[[1]],   file = file.path(out_dir, "mr_funnel_plot.pdf"),      height = 7, width = 7)
ggsave(p4[[1]],   file = file.path(out_dir, "mr_leaveoneout_plot.pdf"), height = 7, width = 7)

write.table(mydata_presso, file.path(out_dir, "mydata_presso.txt"), sep = "\t", quote = FALSE, row.names = FALSE)  # instruments used
write.table(res,           file.path(out_dir, "res.txt"),           sep = "\t", quote = FALSE, row.names = FALSE)  # MR estimates
write.table(res_or,        file.path(out_dir, "res_or.txt"),        sep = "\t", quote = FALSE, row.names = FALSE)  # MR estimates as OR
write.table(het,           file.path(out_dir, "het.txt"),           sep = "\t", quote = FALSE, row.names = FALSE)  # heterogeneity
write.table(plt,           file.path(out_dir, "plt.txt"),           sep = "\t", quote = FALSE, row.names = FALSE)  # pleiotropy (MR-Egger intercept)
write.table(filter_data,   file.path(out_dir, "filter_data.txt"),   sep = "\t", quote = FALSE, row.names = FALSE)  # instrument counts through pipeline
write.table(I2,            file.path(out_dir, "I2.txt"),            sep = "\t", quote = FALSE, row.names = FALSE)  # I^2 statistic

t2 <- Sys.time()
message("MR for ", target_tissue, " ", target_protein, " -> ", eos_phenotype,
        " done in ", round(difftime(t2, t1, units = "mins"), 3), " minutes.")
