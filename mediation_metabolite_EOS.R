# =============================================================================
# Two-step Mendelian randomization: mediation analysis
# (product-of-coefficients method)
#
#   Exposure : EPHX2 (brain pQTL)
#   Mediator : brain / CSF metabolite
#   Outcome  : Early-onset ischemic stroke (EOS)
#
# Effect definitions (standard two-step MR):
#   beta_a = EPHX2    -> metabolite   (step 1)
#   beta_b = metabolite -> EOS        (step 2)
#   beta_c = EPHX2    -> EOS          (total effect)
#   indirect (mediated) effect = beta_a * beta_b
#   mediation proportion       = (beta_a * beta_b) / beta_c
# The 95% CI of the indirect effect is obtained with RMediation::medci.
# =============================================================================

# ---- Libraries --------------------------------------------------------------
library(data.table)
library(openxlsx)
library(RMediation)

# ---- User settings ----------------------------------------------------------
# Each of the three MR result files below should contain, at minimum, columns
# for exposure, outcome, beta and se. The pipeline merges them into one table
# carrying beta_a/se_a (step 1), beta_b/se_b (step 2) and beta_c/se_c (total).
#
# >>> CHECK THESE PATHS AND COLUMN NAMES AGAINST YOUR OWN FILES <<<
# If your column names differ, edit the *_col variables; nothing else needs to
# change.

# --- input files ---
path_c <- "data/EPHX2_EOS_total.csv"        # EPHX2 -> EOS        (total effect, beta_c)
path_a <- "data/EPHX2_metabolite.csv"       # EPHX2 -> metabolite (step 1, beta_a)
path_b <- "data/metabolite_EOS.csv"         # metabolite -> EOS   (step 2, beta_b)

# --- output files ---
out_merged    <- "output/EPHX2_metabolite_EOS_merged.xlsx"
out_mediation <- "output/EPHX2_metabolite_EOS_mediation.xlsx"

# --- labels ---
exposure_name <- "EPHX2"
outcome_name  <- "EOS"

# --- column names in the three input files ---
# total-effect file (EPHX2 -> EOS)
c_exposure_col <- "exposure"
c_outcome_col  <- "outcome"
c_beta_col     <- "b"        # or "beta"
c_se_col       <- "se"

# step-1 file (EPHX2 -> metabolite): exposure = EPHX2, outcome = metabolite
a_exposure_col <- "exposure"
a_outcome_col  <- "outcome"   # the metabolite identifier
a_beta_col     <- "b"
a_se_col       <- "se"

# step-2 file (metabolite -> EOS): exposure = metabolite, outcome = EOS
b_exposure_col <- "exposure"  # the metabolite identifier
b_outcome_col  <- "outcome"
b_beta_col     <- "b"
b_se_col       <- "se"

# ---- Read inputs ------------------------------------------------------------
tot <- as.data.frame(fread(path_c))   # EPHX2 -> EOS   (beta_c)
s1  <- as.data.frame(fread(path_a))   # EPHX2 -> metab (beta_a)
s2  <- as.data.frame(fread(path_b))   # metab -> EOS   (beta_b)

# Standardise column names so the merges below are robust to the source layout
tot <- data.frame(
  exposure = tot[[c_exposure_col]],
  outcome  = tot[[c_outcome_col]],
  beta_c   = tot[[c_beta_col]],
  se_c     = tot[[c_se_col]],
  stringsAsFactors = FALSE
)
s1 <- data.frame(
  exposure  = s1[[a_exposure_col]],
  metabolite = s1[[a_outcome_col]],
  beta_a    = s1[[a_beta_col]],
  se_a      = s1[[a_se_col]],
  stringsAsFactors = FALSE
)
s2 <- data.frame(
  metabolite = s2[[b_exposure_col]],
  outcome   = s2[[b_outcome_col]],
  beta_b    = s2[[b_beta_col]],
  se_b      = s2[[b_se_col]],
  stringsAsFactors = FALSE
)

# ---- Assemble exposure / mediator / outcome table ---------------------------
# Restrict to the target exposure and outcome
tot <- tot[tot$exposure == exposure_name & tot$outcome == outcome_name, ]

# Step 2 (metabolite -> EOS) defines the set of candidate mediators
da1 <- s2[s2$outcome == outcome_name, ]

# Attach step 1 (EPHX2 -> metabolite) on the metabolite key
da1 <- merge(da1, s1[s1$exposure == exposure_name, c("metabolite", "beta_a", "se_a")],
             by = "metabolite")

# Attach the total effect (EPHX2 -> EOS); it is constant across mediators
da1$exposure <- exposure_name
da1 <- merge(da1, tot[, c("exposure", "outcome", "beta_c", "se_c")],
             by = c("exposure", "outcome"))

# Final column order
da1 <- da1[, c("exposure", "metabolite", "outcome",
               "beta_a", "se_a", "beta_b", "se_b", "beta_c", "se_c")]
colnames(da1)[2] <- "mediate"

write.xlsx(da1, out_merged)

# ---- Direction consistency --------------------------------------------------
info <- da1
info$total <- info$beta_c
info$med   <- info$beta_a * info$beta_b
info$dir   <- info$total - info$med

# rows where total, mediated and direct effects share the same sign
pos1 <- which(info$total > 0 & info$med > 0 & info$dir > 0)
pos2 <- which(info$total < 0 & info$med < 0 & info$dir < 0)

dat <- info
dat$direction_consistency <- "No"
dat$direction_consistency[c(pos1, pos2)] <- "Yes"

# ---- Indirect effect, CI, p-value and mediation proportion ------------------
res_m <- data.frame()
for (k in seq_len(nrow(dat))) {
  res_med <- medci(mu.x = dat$beta_a[k], se.x = dat$se_a[k],
                   mu.y = dat$beta_b[k], se.y = dat$se_b[k],
                   type = "asymp")

  lci  <- res_med[[1]][1]
  uci  <- res_med[[1]][2]
  beta <- res_med[[2]]
  se   <- res_med[[3]]

  # two-sided p-value for the indirect effect
  p_mediation <- 2 * pnorm(q = abs(res_med$Estimate / res_med$SE), lower.tail = FALSE)

  # proportion of the total effect explained by mediation, with 95% CI
  med_por <- res_med$Estimate    / dat$beta_c[k]
  med_lci <- res_med$`95% CI`[1] / dat$beta_c[k]
  med_uci <- res_med$`95% CI`[2] / dat$beta_c[k]

  res_m <- rbind(res_m,
                 data.frame(lci, uci, beta, se, p_mediation,
                            med_lci, med_uci, med_por))
}

final <- cbind(dat, res_m)
write.xlsx(final, out_mediation)

message("Mediation analysis complete: ", nrow(final),
        " mediator(s) for ", exposure_name, " -> ", outcome_name, ".")
