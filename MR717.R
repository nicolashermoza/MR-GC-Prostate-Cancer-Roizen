library(TwoSampleMR)
library(R.utils)
library(dplyr)
library(jsonlite)
library(httr)
library(MRPRESSO)
library(RadialMR)

exposure_dat <- read_exposure_data(
  filename = "albinanavitd_clean.tsv",
  sep = "\t",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "standard_error",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  eaf_col = "effect_allele_frequency",
  pval_col = "p_value"
)


# Outcome GWAS with altered names, this is for prostate cancer
outcome_raw <- read.delim("prostateverma.tsv.gz", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
# Remove rows with missing values
outcome_raw <- outcome_raw[!is.na(outcome_raw$odds_ratio) & !is.na(outcome_raw$ci_upper) & !is.na(outcome_raw$ci_lower), ]

# Calculate beta and standard error
outcome_raw$beta <- log(outcome_raw$odds_ratio)
outcome_raw$standard_error <- (log(outcome_raw$ci_upper) - log(outcome_raw$ci_lower)) / (2 * 1.96)

outcome_dat <- format_data(
  outcome_raw,
  type = "outcome",
  snp_col = "rsid",
  beta_col = "beta",
  se_col = "standard_error",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  eaf_col = "effect_allele_frequency",  # or remove this if missing
  pval_col = "p_value"
)

exposure_dat <- exposure_dat[exposure_dat$pval.exposure < 5e-8, ]

harmonised_dat <- harmonise_data(
  exposure_dat = exposure_dat,
  outcome_dat = outcome_dat,
)

mr_results <- mr(harmonised_dat)
print(mr_results)
write.table(mr_results, "MR_results.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
mr_scatter_plot(mr_results, harmonised_dat)

results <- mr(harmonised_dat, method_list = c(
  "mr_ivw", 
  "mr_egger_regression", 
  "mr_weighted_median", 
  "mr_weighted_mode"
))
print(results)

# Harmonise
harmonised_dat <- harmonise_data(exposure_dat, outcome_dat)
harmonised_dat <- harmonised_dat[harmonised_dat$pval.exposure < 5e-8, ]
# Create radial-ready dataset

radial_data <- format_radial(
  BXG = harmonised_dat$beta.exposure,
  BYG = harmonised_dat$beta.outcome,
  seBXG = harmonised_dat$se.exposure,
  seBYG = harmonised_dat$se.outcome,
)

results_radial <- ivw_radial(radial_data, alpha = 0.05)
str(results_radial$outliers)
outlier_indices <- results_radial$outliers$SNP
outlier_snps <- harmonised_dat$SNP[outlier_indices]
cleaned_data <- harmonised_dat[!harmonised_dat$SNP %in% outlier_snps, ]

results = mr(cleaned_data)
mr_scatter_plot(results, cleaned_data)
print(mr_results)
print(mr_pleiotropy_test(cleaned_data)) # see pleiotropy
print(mr_heterogeneity(cleaned_data)) # see heterogeneity
table(cleaned_data$mr_keep) # 3175 snps used