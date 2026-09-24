# SAVE-MORE proteomics preprocessing
# Run chunk-by-chunk from the project root (or from scripts/; it will setwd("..")).
# Each step leaves an object you can inspect (clinical_data, data, data1, data2, ...).
#
# QC order, then mice PMM per assay (not Allocation / srf14):
#   NPX_imputed      = first mice draw
#   NPX_imputed_mean = mean of m = 5 draws
#   NPX_wins         = soft-winsor of NPX_imputed
#   NPX_invnorm      = orderNorm(jitter(NPX_wins)) per Assay
#
# Writes:
#   data/clinical_data.tsv
#   data/Combined_dataset_mice.csv
#   data/qc_unmatched_ids.csv

if (file.exists("data/SAVE_MORE_NPX.xlsx")) {
  # project root
} else if (file.exists("../data/SAVE_MORE_NPX.xlsx")) {
  setwd("..")
} else {
  stop("Open the SAVE_MORE_Proteomics project root, then source this script.")
}

library(haven)
library(readxl)
library(tidyverse)
library(OlinkAnalyze)
library(mice)
library(bestNormalize)

# =============================================================================
# 1. Clinical data
# SPSS = one row per patient. SM-DATA.xlsx has patient outcomes plus
# day-specific labs (crp1/4/7, ferritin1/4/7, il61/64/67). Those labs are
# pivoted to long form so they join to Olink on PatientID + Day.
# =============================================================================
dataset <- read_sav("data/DATA.sav")
dataset_c <- haven::as_factor(dataset) %>%
  dplyr::rename(PatientID = Code)
dataset_c$who28_d <- as.numeric(dataset$who28)

clinical_xlsx <- readxl::read_xlsx("data/SM-DATA.xlsx")

clinical_data2 <- clinical_xlsx %>%
  dplyr::rename(PatientID = Code, who28_v2 = who28) %>%
  dplyr::select(PatientID, sofa14, who14, who28_v2)

clinical_labs <- clinical_xlsx %>%
  dplyr::rename(PatientID = Code) %>%
  dplyr::select(
    PatientID,
    ferritin1, ferritin4, ferritin7,
    crp1, crp4, crp7,
    il61, il64, il67
  ) %>%
  pivot_longer(
    cols = -PatientID,
    names_to = c(".value", "Day"),
    names_pattern = "^(ferritin|crp|il6)([147])$"
  ) %>%
  mutate(
    Day = dplyr::recode(as.character(Day), "1" = "01", "4" = "04", "7" = "07"),
    CRP = crp,
    IL6 = il6
  ) %>%
  dplyr::select(PatientID, Day, CRP, ferritin, IL6)

clinical_data <- dataset_c %>%
  left_join(clinical_data2, by = "PatientID")

clinical_data
clinical_labs
write_tsv(clinical_data, "data/clinical_data.tsv")

# =============================================================================
# 2. NPX, parse SampleID -> PatientID + Day, drop plate controls
# =============================================================================
data <- as_tibble(OlinkAnalyze::read_NPX("data/SAVE_MORE_NPX.xlsx"))
data <- data %>%
  mutate(
    NPX = as.numeric(NPX),
    LOD = as.numeric(LOD),
    MissingFreq = as.numeric(MissingFreq)
  )

data_ids <- data %>%
  mutate(
    SampleID_clean = str_replace_all(as.character(SampleID), "\\s+", ""),
    Day = str_extract(SampleID_clean, "(01|04|07)$"),
    PatientID = if_else(
      !is.na(Day),
      str_remove(SampleID_clean, "(01|04|07)$"),
      NA_character_
    )
  ) %>%
  relocate(PatientID, Day, SampleID_clean, .before = SampleID)

control_samples <- data_ids %>%
  filter(is.na(Day) | str_starts(SampleID_clean, "SC")) %>%
  distinct(SampleID)

data_test <- data_ids %>%
  filter(!SampleID %in% control_samples$SampleID) %>%
  filter(Day %in% c("01", "04", "07"))

nrow(data_ids)
nrow(data_test)
n_distinct(data_test$PatientID)

# =============================================================================
# 3. Join clinical; keep rows with Allocation
# =============================================================================
npx_only <- sort(setdiff(unique(data_test$PatientID), unique(clinical_data$PatientID)))
clinical_only <- sort(setdiff(unique(clinical_data$PatientID), unique(data_test$PatientID)))
npx_only
clinical_only

unmatched_ids <- bind_rows(
  tibble(PatientID = clinical_only, reason = "in_clinical_not_in_npx"),
  tibble(PatientID = npx_only, reason = "in_npx_not_in_clinical")
)
write_csv(unmatched_ids, "data/qc_unmatched_ids.csv")

data1 <- data_test %>%
  left_join(clinical_data, by = "PatientID") %>%
  left_join(clinical_labs, by = c("PatientID", "Day")) %>%
  drop_na(Allocation)

n_distinct(data1$PatientID)
n_distinct(data1$SampleID)

# =============================================================================
# 4. Sample QC_Warning (whole sample-panel) and assays with MissingFreq >= 0.20
# =============================================================================
failed_samples <- data1 %>%
  filter(QC_Warning != "Pass") %>%
  distinct(SampleID)

failed_proteins <- data1 %>%
  distinct(Assay, MissingFreq) %>%
  filter(MissingFreq >= 0.20)

failed_samples
failed_proteins

data2 <- data1 %>%
  filter(!SampleID %in% failed_samples$SampleID) %>%
  filter(!Assay %in% failed_proteins$Assay)

nrow(data2)
n_distinct(data2$Assay)

# =============================================================================
# 5. NPX below LOD -> NA
# =============================================================================
data_lod <- data2 %>%
  mutate(NPX = if_else(!is.na(NPX) & !is.na(LOD) & NPX < LOD, NA_real_, NPX))

sum(is.na(data_lod$NPX)) - sum(is.na(data2$NPX))

# =============================================================================
# 6. Sample-level Olink IQR / median outliers (3 SD)
# Use Olink columns only. Clinical labs (CRP, ferritin, IL6) stay on data_lod
# but are not passed into olink_qc_plot. Fence lines are drawn once per Panel
# (ggplot2 4 draws aes(yintercept=...) once per sample and paints over points).
# =============================================================================
npx_qc_cols <- c(
  "SampleID", "Index", "OlinkID", "UniProt", "Assay", "MissingFreq",
  "Panel_Version", "PlateID", "QC_Warning", "LOD", "NPX", "Normalization",
  "Assay_Warning", "Panel", "PatientID", "Day"
)
data_qc <- data_lod %>%
  dplyr::select(dplyr::any_of(npx_qc_cols))

qc_plot_olink <- olink_qc_plot(
  data_qc,
  color_g = "Day",
  IQR_outlierDef = 3,
  median_outlierDef = 3,
  label_outliers = FALSE,
  outlierLines = FALSE
)

qc_df <- qc_plot_olink$data
outlier_samples <- qc_df %>%
  dplyr::filter(Outlier == 1) %>%
  distinct(SampleID)

outlier_samples
nrow(outlier_samples)

qc_fences <- qc_df %>%
  distinct(Panel, median_low, median_high, iqr_low, iqr_high)

qc_plot <- ggplot(qc_df, aes(x = sample_median, y = IQR, color = factor(Outlier), shape = Day)) +
  geom_point(size = 2.2) +
  geom_vline(data = qc_fences, aes(xintercept = median_low), linetype = "dashed", colour = "grey50") +
  geom_vline(data = qc_fences, aes(xintercept = median_high), linetype = "dashed", colour = "grey50") +
  geom_hline(data = qc_fences, aes(yintercept = iqr_low), linetype = "dashed", colour = "grey50") +
  geom_hline(data = qc_fences, aes(yintercept = iqr_high), linetype = "dashed", colour = "grey50") +
  scale_color_manual(values = c("0" = "grey50", "1" = "#D64541"), labels = c("0" = "in", "1" = "outlier")) +
  facet_wrap(~ Panel, scales = "free") +
  labs(x = "Sample median NPX", y = "IQR", color = NULL, shape = "Day") +
  theme_bw()

qc_plot_olink

data3 <- data_lod %>%
  filter(!SampleID %in% outlier_samples$SampleID)

nrow(data3)
n_distinct(data3$SampleID)
n_distinct(data3$PatientID)

# =============================================================================
# 7. mice PMM per Assay
# Predictors: Day, Age, cbmi, Gender (not Allocation / srf14)
# NPX_imputed = complete(..., 1); NPX_imputed_mean = mean of 5 draws
# =============================================================================
set.seed(500)
assays <- unique(data3$Assay)
mice_pieces <- vector("list", length(assays))

for (i in seq_along(assays)) {
  assay <- assays[[i]]
  dat <- data3 %>% filter(Assay == assay)
  npx <- as.numeric(dat$NPX)
  n_miss <- sum(is.na(npx))
  message(i, "/", length(assays), "  ", assay, "  missing NPX = ", n_miss)

  if (n_miss == 0L || sum(!is.na(npx)) < 5L) {
    dat$NPX_imputed <- npx
    dat$NPX_imputed_mean <- npx
    mice_pieces[[i]] <- dat
    next
  }

  sub <- dat %>%
    transmute(
      NPX = npx,
      Day = factor(Day),
      Age = as.numeric(Age),
      cbmi = as.numeric(cbmi),
      Gender = factor(Gender)
    )

  obs <- !is.na(sub$NPX)
  usable <- c("Day", "Age", "cbmi", "Gender")
  usable <- usable[vapply(usable, function(p) {
    x <- sub[[p]][obs]
    length(unique(x[!is.na(x)])) >= 2L
  }, logical(1))]

  if (length(usable) == 0L) {
    dat$NPX_imputed <- npx
    dat$NPX_imputed_mean <- npx
    mice_pieces[[i]] <- dat
    next
  }

  sub <- sub[, c("NPX", usable), drop = FALSE]

  mice_pieces[[i]] <- tryCatch({
    ini <- mice(sub, maxit = 0, printFlag = FALSE)
    meth <- ini$method
    meth[] <- ""
    meth["NPX"] <- "pmm"
    pred <- ini$predictorMatrix
    pred[] <- 0
    pred["NPX", usable] <- 1
    imp <- mice(
      sub,
      m = 5,
      maxit = 5,
      method = meth,
      predictorMatrix = pred,
      seed = 500 + i,
      printFlag = FALSE
    )
    draws <- sapply(seq_len(5), function(j) complete(imp, j)$NPX)
    dat$NPX_imputed <- draws[, 1]
    dat$NPX_imputed_mean <- rowMeans(draws)
    dat
  }, error = function(e) {
    warning("mice failed for ", assay, ": ", conditionMessage(e))
    dat$NPX_imputed <- npx
    dat$NPX_imputed_mean <- npx
    dat
  })
}

data_imputed <- bind_rows(mice_pieces)
sum(is.na(data3$NPX) & !is.na(data_imputed$NPX_imputed))

# =============================================================================
# 8. Soft-winsorize NPX_imputed per Assay (Tukey fences, log1p pull-in)
# =============================================================================
data_wins <- data_imputed %>%
  group_by(Assay) %>%
  mutate(
    q1 = quantile(NPX_imputed, 0.25, na.rm = TRUE),
    q3 = quantile(NPX_imputed, 0.75, na.rm = TRUE),
    iqr = q3 - q1,
    lower = q1 - 1.5 * iqr,
    upper = q3 + 1.5 * iqr,
    fence = case_when(
      is.na(NPX_imputed) ~ NA_real_,
      NPX_imputed < lower ~ lower,
      NPX_imputed > upper ~ upper,
      TRUE ~ NA_real_
    ),
    NPX_wins = if_else(
      is.na(fence),
      NPX_imputed,
      fence + sign(NPX_imputed - fence) * log1p(abs(NPX_imputed - fence))
    )
  ) %>%
  ungroup() %>%
  select(-q1, -q3, -iqr, -lower, -upper, -fence)

sum(data_wins$NPX_imputed != data_wins$NPX_wins, na.rm = TRUE)

# =============================================================================
# 9. Inverse-rank normalization per Assay
# =============================================================================
set.seed(500)
combined <- data_wins %>%
  group_by(Assay) %>%
  mutate(
    NPX_invnorm = orderNorm(jitter(NPX_imputed, factor = 1e-6))$x.t
  ) %>%
  ungroup() %>%
  relocate(
    PatientID, Day, SampleID, Assay,
    NPX, NPX_imputed, NPX_imputed_mean, NPX_wins, NPX_invnorm
  )

combined
n_distinct(combined$Assay)
n_distinct(combined$PatientID)

count_data <- combined %>% distinct(PatientID, Day, Allocation)
table(count_data$Day, count_data$Allocation)

# =============================================================================
# 10. Write analysis table
# =============================================================================
write.table(
  combined,
  file = "data/Combined_dataset_mice.csv",
  col.names = TRUE,
  row.names = FALSE,
  quote = FALSE,
  sep = ";"
)
