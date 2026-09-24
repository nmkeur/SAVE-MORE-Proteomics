library(dplyr)
library(tidyr)
library(stringr)
library(tibble)

# Parse Olink SampleID values such as "SM AA 027 07" into PatientID + Day.
# Controls (e.g. "SC 1-1") do not match and return NA for PatientID/Day.
parse_olink_sample_id <- function(sample_id) {
  clean <- str_replace_all(as.character(sample_id), "\\s+", "")
  day <- str_extract(clean, "(01|04|07)$")
  tibble(
    SampleID = as.character(sample_id),
    SampleID_clean = clean,
    Day = day,
    PatientID = if_else(!is.na(day), str_remove(clean, "(01|04|07)$"), NA_character_)
  )
}

# Append one QC census row. `df` may be clinical-only (no SampleID/Assay).
qc_census <- function(step, df, note = "") {
  tibble(
    step = step,
    n_rows = nrow(df),
    n_samples = if ("SampleID" %in% names(df)) dplyr::n_distinct(df$SampleID) else NA_integer_,
    n_patients = if ("PatientID" %in% names(df)) dplyr::n_distinct(df$PatientID) else NA_integer_,
    n_assays = if ("Assay" %in% names(df)) dplyr::n_distinct(df$Assay) else NA_integer_,
    note = note
  )
}

# Median-impute NPX within Assay x Day only (no Allocation / srf14).
impute_npx_by_assay_day <- function(df, value = "NPX") {
  df %>%
    group_by(Assay, Day) %>%
    mutate(
      NPX_imputed = if_else(
        is.na(.data[[value]]),
        median(.data[[value]], na.rm = TRUE),
        .data[[value]]
      )
    ) %>%
    ungroup()
}

# PMM-impute NPX separately for each Assay.
# Predictors are covariates only (not Allocation / srf14). Observed NPX is kept.
# NPX_imputed is complete(..., action = 1), matching the old Combined_dataset_v44
# choice of a single imputed draw. NPX_imputed_mean averages the m draws.
impute_npx_mice <- function(df,
                            m = 5,
                            maxit = 5,
                            seed = 500,
                            predictors = c("Day", "Age", "cbmi", "Gender"),
                            value = "NPX") {
  if (!requireNamespace("mice", quietly = TRUE)) {
    stop("Package 'mice' is required for impute_npx_mice().")
  }

  missing_pred <- setdiff(c(value, "Assay", predictors), names(df))
  if (length(missing_pred) > 0) {
    stop("Missing columns for mice: ", paste(missing_pred, collapse = ", "))
  }

  df$.row_id <- seq_len(nrow(df))
  assays <- unique(df$Assay)
  pieces <- vector("list", length(assays))

  for (i in seq_along(assays)) {
    assay <- assays[[i]]
    dat <- df[df$Assay == assay, , drop = FALSE]
    npx <- as.numeric(dat[[value]])
    n_miss <- sum(is.na(npx))

    if (n_miss == 0L || sum(!is.na(npx)) < 5L) {
      dat$NPX_imputed <- npx
      dat$NPX_imputed_mean <- npx
      pieces[[i]] <- dat
      next
    }

    sub <- dat[, c(value, predictors), drop = FALSE]
    names(sub)[names(sub) == value] <- "NPX"
    sub$NPX <- npx
    if ("Day" %in% names(sub)) sub$Day <- factor(sub$Day)
    if ("Gender" %in% names(sub)) sub$Gender <- factor(sub$Gender)
    if ("Age" %in% names(sub)) sub$Age <- as.numeric(sub$Age)
    if ("cbmi" %in% names(sub)) sub$cbmi <- as.numeric(sub$cbmi)

    usable <- predictors[predictors %in% names(sub)]
    # Drop predictors that are all NA or have no variation among observed NPX rows.
    obs <- !is.na(sub$NPX)
    keep <- vapply(usable, function(p) {
      x <- sub[[p]][obs]
      x <- x[!is.na(x)]
      length(unique(x)) >= 2L
    }, logical(1))
    usable <- usable[keep]

    if (length(usable) == 0L) {
      dat$NPX_imputed <- npx
      dat$NPX_imputed_mean <- npx
      pieces[[i]] <- dat
      next
    }

    sub <- sub[, c("NPX", usable), drop = FALSE]
    assay_seed <- seed + i

    pieces[[i]] <- tryCatch({
      ini <- mice::mice(sub, maxit = 0, printFlag = FALSE)
      meth <- ini$method
      meth[] <- ""
      meth["NPX"] <- "pmm"
      pred <- ini$predictorMatrix
      pred[] <- 0
      pred["NPX", usable] <- 1
      imp <- mice::mice(
        sub,
        m = m,
        maxit = maxit,
        method = meth,
        predictorMatrix = pred,
        seed = assay_seed,
        printFlag = FALSE
      )
      draws <- vapply(
        seq_len(m),
        function(j) mice::complete(imp, j)$NPX,
        numeric(nrow(sub))
      )
      dat$NPX_imputed <- draws[, 1]
      dat$NPX_imputed_mean <- rowMeans(draws)
      dat
    }, error = function(e) {
      warning("mice failed for assay ", assay, ": ", conditionMessage(e), call. = FALSE)
      dat$NPX_imputed <- npx
      dat$NPX_imputed_mean <- npx
      dat
    })
  }

  bind_rows(pieces) %>%
    arrange(.row_id) %>%
    select(-.row_id)
}

# Soft-winsorize a numeric column per Assay: values outside Tukey fences
# are pulled toward the fence by log1p(|distance|), not dropped.
soft_winsor <- function(df, value = "NPX_imputed", m = 1.5, group = "Assay") {
  df %>%
    group_by(.data[[group]]) %>%
    mutate(
      .q1 = quantile(.data[[value]], 0.25, na.rm = TRUE),
      .q3 = quantile(.data[[value]], 0.75, na.rm = TRUE),
      .iqr = .q3 - .q1,
      .lower = .q1 - m * .iqr,
      .upper = .q3 + m * .iqr,
      .x = .data[[value]],
      .fence = case_when(
        is.na(.x) ~ NA_real_,
        .x < .lower ~ .lower,
        .x > .upper ~ .upper,
        TRUE ~ NA_real_
      ),
      NPX_wins = if_else(
        is.na(.fence),
        .x,
        .fence + sign(.x - .fence) * log1p(abs(.x - .fence))
      )
    ) %>%
    ungroup() %>%
    select(-.q1, -.q3, -.iqr, -.lower, -.upper, -.x, -.fence)
}

# Per-assay ordered quantile (inverse-rank) normalization of NPX_wins.
# Tiny jitter breaks ties so orderNorm ranks are unique (same as the v4 snippet).
invnorm_npx <- function(df, value = "NPX_wins", jitter_factor = 1e-6) {
  if (!requireNamespace("bestNormalize", quietly = TRUE)) {
    stop("Package 'bestNormalize' is required for invnorm_npx().")
  }
  df %>%
    ungroup() %>%
    group_by(Assay) %>%
    mutate(
      .j = jitter(.data[[value]], factor = jitter_factor),
      NPX_invnorm = bestNormalize::orderNorm(.j)$x.t
    ) %>%
    ungroup() %>%
    select(-.j)
}
