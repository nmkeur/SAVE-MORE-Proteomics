library(tidyverse)
library(survival)
library(survminer)
library(broom)
library(patchwork)
library(pROC)

col_placebo <- "#1f77b4"
col_anakinra <- "#B22222"
time_xlim <- c(0, 14)
time_by <- 2

tidy_cox <- function(fit, model_name = NA_character_) {
  broom::tidy(fit, exponentiate = TRUE, conf.int = TRUE) %>%
    mutate(model = model_name)
}

protein_label <- function(x) {
  dplyr::recode(
    as.character(x),
    "IL6" = "IL-6",
    "PD.L1" = "PD-L1",
    "PD-L1" = "PD-L1",
    "TNFB" = "TNF-B",
    "MCP.1" = "MCP-1",
    "MCP.3" = "MCP-3",
    "MCP.4" = "MCP-4",
    "LIF.R" = "LIF-R",
    "LIF-R" = "LIF-R",
    "IL10" = "IL-10",
    .default = as.character(x)
  )
}

prepare_surv_day1 <- function(npx, prot_raw) {
  day1 <- npx %>%
    mutate(Day = as.character(Day)) %>%
    filter(Day %in% c("1", "01")) %>%
    group_by(Assay) %>%
    mutate(NPX_scaled = as.numeric(scale(NPX_invnorm))) %>%
    ungroup()

  long <- day1 %>%
    dplyr::select(
      PatientID, Assay, NPX_scaled,
      Age, cbmi, Gender, DEXcorrected, Allocation,
      srf14, timesrf14
    ) %>%
    filter(is.na(timesrf14) | timesrf14 != 1) %>%
    mutate(
      status = if_else(srf14 == "Yes", 1L, 0L),
      time = timesrf14,
      Gender = factor(Gender),
      DEXcorrected = factor(DEXcorrected),
      Allocation = factor(Allocation, levels = c("Placebo", "Anakinra"))
    ) %>%
    drop_na(NPX_scaled, time, status, Age, cbmi, Gender, DEXcorrected, Allocation)

  wide <- long %>%
    filter(Assay %in% prot_raw) %>%
    dplyr::select(
      PatientID, Assay, NPX_scaled,
      Age, cbmi, Gender, DEXcorrected, Allocation,
      time, status
    ) %>%
    pivot_wider(names_from = Assay, values_from = NPX_scaled) %>%
    distinct()

  colnames(wide) <- make.names(colnames(wide))
  prot_vars <- make.names(prot_raw)
  wide <- wide %>%
    drop_na(all_of(c(prot_vars, "time", "status", "Age", "cbmi", "Gender", "DEXcorrected", "Allocation")))

  list(long = long, wide = wide, prot_vars = prot_vars)
}

cox_univariate <- function(data, var) {
  fml <- as.formula(paste("Surv(time, status) ~", var))
  coxph(fml, data = data, na.action = na.exclude)
}

cox_proteins_by_arm <- function(data, prot_vars, arm) {
  d <- dplyr::filter(data, Allocation == arm)
  purrr::map_dfr(prot_vars, function(v) {
    tidy_cox(cox_univariate(d, v), arm) %>%
      mutate(protein = protein_label(v), term = v)
  })
}

cox_multivariate_arm <- function(data, prot_vars, arm) {
  d <- dplyr::filter(data, Allocation == arm)
  rhs <- paste(c(prot_vars, "Age", "cbmi", "Gender", "DEXcorrected"), collapse = " + ")
  fml <- as.formula(paste("Surv(time, status) ~", rhs))
  coxph(fml, data = d, na.action = na.exclude)
}

cox_protein_allocation_int <- function(data, prot_var) {
  rhs <- paste0(prot_var, " * Allocation + Age + cbmi + Gender + DEXcorrected")
  coxph(as.formula(paste("Surv(time, status) ~", rhs)), data = data, na.action = na.exclude)
}

cox_joint_allocation_int <- function(data, prot_vars) {
  rhs <- paste0(
    "(", paste(prot_vars, collapse = " + "), ")*Allocation + Age + cbmi + Gender + DEXcorrected"
  )
  coxph(as.formula(paste("Surv(time, status) ~", rhs)), data = data, na.action = na.exclude)
}

cox_joint_allocation_add <- function(data, prot_vars) {
  rhs <- paste(c(prot_vars, "Allocation", "Age", "cbmi", "Gender", "DEXcorrected"), collapse = " + ")
  coxph(as.formula(paste("Surv(time, status) ~", rhs)), data = data, na.action = na.exclude)
}

cox_clinical_allocation <- function(data) {
  coxph(
    Surv(time, status) ~ Allocation + Age + cbmi + Gender + DEXcorrected,
    data = data,
    na.action = na.exclude
  )
}

tidy_allocation_int <- function(fit, prot_vars, model_name) {
  tidy_cox(fit, model_name) %>%
    dplyr::filter(grepl(":Allocation", term)) %>%
    mutate(
      protein_var = sub(":Allocation.*", "", term),
      protein = protein_label(protein_var)
    ) %>%
    dplyr::filter(protein_var %in% prot_vars) %>%
    mutate(
      protein = factor(protein, levels = protein_label(prot_vars)),
      hr_ci = fmt_ci(estimate, conf.low, conf.high)
    ) %>%
    arrange(protein)
}

reference_profile <- function(data, prot_vars) {
  prot_profile <- data %>%
    summarise(across(all_of(prot_vars), ~ median(.x, na.rm = TRUE)))
  num_profile <- tibble(
    Age = mean(data$Age, na.rm = TRUE),
    cbmi = mean(data$cbmi, na.rm = TRUE)
  )
  fac_profile <- tibble(
    Gender = factor(
      names(which.max(table(data$Gender))),
      levels = levels(data$Gender)
    ),
    DEXcorrected = factor(
      names(which.max(table(data$DEXcorrected))),
      levels = levels(data$DEXcorrected)
    )
  )
  bind_cols(prot_profile, num_profile, fac_profile) %>%
    tidyr::crossing(
      Allocation = factor(c("Placebo", "Anakinra"), levels = levels(data$Allocation))
    ) %>%
    relocate(Allocation)
}

km_theme <- function() {
  theme_classic() +
    theme(
      legend.position = "top",
      legend.direction = "horizontal"
    )
}

event_ylim <- function(...) {
  y_max <- max(vapply(list(...), function(fit) {
    s <- fit$surv
    lo <- if (!is.null(fit$lower)) fit$lower else s
    1 - min(c(s, lo), na.rm = TRUE)
  }, numeric(1)), na.rm = TRUE)
  c(0, min(1, ceiling(y_max * 20) / 20))
}

plot_km_unadjusted <- function(fit, data, ylim = c(0, 1)) {
  ggsurvplot(
    fit,
    data = data,
    fun = "event",
    conf.int = TRUE,
    pval = TRUE,
    risk.table = TRUE,
    legend.title = "Treatment",
    legend.labs = c("Placebo", "Anakinra"),
    palette = c(col_placebo, col_anakinra),
    xlab = "Time to SRF (days)",
    ylab = "Cumulative incidence of SRF",
    xlim = time_xlim,
    ylim = ylim,
    break.time.by = time_by,
    ggtheme = km_theme()
  )
}

plot_km_adjusted <- function(fit, newdat, ylim = c(0, 1)) {
  ggsurvplot(
    fit,
    data = newdat,
    fun = "event",
    conf.int = TRUE,
    pval = FALSE,
    risk.table = FALSE,
    legend.title = "Treatment",
    legend.labs = c("Placebo", "Anakinra"),
    palette = c(col_placebo, col_anakinra),
    xlab = "Time to SRF (days)",
    ylab = "Cumulative incidence of SRF",
    xlim = time_xlim,
    ylim = ylim,
    break.time.by = time_by,
    ggtheme = km_theme()
  )
}

allocation_p_label <- function(cox_fit) {
  p <- broom::tidy(cox_fit) %>%
    dplyr::filter(grepl("^Allocation", term)) %>%
    dplyr::pull(p.value)
  paste0("p = ", ifelse(length(p) == 0 || is.na(p), "NA", signif(p[[1]], 2)))
}

allocation_hr_p_label <- function(cox_fit) {
  row <- broom::tidy(cox_fit, exponentiate = TRUE, conf.int = TRUE) %>%
    dplyr::filter(grepl("^Allocation", term))
  if (nrow(row) == 0) {
    return("HR = NA\np = NA")
  }
  paste0(
    "HR = ",
    fmt_ci(row$estimate[[1]], row$conf.low[[1]], row$conf.high[[1]]),
    "\np = ",
    signif(row$p.value[[1]], 2)
  )
}

add_surv_p <- function(ggsurv, label, ylim) {
  ggsurv$plot <- ggsurv$plot +
    annotate(
      "text",
      x = 1.2,
      y = ylim[2] * 0.92,
      label = label,
      hjust = 0,
      vjust = 1,
      size = 4,
      lineheight = 0.95
    )
  ggsurv
}

to_mid_dot <- function(x) gsub("(?<=\\d)\\.(?=\\d)", "·", x, perl = TRUE)

fmt_ci <- function(hr, lo, hi) {
  to_mid_dot(sprintf("%.2f (%.2f–%.2f)", as.numeric(hr), as.numeric(lo), as.numeric(hi)))
}

fmt_p <- function(p) {
  p <- as.numeric(p)
  ifelse(is.na(p), "p=NA", to_mid_dot(sprintf("p=%.2e", p)))
}

plot_arm_forest <- function(forest_df, title = "Forest Plot of Top Proteins Across Models") {
  COL_ANI <- col_anakinra
  COL_PLA <- col_placebo
  BASE_FAMILY <- "Arial"
  BASE_PT <- 12
  TXT_SIZE <- BASE_PT / ggplot2::.pt
  PDGE <- position_dodge(width = 0.70)

  forest_all2 <- forest_df %>%
    mutate(
      model = factor(model, levels = c("Anakinra", "Placebo")),
      HR = estimate,
      CI_low = conf.low,
      CI_high = conf.high
    )

  order_vec <- forest_all2 %>%
    filter(model == "Placebo") %>%
    arrange(desc(HR)) %>%
    pull(protein) %>%
    as.character()

  forest_all2 <- forest_all2 %>%
    mutate(protein = factor(protein, levels = order_vec))

  annot_long <- forest_all2 %>%
    transmute(
      protein = factor(protein, levels = levels(forest_all2$protein)),
      model,
      hrci = fmt_ci(HR, CI_low, CI_high),
      pval = fmt_p(p.value)
    )

  n_prot <- length(levels(forest_all2$protein))

  table_left <- ggplot() +
    geom_text(
      data = dplyr::distinct(annot_long, protein),
      aes(x = 0.02, y = protein, label = protein),
      hjust = 0, family = BASE_FAMILY, size = TXT_SIZE
    ) +
    geom_text(
      data = annot_long,
      aes(x = 0.72, y = protein, label = hrci, color = model),
      position = PDGE, hjust = 0.5, family = BASE_FAMILY, size = TXT_SIZE
    ) +
    scale_y_discrete(limits = rev(levels(forest_all2$protein)), expand = expansion(add = 1)) +
    scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
    scale_color_manual(values = c(Anakinra = COL_ANI, Placebo = COL_PLA)) +
    coord_cartesian(clip = "off") +
    theme_void() +
    theme(
      legend.position = "none",
      plot.margin = margin(18, 4, 5, 8)
    ) +
    annotate("text", x = 0.02, y = n_prot + 0.85, label = "Protein",
             hjust = 0, vjust = 0, fontface = "bold", family = BASE_FAMILY, size = TXT_SIZE) +
    annotate("text", x = 0.72, y = n_prot + 0.85, label = "HR (95% CI)",
             hjust = 0.5, vjust = 0, fontface = "bold", family = BASE_FAMILY, size = TXT_SIZE)

  forest_panel <- ggplot(forest_all2, aes(x = HR, y = protein, color = model)) +
    geom_point(position = PDGE, size = 2.4) +
    geom_errorbarh(aes(xmin = CI_low, xmax = CI_high), position = PDGE, height = 0.18) +
    geom_vline(xintercept = 1, linetype = "dashed") +
    scale_y_discrete(limits = rev(levels(forest_all2$protein)), expand = expansion(add = 1)) +
    scale_x_log10(
      "Hazard Ratio (log scale)",
      labels = function(x) to_mid_dot(format(x, trim = TRUE, scientific = FALSE))
    ) +
    scale_color_manual(values = c(Anakinra = COL_ANI, Placebo = COL_PLA), name = NULL) +
    theme_classic(base_family = BASE_FAMILY, base_size = BASE_PT) +
    theme(
      legend.position = "top",
      axis.title.y = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank()
    )

  table_right <- ggplot(annot_long, aes(y = protein, x = 0.95, label = pval, color = model)) +
    geom_text(position = position_dodge(width = 0.75), hjust = 1, family = BASE_FAMILY, size = TXT_SIZE) +
    scale_y_discrete(limits = rev(levels(annot_long$protein)), expand = expansion(add = 1)) +
    scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
    coord_cartesian(clip = "off") +
    scale_color_manual(values = c(Anakinra = COL_ANI, Placebo = COL_PLA)) +
    theme_void() +
    theme(legend.position = "none", plot.margin = margin(18, 10, 5, 4)) +
    annotate("text", x = 0.95, y = n_prot + 0.85, label = "p-value",
             hjust = 1, vjust = 0, fontface = "bold", family = BASE_FAMILY, size = TXT_SIZE)

  out <- table_left + forest_panel + table_right +
    plot_layout(widths = c(3.2, 3.6, 1.3))
  if (!is.null(title) && !identical(title, "")) {
    out <- out + plot_annotation(
      title = title,
      theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", family = BASE_FAMILY))
    )
  }
  out
}

format_mid <- function(x, digits = 2) {
  gsub("\\.", "·", format(round(x, digits), nsmall = digits))
}

prepare_auc_day1 <- function(npx, lab_raw = c("CRP", "suPAR")) {
  day1 <- npx %>%
    mutate(Day = as.character(Day)) %>%
    dplyr::filter(Day %in% c("1", "01")) %>%
    dplyr::filter(is.na(timesrf14) | timesrf14 != 1)

  prot_all <- sort(unique(as.character(day1$Assay)))
  missing_lab <- setdiff(lab_raw, names(day1))
  if (length(missing_lab) > 0) {
    stop("Clinical lab columns missing: ", paste(missing_lab, collapse = ", "))
  }

  wide <- day1 %>%
    dplyr::select(PatientID, SampleID, Allocation, Assay, NPX_invnorm, srf14, all_of(lab_raw)) %>%
    mutate(across(all_of(lab_raw), as.numeric)) %>%
    pivot_wider(
      id_cols = c(PatientID, SampleID, Allocation, srf14, all_of(lab_raw)),
      names_from = Assay,
      values_from = NPX_invnorm
    ) %>%
    mutate(
      StatusBin = dplyr::case_when(
        srf14 == "Yes" ~ 1,
        srf14 == "No" ~ 0,
        TRUE ~ NA_real_
      ),
      Allocation = factor(Allocation, levels = c("Placebo", "Anakinra"))
    )

  colnames(wide) <- make.names(colnames(wide))
  prot_vars_all <- make.names(prot_all)
  names(prot_vars_all) <- prot_all
  lab_vars <- make.names(lab_raw)
  names(lab_vars) <- lab_raw
  missing_prot <- setdiff(unname(prot_vars_all), names(wide))
  if (length(missing_prot) > 0) {
    stop("Assay columns missing after pivot: ", paste(missing_prot, collapse = ", "))
  }

  list(
    wide = wide,
    prot_raw = prot_all,
    prot_vars = prot_vars_all,
    lab_raw = lab_raw,
    lab_vars = lab_vars
  )
}

youden_coords <- function(roc_obj) {
  opt <- coords(
    roc_obj,
    x = "best",
    best.method = "youden",
    ret = c("threshold", "sensitivity", "specificity"),
    transpose = FALSE
  )
  list(
    threshold = as.numeric(opt$threshold[[1]]),
    sensitivity = as.numeric(opt$sensitivity[[1]]),
    specificity = as.numeric(opt$specificity[[1]])
  )
}

stratified_folds <- function(grp, k) {
  idx <- seq_along(grp)
  folds <- integer(length(grp))
  for (g in unique(grp)) {
    gi <- idx[grp == g]
    gi <- sample(gi)
    folds[gi] <- rep(seq_len(k), length.out = length(gi))
  }
  folds
}

youden_cv <- function(vals, grp, direction, k = 5, seed = 1) {
  ok <- complete.cases(vals, grp)
  vals <- vals[ok]
  grp <- grp[ok]
  n0 <- sum(grp == 0)
  n1 <- sum(grp == 1)
  k_use <- min(k, n0, n1)
  if (k_use < 2) {
    return(list(
      threshold = NA_real_,
      sensitivity = NA_real_,
      specificity = NA_real_,
      n_folds = NA_integer_
    ))
  }

  set.seed(seed)
  folds <- stratified_folds(grp, k_use)
  thr <- sens <- spec <- rep(NA_real_, k_use)

  for (f in seq_len(k_use)) {
    train <- folds != f
    test <- folds == f
    if (length(unique(grp[train])) < 2 || length(unique(grp[test])) < 2) next
    roc_train <- roc(
      response = grp[train],
      predictor = vals[train],
      levels = c(0, 1),
      direction = direction,
      quiet = TRUE
    )
    cut <- youden_coords(roc_train)$threshold
    roc_test <- roc(
      response = grp[test],
      predictor = vals[test],
      levels = c(0, 1),
      direction = direction,
      quiet = TRUE
    )
    held <- coords(
      roc_test,
      x = cut,
      input = "threshold",
      ret = c("sensitivity", "specificity"),
      transpose = FALSE
    )
    thr[f] <- cut
    sens[f] <- as.numeric(held$sensitivity[[1]])
    spec[f] <- as.numeric(held$specificity[[1]])
  }

  list(
    threshold = mean(thr, na.rm = TRUE),
    sensitivity = mean(sens, na.rm = TRUE),
    specificity = mean(spec, na.rm = TRUE),
    n_folds = sum(!is.na(sens))
  )
}

roc_one_protein <- function(dat, prot_var, prot_name, arm, source = "Olink", k = 5, seed = 1) {
  out <- tibble(
    Protein = prot_name,
    Protein_label = protein_label(prot_name),
    Source = source,
    Arm = arm,
    AUC = NA_real_,
    CI_lower = NA_real_,
    CI_upper = NA_real_,
    Direction = NA_character_,
    Threshold_in_sample = NA_real_,
    Sensitivity_in_sample = NA_real_,
    Specificity_in_sample = NA_real_,
    Threshold_cv = NA_real_,
    Sensitivity_cv = NA_real_,
    Specificity_cv = NA_real_,
    n_cv_folds = NA_integer_,
    n = NA_integer_,
    n_srf = NA_integer_
  )

  vals <- dat[[prot_var]]
  grp <- dat$StatusBin
  ok <- complete.cases(vals, grp)
  vals <- vals[ok]
  grp <- grp[ok]
  out$n <- length(grp)
  out$n_srf <- sum(grp == 1)

  if (length(unique(grp)) < 2) return(out)
  if (length(unique(vals)) < 2) return(out)

  roc_obj <- roc(
    response = grp,
    predictor = vals,
    levels = c(0, 1),
    direction = "auto",
    quiet = TRUE
  )
  ci_auc <- ci.auc(roc_obj)
  youden <- youden_coords(roc_obj)
  cv <- youden_cv(vals, grp, direction = roc_obj$direction, k = k, seed = seed)

  out$AUC <- as.numeric(auc(roc_obj))
  out$CI_lower <- as.numeric(ci_auc[1])
  out$CI_upper <- as.numeric(ci_auc[3])
  out$Direction <- as.character(roc_obj$direction)
  out$Threshold_in_sample <- youden$threshold
  out$Sensitivity_in_sample <- youden$sensitivity
  out$Specificity_in_sample <- youden$specificity
  out$Threshold_cv <- cv$threshold
  out$Sensitivity_cv <- cv$sensitivity
  out$Specificity_cv <- cv$specificity
  out$n_cv_folds <- cv$n_folds
  out
}

compute_roc <- function(dat, prot_vars, arm, source = "Olink", k = 5, seed = 1) {
  purrr::map_dfr(seq_along(prot_vars), function(i) {
    roc_one_protein(
      dat,
      prot_var = unname(prot_vars[[i]]),
      prot_name = names(prot_vars)[[i]],
      arm = arm,
      source = source,
      k = k,
      seed = seed
    )
  })
}

plot_auc_top <- function(all_results, n_top = 5) {
  PT10 <- 10 / ggplot2::.pt
  top_by_arm <- all_results %>%
    dplyr::filter(Arm %in% c("Placebo", "Anakinra"), !is.na(AUC)) %>%
    group_by(Arm) %>%
    slice_max(order_by = AUC, n = n_top, with_ties = FALSE) %>%
    ungroup()

  all_results %>%
    dplyr::filter(
      Arm %in% c("Placebo", "Anakinra"),
      Protein %in% top_by_arm$Protein
    ) %>%
    mutate(Protein_label = forcats::fct_reorder(Protein_label, AUC, .fun = max)) %>%
    ggplot(aes(y = Protein_label, x = AUC, fill = Arm)) +
    geom_col(position = position_dodge(width = 0.9), color = "black", alpha = 0.75) +
    geom_text(
      aes(label = format_mid(AUC, 2)),
      position = position_dodge(width = 0.9),
      hjust = 1.2,
      size = PT10
    ) +
    labs(
      title = "Top markers by AUC per arm",
      subtitle = "Day 1, SRF-free at baseline; Olink plus clinical CRP and suPAR; direction auto (AUC >= 0.5)",
      x = "AUC",
      y = NULL
    ) +
    scale_fill_manual(values = c(Anakinra = col_anakinra, Placebo = col_placebo)) +
    coord_cartesian(xlim = c(0, 1)) +
    scale_x_continuous(labels = function(x) format_mid(x, 2)) +
    theme_classic() +
    theme(
      legend.position = "top",
      axis.text = element_text(size = 10)
    )
}
