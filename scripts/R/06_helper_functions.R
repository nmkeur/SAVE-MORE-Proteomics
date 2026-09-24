library(tidyverse)
library(ggplot2)
library(ggpubr)
library(Hmisc)

col_placebo <- "#1f77b4"
col_anakinra <- "#B22222"

clin_raw <- c("Gender", "cbmi", "Age", "suPAR", "ferritin", "IL6", "CRP")
clin_labels <- c("Sex", "BMI", "Age", "suPAR", "Ferritin", "IL6", "CRP")
names(clin_labels) <- clin_raw
clin_plot_order <- c("Age", "BMI", "Sex", "suPAR", "Ferritin", "IL6", "CRP")

add_onset <- function(npx) {
  npx %>%
    mutate(
      Day_num = as.numeric(as.character(Day)),
      timesrf14 = as.numeric(timesrf14),
      Allocation = factor(Allocation, levels = c("Placebo", "Anakinra")),
      onset = dplyr::case_when(
        srf14 == "No" ~ "neverSRF",
        Day_num < timesrf14 ~ "preSRF",
        Day_num >= timesrf14 ~ "activeSRF",
        TRUE ~ NA_character_
      )
    )
}

make_delta_df <- function(npx, npx_col = "NPX_invnorm") {
  df <- add_onset(npx) %>%
    dplyr::filter(
      timesrf14 != 1,
      onset %in% c("preSRF", "activeSRF"),
      !is.na(.data[[npx_col]])
    ) %>%
    mutate(NPX_use = .data[[npx_col]])

  df %>%
    group_by(PatientID, Assay) %>%
    arrange(Day_num, .by_group = TRUE) %>%
    summarise(
      Allocation = dplyr::first(Allocation),
      NPX_pre = dplyr::last(NPX_use[onset == "preSRF"]),
      NPX_act = dplyr::first(NPX_use[onset == "activeSRF"]),
      .groups = "drop"
    ) %>%
    dplyr::filter(!is.na(NPX_pre), !is.na(NPX_act)) %>%
    mutate(delta = NPX_act - NPX_pre)
}

wilcox_delta_by_assay <- function(delta_df) {
  delta_df %>%
    group_by(Assay) %>%
    group_modify(function(d, key) {
      if (n_distinct(d$Allocation) < 2L) {
        return(tibble(
          n_placebo = sum(d$Allocation == "Placebo"),
          n_anakinra = sum(d$Allocation == "Anakinra"),
          median_placebo = median(d$delta[d$Allocation == "Placebo"], na.rm = TRUE),
          median_anakinra = median(d$delta[d$Allocation == "Anakinra"], na.rm = TRUE),
          p.value = NA_real_
        ))
      }
      wt <- wilcox.test(delta ~ Allocation, data = d, exact = FALSE)
      tibble(
        n_placebo = sum(d$Allocation == "Placebo"),
        n_anakinra = sum(d$Allocation == "Anakinra"),
        median_placebo = median(d$delta[d$Allocation == "Placebo"], na.rm = TRUE),
        median_anakinra = median(d$delta[d$Allocation == "Anakinra"], na.rm = TRUE),
        p.value = unname(wt$p.value)
      )
    }) %>%
    ungroup() %>%
    mutate(adj.p.value = p.adjust(p.value, method = "BH")) %>%
    arrange(p.value)
}

fmt_p <- function(p) {
  p <- as.numeric(p)
  lab <- ifelse(
    is.na(p),
    "NA",
    ifelse(p < 0.001, sprintf("%.2e", p), sprintf("%.3f", p))
  )
  gsub("\\.", "·", lab)
}

fmt_num <- function(x, digits = 2) {
  gsub("\\.", "·", sprintf(paste0("%+.", digits, "f"), as.numeric(x)))
}

plot_delta_violin <- function(delta_df, assays, tests) {
  d <- delta_df %>%
    dplyr::filter(Assay %in% assays) %>%
    mutate(Assay = factor(Assay, levels = assays))

  y_top <- d %>%
    group_by(Assay) %>%
    summarise(y = max(delta, na.rm = TRUE), .groups = "drop")

  p_lab <- tests %>%
    dplyr::filter(Assay %in% assays) %>%
    mutate(Assay = factor(Assay, levels = assays)) %>%
    left_join(y_top, by = "Assay") %>%
    mutate(
      diff_med = median_anakinra - median_placebo,
      label = paste0("p=", fmt_p(p.value), "\ndiff=", fmt_num(diff_med)),
      x = 1.5
    )

  ggplot(d, aes(x = Allocation, y = delta, fill = Allocation)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    geom_violin(alpha = 0.3, width = 0.8, color = "black") +
    geom_boxplot(width = 0.2, outlier.size = 0.5) +
    stat_summary(
      fun = mean,
      geom = "point",
      shape = 23,
      size = 3,
      color = "black",
      fill = "white"
    ) +
    geom_text(
      data = p_lab,
      aes(x = x, y = y, label = label),
      inherit.aes = FALSE,
      vjust = 1.15,
      size = 3,
      lineheight = 0.95
    ) +
    scale_fill_manual(values = c(Anakinra = col_anakinra, Placebo = col_placebo)) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.08))) +
    facet_wrap(~Assay, scales = "free_y", ncol = 4) +
    labs(
      x = NULL,
      y = "Delta NPX (first activeSRF - last preSRF)"
    ) +
    theme_classic() +
    theme(
      strip.text = element_text(face = "bold"),
      legend.position = "none"
    )
}

format_lancet <- function(x, digits = 2) {
  gsub("\\.", "·", formatC(x, format = "f", digits = digits))
}

clinical_by_visit <- function(npx) {
  npx %>%
    mutate(Day = as.character(Day)) %>%
    dplyr::filter(Day %in% c("1", "4", "7", "01", "04", "07")) %>%
    mutate(
      Day = dplyr::recode(Day, "01" = "1", "04" = "4", "07" = "7"),
      Allocation = factor(Allocation, levels = c("Placebo", "Anakinra")),
      Gender_num = if_else(Gender == "Female", 0, 1)
    ) %>%
    dplyr::distinct(
      PatientID, Day, Allocation, SampleID,
      Gender_num, cbmi, Age, suPAR, ferritin, IL6, CRP
    ) %>%
    dplyr::rename(Gender = Gender_num)
}

patients_with_three_days <- function(clin) {
  clin %>%
    dplyr::distinct(PatientID, Day) %>%
    count(PatientID) %>%
    dplyr::filter(n == 3) %>%
    pull(PatientID)
}

spearman_r_p05 <- function(df, vars = clin_raw) {
  mat <- df %>%
    dplyr::select(all_of(vars)) %>%
    mutate(across(everything(), as.numeric)) %>%
    as.matrix()
  res <- Hmisc::rcorr(mat, type = "spearman")
  res$P[is.na(res$P)] <- 1
  r <- res$r
  r[res$P >= 0.05] <- 0
  diag(r) <- 1
  r
}

split_triangle_mat <- function(r_anakinra, r_placebo) {
  combined <- r_anakinra
  combined[lower.tri(combined)] <- r_placebo[lower.tri(r_placebo)]
  diag(combined) <- 1
  combined
}

clinical_split_mat <- function(clin, day) {
  day_df <- clin %>% dplyr::filter(Day == as.character(day))
  r_ana <- spearman_r_p05(dplyr::filter(day_df, Allocation == "Anakinra"))
  r_plc <- spearman_r_p05(dplyr::filter(day_df, Allocation == "Placebo"))
  split_triangle_mat(r_ana, r_plc)
}

label_clin_mat <- function(M) {
  dimnames(M) <- list(unname(clin_labels[rownames(M)]), unname(clin_labels[colnames(M)]))
  M
}

plot_tri_corr_from_mat <- function(
    M, title = NULL,
    palette_low = "#2C6BB0", palette_mid = "white", palette_high = "#C73E3A",
    base_family = "sans", base_size = 12,
    var_order = clin_plot_order,
    show_labels = TRUE, digits = 2,
    zero_tol = 1e-8,
    show_diag = TRUE,
    diag_label = FALSE,
    diag_fill = "#C73E3A",
    border_col = "black",
    border_lwd = 0.2
) {
  stopifnot(is.matrix(M))
  vars <- colnames(M)

  if (!is.null(var_order)) {
    stopifnot(all(var_order %in% vars))
    M <- M[var_order, var_order, drop = FALSE]
    vars <- var_order
  }

  df <- as_tibble(as.table(M), .name_repair = "minimal") |>
    set_names(c("Var1", "Var2", "r")) |>
    mutate(
      Var1 = factor(Var1, levels = vars),
      Var2 = factor(Var2, levels = vars),
      i = as.integer(Var1), j = as.integer(Var2),
      is_diag = i == j,
      r_main = ifelse(is_diag, NA_real_, r),
      label_ok = show_labels & !is.na(r_main) & abs(r_main) > zero_tol,
      lab = ifelse(label_ok, gsub("\\.", "·", sprintf(paste0("%.", digits, "f"), r_main)), NA_character_)
    )

  diag_df <- dplyr::filter(df, is_diag)

  ggplot(df, aes(Var2, Var1)) +
    geom_tile(aes(fill = r_main), colour = border_col, linewidth = border_lwd, na.rm = TRUE) +
    { if (show_diag) geom_tile(data = diag_df, fill = diag_fill, colour = border_col, linewidth = border_lwd) } +
    { if (show_labels) geom_text(data = dplyr::filter(df, label_ok),
                                 aes(label = lab), size = base_size / 2.845 * 0.6) } +
    { if (show_diag && diag_label)
      geom_text(
        data = diag_df,
        aes(label = gsub("\\.", "·", sprintf(paste0("%.", digits, "f"), 1))),
        colour = "grey25", size = base_size / 2.845 * 0.6
      ) } +
    scale_fill_gradient2(
      name = NULL,
      low = palette_low, mid = palette_mid, high = palette_high,
      limits = c(-1, 1), breaks = seq(-1, 1, 0.5), labels = format_lancet,
      na.value = NA
    ) +
    scale_y_discrete(limits = rev(vars)) +
    coord_fixed() +
    labs(x = NULL, y = NULL, title = title) +
    theme_minimal(base_family = base_family, base_size = base_size) +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.text.y = element_text(angle = 0, hjust = 1),
      plot.title = element_text(face = "bold", size = 10),
      legend.position = "right",
      axis.ticks = element_blank()
    )
}

olink_il6_by_day <- function(npx, day) {
  npx %>%
    mutate(Day = as.character(Day)) %>%
    dplyr::filter(Assay == "IL6", Day == as.character(day))
}

plot_olink_elisa_panel <- function(df, y, title, y_lab) {
  d <- df %>%
    mutate(
      NPX = as.numeric(NPX),
      ELISA_log2 = log2(as.numeric(.data[[y]]))
    ) %>%
    dplyr::filter(is.finite(NPX), is.finite(ELISA_log2))

  ggplot(d, aes(x = NPX, y = ELISA_log2)) +
    geom_point(alpha = 0.6) +
    geom_smooth(method = "lm", se = TRUE, color = "red") +
    stat_cor(method = "spearman", label.x.npc = "left", label.y.npc = "top") +
    labs(x = "IL-6 (Olink NPX)", y = y_lab, title = title) +
    theme_bw(base_size = 8)
}
