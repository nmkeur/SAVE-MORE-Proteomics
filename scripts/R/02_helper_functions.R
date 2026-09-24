library(tidyverse)
library(lme4)
library(lmerTest)
library(emmeans)
library(broom.mixed)

safe_lmer_threeway <- purrr::safely(function(df) {
  lmer(NPX_invnorm ~ Day * Allocation * srf14 + DEXcorrected + cbmi + Age + Gender + (1|PatientID), data = df)
})

safe_lmer_day_alloc <- purrr::safely(function(df) {
  lmer(NPX_invnorm ~ Day * Allocation + DEXcorrected + cbmi + Age + Gender + (1|PatientID), data = df)
})

safe_lmer_day_srf <- purrr::safely(function(df) {
  lmer(NPX_invnorm ~ Day * srf14 + DEXcorrected + cbmi + Age + Gender + (1|PatientID), data = df)
})

fdr_by_term <- function(anova_df) {
  anova_df %>%
    ungroup() %>%
    group_by(term) %>%
    mutate(adj.p.value = p.adjust(`Pr(>F)`, method = "fdr")) %>%
    ungroup()
}

fdr_by_contrast <- function(contr_df) {
  contr_df %>%
    ungroup() %>%
    group_by(contrast) %>%
    mutate(adj.p.value = p.adjust(p.value, method = "fdr")) %>%
    ungroup()
}

fit_by_assay <- function(data, safe_fun) {
  data %>%
    group_by(Assay) %>%
    nest() %>%
    mutate(
      model_obj = map(data, safe_fun),
      model = map(model_obj, "result"),
      error = map(model_obj, "error")
    )
}

failed_assays <- function(model_tbl) {
  model_tbl %>%
    filter(map_lgl(model, is.null)) %>%
    mutate(error = map_chr(error, ~ if (is.null(.x)) NA_character_ else conditionMessage(.x))) %>%
    dplyr::select(Assay, error)
}

keep_fitted <- function(model_tbl) {
  model_tbl %>%
    ungroup() %>%
    filter(map_lgl(model, ~ !is.null(.x)))
}

anova_by_assay <- function(model_tbl) {
  model_tbl %>%
    ungroup() %>%
    mutate(anova_tab = map(model, ~ as_tibble(anova(.x), rownames = "term"))) %>%
    dplyr::select(Assay, anova_tab) %>%
    unnest(anova_tab) %>%
    fdr_by_term()
}

hits_term <- function(anova_df, term_name, alpha = 0.05) {
  dplyr::filter(anova_df, term == term_name, adj.p.value < alpha)
}

unnest_emm <- function(model_tbl, emm_col) {
  model_tbl %>%
    mutate(emm_data = map(.data[[emm_col]], as.data.frame)) %>%
    dplyr::select(Assay, emm_data) %>%
    unnest(emm_data)
}

unnest_contrasts <- function(model_tbl, tidy_col) {
  model_tbl %>%
    dplyr::select(Assay, all_of(tidy_col)) %>%
    unnest(all_of(tidy_col)) %>%
    fdr_by_contrast()
}

add_emm_pairs <- function(model_tbl, specs) {
  model_tbl <- ungroup(model_tbl)
  for (nm in names(specs)) {
    emm_nm <- nm
    contr_nm <- paste0("contrasts_", nm)
    tidy_nm <- paste0("tidy_", nm)
    spec <- specs[[nm]]
    model_tbl[[emm_nm]] <- map(model_tbl$model, ~ emmeans(.x, spec))
    model_tbl[[contr_nm]] <- map(model_tbl[[emm_nm]], ~ pairs(.x, adjust = "none"))
    model_tbl[[tidy_nm]] <- map(model_tbl[[contr_nm]], broom::tidy)
  }
  model_tbl
}

day_levels <- c("1", "4", "7")

read_lmm_emm <- function(path) {
  read_tsv(path, show_col_types = FALSE) %>%
    mutate(Day = factor(as.character(as.integer(as.character(Day))), levels = day_levels))
}

star_label <- function(p) {
  dplyr::case_when(
    p < 0.001 ~ "***",
    p < 0.01  ~ "**",
    p < 0.05  ~ "*",
    TRUE ~ NA_character_
  )
}

fdr_among <- function(contrasts, assays) {
  contrasts %>%
    dplyr::filter(Assay %in% assays) %>%
    group_by(contrast) %>%
    mutate(adj.p.value = p.adjust(p.value, method = "fdr")) %>%
    ungroup()
}

# One row per Assay: the group with the largest Day7 − Day1 change (direction),
# then assays ordered by direction and desc(NPX_diff) within direction.
order_by_direction <- function(emm_df, dir_col = "group") {
  HH <- emm_df %>%
    group_by(Assay, .data[[dir_col]]) %>%
    mutate(
      Day1_NPX = emmean[Day == "1"][1],
      NPX_diff = emmean - Day1_NPX
    ) %>%
    filter(Day == "7") %>%
    dplyr::arrange(desc(NPX_diff)) %>%
    ungroup() %>%
    group_by(Assay) %>%
    slice_head(n = 1) %>%
    ungroup() %>%
    mutate(direction = .data[[dir_col]]) %>%
    dplyr::select(Assay, direction, NPX_diff)

  HH %>%
    dplyr::select(Assay, direction, NPX_diff) %>%
    group_by(direction) %>%
    dplyr::arrange(desc(NPX_diff), .by_group = TRUE) %>%
    pull(Assay)
}

apply_facet_order <- function(emm_df, brackets_df, dir_col = "group") {
  grp_str <- order_by_direction(emm_df, dir_col)
  emm_df$Assay <- factor(emm_df$Assay, levels = grp_str)
  if (nrow(brackets_df) > 0) {
    brackets_df$Assay <- factor(brackets_df$Assay, levels = grp_str)
  }
  list(emm = emm_df, brackets = brackets_df)
}

attach_offsets <- function(brackets_df, emm_df, offsets) {
  if (nrow(brackets_df) == 0) {
    return(mutate(brackets_df, x_offset = numeric(), x_pos = numeric()))
  }
  brackets_df %>%
    left_join(offsets, by = c("group1", "group2")) %>%
    mutate(
      Day = factor(Day, levels = levels(emm_df$Day)),
      x_pos = as.numeric(Day) + x_offset
    )
}

join_emmeans <- function(brackets_df, emm_df) {
  brackets_df %>%
    left_join(
      emm_df %>% dplyr::select(Assay, Day, group, emmean) %>% dplyr::rename(emmean1 = emmean),
      by = c("Assay", "Day", "group1" = "group")
    ) %>%
    left_join(
      emm_df %>% dplyr::select(Assay, Day, group, emmean) %>% dplyr::rename(emmean2 = emmean),
      by = c("Assay", "Day", "group2" = "group")
    ) %>%
    mutate(
      y_position = pmax(emmean1, emmean2, na.rm = TRUE) + 0.2,
      label = star_label(adj.p.value)
    )
}

brackets_allocation <- function(contrasts, emm_df, offsets) {
  d <- contrasts %>% filter(adj.p.value < 0.05)
  if (nrow(d) == 0 || nrow(emm_df) == 0) {
    return(tibble())
  }
  d <- d %>%
    separate(contrast, into = c("from", "to"), sep = " - ") %>%
    separate(from, into = c("Day1", "TR1"), sep = " ") %>%
    separate(to, into = c("Day2", "TR2"), sep = " ") %>%
    mutate(
      Day = factor(gsub("Day", "", Day1), levels = levels(emm_df$Day)),
      group1 = TR1,
      group2 = TR2
    ) %>%
    filter(Day1 == Day2, group1 == "Anakinra", group2 == "Placebo") %>%
    join_emmeans(emm_df)
  attach_offsets(d, emm_df, offsets)
}

brackets_srf <- function(contrasts, emm_df, offsets) {
  d <- contrasts %>% filter(adj.p.value < 0.05)
  if (nrow(d) == 0 || nrow(emm_df) == 0) {
    return(tibble())
  }
  d <- d %>%
    separate(contrast, into = c("from", "to"), sep = " - ") %>%
    separate(from, into = c("Day1", "SRF1"), sep = " ") %>%
    separate(to, into = c("Day2", "SRF2"), sep = " ") %>%
    mutate(
      Day = factor(gsub("Day", "", Day1), levels = levels(emm_df$Day)),
      group1 = SRF1,
      group2 = SRF2
    ) %>%
    filter(Day1 == Day2, group1 == "No", group2 == "Yes") %>%
    join_emmeans(emm_df)
  attach_offsets(d, emm_df, offsets)
}

brackets_threeway <- function(contrasts, emm_df, offsets) {
  d <- contrasts %>% filter(adj.p.value < 0.05)
  if (nrow(d) == 0 || nrow(emm_df) == 0) {
    return(tibble())
  }
  d <- d %>%
    separate(contrast, into = c("from", "to"), sep = " - ") %>%
    separate(from, into = c("Day1", "TR1", "SRF1"), sep = " ") %>%
    separate(to, into = c("Day2", "TR2", "SRF2"), sep = " ") %>%
    mutate(
      Day = factor(gsub("Day", "", Day1), levels = levels(emm_df$Day)),
      group1 = paste(TR1, SRF1, sep = "."),
      group2 = paste(TR2, SRF2, sep = ".")
    ) %>%
    filter(Day1 == Day2, TR1 == "Anakinra", TR2 == "Placebo", SRF1 == SRF2) %>%
    join_emmeans(emm_df)
  attach_offsets(d, emm_df, offsets)
}

hits_srf_within_arm <- function(contrasts, arm) {
  contrasts %>%
    dplyr::filter(adj.p.value < 0.05) %>%
    separate(contrast, into = c("from", "to"), sep = " - ", remove = FALSE) %>%
    separate(from, into = c("Day1", "TR1", "SRF1"), sep = " ") %>%
    separate(to, into = c("Day2", "TR2", "SRF2"), sep = " ") %>%
    dplyr::filter(
      Day1 == Day2,
      TR1 == arm,
      TR2 == arm,
      SRF1 == "No",
      SRF2 == "Yes"
    )
}

brackets_srf_within_arm <- function(contrasts, emm_df, offsets, arm) {
  d <- hits_srf_within_arm(contrasts, arm)
  if (nrow(d) == 0 || nrow(emm_df) == 0) {
    return(tibble())
  }
  d <- d %>%
    mutate(
      Day = factor(gsub("Day", "", Day1), levels = levels(emm_df$Day)),
      group1 = paste(TR1, SRF1, sep = "."),
      group2 = paste(TR2, SRF2, sep = ".")
    ) %>%
    join_emmeans(emm_df)
  attach_offsets(d, emm_df, offsets)
}

threeway_emm_plot_df <- function(emm_df, assays) {
  emm_df %>%
    mutate(
      group = interaction(Allocation, srf14, sep = "."),
      color_class = case_when(
        group == "Anakinra.Yes" ~ "Anakinra - SRF+",
        group == "Placebo.Yes"  ~ "Placebo - SRF+",
        group == "Anakinra.No"  ~ "Anakinra - SRF-",
        group == "Placebo.No"   ~ "Placebo - SRF-"
      )
    ) %>%
    dplyr::filter(Assay %in% assays)
}

make_pointrange <- function(emm_df, brackets_df, colors, group_levels, ylab = "Estimated NPX", connect = TRUE) {
  if (nrow(emm_df) == 0) {
    return(NULL)
  }
  emm_df$group <- factor(emm_df$group, levels = group_levels)
  dodge <- position_dodge(width = 0.6)
  p <- ggplot(emm_df, aes(x = Day, y = emmean, color = color_class, group = group))
  if (isTRUE(connect)) {
    p <- p + geom_line(position = dodge, linewidth = 0.5)
  }
  p <- p +
    geom_pointrange(
      aes(ymin = lower.CL, ymax = upper.CL),
      position = dodge,
      size = 0.6
    ) +
    facet_wrap(~ Assay, scales = "free_y") +
    scale_color_manual(values = colors, name = NULL) +
    theme_minimal() +
    theme(
      strip.text = element_text(face = "bold", size = 9),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "top"
    ) +
    labs(y = ylab, x = "Day")
  if (nrow(brackets_df) > 0) {
    p <- p + geom_text(
      data = brackets_df,
      aes(x = x_pos + x_offset, y = y_position, label = label),
      inherit.aes = FALSE,
      size = 3
    )
  }
  p
}

save_pointrange <- function(p, emm_df, brackets_df, tsv_emm, tsv_brackets, pdf_file) {
  if (is.null(p)) {
    return(invisible(NULL))
  }
  write_tsv(emm_df, tsv_emm)
  write_tsv(brackets_df, tsv_brackets)
  ggsave(pdf_file, plot = p, width = 16, height = 10, units = "in", dpi = 300)
  invisible(p)
}

day_vs_baseline_contrasts <- function(contrasts) {
  contrasts %>%
    dplyr::filter(
      contrast %in% c(
        "Day1 Anakinra - Day4 Anakinra",
        "Day1 Anakinra - Day7 Anakinra",
        "Day1 Placebo - Day4 Placebo",
        "Day1 Placebo - Day7 Placebo"
      )
    ) %>%
    group_by(contrast) %>%
    mutate(adj.p.value = p.adjust(p.value, method = "fdr")) %>%
    ungroup() %>%
    separate(contrast, into = c("from", "to"), sep = " - ", remove = FALSE) %>%
    separate(from, into = c("Day1", "TR1"), sep = " ") %>%
    separate(to, into = c("Day2", "TR2"), sep = " ") %>%
    mutate(
      Day1 = gsub("Day", "", Day1),
      Day2 = gsub("Day", "", Day2),
      Allocation = TR1
    )
}

brackets_day_vs_baseline <- function(contrasts, emm_df, offsets) {
  d <- contrasts %>% dplyr::filter(adj.p.value < 0.05)
  if (nrow(d) == 0 || nrow(emm_df) == 0) {
    return(tibble())
  }
  d <- d %>%
    mutate(
      Day = factor(Day2, levels = levels(emm_df$Day)),
      group = Allocation
    ) %>%
    left_join(
      emm_df %>% dplyr::select(Assay, Day, group, emmean),
      by = c("Assay", "Day", "group")
    ) %>%
    left_join(offsets, by = "Allocation") %>%
    mutate(
      y_position = emmean + 0.2,
      x_pos = as.numeric(Day),
      label = star_label(adj.p.value)
    )
  d
}
