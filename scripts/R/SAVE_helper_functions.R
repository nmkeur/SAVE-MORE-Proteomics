library(tidyverse)
library(lmerTest)
library(broom.mixed)

col_placebo <- "#1f77b4"
col_anakinra <- "#B22222"
col_srf_no <- "steelblue"
col_srf_yes <- "firebrick"

read_save_cohort <- function(complete_inflammation_d1 = TRUE) {
  clin <- readr::read_tsv("data/SAVE/clinical_data.tsv", show_col_types = FALSE) %>%
    dplyr::filter(!is.na(AGE), !is.na(GENDER))
  npx <- readr::read_delim("data/SAVE/NPX_long.txt", delim = ";", show_col_types = FALSE)
  out <- npx %>%
    mutate(
      NPX = as.numeric(NPX),
      Day = dplyr::recode(as.character(Day), "D1" = "1", "D7" = "7")
    ) %>%
    inner_join(clin, by = c("SampleID" = "SAVE_CODE")) %>%
    dplyr::rename(srf14 = DEVELOPED_SEVERE_RESPIRATORY_FAILURE) %>%
    dplyr::filter(srf14 %in% c("No", "Yes"), Day %in% c("1", "7")) %>%
    mutate(srf14 = factor(srf14, levels = c("No", "Yes")))
  # 153 = 135 SRF No + 18 SRF Yes: Day-1 NPX for ADA (typical inflammation assay).
  # Without this, unique AGE+GENDER+any-NPX patients are 155 (adds C052, N020).
  if (isTRUE(complete_inflammation_d1)) {
    keep <- out %>%
      dplyr::filter(Day == "1", Assay == "ADA", !is.na(NPX)) %>%
      dplyr::distinct(SampleID)
    out <- dplyr::semi_join(out, keep, by = "SampleID")
  }
  out
}

shared_uniprot_map <- function(savemore, save_df) {
  more_map <- savemore %>%
    dplyr::distinct(Assay, UniProt)
  shared <- intersect(more_map$UniProt, unique(save_df$UniProt))
  missing_more <- setdiff(more_map$UniProt, unique(save_df$UniProt))
  list(
    map = more_map %>% dplyr::filter(UniProt %in% shared),
    shared = shared,
    missing_in_save = missing_more
  )
}

more_anakinra_srf_hits <- function(brackets) {
  brackets %>%
    dplyr::filter(
      TR1 == "Anakinra",
      TR2 == "Anakinra",
      SRF1 == "No",
      SRF2 == "Yes",
      Day %in% c(1, 7)
    ) %>%
    mutate(
      Day = as.character(Day),
      more_diff = emmean2 - emmean1,
      more_p = p.value,
      more_p_adj = adj.p.value,
      more_sig = adj.p.value < 0.05
    ) %>%
    dplyr::select(Assay, Day, more_diff, more_p, more_p_adj, more_sig)
}

wilcox_srf_day <- function(save_df, assay_map) {
  d <- save_df %>%
    inner_join(assay_map, by = "UniProt", suffix = c("_save", "")) %>%
    dplyr::filter(!is.na(NPX), !is.na(srf14))

  d %>%
    group_by(Assay, UniProt, Day) %>%
    group_modify(function(x, key) {
      if (n_distinct(x$srf14) < 2L || length(unique(na.omit(x$NPX))) < 2L) {
        return(tibble(
          n_no = sum(x$srf14 == "No"),
          n_yes = sum(x$srf14 == "Yes"),
          median_no = median(x$NPX[x$srf14 == "No"], na.rm = TRUE),
          median_yes = median(x$NPX[x$srf14 == "Yes"], na.rm = TRUE),
          save_diff = NA_real_,
          save_p = NA_real_
        ))
      }
      wt <- wilcox.test(NPX ~ srf14, data = x, exact = FALSE)
      med_no <- median(x$NPX[x$srf14 == "No"], na.rm = TRUE)
      med_yes <- median(x$NPX[x$srf14 == "Yes"], na.rm = TRUE)
      tibble(
        n_no = sum(x$srf14 == "No"),
        n_yes = sum(x$srf14 == "Yes"),
        median_no = med_no,
        median_yes = med_yes,
        save_diff = med_yes - med_no,
        save_p = unname(wt$p.value)
      )
    }) %>%
    ungroup() %>%
    group_by(Day) %>%
    mutate(save_p_adj = p.adjust(save_p, method = "BH")) %>%
    ungroup()
}

join_overlap <- function(more_hits, save_tests) {
  more_hits %>%
    full_join(save_tests, by = c("Assay", "Day")) %>%
    mutate(
      save_sig_fdr = !is.na(save_p_adj) & save_p_adj < 0.05,
      save_sig_p = !is.na(save_p) & save_p < 0.05,
      same_direction = sign(more_diff) == sign(save_diff) &
        !is.na(more_diff) & !is.na(save_diff),
      overlap_fdr = more_sig & save_sig_fdr,
      overlap_p = more_sig & save_sig_p
    ) %>%
    arrange(Day, save_p)
}

plot_overlap_scatter <- function(overlap_df) {
  d <- overlap_df %>%
    dplyr::filter(!is.na(more_diff), !is.na(save_diff)) %>%
    mutate(
      status = dplyr::case_when(
        overlap_p ~ "Hit in both (p<0.05 SAVE)",
        more_sig ~ "SAVE-MORE hit only",
        TRUE ~ "Not a SAVE-MORE hit"
      ),
      Day = paste("Day", Day)
    )

  ggplot(d, aes(x = more_diff, y = save_diff, color = status)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
    geom_point(alpha = 0.85, size = 2) +
    ggrepel::geom_text_repel(
      data = dplyr::filter(d, overlap_p | (more_sig & save_p < 0.1)),
      aes(label = Assay),
      size = 3,
      max.overlaps = 40,
      show.legend = FALSE
    ) +
    scale_color_manual(
      values = c(
        "Hit in both (p<0.05 SAVE)" = col_anakinra,
        "SAVE-MORE hit only" = col_placebo,
        "Not a SAVE-MORE hit" = "grey70"
      )
    ) +
    facet_wrap(~Day) +
    labs(
      x = "SAVE-MORE Anakinra: EMM Yes - No",
      y = "SAVE: median NPX Yes - No",
      color = NULL
    ) +
    theme_classic() +
    theme(legend.position = "top")
}

plot_save_violins <- function(save_df, assay_map, assays, day, tests = NULL) {
  d <- save_df %>%
    inner_join(assay_map, by = "UniProt", suffix = c("_save", "")) %>%
    dplyr::filter(Assay %in% assays, Day == as.character(day)) %>%
    mutate(Assay = factor(Assay, levels = assays))

  p <- ggplot(d, aes(x = srf14, y = NPX, fill = srf14)) +
    geom_violin(alpha = 0.35, color = "black", width = 0.8) +
    geom_boxplot(width = 0.2, outlier.size = 0.4) +
    scale_fill_manual(values = c(No = col_srf_no, Yes = col_srf_yes)) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.22))) +
    facet_wrap(~Assay, scales = "free_y", ncol = 4) +
    labs(x = "SRF", y = "SAVE NPX") +
    theme_classic() +
    theme(
      strip.text = element_text(face = "bold"),
      legend.position = "none"
    )

  if (!is.null(tests) && length(assays) > 0L) {
    y_top <- d %>%
      group_by(Assay) %>%
      summarise(y = max(NPX, na.rm = TRUE), .groups = "drop")
    p_lab <- tests %>%
      mutate(Day = as.character(Day)) %>%
      dplyr::filter(Day == as.character(day), Assay %in% assays) %>%
      mutate(Assay = factor(Assay, levels = assays)) %>%
      left_join(y_top, by = "Assay") %>%
      mutate(
        label = paste0(
          "p=",
          ifelse(is.na(save_p), "NA", signif(save_p, 2)),
          "; diff=",
          ifelse(is.na(save_diff), "NA", signif(save_diff, 2))
        )
      )
    p <- p + geom_text(
      data = p_lab,
      aes(x = 1.5, y = y, label = gsub("; ", "\n", label)),
      inherit.aes = FALSE,
      vjust = -0.15,
      size = 3
    )
  }
  p
}

prepare_save_paired <- function(save_df, assay_map, assays) {
  save_df %>%
    inner_join(assay_map, by = "UniProt", suffix = c("_save", "")) %>%
    dplyr::filter(Assay %in% assays, !is.na(NPX)) %>%
    mutate(
      Day = factor(Day, levels = c("1", "7")),
      srf14 = factor(srf14, levels = c("No", "Yes"))
    ) %>%
    group_by(SampleID, Assay) %>%
    dplyr::filter(n_distinct(Day) == 2) %>%
    ungroup()
}

fit_save_day_srf <- function(paired_df) {
  paired_df %>%
    group_by(Assay) %>%
    nest() %>%
    mutate(
      model = purrr::map(data, function(d) {
        tryCatch(
          lmerTest::lmer(NPX ~ Day * srf14 + AGE + GENDER + DEXAMETHAZONE +(1 | SampleID), data = d),
          error = function(e) NULL
        )
      }),
      tidied = purrr::map(model, function(m) {
        if (is.null(m)) {
          return(tibble(term = NA_character_, estimate = NA_real_, p.value = NA_real_))
        }
        broom.mixed::tidy(m)
      }),
      y.position = purrr::map_dbl(data, ~ mean(.x$NPX[.x$Day == "1"], na.rm = TRUE))
    ) %>%
    tidyr::unnest(tidied) %>%
    dplyr::filter(term == "Day7:srf14Yes") %>%
    ungroup() %>%
    mutate(
      p.adj = p.adjust(p.value, method = "holm"),
      label = paste0(
        "p=",
        ifelse(is.na(p.adj), "NA", ifelse(p.adj >= 0.995, "1", signif(p.adj, 2))),
        "; Δ=",
        ifelse(is.na(estimate), "NA", signif(estimate, 2))
      )
    )
}

plot_save_slopes <- function(paired_df, tests) {
  means <- paired_df %>%
    group_by(Assay, Day, srf14) %>%
    summarise(NPX = mean(NPX, na.rm = TRUE), .groups = "drop")

  y_top <- paired_df %>%
    group_by(Assay) %>%
    summarise(y = max(NPX, na.rm = TRUE), .groups = "drop")
  p_lab <- tests %>%
    dplyr::select(Assay, label) %>%
    left_join(y_top, by = "Assay") %>%
    mutate(label = gsub("; ", "\n", label))

  ggplot(paired_df, aes(x = Day, y = NPX, color = srf14)) +
    geom_line(aes(group = SampleID), alpha = 0.22, linewidth = 0.35) +
    geom_point(aes(group = SampleID), alpha = 0.22, size = 0.9) +
    geom_line(data = means, aes(group = srf14), linewidth = 1.1) +
    geom_point(data = means, aes(group = srf14), size = 2.4) +
    geom_text(
      data = p_lab,
      aes(x = 1.5, y = y, label = label),
      inherit.aes = FALSE,
      vjust = -0.15,
      size = 2.6
    ) +
    scale_x_discrete(labels = c("1" = "Day 1", "7" = "Day 7")) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.22))) +
    scale_color_manual(values = c(No = "#4C78A8", Yes = "#E07A3D"), name = "SRF") +
    facet_wrap(~Assay, scales = "free_y", ncol = 6) +
    labs(x = NULL, y = "NPX") +
    theme_classic() +
    theme(
      legend.position = "top",
      strip.text = element_text(face = "bold", size = 9)
    )
}
