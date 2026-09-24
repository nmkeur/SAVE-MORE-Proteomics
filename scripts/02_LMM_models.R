library(tidyverse)
source("scripts/R/02_helper_functions.R")

dir.create("results", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)

data <- read.table("data/Combined_dataset_mice.csv", header = TRUE, sep = ";")
data$Day <- as.factor(data$Day)



colors_alloc <- c("Anakinra" = "#B22222", "Placebo" = "#1f77b4")
colors_srf <- c("SRF+" = "orange", "SRF-" = "#666666")
colors_threeway <- c(
  "Anakinra - SRF+" = "#B22222",
  "Placebo - SRF+"  = "#1f77b4",
  "Anakinra - SRF-" = "#F08080",
  "Placebo - SRF-"  = "#aec7e8"
)
offsets_alloc <- tribble(
  ~group1, ~group2, ~x_offset,
  "Anakinra", "Placebo", -0.15
)
offsets_srf <- tribble(
  ~group1, ~group2, ~x_offset,
  "No", "Yes", 0.15
)
offsets_threeway <- tribble(
  ~group1, ~group2, ~x_offset,
  "Anakinra.No", "Placebo.No", -0.15,
  "Anakinra.Yes", "Placebo.Yes", 0.15
)
offsets_srf_arm <- tribble(
  ~group1, ~group2, ~x_offset,
  "Anakinra.No", "Anakinra.Yes", -0.15,
  "Placebo.No", "Placebo.Yes", 0.15
)


# -----------------------------------------------------------------------------
# Three-way: Day * Allocation * srf14
# -----------------------------------------------------------------------------
model_threeway <- fit_by_assay(data, safe_lmer_threeway)
write_tsv(failed_assays(model_threeway), "results/lmm_threeway_failed.tsv")

model_threeway <- model_threeway %>%
  keep_fitted() %>%
  add_emm_pairs(list(
    emm = ~ Day * Allocation * srf14,
    emm_day_alloc = ~ Day * Allocation,
    emm_day_srf = ~ Day * srf14
  ))

omnibus_anova <- anova_by_assay(model_threeway)
protein_omni_day_alloc <- hits_term(omnibus_anova, "Day:Allocation")
protein_omni_day_srf <- hits_term(omnibus_anova, "Day:srf14")
protein_omni_threeway <- hits_term(omnibus_anova, "Day:Allocation:srf14")

write_tsv(omnibus_anova, "results/omnibus_interaction_anova.tsv")
write_tsv(protein_omni_day_alloc, "results/omnibus_hits_DayAllocation.tsv")
write_tsv(protein_omni_day_srf, "results/omnibus_hits_DaySRF14.tsv")
write_tsv(protein_omni_threeway, "results/omnibus_hits_DayAllocationSRF.tsv")

lmm_threeway_fixed <- model_threeway %>%
  mutate(fixed = map(model, ~ broom.mixed::tidy(.x, effects = "fixed"))) %>%
  dplyr::select(Assay, fixed) %>%
  unnest(fixed)
write_tsv(lmm_threeway_fixed, "results/lmm_threeway_fixed.tsv")

lmm_threeway_emm <- unnest_emm(model_threeway, "emm") %>%
  mutate(Day = factor(as.character(as.integer(as.character(Day))), levels = day_levels))
lmm_threeway_emm_day_alloc <- unnest_emm(model_threeway, "emm_day_alloc") %>%
  mutate(Day = factor(as.character(as.integer(as.character(Day))), levels = day_levels))
lmm_threeway_emm_day_srf <- unnest_emm(model_threeway, "emm_day_srf") %>%
  mutate(Day = factor(as.character(as.integer(as.character(Day))), levels = day_levels))
lmm_threeway_contrasts <- unnest_contrasts(model_threeway, "tidy_emm")
lmm_threeway_contrasts_day_alloc <- unnest_contrasts(model_threeway, "tidy_emm_day_alloc")
lmm_threeway_contrasts_day_srf <- unnest_contrasts(model_threeway, "tidy_emm_day_srf")

write_tsv(lmm_threeway_emm, "results/lmm_threeway_emm.tsv")
write_tsv(lmm_threeway_emm_day_alloc, "results/lmm_threeway_emm_day_alloc.tsv")
write_tsv(lmm_threeway_emm_day_srf, "results/lmm_threeway_emm_day_srf.tsv")
write_tsv(lmm_threeway_contrasts, "results/lmm_threeway_contrasts.tsv")
write_tsv(lmm_threeway_contrasts_day_alloc, "results/lmm_threeway_contrasts_day_alloc.tsv")
write_tsv(lmm_threeway_contrasts_day_srf, "results/lmm_threeway_contrasts_day_srf.tsv")

rm(model_threeway)

protein_omni_day_alloc
protein_omni_day_srf
protein_omni_threeway

day_alloc_omni_df <- lmm_threeway_emm_day_alloc %>%
  mutate(group = Allocation, color_class = Allocation) %>%
  dplyr::filter(Assay %in% protein_omni_day_alloc$Assay)
posthoc_omni_day_alloc <- fdr_among(lmm_threeway_contrasts_day_alloc, protein_omni_day_alloc$Assay)
brackets_omni_day_alloc_df <- brackets_allocation(posthoc_omni_day_alloc, day_alloc_omni_df, offsets_alloc)
if (nrow(day_alloc_omni_df) > 0) {
  ord <- apply_facet_order(day_alloc_omni_df, brackets_omni_day_alloc_df)
  day_alloc_omni_df <- ord$emm
  brackets_omni_day_alloc_df <- ord$brackets
}
p_omni_day_alloc <- make_pointrange(
  day_alloc_omni_df, brackets_omni_day_alloc_df, colors_alloc, c("Anakinra", "Placebo")
)
save_pointrange(
  p_omni_day_alloc, day_alloc_omni_df, brackets_omni_day_alloc_df,
  "results/Results_Assay_DayAllocation_omnibus_NPXinvnorm.tsv",
  "results/Results_Assay_DayAllocation_omnibus_brackets.tsv",
  "figures/Assay_Pointrange_DayAllocation_omnibus.pdf"
)
p_omni_day_alloc

day_srf_omni_df <- lmm_threeway_emm_day_srf %>%
  mutate(
    group = srf14,
    color_class = case_when(srf14 == "Yes" ~ "SRF+", srf14 == "No" ~ "SRF-")
  ) %>%
  dplyr::filter(Assay %in% protein_omni_day_srf$Assay)
posthoc_omni_day_srf <- fdr_among(lmm_threeway_contrasts_day_srf, protein_omni_day_srf$Assay)
brackets_omni_day_srf_df <- brackets_srf(posthoc_omni_day_srf, day_srf_omni_df, offsets_srf)
if (nrow(day_srf_omni_df) > 0) {
  ord <- apply_facet_order(day_srf_omni_df, brackets_omni_day_srf_df)
  day_srf_omni_df <- ord$emm
  brackets_omni_day_srf_df <- ord$brackets
}
p_omni_day_srf <- make_pointrange(
  day_srf_omni_df, brackets_omni_day_srf_df, colors_srf, c("No", "Yes"),
  ylab = "Estimated INT(NPX)"
)
save_pointrange(
  p_omni_day_srf, day_srf_omni_df, brackets_omni_day_srf_df,
  "results/Results_Assay_DaySRF14_omnibus_NPXinvnorm.tsv",
  "results/Results_Assay_DaySRF14_omnibus_brackets.tsv",
  "figures/Assay_Pointrange_DaySRF14_omnibus.pdf"
)

p_omni_day_srf

threeway_omni_df <- lmm_threeway_emm %>%
  mutate(
    group = interaction(Allocation, srf14, sep = "."),
    color_class = case_when(
      group == "Anakinra.Yes" ~ "Anakinra - SRF+",
      group == "Placebo.Yes"  ~ "Placebo - SRF+",
      group == "Anakinra.No"  ~ "Anakinra - SRF-",
      group == "Placebo.No"   ~ "Placebo - SRF-"
    )
  ) %>%
  dplyr::filter(Assay %in% protein_omni_threeway$Assay)
posthoc_omni_threeway <- fdr_among(lmm_threeway_contrasts, protein_omni_threeway$Assay)
brackets_omni_threeway_df <- brackets_threeway(posthoc_omni_threeway, threeway_omni_df, offsets_threeway)
if (nrow(threeway_omni_df) > 0) {
  ord <- apply_facet_order(threeway_omni_df, brackets_omni_threeway_df)
  threeway_omni_df <- ord$emm
  brackets_omni_threeway_df <- ord$brackets
}
p_omni_threeway <- make_pointrange(
  threeway_omni_df, brackets_omni_threeway_df, colors_threeway,
  c("Anakinra.No", "Anakinra.Yes", "Placebo.No", "Placebo.Yes")
)
save_pointrange(
  p_omni_threeway, threeway_omni_df, brackets_omni_threeway_df,
  "results/Results_Assay_DayAllocationSRF_omnibus_NPXinvnorm.tsv",
  "results/Results_Assay_DayAllocationSRF_omnibus_brackets.tsv",
  "figures/Assay_Pointrange_DayAllocationSRF_omnibus.pdf"
)


# -----------------------------------------------------------------------------
# Same-day SRF No vs Yes within each allocation (all FDR hits, not omnibus-only)
# -----------------------------------------------------------------------------
protein_srf_anakinra <- hits_srf_within_arm(lmm_threeway_contrasts, "Anakinra")
protein_srf_placebo <- hits_srf_within_arm(lmm_threeway_contrasts, "Placebo")
write_tsv(protein_srf_anakinra, "results/lmm_threeway_srf_hits_anakinra.tsv")
write_tsv(protein_srf_placebo, "results/lmm_threeway_srf_hits_placebo.tsv")
protein_srf_anakinra %>% distinct(Assay, Day1)
protein_srf_placebo %>% distinct(Assay, Day1)

assays_srf_both <- union(unique(protein_srf_anakinra$Assay), unique(protein_srf_placebo$Assay))
srf_both_df <- threeway_emm_plot_df(lmm_threeway_emm, assays_srf_both)
brackets_anakinra_srf_df <- brackets_srf_within_arm(
  lmm_threeway_contrasts, srf_both_df, offsets_srf_arm, "Anakinra"
)
brackets_placebo_srf_df <- brackets_srf_within_arm(
  lmm_threeway_contrasts, srf_both_df, offsets_srf_arm, "Placebo"
)
brackets_srf_both_df <- bind_rows(brackets_anakinra_srf_df, brackets_placebo_srf_df)
if (nrow(srf_both_df) > 0) {
  ord <- apply_facet_order(srf_both_df, brackets_srf_both_df)
  srf_both_df <- ord$emm
  brackets_srf_both_df <- ord$brackets
}
p_srf_both <- make_pointrange(
  srf_both_df, brackets_srf_both_df, colors_threeway,
  c("Anakinra.No", "Anakinra.Yes", "Placebo.No", "Placebo.Yes"),
  connect = FALSE
)
save_pointrange(
  p_srf_both, srf_both_df, brackets_srf_both_df,
  "results/Results_Assay_DayAllocationSRF_both_NPXinvnorm.tsv",
  "results/Results_Assay_DayAllocationSRF_both_brackets.tsv",
  "figures/Assay_Pointrange_Allocation_srf14.pdf"
)
write_tsv(
  dplyr::filter(srf_both_df, Assay %in% unique(protein_srf_anakinra$Assay)),
  "results/Results_Assay_DayAllocationSRF_Anakinra_NPXinvnorm.tsv"
)
write_tsv(brackets_anakinra_srf_df, "results/Results_Assay_DayAllocationSRF_Anakinra_brackets_NPXinvnorm.tsv")
write_tsv(
  dplyr::filter(srf_both_df, Assay %in% unique(protein_srf_placebo$Assay)),
  "results/Results_Assay_DayAllocationSRF_Placebo_NPXinvnorm.tsv"
)
write_tsv(brackets_placebo_srf_df, "results/Results_Assay_DayAllocationSRF_Placebo_brackets_NPXinvnorm.tsv")
p_srf_both


# -----------------------------------------------------------------------------
# Two-way: Day * Allocation (no srf14)
# -----------------------------------------------------------------------------
model_twoway_day_alloc <- fit_by_assay(data, safe_lmer_day_alloc)
write_tsv(failed_assays(model_twoway_day_alloc), "results/lmm_twoway_day_alloc_failed.tsv")

model_twoway_day_alloc <- model_twoway_day_alloc %>%
  keep_fitted() %>%
  add_emm_pairs(list(emm_day_alloc = ~ Day * Allocation))

anova_twoway_day_alloc <- anova_by_assay(model_twoway_day_alloc)
protein_twoway_day_alloc <- hits_term(anova_twoway_day_alloc, "Day:Allocation")
lmm_twoway_day_alloc_emm <- unnest_emm(model_twoway_day_alloc, "emm_day_alloc") %>%
  mutate(Day = factor(as.character(as.integer(as.character(Day))), levels = day_levels))
lmm_twoway_day_alloc_contrasts <- unnest_contrasts(model_twoway_day_alloc, "tidy_emm_day_alloc")

write_tsv(anova_twoway_day_alloc, "results/twoway_DayAllocation_anova.tsv")
write_tsv(protein_twoway_day_alloc, "results/twoway_hits_DayAllocation.tsv")
write_tsv(lmm_twoway_day_alloc_emm, "results/lmm_twoway_day_alloc_emm.tsv")
write_tsv(lmm_twoway_day_alloc_contrasts, "results/lmm_twoway_day_alloc_contrasts.tsv")

rm(model_twoway_day_alloc)

protein_twoway_day_alloc

day_alloc_twoway_df <- lmm_twoway_day_alloc_emm %>%
  mutate(group = Allocation, color_class = Allocation) %>%
  dplyr::filter(Assay %in% protein_twoway_day_alloc$Assay)
posthoc_twoway_day_alloc <- fdr_among(lmm_twoway_day_alloc_contrasts, protein_twoway_day_alloc$Assay)
brackets_twoway_day_alloc_df <- brackets_allocation(posthoc_twoway_day_alloc, day_alloc_twoway_df, offsets_alloc)
if (nrow(day_alloc_twoway_df) > 0) {
  ord <- apply_facet_order(day_alloc_twoway_df, brackets_twoway_day_alloc_df)
  day_alloc_twoway_df <- ord$emm
  brackets_twoway_day_alloc_df <- ord$brackets
}
p_twoway_day_alloc <- make_pointrange(
  day_alloc_twoway_df, brackets_twoway_day_alloc_df, colors_alloc, c("Anakinra", "Placebo")
)
save_pointrange(
  p_twoway_day_alloc, day_alloc_twoway_df, brackets_twoway_day_alloc_df,
  "results/Results_Assay_DayAllocation_twoway_NPXinvnorm.tsv",
  "results/Results_Assay_DayAllocation_twoway_brackets.tsv",
  "figures/Assay_Pointrange_DayAllocation_twoway.pdf"
)
p_twoway_day_alloc


# -----------------------------------------------------------------------------
# Within-arm Day 1 vs Day 4 and Day 1 vs Day 7 (two-way Day * Allocation)
# -----------------------------------------------------------------------------
offsets_day_arm <- tribble(
  ~Allocation, ~x_offset,
  "Anakinra", -0.15,
  "Placebo", 0.15
)

day_vs_bl <- day_vs_baseline_contrasts(lmm_twoway_day_alloc_contrasts)
write_tsv(day_vs_bl, "results/lmm_twoway_day_alloc_day_vs_baseline.tsv")
day_vs_bl

protein_day_vs_bl <- day_vs_bl %>%
  dplyr::filter(adj.p.value < 0.05) %>%
  distinct(Assay)
protein_day_vs_bl

day_vs_bl_df <- lmm_twoway_day_alloc_emm %>%
  mutate(group = Allocation, color_class = Allocation) %>%
  dplyr::filter(Assay %in% protein_day_vs_bl$Assay)
brackets_day_vs_bl_df <- brackets_day_vs_baseline(day_vs_bl, day_vs_bl_df, offsets_day_arm)
if (nrow(day_vs_bl_df) > 0) {
  ord <- apply_facet_order(day_vs_bl_df, brackets_day_vs_bl_df)
  day_vs_bl_df <- ord$emm
  brackets_day_vs_bl_df <- ord$brackets
}
p_day_vs_bl <- make_pointrange(
  day_vs_bl_df, brackets_day_vs_bl_df, colors_alloc, c("Anakinra", "Placebo")
)
save_pointrange(
  p_day_vs_bl, day_vs_bl_df, brackets_day_vs_bl_df,
  "results/Results_Assay_DayAllocation_day_vs_baseline_NPXinvnorm.tsv",
  "results/Results_Assay_DayAllocation_day_vs_baseline_brackets.tsv",
  "figures/Assay_Pointrange_DayAllocation_day_vs_baseline.pdf"
)
p_day_vs_bl


# -----------------------------------------------------------------------------
# Two-way: Day * srf14 (no Allocation)
# -----------------------------------------------------------------------------
model_twoway_day_srf <- fit_by_assay(data, safe_lmer_day_srf)
write_tsv(failed_assays(model_twoway_day_srf), "results/lmm_twoway_day_srf_failed.tsv")

model_twoway_day_srf <- model_twoway_day_srf %>%
  keep_fitted() %>%
  add_emm_pairs(list(emm_day_srf = ~ Day * srf14))

anova_twoway_day_srf <- anova_by_assay(model_twoway_day_srf)
protein_twoway_day_srf <- hits_term(anova_twoway_day_srf, "Day:srf14")
lmm_twoway_day_srf_emm <- unnest_emm(model_twoway_day_srf, "emm_day_srf") %>%
  mutate(Day = factor(as.character(as.integer(as.character(Day))), levels = day_levels))
lmm_twoway_day_srf_contrasts <- unnest_contrasts(model_twoway_day_srf, "tidy_emm_day_srf")

write_tsv(anova_twoway_day_srf, "results/twoway_DaySRF14_anova.tsv")
write_tsv(protein_twoway_day_srf, "results/twoway_hits_DaySRF14.tsv")
write_tsv(lmm_twoway_day_srf_emm, "results/lmm_twoway_day_srf_emm.tsv")
write_tsv(lmm_twoway_day_srf_contrasts, "results/lmm_twoway_day_srf_contrasts.tsv")

rm(model_twoway_day_srf)

protein_twoway_day_srf

day_srf_twoway_df <- lmm_twoway_day_srf_emm %>%
  mutate(
    group = srf14,
    color_class = case_when(srf14 == "Yes" ~ "SRF+", srf14 == "No" ~ "SRF-")
  ) %>%
  dplyr::filter(Assay %in% protein_twoway_day_srf$Assay)
posthoc_twoway_day_srf <- fdr_among(lmm_twoway_day_srf_contrasts, protein_twoway_day_srf$Assay)
brackets_twoway_day_srf_df <- brackets_srf(posthoc_twoway_day_srf, day_srf_twoway_df, offsets_srf)
if (nrow(day_srf_twoway_df) > 0) {
  ord <- apply_facet_order(day_srf_twoway_df, brackets_twoway_day_srf_df)
  day_srf_twoway_df <- ord$emm
  brackets_twoway_day_srf_df <- ord$brackets
}
p_twoway_day_srf <- make_pointrange(
  day_srf_twoway_df, brackets_twoway_day_srf_df, colors_srf, c("No", "Yes"),
  ylab = "Estimated INT(NPX)"
)
save_pointrange(
  p_twoway_day_srf, day_srf_twoway_df, brackets_twoway_day_srf_df,
  "results/Results_Assay_DaySRF14_twoway_NPXinvnorm.tsv",
  "results/Results_Assay_DaySRF14_twoway_brackets.tsv",
  "figures/Assay_Pointrange_DaySRF14_twoway.pdf"
)
p_twoway_day_srf
