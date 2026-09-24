library(tidyverse)
source("scripts/R/SAVE_helper_functions.R")

dir.create("results", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)

savemore <- read.table("data/Combined_dataset_mice.csv", header = TRUE, sep = ";")
save_df <- read_save_cohort()
save_n <- save_df %>%
  dplyr::distinct(SampleID, srf14) %>%
  dplyr::count(srf14)
write_tsv(save_n, "results/save_cohort_n.tsv")
save_n
n_distinct(save_df$SampleID)
uniprot <- shared_uniprot_map(savemore, save_df)
uniprot$missing_in_save
length(uniprot$shared)

brackets_anakinra <- read_tsv(
  "results/Results_Assay_DayAllocationSRF_Anakinra_brackets_NPXinvnorm.tsv",
  show_col_types = FALSE
)
more_hits <- more_anakinra_srf_hits(brackets_anakinra)
write_tsv(more_hits, "results/save_more_anakinra_srf_day17.tsv")
more_hits %>% dplyr::count(Day, more_sig)

save_tests <- wilcox_srf_day(save_df, uniprot$map)
write_tsv(save_tests, "results/save_srf_wilcox_day17.tsv")
save_tests

overlap <- join_overlap(more_hits, save_tests)
write_tsv(overlap, "results/save_replication_overlap.tsv")
overlap

overlap_summary <- overlap %>%
  dplyr::filter(Assay %in% uniprot$map$Assay) %>%
  group_by(Day) %>%
  summarise(
    n_shared = n(),
    n_more_sig = sum(more_sig, na.rm = TRUE),
    n_save_fdr = sum(save_sig_fdr, na.rm = TRUE),
    n_save_p = sum(save_sig_p, na.rm = TRUE),
    n_overlap_fdr = sum(overlap_fdr, na.rm = TRUE),
    n_overlap_p = sum(overlap_p, na.rm = TRUE),
    n_overlap_p_same_dir = sum(overlap_p & same_direction, na.rm = TRUE),
    .groups = "drop"
  )
write_tsv(overlap_summary, "results/save_replication_overlap_summary.tsv")
overlap_summary

overlap_hits <- overlap %>%
  dplyr::filter(overlap_p) %>%
  arrange(Day, save_p)
overlap_hits

slope_assays <- more_hits %>%
  dplyr::filter(more_sig, Assay %in% uniprot$map$Assay) %>%
  distinct(Assay) %>%
  pull(Assay)

save_paired <- prepare_save_paired(save_df, uniprot$map, slope_assays)
save_day_srf <- fit_save_day_srf(save_paired)
write_tsv(save_day_srf %>% dplyr::select(Assay, estimate, p.value, p.adj, label), "results/save_day_srf_interaction.tsv")
save_day_srf

p_slopes <- plot_save_slopes(save_paired, save_day_srf)
p_slopes
ggsave(
  "figures/SAVE_replication_slopes.pdf",
  plot = p_slopes,
  width = 16,
  height = max(8, 2.4 * ceiling(length(slope_assays) / 6)),
  units = "in"
)

p_scatter <- plot_overlap_scatter(overlap)
p_scatter
ggsave(
  "figures/SAVE_replication_overlap_scatter.pdf",
  plot = p_scatter,
  width = 10,
  height = 5.5,
  units = "in"
)

assays_d1 <- overlap_hits %>% dplyr::filter(Day == "1") %>% pull(Assay)
assays_d7 <- overlap_hits %>% dplyr::filter(Day == "7") %>% pull(Assay)

p_save_d1 <- plot_save_violins(save_df, uniprot$map, assays_d1, 1, overlap)
p_save_d1
ggsave(
  "figures/SAVE_replication_violins_Day1.pdf",
  plot = p_save_d1,
  width = 12,
  height = max(4, 2.2 * ceiling(length(assays_d1) / 4)),
  units = "in"
)

p_save_d7 <- plot_save_violins(save_df, uniprot$map, assays_d7, 7, overlap)
p_save_d7
ggsave(
  "figures/SAVE_replication_violins_Day7.pdf",
  plot = p_save_d7,
  width = 12,
  height = max(4, 2.2 * ceiling(length(assays_d7) / 4)),
  units = "in"
)
