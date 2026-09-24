library(tidyverse)
library(patchwork)
source("scripts/R/06_helper_functions.R")

dir.create("results", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)

npx <- read.table("data/Combined_dataset_mice.csv", header = TRUE, sep = ";")

# Last pre-SRF vs first active-SRF NPX_invnorm (drop Day-1 SRF and neverSRF)
delta_df <- make_delta_df(npx, npx_col = "NPX_invnorm")
write_tsv(delta_df, "results/delta_pre_active.tsv")
delta_df

delta_wilcox <- wilcox_delta_by_assay(delta_df)
write_tsv(delta_wilcox, "results/delta_pre_active_wilcox.tsv")
delta_wilcox

delta_sig <- delta_wilcox %>%
  dplyr::filter(!is.na(p.value), p.value < 0.05) %>%
  pull(Assay)

p_delta <- plot_delta_violin(delta_df, delta_sig, delta_wilcox)
p_delta
ggsave(
  "figures/Delta_pre_active_violin.pdf",
  plot = p_delta,
  width = 16,
  height = max(6, 2.4 * ceiling(length(delta_sig) / 4)),
  units = "in"
)

# Clinical Spearman split heatmaps (upper Anakinra, lower Placebo; p >= 0.05 set to 0)
clin <- clinical_by_visit(npx)
complete_ids <- patients_with_three_days(clin)
clin_complete <- clin %>% dplyr::filter(PatientID %in% complete_ids)
length(complete_ids)
clin_complete %>% dplyr::count(Day, Allocation)

mat_D1 <- clinical_split_mat(clin_complete, 1)
mat_D4 <- clinical_split_mat(clin_complete, 4)
mat_D7 <- clinical_split_mat(clin_complete, 7)
mat_D1
mat_D4
mat_D7

write.csv(mat_D1, "results/combined_mat_D1.csv", row.names = TRUE, quote = FALSE)
write.csv(mat_D4, "results/combined_mat_D4.csv", row.names = TRUE, quote = FALSE)
write.csv(mat_D7, "results/combined_mat_D7.csv", row.names = TRUE, quote = FALSE)

mat_D1_lab <- label_clin_mat(mat_D1)
mat_D4_lab <- label_clin_mat(mat_D4)
mat_D7_lab <- label_clin_mat(mat_D7)

pD1 <- plot_tri_corr_from_mat(mat_D1_lab, title = "Correlations (D1)")
pD4 <- plot_tri_corr_from_mat(mat_D4_lab, title = "Correlations (D4)")
pD7 <- plot_tri_corr_from_mat(mat_D7_lab, title = "Correlations (D7)")

corr_block <- (pD1 + pD4 + pD7) + plot_layout(guides = "collect") &
  theme(legend.position = "right", plot.margin = margin(4, 4, 2, 4))
corr_block

ggsave(
  "figures/Clinical_corr_heatmaps_D147.pdf",
  plot = corr_block,
  width = 12,
  height = 4.2,
  units = "in"
)

# Olink IL-6 NPX vs day-matched ELISA (Spearman; log2 pg/mL)
Data_D1 <- olink_il6_by_day(npx, 1)
Data_D4 <- olink_il6_by_day(npx, 4)
Data_D7 <- olink_il6_by_day(npx, 7)
Data_D1
Data_D4
Data_D7

p1_il6 <- plot_olink_elisa_panel(Data_D1, "IL6", "IL-6: Olink vs ELISA (Day 1)", "IL-6 ELISA (log2 pg/mL)")
p2_il6 <- plot_olink_elisa_panel(Data_D4, "IL6", "IL-6: Olink vs ELISA (Day 4)", "IL-6 ELISA (log2 pg/mL)")
p3_il6 <- plot_olink_elisa_panel(Data_D7, "IL6", "IL-6: Olink vs ELISA (Day 7)", "IL-6 ELISA (log2 pg/mL)")

p1_crp <- plot_olink_elisa_panel(Data_D1, "CRP", "IL-6 (Olink) vs CRP (ELISA) (Day 1)", "CRP ELISA (log2 pg/mL)")
p2_crp <- plot_olink_elisa_panel(Data_D4, "CRP", "IL-6 (Olink) vs CRP (ELISA) (Day 4)", "CRP ELISA (log2 pg/mL)")
p3_crp <- plot_olink_elisa_panel(Data_D7, "CRP", "IL-6 (Olink) vs CRP (ELISA) (Day 7)", "CRP ELISA (log2 pg/mL)")

p1_fer <- plot_olink_elisa_panel(Data_D1, "ferritin", "IL-6 (Olink) vs Ferritin (ELISA) (Day 1)", "Ferritin ELISA (log2 pg/mL)")
p2_fer <- plot_olink_elisa_panel(Data_D4, "ferritin", "IL-6 (Olink) vs Ferritin (ELISA) (Day 4)", "Ferritin ELISA (log2 pg/mL)")
p3_fer <- plot_olink_elisa_panel(Data_D7, "ferritin", "IL-6 (Olink) vs Ferritin (ELISA) (Day 7)", "Ferritin ELISA (log2 pg/mL)")

p_il6_olink_elisa <- ggpubr::ggarrange(
  p1_il6, p2_il6, p3_il6,
  p1_crp, p2_crp, p3_crp,
  p1_fer, p2_fer, p3_fer,
  nrow = 3, ncol = 3, labels = "AUTO"
)
p_il6_olink_elisa

ggsave(
  "figures/IL6_Olink_vs_ELISA_correlation.pdf",
  plot = p_il6_olink_elisa,
  width = 10,
  height = 10,
  units = "in"
)
