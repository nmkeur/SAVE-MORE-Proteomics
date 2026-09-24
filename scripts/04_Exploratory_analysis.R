library(tidyverse)
library(patchwork)
source("scripts/R/04_helper_functions.R")

dir.create("results", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)

# Frozen STRING layout. Do not overwrite ppi_scaffold_*_xy.tsv.
# Mice brackets only recolor proteins already on the scaffold.
data_anakinra <- read_tsv(
  "results/Results_Assay_DayAllocationSRF_Anakinra_brackets_NPXinvnorm.tsv",
  show_col_types = FALSE
)
data_placebo <- read_tsv(
  "results/Results_Assay_DayAllocationSRF_Placebo_brackets_NPXinvnorm.tsv",
  show_col_types = FALSE
)
scaffold_anakinra <- read_tsv("results/ppi_scaffold_anakinra_xy.tsv", show_col_types = FALSE)
scaffold_placebo <- read_tsv("results/ppi_scaffold_placebo_xy.tsv", show_col_types = FALSE)

arm_anakinra <- build_arm_network(data_anakinra, scaffold_anakinra, "Anakinra")
arm_placebo <- build_arm_network(data_placebo, scaffold_placebo, "Placebo")

table(arm_anakinra$nodes$node_status)
table(arm_placebo$nodes$node_status)
arm_anakinra$mice_only_omitted
arm_placebo$mice_only_omitted

write_tsv(arm_anakinra$nodes, "results/ppi_nodes_anakinra_scaffold.tsv")
write_tsv(arm_placebo$nodes, "results/ppi_nodes_placebo_scaffold.tsv")
write_tsv(arm_anakinra$edges, "results/ppi_edges_anakinra_scaffold.tsv")
write_tsv(arm_placebo$edges, "results/ppi_edges_placebo_scaffold.tsv")

p_net_anakinra <- arm_anakinra$plot
p_net_placebo <- arm_placebo$plot

p_net_placebo
p_net_anakinra

combined <- (p_net_placebo | p_net_anakinra) +
  plot_annotation(tag_levels = "A") &
  theme(legend.position = "none")
combined

ggsave(
  "figures/PPI_Allocation_scaffold.pdf",
  plot = combined,
  width = 14,
  height = 7,
  units = "in"
)


# -----------------------------------------------------------------------------
# Heatmap: row-scaled 3-way LMM EMMs (proteins as columns)
# -----------------------------------------------------------------------------
emm_threeway <- read_tsv("results/lmm_threeway_emm.tsv", show_col_types = FALSE)
emm_matrix <- prepare_emm_matrix(emm_threeway)
scaled_traj <- scale_protein_trajectories(emm_matrix)

heat_emm <- cluster_emm_heatmap(scaled_traj)
protein_clusters <- heat_emm$protein_clusters
condition_clusters <- heat_emm$condition_clusters

write_tsv(protein_clusters, "results/emm_heatmap_protein_clusters.tsv")
write_tsv(condition_clusters, "results/emm_heatmap_condition_clusters.tsv")
write_csv(
  protein_clusters %>% dplyr::rename(Protein = Assay),
  "results/Protein_Clusters.csv"
)

p_heat <- heat_emm$plot
p_heat
ggsave(
  "figures/EMM_heatmap.pdf",
  plot = p_heat,
  width = 18,
  height = 6,
  units = "in"
)

