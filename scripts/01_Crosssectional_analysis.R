library(tidyverse)
library(patchwork)
source("scripts/R/01_helper_functions.R")

data <- read.table("data/Combined_dataset_mice.csv", header = TRUE, sep = ";")
data$Day <- factor(sprintf("%02d", as.integer(data$Day)), levels = c("01", "04", "07"))
data$Allocation <- factor(data$Allocation, levels = c("Placebo", "Anakinra"))

cross_sectional_n <- data %>%
  dplyr::filter(
    !is.na(NPX_wins),
    !is.na(Allocation),
    !is.na(cbmi),
    !is.na(Age),
    !is.na(Gender),
    !is.na(DEXcorrected)
  ) %>%
  dplyr::distinct(SampleID, Day, Allocation) %>%
  dplyr::count(Day, Allocation, name = "n") %>%
  tidyr::pivot_wider(names_from = Allocation, values_from = n, values_fill = 0) %>%
  mutate(
    Day = recode(as.character(Day), "01" = "1", "04" = "4", "07" = "7"),
    Total = Placebo + Anakinra
  ) %>%
  dplyr::select(Day, Placebo, Anakinra, Total)
write_tsv(cross_sectional_n, "results/Results_Assay_crossectional_n.tsv")
cross_sectional_n

# Within-day Anakinra vs Placebo (not the LMM)
cross_sectional <- data %>%
  group_by(Assay, Day) %>%
  nest() %>%
  dplyr::mutate(
    model = purrr::map(data, ~ lm(NPX_wins ~ Allocation + cbmi + Age + Gender + DEXcorrected, data = .x)),
    results = purrr::map(model, broom::tidy)
  ) %>%
  dplyr::select(-data, -model) %>%
  unnest(results) %>%
  dplyr::filter(term == "AllocationAnakinra") %>%
  mutate(
    log2FC = estimate,
    p_value = p.value
  ) %>%
  arrange(Assay, Day) %>%
  ungroup() %>%
  group_by(Day) %>%
  dplyr::mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
  ungroup()

cross_sectional_D1 <- subset(cross_sectional, Day == "01")
cross_sectional_D4 <- subset(cross_sectional, Day == "04")
cross_sectional_D7 <- subset(cross_sectional, Day == "07")

write_tsv(cross_sectional, "results/Results_Assay_crossectional_NPXinvnorm.tsv")
write_tsv(cross_sectional_D1, "results/Results_Assay_crossectional_Day1_NPXinvnorm.tsv")
write_tsv(cross_sectional_D4, "results/Results_Assay_crossectional_Day4_NPXinvnorm.tsv")
write_tsv(cross_sectional_D7, "results/Results_Assay_crossectional_Day7_NPXinvnorm.tsv")

Vol_D1 <- make_volcano(cross_sectional_D1, "p.value", "Baseline - Anakinra vs Placebo", ylab_raw)
Vol_D4 <- make_volcano(cross_sectional_D4, "p.value", "Day 4 - Anakinra vs Placebo", ylab_raw)
Vol_D7 <- make_volcano(cross_sectional_D7, "p.value", "Day 7 - Anakinra vs Placebo", ylab_raw)

Vol_D1a <- make_volcano(cross_sectional_D1, "p_adj", "Baseline - Anakinra vs Placebo", ylab_adj)
Vol_D4a <- make_volcano(cross_sectional_D4, "p_adj", "Day 4 - Anakinra vs Placebo", ylab_adj)
Vol_D7a <- make_volcano(cross_sectional_D7, "p_adj", "Day 7 - Anakinra vs Placebo", ylab_adj)

volcano_grid <- (Vol_D1a | Vol_D4a | Vol_D7a) / (Vol_D1 | Vol_D4 | Vol_D7)
volcano_grid <- volcano_grid + plot_layout(guides = "collect") &
  theme(legend.position = "right")

volcano_grid

ggsave(
  "figures/Volcano_Allocation_withinDay.pdf",
  plot = volcano_grid,
  width = 16,
  height = 10,
  units = "in",
  dpi = 300
)
