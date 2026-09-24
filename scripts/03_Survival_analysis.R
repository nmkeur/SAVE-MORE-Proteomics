library(tidyverse)
library(survival)
library(survminer)
library(patchwork)
source("scripts/R/03_helper_functions.R")

dir.create("results", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)

npx <- read.table("data/Combined_dataset_mice.csv", header = TRUE, sep = ";")

# Day × SRF proteins that already differ at Day 1 (two-way LMM same-day contrast)
day_srf_day1 <- read_tsv("results/Results_Assay_DaySRF14_twoway_brackets.tsv", show_col_types = FALSE) %>%
  dplyr::filter(Day1 == "Day1", group1 == "No", group2 == "Yes", adj.p.value < 0.05)

# Working set for univariate / 6-protein multivariate (not the full Day-1 SRF list)
prot_raw <- c("PD-L1", "IL6", "MCP-3", "MCP-1", "CD244", "TNFB")
prot_raw_selected <- c("PD-L1", "IL6", "TNFB")

surv <- prepare_surv_day1(npx, prot_raw)
surv_data <- surv$wide
prot_vars <- surv$prot_vars
prot_vars_selected <- make.names(prot_raw_selected)
surv_anakinra <- dplyr::filter(surv_data, Allocation == "Anakinra")
surv_placebo <- dplyr::filter(surv_data, Allocation == "Placebo")



# -----------------------------------------------------------------------------
# Univariate Cox per protein, separately by arm
# -----------------------------------------------------------------------------
uni_anakinra <- cox_proteins_by_arm(surv_data, prot_vars, "Anakinra")
uni_placebo <- cox_proteins_by_arm(surv_data, prot_vars, "Placebo")
uni_by_arm <- bind_rows(uni_anakinra, uni_placebo)

write_tsv(uni_by_arm, "results/survival_univariate_by_arm.tsv")
uni_by_arm

p_forest_uni <- plot_arm_forest(
  uni_by_arm,
  title = "Univariate Cox: Day 1 Day × SRF proteins"
)

p_forest_uni
ggsave(
  "figures/Survival_forest_univariate_by_arm.pdf",
  plot = p_forest_uni,
  width = 12,
  height = max(4.5, 0.5 * length(prot_raw) + 2.2),
  units = "in"
)


# -----------------------------------------------------------------------------
# Multivariate Cox per arm (6 proteins + Age + cbmi + Gender + DEXcorrected)
# -----------------------------------------------------------------------------
cox_multi_anakinra <- cox_multivariate_arm(surv_data, prot_vars, "Anakinra")
cox_multi_placebo <- cox_multivariate_arm(surv_data, prot_vars, "Placebo")

multi_anakinra <- tidy_cox(cox_multi_anakinra, "Anakinra") %>%
  dplyr::filter(term %in% prot_vars) %>%
  mutate(protein = protein_label(term))
multi_placebo <- tidy_cox(cox_multi_placebo, "Placebo") %>%
  dplyr::filter(term %in% prot_vars) %>%
  mutate(protein = protein_label(term))
multi_by_arm <- bind_rows(multi_anakinra, multi_placebo)

write_tsv(tidy_cox(cox_multi_anakinra, "Anakinra"), "results/survival_multivariate_anakinra.tsv")
write_tsv(tidy_cox(cox_multi_placebo, "Placebo"), "results/survival_multivariate_placebo.tsv")
write_tsv(multi_by_arm, "results/survival_multivariate_proteins_by_arm.tsv")

summary(cox_multi_anakinra)
summary(cox_multi_placebo)
multi_by_arm

p_forest_multi <- plot_arm_forest(
  multi_by_arm,
  title = "Multivariate Cox: Day 1 Day × SRF proteins"
)
p_forest_multi
ggsave(
  "figures/Survival_forest_multivariate_by_arm.pdf",
  plot = p_forest_multi,
  width = 12,
  height = max(4.5, 0.5 * length(prot_raw) + 2.2),
  units = "in"
)


# -----------------------------------------------------------------------------
# Multivariate Cox per arm, proteins kept after the 6-protein model
# (PD-L1, IL6, TNFB + Age + cbmi + Gender + DEXcorrected)
# -----------------------------------------------------------------------------
cox_multiSelected_anakinra <- cox_multivariate_arm(surv_data, prot_vars_selected, "Anakinra")
cox_multiSelected_placebo <- cox_multivariate_arm(surv_data, prot_vars_selected, "Placebo")

multiSelected_anakinra <- tidy_cox(cox_multiSelected_anakinra, "Anakinra") %>%
  dplyr::filter(term %in% prot_vars_selected) %>%
  mutate(protein = protein_label(term))
multiSelected_placebo <- tidy_cox(cox_multiSelected_placebo, "Placebo") %>%
  dplyr::filter(term %in% prot_vars_selected) %>%
  mutate(protein = protein_label(term))
multiSelected_by_arm <- bind_rows(multiSelected_anakinra, multiSelected_placebo)

write_tsv(tidy_cox(cox_multiSelected_anakinra, "Anakinra"), "results/survival_multivariateSelected_anakinra.tsv")
write_tsv(tidy_cox(cox_multiSelected_placebo, "Placebo"), "results/survival_multivariateSelected_placebo.tsv")
write_tsv(multiSelected_by_arm, "results/survival_multivariateSelected_proteins_by_arm.tsv")

summary(cox_multiSelected_anakinra)
summary(cox_multiSelected_placebo)
multiSelected_by_arm

p_forest_multiSelected <- plot_arm_forest(
  multiSelected_by_arm,
  title = "Multivariate Cox: selected Day 1 proteins"
)
p_forest_multiSelected
ggsave(
  "figures/Survival_forest_multivariate_selected_by_arm.pdf",
  plot = p_forest_multiSelected,
  width = 12,
  height = max(4.5, 0.5 * length(prot_raw_selected) + 2.2),
  units = "in"
)


# -----------------------------------------------------------------------------
# Formal Allocation × protein interaction (Day 1; IL6, PD-L1, TNFB)
# Per-protein models, then one joint model (proteins)*Allocation + covariates.
# Supplementary Table 2 is the joint interaction terms.
# -----------------------------------------------------------------------------
cox_int_IL6 <- cox_protein_allocation_int(surv_data, "IL6")
cox_int_PDL1 <- cox_protein_allocation_int(surv_data, "PD.L1")
cox_int_TNFB <- cox_protein_allocation_int(surv_data, "TNFB")
summary(cox_int_IL6)
summary(cox_int_PDL1)
summary(cox_int_TNFB)

cox_int_by_protein <- bind_rows(
  tidy_allocation_int(cox_int_IL6, "IL6", "per_protein"),
  tidy_allocation_int(cox_int_PDL1, "PD.L1", "per_protein"),
  tidy_allocation_int(cox_int_TNFB, "TNFB", "per_protein")
)
write_tsv(cox_int_by_protein, "results/survival_allocation_interaction_by_protein.tsv")
cox_int_by_protein

cox_int_add <- cox_joint_allocation_add(surv_data, prot_vars_selected)
cox_int_joint <- cox_joint_allocation_int(surv_data, prot_vars_selected)
anova_int <- anova(cox_int_add, cox_int_joint, test = "LRT")
anova_int
summary(cox_int_joint)

cox_int_joint_terms <- tidy_allocation_int(
  cox_int_joint,
  c("IL6", "PD.L1", "TNFB"),
  "joint"
)
write_tsv(cox_int_joint_terms, "results/survival_allocation_interaction_joint.tsv")
cox_int_joint_terms


# -----------------------------------------------------------------------------
# Cox: protein-unadjusted vs protein-adjusted (same clinical covariates)
# Unadjusted = Allocation + Age + cbmi + Gender + DEXcorrected
# Adjusted   = those terms + PD-L1 + IL6 + TNFB
# Both curves are Cox predictions at the same reference profile.
# -----------------------------------------------------------------------------
cox_clin <- cox_clinical_allocation(surv_data)
write_tsv(tidy_cox(cox_clin, "pooled_clinical"), "results/survival_clinical_pooled.tsv")
summary(cox_clin)

cox_multi <- coxph(
  as.formula(paste(
    "Surv(time, status) ~",
    paste(c(prot_vars_selected, "Allocation", "Age", "cbmi", "Gender", "DEXcorrected"), collapse = " + ")
  )),
  data = surv_data,
  na.action = na.exclude
)
write_tsv(tidy_cox(cox_multi, "pooled_multivariate"), "results/survival_multivariate_pooled.tsv")
summary(cox_multi)

newdat <- reference_profile(surv_data, prot_vars_selected)
fit_clin <- survfit(cox_clin, newdata = newdat)
fit_adj <- survfit(cox_multi, newdata = newdat)

km_ylim <- event_ylim(fit_clin, fit_adj)
km_obs <- survfit(Surv(time, status) ~ Allocation, data = surv_data)
p_tab <- plot_km_unadjusted(km_obs, surv_data, ylim = km_ylim)

p_km <- add_surv_p(
  plot_km_adjusted(fit_clin, newdat, ylim = km_ylim),
  allocation_p_label(cox_clin),
  km_ylim
)
p_km_print <- p_km$plot /
  (p_tab$table + labs(title = NULL, x = "Time to SRF (days)")) +
  plot_layout(heights = c(2.4, 1))
p_km_print

pdf("figures/Survival_KM_unadjusted_Allocation.pdf", width = 7, height = 8)
print(p_km_print)
dev.off()

p_adj <- add_surv_p(
  plot_km_adjusted(fit_adj, newdat, ylim = km_ylim),
  allocation_p_label(cox_multi),
  km_ylim
)
p_adj

pdf("figures/Survival_KM_adjusted_Allocation.pdf", width = 7, height = 6)
print(p_adj)
dev.off()

p_side <- wrap_elements(p_km_print) +
  wrap_elements(p_adj$plot + labs(title = "With proteins (PD-L1 + IL-6 + TNF-B)")) +
  plot_layout(widths = c(1, 1))

p_side
ggsave("figures/Survival_KM_unadjusted_vs_adjusted.pdf", plot = p_side, width = 12, height = 5.5, units = "in")


# -----------------------------------------------------------------------------
# AUC: Day 1, SRF-free at baseline; Olink assays, then clinical CRP / suPAR
# -----------------------------------------------------------------------------
auc_day1 <- prepare_auc_day1(npx)
auc_wide <- auc_day1$wide
auc_prot_vars <- auc_day1$prot_vars
auc_lab_vars <- auc_day1$lab_vars

auc_placebo <- dplyr::filter(auc_wide, Allocation == "Placebo")
auc_anakinra <- dplyr::filter(auc_wide, Allocation == "Anakinra")

auc_res_placebo <- compute_roc(auc_placebo, auc_prot_vars, "Placebo", source = "Olink")
auc_res_anakinra <- compute_roc(auc_anakinra, auc_prot_vars, "Anakinra", source = "Olink")
auc_res_global <- compute_roc(auc_wide, auc_prot_vars, "Global", source = "Olink")

auc_lab_placebo <- compute_roc(auc_placebo, auc_lab_vars, "Placebo", source = "Clinical")
auc_lab_anakinra <- compute_roc(auc_anakinra, auc_lab_vars, "Anakinra", source = "Clinical")
auc_lab_global <- compute_roc(auc_wide, auc_lab_vars, "Global", source = "Clinical")

auc_olink <- bind_rows(auc_res_placebo, auc_res_anakinra, auc_res_global)
auc_lab <- bind_rows(auc_lab_placebo, auc_lab_anakinra, auc_lab_global)
auc_all <- bind_rows(auc_olink, auc_lab)
write_tsv(auc_all, "results/survival_auc_day1_by_arm.tsv")
write_tsv(auc_lab, "results/survival_auc_clinical_labs.tsv")
auc_all
auc_lab

auc_cox <- auc_all %>%
  dplyr::filter(Protein %in% prot_raw)
write_tsv(auc_cox, "results/survival_auc_cox_proteins.tsv")
auc_cox

auc_top5 <- auc_all %>%
  dplyr::filter(!is.na(AUC)) %>%
  group_by(Arm) %>%
  slice_max(order_by = AUC, n = 5, with_ties = FALSE) %>%
  ungroup()
auc_top5

AUC_plot <- plot_auc_top(auc_all)
AUC_plot
ggsave(
  "figures/Survival_AUC_day1_top_proteins.pdf",
  plot = AUC_plot,
  width = 10,
  height = 7,
  units = "in"
)
