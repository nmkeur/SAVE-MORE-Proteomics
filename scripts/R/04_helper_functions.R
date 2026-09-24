library(rbioapi)
library(tidyverse)
library(igraph)
library(tidygraph)
library(ggraph)
library(patchwork)
library(pheatmap)
library(ggplotify)

string_required_score <- 700

one_estimate <- function(brackets, map_df) {
  brackets %>%
    dplyr::select(Assay, estimate) %>%
    group_by(Assay) %>%
    arrange(estimate) %>%
    slice_head(n = 1) %>%
    ungroup() %>%
    left_join(
      map_df %>% dplyr::select(queryItem, preferredName, stringId),
      by = c("Assay" = "queryItem")
    )
}

build_arm_network <- function(brackets, scaffold, title, required_score = string_required_score) {
  assays_mice <- unique(brackets$Assay)
  assays_scaf <- unique(na.omit(scaffold$Assay))

  mapped <- rba_string_map_ids(ids = assays_scaf, species = 9606)
  est <- one_estimate(brackets, mapped)

  dropped_mice <- setdiff(assays_mice, assays_scaf)
  if (length(dropped_mice) > 0) {
    message(
      title, ": mice-only proteins not added to the map: ",
      paste(dropped_mice, collapse = ", ")
    )
  }

  int_net <- rba_string_interactions_network(
    ids = assays_scaf,
    species = 9606,
    required_score = required_score
  )

  if (is.null(int_net) || nrow(int_net) == 0L) {
    edges <- tibble(from = character(), to = character(), score = numeric())
  } else {
    edges <- int_net %>%
      transmute(from = preferredName_A, to = preferredName_B, score = score) %>%
      dplyr::filter(from %in% scaffold$name, to %in% scaffold$name)
  }

  nodes <- scaffold %>%
    dplyr::select(name, x, y, Assay) %>%
    left_join(est %>% dplyr::select(Assay, estimate), by = "Assay") %>%
    mutate(
      in_mice = Assay %in% assays_mice,
      node_status = if_else(in_mice, "hit", "scaffold_only"),
      regulation = case_when(
        !in_mice ~ NA_character_,
        estimate > 0 ~ "Up",
        TRUE ~ "Down"
      ),
      plot_col = case_when(
        !in_mice ~ "Not in mice hits",
        regulation == "Up" ~ "Up",
        TRUE ~ "Down"
      )
    )

  g <- graph_from_data_frame(
    d = as.data.frame(edges),
    vertices = as.data.frame(nodes),
    directed = FALSE
  )
  tg <- as_tbl_graph(g) %>%
    activate(nodes) %>%
    mutate(deg = centrality_degree())

  list(
    mapped = mapped,
    est = est,
    nodes = nodes,
    edges = edges,
    mice_only_omitted = dropped_mice,
    tg = tg,
    plot = plot_ppi_ggraph(tg, title)
  )
}

plot_ppi_ggraph <- function(tg, title) {
  hub_nodes <- tg %>%
    activate(nodes) %>%
    as_tibble() %>%
    dplyr::filter(in_mice) %>%
    slice_max(deg, n = 50, with_ties = FALSE) %>%
    pull(name)

  ggraph(tg, layout = "manual", x = x, y = y) +
    geom_edge_link(
      aes(width = score),
      colour = "grey40",
      alpha = 0.25,
      show.legend = FALSE
    ) +
    scale_edge_width(range = c(0.2, 1.2)) +
    geom_node_point(aes(color = plot_col, shape = node_status), size = 4) +
    geom_node_text(
      data = function(d) dplyr::filter(d, name %in% hub_nodes),
      aes(label = Assay),
      family = "sans",
      repel = TRUE,
      size = 3
    ) +
    scale_color_manual(
      values = c(
        "Down" = "skyblue",
        "Up" = "tomato",
        "Not in mice hits" = "grey75"
      ),
      breaks = c("Up", "Down"),
      guide = "none"
    ) +
    scale_shape_manual(values = c(hit = 16, scaffold_only = 1), guide = "none") +
    labs(title = title, color = NULL, shape = NULL) +
    theme_graph(base_family = "sans") +
    theme(legend.position = "none")
}

ppi_graph_from_tables <- function(nodes, edges, title) {
  g <- graph_from_data_frame(
    d = as.data.frame(edges),
    vertices = as.data.frame(nodes),
    directed = FALSE
  )
  tg <- as_tbl_graph(g) %>%
    activate(nodes) %>%
    mutate(deg = centrality_degree())
  plot_ppi_ggraph(tg, title)
}

k_proteins <- 4L
k_conditions <- 3L

prepare_emm_matrix <- function(emm_df) {
  wide <- emm_df %>%
    mutate(
      Day = factor(as.integer(Day), levels = c(1, 4, 7)),
      Allocation = factor(Allocation, levels = c("Placebo", "Anakinra")),
      srf14 = factor(srf14, levels = c("No", "Yes")),
      group_time = paste(Day, Allocation, srf14, sep = "_")
    ) %>%
    dplyr::select(Assay, group_time, emmean) %>%
    pivot_wider(names_from = group_time, values_from = emmean) %>%
    column_to_rownames("Assay")
  as.matrix(wide)
}

scale_protein_trajectories <- function(emm_mat) {
  scaled <- t(scale(t(emm_mat), center = TRUE, scale = TRUE))
  bad <- !is.finite(scaled)
  if (any(bad)) {
    warning("Non-finite values after row-scaling; setting them to 0")
    scaled[bad] <- 0
  }
  scaled
}

cluster_emm_heatmap <- function(scaled_traj, k_prot = k_proteins, k_cond = k_conditions) {
  heat_mat <- t(scaled_traj)
  heat <- pheatmap(
    heat_mat,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    cellwidth = 10,
    cellheight = 12,
    color = colorRampPalette(c("navyblue", "white", "firebrick"))(100),
    clustering_distance_rows = "euclidean",
    clustering_distance_cols = "euclidean",
    clustering_method = "ward.D2",
    cutree_rows = k_cond,
    cutree_cols = k_prot,
    border_color = "black",
    fontsize = 10,
    fontsize_row = 9,
    fontsize_col = 8,
    angle_col = 90,
    main = "Protein trajectories (row-scaled EMMs: Day x Allocation x SRF)",
    silent = TRUE
  )

  protein_clusters <- tibble(
    Assay = names(cutree(heat$tree_col, k = k_prot)),
    Cluster = as.integer(cutree(heat$tree_col, k = k_prot))
  )
  condition_clusters <- tibble(
    Condition = names(cutree(heat$tree_row, k = k_cond)),
    Cluster = as.integer(cutree(heat$tree_row, k = k_cond))
  )

  list(
    heat = heat,
    plot = as.ggplot(heat$gtable),
    protein_clusters = protein_clusters,
    condition_clusters = condition_clusters
  )
}
