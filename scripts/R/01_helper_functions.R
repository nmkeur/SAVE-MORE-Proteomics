library(ggplot2)
library(EnhancedVolcano)

ylab_raw <- bquote(~-Log[10]~("p·value"))
ylab_adj <- bquote(~-Log[10]~("Adjusted p·value"))

BASE_SIZE <- 10
BASE_FAMILY <- "arial"

format_lancet <- function(x, digits = 1) {
  gsub("\\.", "·", formatC(x, format = "f", digits = digits))
}

make_volcano <- function(df, y, title, ylab, p_cut = 0.05) {
  is_dep <- !is.na(df[[y]]) & df[[y]] < p_cut
  col_dep <- "#D64541"
  col_ns  <- "#3A7DCE"
  keyvals <- ifelse(is_dep, col_dep, col_ns)
  names(keyvals) <- ifelse(is_dep, "DEP", "non-DEP")

  EnhancedVolcano(
    toptable        = df,
    x               = "log2FC",
    y               = y,
    lab             = df$Assay,
    xlab            = bquote(~Log[2]~fold~change),
    ylab            = ylab,
    pCutoff         = p_cut,
    FCcutoff        = 0,
    xlim            = c(-0.6, 0.6),
    ylim            = c(0, 5),
    title           = title,
    subtitle        = NULL,
    caption         = NULL,
    titleLabSize    = 10,
    boxedLabels     = TRUE,
    drawConnectors  = TRUE,
    widthConnectors = 0.4,
    colConnectors   = "grey40",
    labSize         = 3.2,
    pointSize       = 1.6,
    colAlpha        = 0.8,
    shape           = 16,
    colCustom       = keyvals,
    legendPosition  = "right",
    legendLabSize   = BASE_SIZE * 0.75,
    legendIconSize  = 1,
    gridlines.major = TRUE,
    gridlines.minor = FALSE,
    border          = "partial",
    raster          = FALSE
  ) +
    theme(
      legend.text  = element_text(size = BASE_SIZE * 0.8, family = BASE_FAMILY),
      legend.title = element_blank(),
      axis.title   = element_text(size = BASE_SIZE * 0.85),
      axis.text    = element_text(size = BASE_SIZE * 0.8),
      axis.text.x  = element_text(size = BASE_SIZE * 0.8),
      axis.text.y  = element_text(size = BASE_SIZE * 0.8),
      legend.spacing.x = unit(2, "mm"),
      legend.key.size  = unit(3, "mm")
    ) +
    scale_x_continuous(
      limits = c(-0.6, 0.6),
      breaks = seq(-0.6, 0.6, by = 0.2),
      labels = function(x) format_lancet(x, 1),
      expand = expansion(mult = 0)
    ) +
    scale_colour_manual(
      values = c("non-DEP" = col_ns, "DEP" = col_dep),
      breaks = c("non-DEP", "DEP"),
      drop = FALSE,
      name = NULL
    )
}
