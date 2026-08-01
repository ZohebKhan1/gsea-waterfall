# ----
# author:
# - Zoheb Khan
#
# script path:
# - functions/gsea_volcano_plot.R
#
# functions:
# - functions/gsea_plot_utils.R
# ----

# 1.0 create symmetric volcano plot -----------------

#' Make a symmetric GSEA volcano plot
#'
#' @param gsea_results GSEA result table.
#' @param term_col Column containing GO term names.
#' @param nes_col Column containing normalized enrichment scores.
#' @param pvalue_col Column containing nominal GSEA p-values.
#' @param padj_col Column containing adjusted p-values.
#' @param id_col Optional column containing GO IDs.
#' @param padj_cutoff Adjusted p-value cutoff for coloring and shaded regions.
#' @param label_n Number of significant terms to label.
#' @param label_terms Optional exact term descriptions or GO IDs to label.
#' @param y_max Upper y-axis limit for capped -log10 adjusted p-values.
#' @param y_min Lower y-axis limit for capped -log10 adjusted p-values.
#' @param point_size Point size.
#' @param point_colors Optional named colors for `Significantly up`,
#'   `Significantly down`, and `Not significant`.
#' @param label_size Label font size in points.
#' @param label_color Text color for GO term labels.
#' @param label_nudge_x Horizontal starting offset for labeled terms.
#' @param label_nudge_y Vertical starting offset for labeled terms.
#' @param count_label_size Font size in points for significant-count labels.
#' @param contrast_label Optional concise contrast label for the NES axis.
#' @param label_words_per_line Number of GO term words shown on each label line.
#' @param legend_position Legend position passed to `ggplot2::theme()`.
#' @param font_family Figure font family.
#'
#' @return A ggplot object without file-system side effects.
#'
#' @export
plot_gsea_volcano <- function(gsea_results,
                              term_col = 'go_description',
                              nes_col = 'NES',
                              pvalue_col = 'pval',
                              padj_col = 'padj',
                              id_col = NULL,
                              padj_cutoff = 0.05,
                              label_n = 14L,
                              label_terms = NULL,
                              y_max = 25,
                              y_min = 0,
                              point_size = 1.0,
                              point_colors = NULL,
                              label_size = 7.5,
                              label_color = 'black',
                              label_nudge_x = 0.18,
                              label_nudge_y = 0,
                              count_label_size = 9,
                              contrast_label = NULL,
                              label_words_per_line = 3L,
                              legend_position = 'top',
                              font_family = 'Nimbus Sans') {
  padj_cutoff <- .gsea_validate_number(
    padj_cutoff,
    'padj_cutoff',
    minimum = 0,
    maximum = 1)
  label_n <- .gsea_validate_count(label_n, 'label_n')
  point_size <- .gsea_validate_number(
    point_size,
    'point_size',
    minimum = 0,
    minimum_inclusive = FALSE)
  label_size <- .gsea_validate_number(
    label_size,
    'label_size',
    minimum = 0,
    minimum_inclusive = FALSE)
  count_label_size <- .gsea_validate_number(
    count_label_size,
    'count_label_size',
    minimum = 0,
    minimum_inclusive = FALSE)
  label_nudge_x <- .gsea_validate_number(label_nudge_x, 'label_nudge_x')
  label_nudge_y <- .gsea_validate_number(label_nudge_y, 'label_nudge_y')
  .gsea_validate_limits(y_min, y_max, 'y_min', 'y_max')
  plot_tbl <- .gsea_standardize_results(
    gsea_results = gsea_results,
    term_col = term_col,
    nes_col = nes_col,
    pvalue_col = pvalue_col,
    padj_col = padj_col,
    id_col = id_col)
  plot_tbl$neg_log10_padj <- .gsea_neg_log10(plot_tbl$padj)
  plot_tbl$plot_neg_log10_padj <- pmin(plot_tbl$neg_log10_padj, y_max)
  plot_tbl$significant <- !is.na(plot_tbl$padj) & plot_tbl$padj < padj_cutoff
  plot_tbl$point_group <- ifelse(
    plot_tbl$significant & plot_tbl$NES >= 0,
    'Significantly up',
    ifelse(plot_tbl$significant & plot_tbl$NES < 0, 'Significantly down', 'Not significant'))
  plot_tbl$point_group <- factor(
    plot_tbl$point_group,
    levels = c('Significantly up', 'Significantly down', 'Not significant'))
  color_values <- .gsea_default_colors()[levels(plot_tbl$point_group)]
  color_values[['Not significant']] <- .gsea_light_gray()
  custom_colors <- .gsea_resolve_named_colors(
    color_values = point_colors,
    required_names = levels(plot_tbl$point_group),
    parameter_name = 'point_colors')
  if (!is.null(custom_colors)) {
    color_values <- custom_colors
  }
  plot_tbl$label_score <- abs(plot_tbl$NES) * plot_tbl$plot_neg_log10_padj
  label_requests <- .gsea_combine_label_requests(label_terms = label_terms)
  if (!is.null(label_requests) && length(label_requests) > 0L) {
    label_tbl <- .gsea_select_labels(plot_tbl, label_terms = label_requests)
  } else {
    label_tbl <- .gsea_select_volcano_labels(plot_tbl, label_n = label_n)
  }
  x_limits <- c(-3.5, 3.5)
  y_limits <- c(y_min, y_max + 3.5)
  label_tbl$plot_neg_log10_padj <- pmin(label_tbl$plot_neg_log10_padj, y_max - 2.0)
  label_tbl$label_text <- .gsea_wrap_label(label_tbl$go_description, words_per_line = label_words_per_line)
  repel_tbl <- plot_tbl
  repel_tbl$label_text <- ''
  repel_tbl$label_nudge_x <- 0
  repel_tbl$label_nudge_y <- 0
  label_match <- match(label_tbl$go_term_id, repel_tbl$go_term_id)
  repel_tbl$label_text[label_match] <- label_tbl$label_text
  repel_tbl$plot_neg_log10_padj[label_match] <- label_tbl$plot_neg_log10_padj
  # inward nudges use central open space and avoid clipped edge labels
  label_nudge_x_values <- ifelse(label_tbl$NES < 0, label_nudge_x, -label_nudge_x)
  repel_tbl$label_nudge_x[label_match] <- label_nudge_x_values
  repel_tbl$label_nudge_y[label_match] <- label_nudge_y +
    rep(c(0.18, -0.12, 0.08, -0.18), length.out = length(label_match))
  negative_count <- sum(plot_tbl$significant & plot_tbl$NES < 0, na.rm = TRUE)
  positive_count <- sum(plot_tbl$significant & plot_tbl$NES > 0, na.rm = TRUE)
  # place count boxes away from upper GO labels
  count_y <- y_min + 0.35
  x_axis_label <- if (is.null(contrast_label)) {
    'Normalized enrichment score (NES)'
  } else {
    paste0(contrast_label, ' NES')
  }
  point_layers <- .gsea_split_point_layers(
    plot_tbl = plot_tbl,
    group_col = 'point_group',
    background_groups = 'Not significant')

  volcano_plot <- ggplot2::ggplot(plot_tbl, ggplot2::aes(x = .data$NES, y = .data$plot_neg_log10_padj))
  volcano_plot <- volcano_plot +
    ggplot2::geom_point(
      data = point_layers$background,
      ggplot2::aes(color = .data$point_group),
      size = point_size,
      alpha = 1,
      stroke = 0) +
    ggplot2::geom_point(
      data = point_layers$foreground,
      ggplot2::aes(color = .data$point_group),
      size = point_size,
      alpha = 1,
      stroke = 0) +
    ggrepel::geom_text_repel(
      data = repel_tbl,
      ggplot2::aes(label = .data$label_text),
      nudge_x = repel_tbl$label_nudge_x,
      nudge_y = repel_tbl$label_nudge_y,
      family = font_family,
      color = label_color,
      size = label_size / ggplot2::.pt,
      lineheight = 0.9,
      min.segment.length = 0,
      segment.color = '#4D4D4D',
      segment.size = 0.13,
      box.padding = 0.45,
      point.padding = 0.16,
      xlim = c(x_limits[[1L]] + 0.18, x_limits[[2L]] - 0.18),
      ylim = c(y_min, y_max),
      direction = 'both',
      force = 1.8,
      force_pull = 3.4,
      max.time = 8,
      max.iter = 20000,
      max.overlaps = Inf,
      seed = 42,
      show.legend = FALSE) +
    ggplot2::annotate(
      'label',
      x = x_limits[[1L]],
      y = count_y,
      label = paste0('Significantly down: ', negative_count),
      hjust = 0,
      vjust = 0,
      family = font_family,
      size = count_label_size / ggplot2::.pt,
      fontface = 'bold',
      color = 'white',
      fill = color_values[['Significantly down']],
      linewidth = 0,
      label.padding = grid::unit(0.12, 'lines'),
      label.r = grid::unit(0.08, 'lines')) +
    ggplot2::annotate(
      'label',
      x = x_limits[[2L]],
      y = count_y,
      label = paste0('Significantly up: ', positive_count),
      hjust = 1,
      vjust = 0,
      family = font_family,
      size = count_label_size / ggplot2::.pt,
      fontface = 'bold',
      color = 'white',
      fill = color_values[['Significantly up']],
      linewidth = 0,
      label.padding = grid::unit(0.12, 'lines'),
      label.r = grid::unit(0.08, 'lines')) +
    ggplot2::scale_color_manual(values = color_values, drop = FALSE) +
    ggplot2::scale_x_continuous(
      trans = scales::pseudo_log_trans(sigma = 1),
      breaks = c(-3, -2, -1, 0, 1, 2, 3),
      expand = ggplot2::expansion(mult = 0.02)) +
    ggplot2::scale_y_continuous(
      trans = scales::pseudo_log_trans(sigma = 1),
      breaks = c(0, 1, 2, 5, 10, 25),
      expand = c(0, 0)) +
    ggplot2::coord_cartesian(xlim = x_limits, ylim = y_limits, clip = 'off') +
    ggplot2::labs(x = x_axis_label, y = '\u2212log10(Adjusted p-value)') +
    ggplot2::guides(color = ggplot2::guide_legend(
      title = 'GSEA direction',
      nrow = 1,
      byrow = TRUE,
      title.position = 'top',
      title.hjust = 0.5,
      override.aes = list(size = point_size + 0.8))) +
    .gsea_theme(font_family) +
    ggplot2::theme(
      axis.title.y = ggplot2::element_text(margin = ggplot2::margin(r = 5.42)),
      legend.position = legend_position,
      legend.direction = 'horizontal',
      legend.box = 'horizontal',
      legend.justification = 'center',
      plot.margin = ggplot2::margin(10, 18, 8, 12))

  volcano_plot
}
