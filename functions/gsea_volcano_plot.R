if (!base::exists('standardize_gsea_results', mode = 'function')) {
  base::stop('Source gsea_plot_utils.R before this plotting function file.', call. = FALSE)
}

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
#' @return A ggplot object.
make_gsea_volcano <- function(gsea_results,
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
  plot_tbl <- standardize_gsea_results(
    gsea_results = gsea_results,
    term_col = term_col,
    nes_col = nes_col,
    pvalue_col = pvalue_col,
    padj_col = padj_col,
    id_col = id_col)
  plot_tbl$neg_log10_padj <- safe_neg_log10(plot_tbl$padj)
  plot_tbl$plot_neg_log10_padj <- base::pmin(plot_tbl$neg_log10_padj, y_max)
  plot_tbl$significant <- !base::is.na(plot_tbl$padj) & plot_tbl$padj < padj_cutoff
  plot_tbl$point_group <- base::ifelse(
    plot_tbl$significant & plot_tbl$NES >= 0,
    'Significantly up',
    base::ifelse(plot_tbl$significant & plot_tbl$NES < 0, 'Significantly down', 'Not significant'))
  plot_tbl$point_group <- base::factor(
    plot_tbl$point_group,
    levels = base::c('Significantly up', 'Significantly down', 'Not significant'))
  color_values <- gsea_default_colors()[base::levels(plot_tbl$point_group)]
  color_values[['Not significant']] <- gsea_light_gray()
  custom_colors <- resolve_named_colors(
    color_values = point_colors,
    required_names = base::levels(plot_tbl$point_group),
    parameter_name = 'point_colors')
  if (!base::is.null(custom_colors)) {
    color_values <- custom_colors
  }
  plot_tbl$label_score <- base::abs(plot_tbl$NES) * plot_tbl$plot_neg_log10_padj
  label_requests <- combine_label_requests(label_terms = label_terms)
  if (!base::is.null(label_requests) && base::length(label_requests) > 0L) {
    label_tbl <- select_gsea_labels(plot_tbl, label_terms = label_requests)
  } else {
    label_tbl <- select_volcano_labels(plot_tbl, label_n = label_n)
  }
  x_limits <- base::c(-3.5, 3.5)
  y_limits <- base::c(y_min, y_max + 3.5)
  label_tbl$plot_neg_log10_padj <- base::pmin(label_tbl$plot_neg_log10_padj, y_max - 2.0)
  label_tbl$label_text <- wrap_gsea_label(label_tbl$go_description, words_per_line = label_words_per_line)
  repel_tbl <- plot_tbl
  repel_tbl$label_text <- ''
  repel_tbl$label_nudge_x <- 0
  repel_tbl$label_nudge_y <- 0
  label_match <- base::match(label_tbl$go_term_id, repel_tbl$go_term_id)
  repel_tbl$label_text[label_match] <- label_tbl$label_text
  repel_tbl$plot_neg_log10_padj[label_match] <- label_tbl$plot_neg_log10_padj
  # inward nudges use central open space and avoid clipped edge labels
  label_nudge_x_values <- base::ifelse(label_tbl$NES < 0, label_nudge_x, -label_nudge_x)
  repel_tbl$label_nudge_x[label_match] <- label_nudge_x_values
  repel_tbl$label_nudge_y[label_match] <- label_nudge_y +
    base::rep(base::c(0.18, -0.12, 0.08, -0.18), length.out = base::length(label_match))
  negative_count <- base::sum(plot_tbl$significant & plot_tbl$NES < 0, na.rm = TRUE)
  positive_count <- base::sum(plot_tbl$significant & plot_tbl$NES > 0, na.rm = TRUE)
  # place count boxes away from upper go labels
  count_y <- y_min + 0.35
  x_axis_label <- if (base::is.null(contrast_label)) {
    'Normalized enrichment score (NES)'
  } else {
    base::paste0(contrast_label, ' NES')
  }
  point_layers <- split_point_layers(
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
      xlim = base::c(x_limits[[1L]] + 0.18, x_limits[[2L]] - 0.18),
      ylim = base::c(y_min, y_max),
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
      label = base::paste0('Significantly down: ', negative_count),
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
      label = base::paste0('Significantly up: ', positive_count),
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
      breaks = base::c(-3, -2, -1, 0, 1, 2, 3),
      expand = ggplot2::expansion(mult = 0.02)) +
    ggplot2::scale_y_continuous(
      trans = scales::pseudo_log_trans(sigma = 1),
      breaks = base::c(0, 1, 2, 5, 10, 25),
      expand = base::c(0, 0)) +
    ggplot2::coord_cartesian(xlim = x_limits, ylim = y_limits, clip = 'off') +
    ggplot2::labs(x = x_axis_label, y = base::quote(-log[10]('Adjusted p-value'))) +
    ggplot2::guides(color = ggplot2::guide_legend(
      title = 'GSEA direction',
      nrow = 1,
      byrow = TRUE,
      title.position = 'top',
      title.hjust = 0.5,
      override.aes = base::list(size = point_size + 0.8))) +
    ggplot2::theme_classic(base_family = font_family, base_size = 9) +
    ggplot2::theme(
      axis.title = ggplot2::element_text(size = 9, color = 'black'),
      axis.text = ggplot2::element_text(size = 8, color = 'black'),
      axis.line = ggplot2::element_line(color = 'black', linewidth = 0.24),
      axis.ticks = ggplot2::element_line(color = 'black', linewidth = 0.20),
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = legend_position,
      legend.direction = 'horizontal',
      legend.box = 'horizontal',
      legend.justification = 'center',
      legend.title = ggplot2::element_text(size = 8, color = 'black', face = 'bold'),
      legend.text = ggplot2::element_text(size = 7.5, color = 'black'),
      legend.key.size = grid::unit(0.25, 'cm'),
      plot.margin = ggplot2::margin(10, 18, 8, 12))

  volcano_plot
}

plot_gsea_volcano <- make_gsea_volcano
