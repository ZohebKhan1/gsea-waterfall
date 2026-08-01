# ----
# author:
# - Zoheb Khan
#
# script path:
# - functions/gsea_half_volcano_plot.R
#
# functions:
# - functions/gsea_plot_utils.R
# ----

# 1.0 create directional volcano plot -----------------

#' Make a directional half-volcano plot for GSEA results
#'
#' @param gsea_results GSEA result table.
#' @param term_col Column containing GO term names.
#' @param nes_col Column containing normalized enrichment scores.
#' @param pvalue_col Column containing nominal GSEA p-values.
#' @param padj_col Column containing adjusted p-values.
#' @param id_col Optional column containing GO IDs.
#' @param direction NES direction to plot.
#' @param p_col P-value column for the y-axis. Use `pvalue` or `padj`.
#' @param padj_cutoff Adjusted p-value cutoff for significance coloring.
#' @param label_n Number of significant terms to label.
#' @param label_terms Optional exact term descriptions or GO IDs to label.
#' @param term_groups Optional named GO category seeds for coloring.
#' @param term_group_colors Optional named colors for `term_groups`.
#' @param inner_nes_limit Inner absolute NES cutoff shown on the x-axis.
#' @param point_size Point size.
#' @param point_colors Optional named colors for the direction label and
#'   `Not significant`.
#' @param label_size Label font size in points.
#' @param label_color Text color for GO term labels.
#' @param contrast_label Optional concise contrast label for the NES axis.
#' @param label_words_per_line Number of GO term words shown on each label line.
#' @param y_min Optional lower y-axis limit.
#' @param y_max Optional upper y-axis limit.
#' @param legend_position Legend position passed to `ggplot2::theme()`.
#' @param font_family Figure font family.
#'
#' @return A ggplot object without file-system side effects.
#'
#' @export
plot_gsea_half_volcano <- function(gsea_results,
                                   term_col = 'go_description',
                                   nes_col = 'NES',
                                   pvalue_col = 'pval',
                                   padj_col = 'padj',
                                   id_col = NULL,
                                   direction = c('positive', 'negative'),
                                   p_col = 'padj',
                                   padj_cutoff = 0.05,
                                   label_n = 12L,
                                   label_terms = NULL,
                                   term_groups = NULL,
                                   term_group_colors = NULL,
                                   inner_nes_limit = 1,
                                   point_size = 1.35,
                                   point_colors = NULL,
                                   label_size = 7.5,
                                   label_color = 'black',
                                   contrast_label = NULL,
                                   label_words_per_line = 3L,
                                   y_min = NULL,
                                   y_max = NULL,
                                   legend_position = NULL,
                                   font_family = 'Nimbus Sans') {
  direction <- match.arg(direction)
  padj_cutoff <- .gsea_validate_number(
    padj_cutoff,
    'padj_cutoff',
    minimum = 0,
    maximum = 1)
  label_n <- .gsea_validate_count(label_n, 'label_n')
  inner_nes_limit <- .gsea_validate_number(
    inner_nes_limit,
    'inner_nes_limit',
    minimum = 0)
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
  .gsea_validate_limits(y_min, y_max, 'y_min', 'y_max')
  .gsea_validate_term_groups(term_groups)
  plot_tbl <- .gsea_standardize_results(
    gsea_results = gsea_results,
    term_col = term_col,
    nes_col = nes_col,
    pvalue_col = pvalue_col,
    padj_col = padj_col,
    id_col = id_col)
  if (!p_col %in% c('pvalue', 'padj')) {
    stop("`p_col` must be either 'pvalue' or 'padj'.", call. = FALSE)
  }

  plot_tbl <- plot_tbl[if (direction == 'positive') plot_tbl$NES > 0 else plot_tbl$NES < 0, , drop = FALSE]
  if (inner_nes_limit > 0) {
    plot_tbl <- plot_tbl[if (direction == 'positive') plot_tbl$NES >= inner_nes_limit else plot_tbl$NES <= -inner_nes_limit, , drop = FALSE]
  }
  if (nrow(plot_tbl) == 0L) {
    stop('No GSEA terms remain after applying the half-volcano NES boundary.', call. = FALSE)
  }
  plot_tbl$neg_log10_p <- .gsea_neg_log10(plot_tbl[[p_col]])
  plot_tbl$significant <- !is.na(plot_tbl$padj) & plot_tbl$padj < padj_cutoff
  if (is.null(term_groups)) {
    direction_label <- if (direction == 'positive') 'Significantly up' else 'Significantly down'
    direction_color <- if (direction == 'positive') '#CC79A7' else '#0072B2'
    plot_tbl$point_group <- ifelse(plot_tbl$significant, direction_label, 'Not significant')
    plot_tbl$point_group <- factor(plot_tbl$point_group, levels = c(direction_label, 'Not significant'))
    color_values <- c(stats::setNames(direction_color, direction_label), 'Not significant' = .gsea_light_gray())
    custom_colors <- .gsea_resolve_named_colors(
      color_values = point_colors,
      required_names = levels(plot_tbl$point_group),
      parameter_name = 'point_colors')
    if (!is.null(custom_colors)) {
      color_values <- custom_colors
    }
    legend_title <- 'GSEA direction'
    background_groups <- 'Not significant'
  } else {
    plot_tbl$point_group <- .gsea_assign_term_groups(
      gsea_results = plot_tbl,
      term_groups = term_groups)
    plot_tbl$point_group <- as.character(plot_tbl$point_group)
    plot_tbl$point_group[!plot_tbl$significant] <- 'Not significant'
    plot_tbl$point_group <- factor(
      plot_tbl$point_group,
      levels = c(names(term_groups), 'Other', 'Not significant'))
    color_values <- .gsea_resolve_term_group_colors(
      term_groups = term_groups,
      term_group_colors = term_group_colors)
    color_values <- c(color_values, 'Not significant' = .gsea_light_gray())
    legend_title <- 'GO:BP category'
    background_groups <- c('Other', 'Not significant')
  }
  plot_tbl$label_score <- abs(plot_tbl$NES) * plot_tbl$neg_log10_p
  significant_tbl <- plot_tbl[plot_tbl$significant, , drop = FALSE]
  significant_tbl <- significant_tbl[order(abs(significant_tbl$NES), significant_tbl$padj), , drop = FALSE]
  significant_tbl$label_rank <- seq_len(nrow(significant_tbl))
  label_requests <- .gsea_combine_label_requests(label_terms = label_terms)
  if (!is.null(label_requests) && length(label_requests) > 0L) {
    label_tbl <- .gsea_select_labels(plot_tbl, label_terms = label_requests)
  } else {
    label_tbl <- .gsea_select_labels(significant_tbl, label_n = label_n, score_col = 'label_score', rank_col = 'label_rank')
  }
  label_tbl$label_text <- .gsea_wrap_label(label_tbl$go_description, words_per_line = label_words_per_line)
  repel_tbl <- plot_tbl
  repel_tbl$label_text <- ''
  repel_tbl$label_nudge_x <- 0
  repel_tbl$label_nudge_y <- 0
  label_match <- match(label_tbl$go_term_id, repel_tbl$go_term_id)
  repel_tbl$label_text[label_match] <- label_tbl$label_text
  label_direction <- if (direction == 'positive') 1 else -1
  repel_tbl$label_nudge_x[label_match] <- label_direction * 0.12
  y_limits <- .gsea_waterfall_y_limits(plot_tbl$neg_log10_p, y_min = y_min, y_max = y_max, buffer_fraction = 0.035)
  x_limits <- .gsea_axis_limits(plot_tbl$NES, buffer_fraction = 0.035)
  if (direction == 'positive') {
    x_limits[[1L]] <- inner_nes_limit
  } else {
    x_limits[[2L]] <- -inner_nes_limit
  }
  x_breaks <- scales::breaks_pretty(n = 4)(x_limits)
  x_breaks <- x_breaks[x_breaks >= x_limits[[1L]] & x_breaks <= x_limits[[2L]]]
  use_default_legend <- is.null(legend_position)
  if (use_default_legend) {
    legend_position <- if (direction == 'negative') c(0.97, 0.97) else c(0.03, 0.97)
  }
  legend_anchor <- if (use_default_legend && direction == 'negative') {
    c(1, 1)
  } else {
    c(0, 1)
  }
  x_axis_label <- if (is.null(contrast_label)) {
    'Normalized enrichment score (NES)'
  } else {
    paste0(contrast_label, ' NES')
  }
  point_layers <- .gsea_split_point_layers(
    plot_tbl = plot_tbl,
    group_col = 'point_group',
    background_groups = background_groups)

  half_plot <- ggplot2::ggplot(plot_tbl, ggplot2::aes(x = .data$NES, y = .data$neg_log10_p)) +
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
    ggplot2::geom_hline(
      yintercept = -log10(padj_cutoff),
      color = 'black',
      linewidth = 0.22,
      linetype = 'dashed') +
    ggrepel::geom_text_repel(
      data = repel_tbl,
      ggplot2::aes(label = .data$label_text),
      nudge_x = repel_tbl$label_nudge_x,
      nudge_y = repel_tbl$label_nudge_y,
      family = font_family,
      color = label_color,
      size = label_size / ggplot2::.pt,
      lineheight = 0.9,
      hjust = 0.5,
      min.segment.length = 0,
      segment.color = 'black',
      segment.size = 0.16,
      box.padding = 0.82,
      point.padding = 0.08,
      force = 3.0,
      force_pull = 0.7,
      max.time = 8,
      max.iter = 22000,
      max.overlaps = Inf,
      seed = 42,
      show.legend = FALSE) +
    ggplot2::scale_color_manual(values = color_values, drop = FALSE) +
    ggplot2::scale_x_continuous(breaks = x_breaks, expand = ggplot2::expansion(mult = 0.015)) +
    ggplot2::scale_y_continuous(expand = c(0, 0), position = if (direction == 'negative') 'right' else 'left') +
    ggplot2::coord_cartesian(xlim = x_limits, ylim = y_limits, clip = 'off') +
    ggplot2::labs(x = x_axis_label, y = paste0('-log10(', p_col, ')'), color = legend_title) +
    ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(size = point_size + 0.5))) +
    .gsea_theme(font_family) +
    ggplot2::theme(
      legend.position = legend_position,
      legend.justification = legend_anchor,
      legend.background = ggplot2::element_rect(fill = 'white', color = NA),
      plot.margin = ggplot2::margin(8, 10, 8, 10))

  half_plot
}
