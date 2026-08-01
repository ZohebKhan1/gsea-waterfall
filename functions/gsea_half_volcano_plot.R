if (!base::exists('standardize_gsea_results', mode = 'function')) {
  base::stop('Source gsea_plot_utils.R before this plotting function file.', call. = FALSE)
}

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
#' @param category_n Optional expected number of categories in `term_groups`.
#' @param max_categories Maximum allowed number of custom categories.
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
#' @return A ggplot object.
make_gsea_half_volcano <- function(gsea_results,
                                   term_col = 'go_description',
                                   nes_col = 'NES',
                                   pvalue_col = 'pval',
                                   padj_col = 'padj',
                                   id_col = NULL,
                                   direction = base::c('positive', 'negative'),
                                   p_col = 'padj',
                                   padj_cutoff = 0.05,
                                   label_n = 12L,
                                   label_terms = NULL,
                                   term_groups = NULL,
                                   category_n = NULL,
                                   max_categories = 6L,
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
  direction <- base::match.arg(direction)
  validate_term_groups(term_groups = term_groups, category_n = category_n)
  validate_max_categories(term_groups = term_groups, max_categories = max_categories)
  plot_tbl <- standardize_gsea_results(
    gsea_results = gsea_results,
    term_col = term_col,
    nes_col = nes_col,
    pvalue_col = pvalue_col,
    padj_col = padj_col,
    id_col = id_col)
  if (!p_col %in% base::c('pvalue', 'padj')) {
    base::stop("`p_col` must be either 'pvalue' or 'padj'.", call. = FALSE)
  }

  plot_tbl <- plot_tbl[if (direction == 'positive') plot_tbl$NES > 0 else plot_tbl$NES < 0, , drop = FALSE]
  inner_nes_limit <- base::as.numeric(inner_nes_limit)
  if (base::is.na(inner_nes_limit) || inner_nes_limit < 0) {
    base::stop('`inner_nes_limit` must be a non-negative number.', call. = FALSE)
  }
  if (inner_nes_limit > 0) {
    plot_tbl <- plot_tbl[if (direction == 'positive') plot_tbl$NES >= inner_nes_limit else plot_tbl$NES <= -inner_nes_limit, , drop = FALSE]
  }
  if (base::nrow(plot_tbl) == 0L) {
    base::stop('No GSEA terms remain after applying the half-volcano NES boundary.', call. = FALSE)
  }
  plot_tbl$neg_log10_p <- safe_neg_log10(plot_tbl[[p_col]])
  plot_tbl$significant <- !base::is.na(plot_tbl$padj) & plot_tbl$padj < padj_cutoff
  if (base::is.null(term_groups)) {
    direction_label <- if (direction == 'positive') 'Significantly up' else 'Significantly down'
    direction_color <- if (direction == 'positive') '#CC79A7' else '#0072B2'
    plot_tbl$point_group <- base::ifelse(plot_tbl$significant, direction_label, 'Not significant')
    plot_tbl$point_group <- base::factor(plot_tbl$point_group, levels = base::c(direction_label, 'Not significant'))
    color_values <- base::c(stats::setNames(direction_color, direction_label), 'Not significant' = gsea_light_gray())
    custom_colors <- resolve_named_colors(
      color_values = point_colors,
      required_names = base::levels(plot_tbl$point_group),
      parameter_name = 'point_colors')
    if (!base::is.null(custom_colors)) {
      color_values <- custom_colors
    }
    legend_title <- 'GSEA direction'
    background_groups <- 'Not significant'
  } else {
    plot_tbl$point_group <- assign_plot_term_groups(
      gsea_results = plot_tbl,
      term_groups = term_groups)
    plot_tbl$point_group <- base::as.character(plot_tbl$point_group)
    plot_tbl$point_group[!plot_tbl$significant] <- 'Not significant'
    plot_tbl$point_group <- base::factor(
      plot_tbl$point_group,
      levels = base::c(base::names(term_groups), 'Other', 'Not significant'))
    color_values <- resolve_term_group_colors(
      term_groups = term_groups,
      term_group_colors = term_group_colors)
    color_values <- base::c(color_values, 'Not significant' = gsea_light_gray())
    legend_title <- 'GO:BP category'
    background_groups <- base::c('Other', 'Not significant')
  }
  plot_tbl$label_score <- base::abs(plot_tbl$NES) * plot_tbl$neg_log10_p
  significant_tbl <- plot_tbl[plot_tbl$significant, , drop = FALSE]
  significant_tbl <- significant_tbl[base::order(base::abs(significant_tbl$NES), significant_tbl$padj), , drop = FALSE]
  significant_tbl$label_rank <- base::seq_len(base::nrow(significant_tbl))
  label_requests <- combine_label_requests(label_terms = label_terms)
  if (!base::is.null(label_requests) && base::length(label_requests) > 0L) {
    label_tbl <- select_gsea_labels(plot_tbl, label_terms = label_requests)
  } else {
    label_tbl <- select_gsea_labels(significant_tbl, label_n = label_n, score_col = 'label_score', rank_col = 'label_rank')
  }
  label_tbl$label_text <- wrap_gsea_label(label_tbl$go_description, words_per_line = label_words_per_line)
  repel_tbl <- plot_tbl
  repel_tbl$label_text <- ''
  repel_tbl$label_nudge_x <- 0
  repel_tbl$label_nudge_y <- 0
  label_match <- base::match(label_tbl$go_term_id, repel_tbl$go_term_id)
  repel_tbl$label_text[label_match] <- label_tbl$label_text
  label_direction <- if (direction == 'positive') 1 else -1
  repel_tbl$label_nudge_x[label_match] <- label_direction * 0.12
  y_limits <- calculate_waterfall_y_limits(plot_tbl$neg_log10_p, y_min = y_min, y_max = y_max, buffer_fraction = 0.035)
  x_limits <- calculate_axis_limits(plot_tbl$NES, buffer_fraction = 0.035)
  if (direction == 'positive') {
    x_limits[[1L]] <- inner_nes_limit
  } else {
    x_limits[[2L]] <- -inner_nes_limit
  }
  x_breaks <- scales::breaks_pretty(n = 4)(x_limits)
  x_breaks <- x_breaks[x_breaks >= x_limits[[1L]] & x_breaks <= x_limits[[2L]]]
  use_default_legend <- base::is.null(legend_position)
  if (use_default_legend) {
    legend_position <- if (direction == 'negative') base::c(0.97, 0.97) else base::c(0.03, 0.97)
  }
  legend_anchor <- if (use_default_legend && direction == 'negative') {
    base::c(1, 1)
  } else {
    base::c(0, 1)
  }
  x_axis_label <- if (base::is.null(contrast_label)) {
    'Normalized enrichment score (NES)'
  } else {
    base::paste0(contrast_label, ' NES')
  }
  point_layers <- split_point_layers(
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
      yintercept = -base::log10(padj_cutoff),
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
    ggplot2::scale_y_continuous(expand = base::c(0, 0), position = if (direction == 'negative') 'right' else 'left') +
    ggplot2::coord_cartesian(xlim = x_limits, ylim = y_limits, clip = 'off') +
    ggplot2::labs(x = x_axis_label, y = base::paste0('-log10(', p_col, ')'), color = legend_title) +
    ggplot2::guides(color = ggplot2::guide_legend(override.aes = base::list(size = point_size + 0.5))) +
    ggplot2::theme_classic(base_family = font_family, base_size = 9) +
    ggplot2::theme(
      axis.title = ggplot2::element_text(size = 9, color = 'black'),
      axis.text = ggplot2::element_text(size = 8, color = 'black'),
      axis.line = ggplot2::element_line(color = 'black', linewidth = 0.24),
      axis.ticks = ggplot2::element_line(color = 'black', linewidth = 0.20),
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = legend_position,
      legend.justification = legend_anchor,
      legend.background = ggplot2::element_rect(fill = 'white', color = NA),
      legend.title = ggplot2::element_text(size = 8, color = 'black', face = 'bold'),
      legend.text = ggplot2::element_text(size = 7.5, color = 'black'),
      legend.key.size = grid::unit(0.25, 'cm'),
      plot.margin = ggplot2::margin(8, 10, 8, 10))

  half_plot
}

plot_gsea_half_volcano <- make_gsea_half_volcano
