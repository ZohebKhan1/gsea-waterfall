if (!base::exists('standardize_gsea_results', mode = 'function')) {
  base::stop('Source gsea_plot_utils.R before this plotting function file.', call. = FALSE)
}

#' Make a GSEA waterfall plot
#'
#' @param gsea_results GSEA result table.
#' @param term_col Column containing GO term names.
#' @param nes_col Column containing normalized enrichment scores.
#' @param pvalue_col Column containing nominal GSEA p-values.
#' @param padj_col Column containing adjusted p-values.
#' @param id_col Optional column containing GO IDs.
#' @param direction Direction to plot.
#' @param top_n Number of terms to show.
#' @param rank_by Column used to rank terms before taking `top_n`. Use `NES`,
#'   `padj`, `pvalue`, or `pval`.
#' @param x_label Optional x-axis title. Defaults to a rank-by label.
#' @param label_terms Optional exact term descriptions or GO IDs to label.
#' @param label_ranks Optional plotted ranks to label after ranking and
#'   filtering.
#' @param label_n Number of labels to automatically select across the ranked
#'   terms when `label_terms` is not supplied.
#' @param term_groups Optional named GO category seeds for coloring.
#' @param category_n Optional expected number of categories in `term_groups`;
#'   automatic term-group colors are generated for this category count.
#' @param max_categories Maximum allowed number of custom categories.
#' @param term_group_colors Optional named colors for `term_groups`.
#' @param label_size Label font size in points.
#' @param label_fontface Label font face.
#' @param label_words_per_line Number of GO term words shown on each label line.
#' @param label_repel_direction Direction passed to `ggrepel::geom_text_repel`.
#'   Use `y` for mostly vertical leader lines.
#' @param label_x_nudge Horizontal label nudge in rank units. The default
#'   `NULL` adds small alternating nudges so leader lines are not all vertical.
#' @param label_y_nudge_fraction Vertical label nudge as a fraction of NES span.
#' @param y_min Optional lower y-axis limit.
#' @param y_max Optional upper y-axis limit.
#' @param legend_position Legend position passed to `ggplot2::theme()`. Use
#'   strings such as `top`, `bottom`, `left`, `right`, `none`, or numeric
#'   coordinates such as `c(0.98, 0.98)`.
#' @param font_family Figure font family.
#'
#' @return A ggplot object.
make_gsea_waterfall <- function(gsea_results,
                                term_col = 'go_description',
                                nes_col = 'NES',
                                pvalue_col = 'pval',
                                padj_col = 'padj',
                                id_col = NULL,
                                direction = base::c('positive', 'negative'),
                                top_n = 100L,
                                rank_by = base::c('NES', 'padj', 'pvalue', 'pval'),
                                x_label = NULL,
                                label_terms = NULL,
                                label_ranks = NULL,
                                label_n = 12L,
                                term_groups = NULL,
                                category_n = NULL,
                                max_categories = 6L,
                                term_group_colors = NULL,
                                label_size = 7.5,
                                label_fontface = 'plain',
                                label_words_per_line = 2L,
                                label_repel_direction = 'y',
                                label_x_nudge = NULL,
                                label_y_nudge_fraction = 0.11,
                                y_min = NULL,
                                y_max = NULL,
                                legend_position = NULL,
                                font_family = 'Nimbus Sans') {
  direction <- base::match.arg(direction)
  rank_by <- normalize_rank_by(rank_by)
  gsea_results <- standardize_gsea_results(
    gsea_results = gsea_results,
    term_col = term_col,
    nes_col = nes_col,
    pvalue_col = pvalue_col,
    padj_col = padj_col,
    id_col = id_col)
  validate_term_groups(term_groups = term_groups, category_n = category_n)
  validate_max_categories(term_groups = term_groups, max_categories = max_categories)

  if (direction == 'positive') {
    plot_tbl <- gsea_results[gsea_results$NES > 0, , drop = FALSE]
  } else {
    plot_tbl <- gsea_results[gsea_results$NES < 0, , drop = FALSE]
  }

  if (base::nrow(plot_tbl) == 0L) {
    base::stop('No GSEA terms remain for the requested direction.', call. = FALSE)
  }

  plot_tbl <- rank_gsea_direction_terms(plot_tbl, direction = direction, rank_by = rank_by)
  plot_tbl <- utils::head(plot_tbl, top_n)
  plot_tbl$waterfall_rank <- base::seq_len(base::nrow(plot_tbl))
  plot_tbl$label_score <- base::abs(plot_tbl$NES) * safe_neg_log10(plot_tbl$padj)
  plot_tbl$term_category <- assign_plot_term_groups(
    gsea_results = plot_tbl,
    term_groups = term_groups)
  color_values <- resolve_term_group_colors(
    term_groups = term_groups,
    term_group_colors = term_group_colors)
  label_requests <- combine_label_requests(
    label_terms = label_terms,
    label_ranks = label_ranks)
  label_tbl <- select_gsea_labels(
    plot_tbl = plot_tbl,
    label_terms = label_requests,
    label_n = label_n,
    score_col = 'label_score',
    rank_col = 'waterfall_rank')
  label_tbl$label_text <- wrap_gsea_label(label_tbl$go_description, words_per_line = label_words_per_line)
  y_limits <- calculate_waterfall_y_limits(plot_tbl$NES, y_min = y_min, y_max = y_max)
  y_breaks <- scales::breaks_pretty(n = 4)(y_limits)
  label_tbl$label_nudge_x <- numeric(base::nrow(label_tbl))
  label_tbl$label_nudge_y <- numeric(base::nrow(label_tbl))
  if (base::nrow(label_tbl) > 0L) {
    label_tbl <- label_tbl[base::order(label_tbl$waterfall_rank), , drop = FALSE]
    nudge_step <- base::max(base::diff(base::range(plot_tbl$NES, na.rm = TRUE)) * label_y_nudge_fraction, 0.08)
    if (base::is.null(label_x_nudge)) {
      label_tbl$label_nudge_x <- base::rep(base::c(-0.65, 0.35, 0, 0.6, -0.3), length.out = base::nrow(label_tbl))
    } else {
      label_tbl$label_nudge_x <- base::rep(label_x_nudge, length.out = base::nrow(label_tbl))
    }
    label_tbl$label_nudge_y <- base::rep(base::c(nudge_step, -nudge_step), length.out = base::nrow(label_tbl))
    if (direction == 'negative') {
      label_tbl$label_nudge_y <- -label_tbl$label_nudge_y
    }
  }
  x_breaks <- base::unique(base::c(0, scales::breaks_pretty(n = 4)(base::c(1, base::nrow(plot_tbl)))))
  if (base::is.null(x_label)) {
    x_label <- base::paste0('Ranked GO terms (by ', rank_by, ')')
  }
  if (base::is.null(legend_position)) {
    legend_position <- if (direction == 'positive') base::c(0.98, 0.98) else base::c(0.03, 0.98)
  }
  legend_anchor <- if (direction == 'positive') base::c(1, 1) else base::c(0, 1)
  point_layers <- split_point_layers(
    plot_tbl = plot_tbl,
    group_col = 'term_category',
    background_groups = 'Other')

  ggplot2::ggplot(plot_tbl, ggplot2::aes(x = .data$waterfall_rank, y = .data$NES)) +
    ggplot2::geom_point(
      data = point_layers$background,
      ggplot2::aes(color = .data$term_category),
      size = 0.95,
      alpha = 1,
      stroke = 0) +
    ggplot2::geom_point(
      data = point_layers$foreground,
      ggplot2::aes(color = .data$term_category),
      size = 0.95,
      alpha = 1,
      stroke = 0) +
    ggrepel::geom_text_repel(
      data = label_tbl,
      ggplot2::aes(label = .data$label_text, color = .data$term_category),
      family = font_family,
      fontface = label_fontface,
      size = label_size / ggplot2::.pt,
      lineheight = 0.9,
      hjust = 0.5,
      nudge_x = label_tbl$label_nudge_x,
      nudge_y = label_tbl$label_nudge_y,
      min.segment.length = 0,
      segment.color = 'black',
      segment.size = 0.16,
      box.padding = 0.76,
      point.padding = 0.06,
      direction = label_repel_direction,
      force = 2.4,
      force_pull = 0.85,
      max.time = 8,
      max.iter = 18000,
      max.overlaps = Inf,
      seed = 42,
      show.legend = FALSE) +
    ggplot2::scale_color_manual(values = color_values, drop = FALSE) +
    ggplot2::scale_x_continuous(breaks = x_breaks, expand = ggplot2::expansion(mult = base::c(0, 0.02))) +
    ggplot2::scale_y_continuous(breaks = y_breaks, expand = base::c(0, 0)) +
    ggplot2::coord_cartesian(xlim = base::c(0, base::nrow(plot_tbl) + 1), ylim = y_limits, clip = 'off') +
    ggplot2::labs(x = x_label, y = 'Normalized enrichment score (NES)', color = 'GO:BP category') +
    ggplot2::guides(color = ggplot2::guide_legend(override.aes = base::list(size = 1.5))) +
    ggplot2::theme_classic(base_family = font_family, base_size = 9) +
    ggplot2::theme(
      axis.title = ggplot2::element_text(size = 9, color = 'black'),
      axis.text = ggplot2::element_text(size = 8, color = 'black'),
      axis.line = ggplot2::element_line(color = 'black', linewidth = 0.24),
      axis.ticks = ggplot2::element_line(color = 'black', linewidth = 0.20),
      axis.ticks.length = grid::unit(1.4, 'pt'),
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = legend_position,
      legend.justification = legend_anchor,
      legend.background = ggplot2::element_blank(),
      legend.box.background = ggplot2::element_blank(),
      legend.title = ggplot2::element_text(size = 8, color = 'black', face = 'bold', hjust = 0),
      legend.box.just = 'left',
      legend.text = ggplot2::element_text(size = 7.5, color = 'black'),
      legend.key.size = grid::unit(0.25, 'cm'),
      plot.margin = ggplot2::margin(8, 10, 8, 8))
}

plot_gsea_waterfall <- make_gsea_waterfall
