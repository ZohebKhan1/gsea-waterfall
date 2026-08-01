# ----
# author:
# - Zoheb Khan
#
# script path:
# - functions/gsea_waterfall_plot.R
#
# functions:
# - functions/gsea_plot_utils.R
# ----

# 1.0 create waterfall plot -----------------

#' Make a GSEA waterfall plot
#'
#' Positive and negative plots exclude zero NES values. Terms are ranked before
#' `top_n` is applied. Term descriptions act as unique keys unless `id_col`
#' supplies a separate complete, unique identifier.
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
#' @return A ggplot object without file-system side effects
#'
#' @export
plot_gsea_waterfall <- function(gsea_results,
                                term_col = 'go_description',
                                nes_col = 'NES',
                                pvalue_col = 'pval',
                                padj_col = 'padj',
                                id_col = NULL,
                                direction = c('positive', 'negative'),
                                top_n = 100L,
                                rank_by = c('NES', 'padj', 'pvalue', 'pval'),
                                x_label = NULL,
                                label_terms = NULL,
                                label_ranks = NULL,
                                label_n = 12L,
                                term_groups = NULL,
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
  direction <- match.arg(direction)
  rank_by <- .gsea_normalize_rank_by(rank_by)
  top_n <- .gsea_validate_count(top_n, 'top_n', minimum = 1L)
  if (length(label_terms) == 0L && length(label_ranks) == 0L) {
    label_n <- .gsea_validate_count(label_n, 'label_n')
  }
  .gsea_validate_limits(y_min, y_max, 'y_min', 'y_max')
  .gsea_validate_positive_integer_vector(label_ranks, 'label_ranks')
  gsea_results <- .gsea_standardize_results(
    gsea_results = gsea_results,
    term_col = term_col,
    nes_col = nes_col,
    pvalue_col = pvalue_col,
    padj_col = padj_col,
    id_col = id_col)

  if (direction == 'positive') {
    plot_tbl <- gsea_results[gsea_results$NES > 0, , drop = FALSE]
  } else {
    plot_tbl <- gsea_results[gsea_results$NES < 0, , drop = FALSE]
  }

  if (nrow(plot_tbl) == 0L) {
    stop('No GSEA terms remain for the requested direction.', call. = FALSE)
  }

  plot_tbl <- .gsea_rank_direction_terms(plot_tbl, direction = direction, rank_by = rank_by)
  plot_tbl <- utils::head(plot_tbl, top_n)
  plot_tbl$waterfall_rank <- seq_len(nrow(plot_tbl))
  plot_tbl$label_score <- abs(plot_tbl$NES) * .gsea_neg_log10(plot_tbl$padj)
  plot_tbl$term_group <- .gsea_assign_term_groups(
    gsea_results = plot_tbl,
    term_groups = term_groups)
  color_values <- .gsea_resolve_term_group_colors(
    term_groups = term_groups,
    term_group_colors = term_group_colors)
  label_requests <- .gsea_combine_label_requests(
    label_terms = label_terms,
    label_ranks = label_ranks)
  label_tbl <- .gsea_select_labels(
    plot_tbl = plot_tbl,
    label_terms = label_requests,
    label_n = label_n,
    score_col = 'label_score',
    rank_col = 'waterfall_rank')
  label_tbl$label_text <- .gsea_wrap_label(label_tbl$go_description, words_per_line = label_words_per_line)
  y_limits <- .gsea_resolve_axis_limits(plot_tbl$NES, lower = y_min, upper = y_max)
  y_breaks <- scales::breaks_pretty(n = 4)(y_limits)
  label_tbl$label_nudge_x <- numeric(nrow(label_tbl))
  label_tbl$label_nudge_y <- numeric(nrow(label_tbl))
  if (nrow(label_tbl) > 0L) {
    label_tbl <- label_tbl[order(label_tbl$waterfall_rank), , drop = FALSE]
    nudge_step <- max(diff(range(plot_tbl$NES, na.rm = TRUE)) * label_y_nudge_fraction, 0.08)
    if (is.null(label_x_nudge)) {
      label_tbl$label_nudge_x <- rep(c(-0.65, 0.35, 0, 0.6, -0.3), length.out = nrow(label_tbl))
    } else {
      label_tbl$label_nudge_x <- rep(label_x_nudge, length.out = nrow(label_tbl))
    }
    label_tbl$label_nudge_y <- rep(c(nudge_step, -nudge_step), length.out = nrow(label_tbl))
    if (direction == 'negative') {
      label_tbl$label_nudge_y <- -label_tbl$label_nudge_y
    }
  }
  x_breaks <- unique(c(0, scales::breaks_pretty(n = 4)(c(1, nrow(plot_tbl)))))
  if (is.null(x_label)) {
    x_label <- paste0('Ranked GO terms (by ', rank_by, ')')
  }
  if (is.null(legend_position)) {
    legend_position <- if (direction == 'positive') c(0.98, 0.98) else c(0.03, 0.98)
  }
  legend_anchor <- if (direction == 'positive') c(1, 1) else c(0, 1)
  point_layers <- .gsea_split_point_layers(
    plot_tbl = plot_tbl,
    group_col = 'term_group',
    background_groups = 'Other')

  ggplot2::ggplot(plot_tbl, ggplot2::aes(x = .data$waterfall_rank, y = .data$NES)) +
    ggplot2::geom_point(
      data = point_layers$background,
      ggplot2::aes(color = .data$term_group),
      size = 0.95,
      alpha = 1,
      stroke = 0) +
    ggplot2::geom_point(
      data = point_layers$foreground,
      ggplot2::aes(color = .data$term_group),
      size = 0.95,
      alpha = 1,
      stroke = 0) +
    ggrepel::geom_text_repel(
      data = label_tbl,
      ggplot2::aes(label = .data$label_text, color = .data$term_group),
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
    ggplot2::scale_x_continuous(breaks = x_breaks, expand = ggplot2::expansion(mult = c(0, 0.02))) +
    ggplot2::scale_y_continuous(breaks = y_breaks, expand = c(0, 0)) +
    ggplot2::coord_cartesian(xlim = c(0, nrow(plot_tbl) + 1), ylim = y_limits, clip = 'off') +
    ggplot2::labs(x = x_label, y = 'Normalized enrichment score (NES)', color = 'GO:BP category') +
    ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(size = 1.5))) +
    .gsea_theme(font_family) +
    ggplot2::theme(
      axis.ticks.length = grid::unit(1.4, 'pt'),
      legend.position = legend_position,
      legend.justification = legend_anchor,
      legend.background = ggplot2::element_blank(),
      legend.box.background = ggplot2::element_blank(),
      legend.title = ggplot2::element_text(hjust = 0),
      legend.box.just = 'left',
      plot.margin = ggplot2::margin(8, 10, 8, 8))
}
