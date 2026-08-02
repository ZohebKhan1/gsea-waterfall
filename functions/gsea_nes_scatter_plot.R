# ----
# author:
# - Zoheb Khan
#
# script path:
# - functions/gsea_nes_scatter_plot.R
#
# functions:
# - functions/gsea_plot_utils.R
# ----

# 1.0 create comparative NES scatterplot -----------------

#' Make a comparative GSEA NES scatterplot
#'
#' Both inputs must contain complete, unique term keys. The plot uses their
#' exact inner intersection, retains descriptions from `gsea_x`, and does not
#' impute unmatched terms. Significance is defined as `padj < padj_cutoff`.
#'
#' @param gsea_x First GSEA result table.
#' @param gsea_y Second GSEA result table.
#' @param x_label X-axis label.
#' @param y_label Y-axis label.
#' @param term_col Column containing term descriptions and the fallback key in
#'   both tables.
#' @param nes_col Column containing normalized enrichment scores in both tables.
#' @param pvalue_col Column containing nominal GSEA p-values in both tables.
#' @param padj_col Column containing adjusted p-values in both tables.
#' @param id_col Optional column containing unique stable term identifiers in
#'   both tables; overrides existing `go_term_id` columns.
#' @param color_by Color points by significance or matched GO term group.
#' @param term_groups Named list used when `color_by = 'term_group'`.
#' @param term_group_colors Optional named colors for `term_groups`.
#' @param point_colors Optional named colors for significance groups when
#'   `color_by = 'significance'`.
#' @param quadrant Terms to include: `all`, `q1`, or `q3`.
#' @param x_name Short name used for the first contrast in the legend.
#' @param y_name Short name used for the second contrast in the legend.
#' @param label_n_x_only Number of terms to label that are significant only in `gsea_x`.
#' @param label_n_y_only Number of terms to label that are significant only in `gsea_y`.
#' @param label_n_both Number of terms to label that are significant in both contrasts.
#' @param label_terms Optional exact term descriptions or GO IDs to label.
#' @param point_size Point size.
#' @param label_size Label font size in points.
#' @param label_color Text color for GO term labels.
#' @param label_words_per_line Number of GO term words shown on each label line.
#' @param padj_cutoff Adjusted p-value cutoff for significance coloring.
#' @param include_nonsignificant Keep terms that are not significant in either contrast.
#' @param equal_axis_limits Use matching x/y limits. For `q1` and `q3` this
#'   zooms to the selected quadrant while keeping the x and y ranges equal.
#'   Cannot be combined with explicit axis limits.
#' @param quadrant_min_abs_nes Inner absolute NES boundary for Q1 and Q3 plots
#'   when `equal_axis_limits = TRUE`.
#' @param x_min Optional lower x-axis limit.
#' @param x_max Optional upper x-axis limit.
#' @param y_min Optional lower y-axis limit.
#' @param y_max Optional upper y-axis limit.
#' @param show_fit_line Draw a linear best-fit line through plotted terms.
#' @param legend_nrow Optional number of legend rows. If `NULL`, significance
#'   legends use two rows and term-group legends use one row.
#' @param legend_position Legend position passed to `ggplot2::theme()`.
#' @param font_family Figure font family.
#'
#' @return A ggplot object built from the exact shared term-key intersection
#'   without file-system side effects
#'
plot_gsea_nes_scatter <- function(gsea_x,
                                  gsea_y,
                                  x_label,
                                  y_label,
                                  term_col = 'go_description',
                                  nes_col = 'NES',
                                  pvalue_col = 'pval',
                                  padj_col = 'padj',
                                  id_col = NULL,
                                  color_by = c('significance', 'term_group'),
                                  term_groups = NULL,
                                  term_group_colors = NULL,
                                  point_colors = NULL,
                                  quadrant = c('all', 'q1', 'q3'),
                                  x_name = 'x contrast',
                                  y_name = 'y contrast',
                                  label_n_x_only = 0L,
                                  label_n_y_only = 0L,
                                  label_n_both = 0L,
                                  label_terms = NULL,
                                  point_size = 1.3,
                                  label_size = 7.5,
                                  label_color = 'black',
                                  label_words_per_line = 2L,
                                  padj_cutoff = 0.05,
                                  include_nonsignificant = FALSE,
                                  equal_axis_limits = FALSE,
                                  quadrant_min_abs_nes = 0,
                                  x_min = NULL,
                                  x_max = NULL,
                                  y_min = NULL,
                                  y_max = NULL,
                                  show_fit_line = TRUE,
                                  legend_nrow = NULL,
                                  legend_position = 'top',
                                  font_family = 'Nimbus Sans') {
  color_by <- match.arg(color_by)
  quadrant <- match.arg(quadrant)
  include_nonsignificant <- .gsea_validate_flag(
    include_nonsignificant,
    'include_nonsignificant')
  equal_axis_limits <- .gsea_validate_flag(equal_axis_limits, 'equal_axis_limits')
  show_fit_line <- .gsea_validate_flag(show_fit_line, 'show_fit_line')
  padj_cutoff <- .gsea_validate_number(
    padj_cutoff,
    'padj_cutoff',
    minimum = 0,
    maximum = 1)
  point_size <- .gsea_validate_number(point_size, 'point_size', minimum = 0)
  if (length(label_terms) == 0L) {
    label_n_x_only <- .gsea_validate_count(label_n_x_only, 'label_n_x_only')
    label_n_y_only <- .gsea_validate_count(label_n_y_only, 'label_n_y_only')
    label_n_both <- .gsea_validate_count(label_n_both, 'label_n_both')
  }
  if (equal_axis_limits && quadrant != 'all') {
    quadrant_min_abs_nes <- .gsea_validate_number(
      quadrant_min_abs_nes,
      'quadrant_min_abs_nes',
      minimum = 0)
  }
  .gsea_validate_limits(x_min, x_max, 'x_min', 'x_max')
  .gsea_validate_limits(y_min, y_max, 'y_min', 'y_max')
  if (equal_axis_limits && any(!vapply(
    list(x_min, x_max, y_min, y_max),
    is.null,
    logical(1)))) {
    stop('`equal_axis_limits` cannot be combined with explicit axis limits.', call. = FALSE)
  }
  if (is.null(legend_nrow)) {
    legend_nrow <- if (color_by == 'significance') 2L else 1L
  } else {
    legend_nrow <- .gsea_validate_count(legend_nrow, 'legend_nrow', minimum = 1L)
  }
  x_tbl <- .gsea_standardize_results(
    gsea_results = gsea_x,
    term_col = term_col,
    nes_col = nes_col,
    pvalue_col = pvalue_col,
    padj_col = padj_col,
    id_col = id_col)
  y_tbl <- .gsea_standardize_results(
    gsea_results = gsea_y,
    term_col = term_col,
    nes_col = nes_col,
    pvalue_col = pvalue_col,
    padj_col = padj_col,
    id_col = id_col)
  plot_tbl <- merge(x_tbl, y_tbl, by = 'go_term_id', suffixes = c('_x', '_y'))
  if (nrow(plot_tbl) == 0L) {
    stop('The GSEA tables have no overlapping term keys.', call. = FALSE)
  }
  plot_tbl$go_description <- plot_tbl$go_description_x
  plot_tbl$go_description_x <- NULL
  plot_tbl$go_description_y <- NULL

  plot_tbl$significant_x <- plot_tbl$padj_x < padj_cutoff
  plot_tbl$significant_y <- plot_tbl$padj_y < padj_cutoff
  if (!include_nonsignificant) {
    plot_tbl <- plot_tbl[plot_tbl$significant_x | plot_tbl$significant_y, , drop = FALSE]
  }
  plot_tbl$label_score <- pmax(abs(plot_tbl$NES_x), abs(plot_tbl$NES_y)) *
    pmax(.gsea_neg_log10(plot_tbl$padj_x), .gsea_neg_log10(plot_tbl$padj_y))
  plot_tbl$scatter_quadrant <- ifelse(
    plot_tbl$NES_x >= 0 & plot_tbl$NES_y >= 0,
    'Q1',
    ifelse(plot_tbl$NES_x < 0 & plot_tbl$NES_y >= 0, 'Q2', ifelse(plot_tbl$NES_x < 0 & plot_tbl$NES_y < 0, 'Q3', 'Q4')))

  if (color_by == 'significance') {
    plot_tbl$point_group <- ifelse(
      plot_tbl$significant_x & plot_tbl$significant_y,
      'Significant in both',
      ifelse(
        plot_tbl$significant_x & !plot_tbl$significant_y,
        paste0('Significant in ', x_name, ' only'),
        ifelse(!plot_tbl$significant_x & plot_tbl$significant_y, paste0('Significant in ', y_name, ' only'), 'Not significant in either')))
    group_names <- c(
      'Significant in both',
      paste0('Significant in ', x_name, ' only'),
      paste0('Significant in ', y_name, ' only'),
      'Not significant in either')
    if (is.null(point_colors)) {
      brewer_colors <- RColorBrewer::brewer.pal(4L, 'Set1')
      color_values <- stats::setNames(
        c(brewer_colors[[4L]], brewer_colors[[1L]], brewer_colors[[2L]], .gsea_nonsignificant_color),
        group_names)
    } else {
      color_values <- .gsea_resolve_named_colors(
        color_values = point_colors,
        required_names = group_names,
        parameter_name = 'point_colors')
    }
  } else {
    if (is.null(term_groups)) {
      stop("Provide `term_groups` when `color_by = 'term_group'.", call. = FALSE)
    }
    plot_tbl$point_group <- .gsea_assign_term_groups(
      plot_tbl,
      term_groups = term_groups)
    color_values <- .gsea_resolve_term_group_colors(
      term_groups = term_groups,
      term_group_colors = term_group_colors)
  }
  if (quadrant == 'q1') {
    plot_tbl <- plot_tbl[plot_tbl$NES_x > 0 & plot_tbl$NES_y > 0, , drop = FALSE]
  } else if (quadrant == 'q3') {
    plot_tbl <- plot_tbl[plot_tbl$NES_x < 0 & plot_tbl$NES_y < 0, , drop = FALSE]
  }
  if (nrow(plot_tbl) == 0L) {
    stop('No matched GSEA terms remain after the requested filters.', call. = FALSE)
  }
  group_counts <- table(plot_tbl$point_group)
  color_values <- color_values[names(color_values) %in% names(group_counts)]
  count_values <- as.integer(group_counts[names(color_values)])
  legend_labels <- stats::setNames(
    paste0(names(color_values), ' (n=', count_values, ')'),
    names(color_values))

  select_scatter_label_rows <- function(label_tbl, label_n) {
    if (label_n <= 0L || nrow(label_tbl) == 0L) {
      return(label_tbl[0L, , drop = FALSE])
    }
    label_tbl <- label_tbl[order(-label_tbl$label_score, label_tbl$padj_x, label_tbl$padj_y, label_tbl$go_description), , drop = FALSE]
    if (quadrant != 'all') {
      return(utils::head(label_tbl, label_n))
    }
    label_tbl$label_group_rank <- stats::ave(
      -label_tbl$label_score,
      label_tbl$scatter_quadrant,
      FUN = rank,
      ties.method = 'first')
    label_tbl <- label_tbl[order(label_tbl$label_group_rank, label_tbl$scatter_quadrant), , drop = FALSE]
    label_tbl$label_group_rank <- NULL
    utils::head(label_tbl, label_n)
  }

  if (length(label_terms) > 0L) {
    label_tbl <- .gsea_select_labels(plot_tbl, label_terms = as.list(label_terms))
  } else {
    label_x_tbl <- select_scatter_label_rows(
      plot_tbl[plot_tbl$significant_x & !plot_tbl$significant_y, , drop = FALSE],
      label_n_x_only)
    label_y_tbl <- select_scatter_label_rows(
      plot_tbl[!plot_tbl$significant_x & plot_tbl$significant_y, , drop = FALSE],
      label_n_y_only)

    label_both_tbl <- plot_tbl[
      plot_tbl$significant_x &
        plot_tbl$significant_y,
      ,
      drop = FALSE]
    label_both_tbl <- select_scatter_label_rows(label_both_tbl, label_n_both)

    label_tbl <- rbind(label_x_tbl, label_y_tbl, label_both_tbl)
  }
  label_tbl$label_text <- .gsea_wrap_label(label_tbl$go_description, words_per_line = label_words_per_line)
  label_tbl$label_nudge_x <- ifelse(label_tbl$NES_x >= 0, 0.08, -0.08)
  label_tbl$label_nudge_y <- ifelse(label_tbl$NES_y >= 0, 0.08, -0.08)

  axis_x_limits <- NULL
  axis_y_limits <- NULL
  if (equal_axis_limits) {
    axis_limit <- max(abs(c(plot_tbl$NES_x, plot_tbl$NES_y)), na.rm = TRUE)
    if (!is.finite(axis_limit) || axis_limit <= 0) {
      axis_limit <- 1
    }
    axis_limit <- axis_limit * 1.08
    quadrant_min_abs_nes <- max(quadrant_min_abs_nes, 0, na.rm = TRUE)
    if (axis_limit <= quadrant_min_abs_nes) {
      axis_limit <- quadrant_min_abs_nes + 0.5
    }
    if (quadrant == 'q1') {
      axis_x_limits <- c(quadrant_min_abs_nes, axis_limit)
      axis_y_limits <- c(quadrant_min_abs_nes, axis_limit)
    } else if (quadrant == 'q3') {
      axis_x_limits <- c(-axis_limit, -quadrant_min_abs_nes)
      axis_y_limits <- c(-axis_limit, -quadrant_min_abs_nes)
    } else {
      axis_x_limits <- c(-axis_limit, axis_limit)
      axis_y_limits <- c(-axis_limit, axis_limit)
    }
  }
  x_breaks <- NULL
  y_breaks <- NULL
  if (is.null(axis_x_limits) && (!is.null(x_min) || !is.null(x_max))) {
    axis_x_limits <- .gsea_axis_limits(plot_tbl$NES_x, buffer_fraction = 0.06)
  }
  if (is.null(axis_y_limits) && (!is.null(y_min) || !is.null(y_max))) {
    axis_y_limits <- .gsea_axis_limits(plot_tbl$NES_y, buffer_fraction = 0.06)
  }
  if (!is.null(x_min)) {
    axis_x_limits[[1L]] <- x_min
  }
  if (!is.null(x_max)) {
    axis_x_limits[[2L]] <- x_max
  }
  if (!is.null(y_min)) {
    axis_y_limits[[1L]] <- y_min
  }
  if (!is.null(y_max)) {
    axis_y_limits[[2L]] <- y_max
  }
  if (!is.null(axis_x_limits)) {
    x_breaks <- scales::breaks_pretty(n = 5)(axis_x_limits)
    if (axis_x_limits[[1L]] <= 0 && axis_x_limits[[2L]] >= 0) {
      x_breaks <- sort(unique(c(x_breaks, 0)))
    }
  }
  if (!is.null(axis_y_limits)) {
    y_breaks <- scales::breaks_pretty(n = 5)(axis_y_limits)
    if (axis_y_limits[[1L]] <= 0 && axis_y_limits[[2L]] >= 0) {
      y_breaks <- sort(unique(c(y_breaks, 0)))
    }
  }

  fit_line_layer <- if (show_fit_line && nrow(plot_tbl) > 1L) {
    ggplot2::geom_smooth(
      method = 'lm',
      formula = y ~ x,
      se = FALSE,
      color = 'black',
      linewidth = 0.32,
      linetype = 'solid')
  } else {
    NULL
  }

  coord_layer <- if (equal_axis_limits) {
    ggplot2::coord_fixed(ratio = 1, xlim = axis_x_limits, ylim = axis_y_limits, clip = 'off')
  } else {
    ggplot2::coord_cartesian(xlim = axis_x_limits, ylim = axis_y_limits, clip = 'off')
  }
  point_layers <- .gsea_split_point_layers(
    plot_tbl = plot_tbl,
    group_col = 'point_group',
    background_groups = c('Not significant in either', 'Not significant', 'Other'))
  label_layer <- ggrepel::geom_text_repel(
    data = label_tbl,
    ggplot2::aes(label = .data$label_text),
    nudge_x = label_tbl$label_nudge_x,
    nudge_y = label_tbl$label_nudge_y,
    family = font_family,
    color = label_color,
    size = label_size / ggplot2::.pt,
    lineheight = 0.95,
    hjust = 0.5,
    min.segment.length = 0,
    segment.color = 'black',
    segment.size = 0.16,
    box.padding = 0.62,
    point.padding = 0.12,
    force = 1.9,
    force_pull = 3.0,
    max.time = 8,
    max.iter = 18000,
    max.overlaps = Inf,
    seed = 42,
    show.legend = FALSE)

  ggplot2::ggplot(plot_tbl, ggplot2::aes(x = .data$NES_x, y = .data$NES_y)) +
    fit_line_layer +
    ggplot2::geom_hline(yintercept = 0, color = 'black', linewidth = 0.22, linetype = 'dashed') +
    ggplot2::geom_vline(xintercept = 0, color = 'black', linewidth = 0.22, linetype = 'dashed') +
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
    label_layer +
    ggplot2::scale_color_manual(values = color_values, labels = legend_labels, drop = FALSE) +
    ggplot2::scale_x_continuous(
      breaks = x_breaks,
      position = if (quadrant == 'q3') 'top' else 'bottom',
      expand = if (is.null(axis_x_limits)) ggplot2::expansion(mult = 0.06) else ggplot2::expansion(mult = 0)) +
    ggplot2::scale_y_continuous(
      breaks = y_breaks,
      position = if (quadrant == 'q3') 'right' else 'left',
      expand = if (is.null(axis_y_limits)) ggplot2::expansion(mult = 0.06) else ggplot2::expansion(mult = 0)) +
    coord_layer +
    ggplot2::labs(x = x_label, y = y_label, color = if (color_by == 'significance') 'GSEA enrichment' else 'GO:BP category') +
    ggplot2::guides(color = ggplot2::guide_legend(
      nrow = legend_nrow,
      byrow = TRUE,
      title.position = 'top',
      title.hjust = 0.5,
      override.aes = list(size = point_size + 0.5))) +
    .gsea_theme(font_family) +
    ggplot2::theme(
      legend.position = legend_position,
      legend.direction = 'horizontal',
      legend.box = 'horizontal',
      legend.justification = 'center',
      legend.title = ggplot2::element_text(hjust = 0.5),
      legend.margin = ggplot2::margin(0, 0, 8, 0),
      plot.margin = ggplot2::margin(10, 10, 8, 10))
}
