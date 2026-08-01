if (!base::exists('standardize_gsea_results', mode = 'function')) {
  base::stop('Source gsea_plot_utils.R before this plotting function file.', call. = FALSE)
}

#' Make a comparative GSEA NES scatterplot
#'
#' @param gsea_x First GSEA result table.
#' @param gsea_y Second GSEA result table.
#' @param x_label X-axis label.
#' @param y_label Y-axis label.
#' @param term_col Column containing GO term names in both tables.
#' @param nes_col Column containing normalized enrichment scores in both tables.
#' @param pvalue_col Column containing nominal GSEA p-values in both tables.
#' @param padj_col Column containing adjusted p-values in both tables.
#' @param id_col Optional column containing GO IDs in both tables.
#' @param color_by Color points by significance or matched GO term group.
#' @param term_groups Named list used when `color_by = 'term_group'`.
#' @param category_n Optional expected number of categories in `term_groups`;
#'   automatic term-group colors are generated for this category count.
#' @param max_categories Maximum allowed number of custom categories.
#' @param term_group_colors Optional named colors for `term_groups`.
#' @param point_colors Optional named colors for significance groups when
#'   `color_by = 'significance'`.
#' @param quadrant Terms to include: all, Q1, or Q3.
#' @param label_n_one_contrast Number of terms to label that are significant in only one contrast.
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
#' @param equal_axis_limits Use matching x/y limits. For Q1 and Q3 this zooms
#'   to the selected quadrant while keeping the x and y ranges equal.
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
#' @return A ggplot object.
make_gsea_scatterplot <- function(gsea_x,
                                  gsea_y,
                                  x_label,
                                  y_label,
                                  term_col = 'go_description',
                                  nes_col = 'NES',
                                  pvalue_col = 'pval',
                                  padj_col = 'padj',
                                  id_col = NULL,
                                  color_by = base::c('significance', 'term_group'),
                                  term_groups = NULL,
                                  category_n = NULL,
                                  max_categories = 6L,
                                  term_group_colors = NULL,
                                  point_colors = NULL,
                                  quadrant = base::c('all', 'q1', 'q3'),
                                  x_name = 'x contrast',
                                  y_name = 'y contrast',
                                  label_n_one_contrast = 0L,
                                  label_n_x_only = NULL,
                                  label_n_y_only = NULL,
                                  label_n_both = 0L,
                                  label_terms = NULL,
                                  point_size = 0.95,
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
  color_by <- base::match.arg(color_by)
  quadrant <- base::match.arg(quadrant)
  validate_term_groups(term_groups = term_groups, category_n = category_n)
  validate_max_categories(term_groups = term_groups, max_categories = max_categories)
  if (base::is.null(legend_nrow)) {
    legend_nrow <- if (color_by == 'significance') 2L else 1L
  }
  if (base::is.null(label_n_x_only)) {
    label_n_x_only <- label_n_one_contrast
  }
  if (base::is.null(label_n_y_only)) {
    label_n_y_only <- label_n_one_contrast
  }
  x_tbl <- standardize_gsea_results(
    gsea_results = gsea_x,
    term_col = term_col,
    nes_col = nes_col,
    pvalue_col = pvalue_col,
    padj_col = padj_col,
    id_col = id_col)
  y_tbl <- standardize_gsea_results(
    gsea_results = gsea_y,
    term_col = term_col,
    nes_col = nes_col,
    pvalue_col = pvalue_col,
    padj_col = padj_col,
    id_col = id_col)
  plot_tbl <- base::merge(
    x_tbl,
    y_tbl,
    by = base::c('go_term_id', 'go_term_class', 'go_description'),
    suffixes = base::c('_x', '_y'))
  plot_tbl <- plot_tbl[stats::complete.cases(plot_tbl[, base::c('NES_x', 'NES_y')]), , drop = FALSE]

  sig_x <- !base::is.na(plot_tbl$padj_x) & plot_tbl$padj_x < padj_cutoff
  sig_y <- !base::is.na(plot_tbl$padj_y) & plot_tbl$padj_y < padj_cutoff
  plot_tbl$significant_x <- sig_x
  plot_tbl$significant_y <- sig_y
  if (!base::isTRUE(include_nonsignificant)) {
    plot_tbl <- plot_tbl[plot_tbl$significant_x | plot_tbl$significant_y, , drop = FALSE]
  }
  plot_tbl$label_score <- base::pmax(base::abs(plot_tbl$NES_x), base::abs(plot_tbl$NES_y)) *
    base::pmax(safe_neg_log10(plot_tbl$padj_x), safe_neg_log10(plot_tbl$padj_y))
  plot_tbl$scatter_quadrant <- base::ifelse(
    plot_tbl$NES_x >= 0 & plot_tbl$NES_y >= 0,
    'Q1',
    base::ifelse(plot_tbl$NES_x < 0 & plot_tbl$NES_y >= 0, 'Q2', base::ifelse(plot_tbl$NES_x < 0 & plot_tbl$NES_y < 0, 'Q3', 'Q4')))

  if (color_by == 'significance') {
    plot_tbl$point_group <- base::ifelse(
      plot_tbl$significant_x & plot_tbl$significant_y,
      'Significant in both',
      base::ifelse(
        plot_tbl$significant_x & !plot_tbl$significant_y,
        base::paste0('Significant in ', x_name, ' only'),
        base::ifelse(!plot_tbl$significant_x & plot_tbl$significant_y, base::paste0('Significant in ', y_name, ' only'), 'Not significant in either')))
    if (!base::requireNamespace('RColorBrewer', quietly = TRUE)) {
      base::stop('The RColorBrewer package is required for default scatterplot significance colors.', call. = FALSE)
    }
    brewer_colors <- RColorBrewer::brewer.pal(4L, 'Set1')
    color_values <- base::c(
      'Significant in both' = brewer_colors[[4L]],
      stats::setNames(brewer_colors[[1L]], base::paste0('Significant in ', x_name, ' only')),
      stats::setNames(brewer_colors[[2L]], base::paste0('Significant in ', y_name, ' only')),
      'Not significant in either' = gsea_light_gray())
    custom_colors <- resolve_named_colors(
      color_values = point_colors,
      required_names = base::names(color_values),
      parameter_name = 'point_colors')
    if (!base::is.null(custom_colors)) {
      color_values <- custom_colors
    }
  } else {
    if (base::is.null(term_groups)) {
      base::stop("Provide `term_groups` when `color_by = 'term_group'.", call. = FALSE)
    }
    term_group_mapping <- assign_plot_term_groups(
      x_tbl,
      term_groups = term_groups)
    plot_tbl$point_group <- term_group_mapping[base::match(plot_tbl$go_term_id, x_tbl$go_term_id)]
    color_values <- resolve_term_group_colors(
      term_groups = term_groups,
      term_group_colors = term_group_colors)
  }
  if (quadrant == 'q1') {
    plot_tbl <- plot_tbl[plot_tbl$NES_x > 0 & plot_tbl$NES_y > 0, , drop = FALSE]
  } else if (quadrant == 'q3') {
    plot_tbl <- plot_tbl[plot_tbl$NES_x < 0 & plot_tbl$NES_y < 0, , drop = FALSE]
  }
  group_counts <- base::table(plot_tbl$point_group)
  color_values <- color_values[base::names(color_values) %in% base::names(group_counts)]
  count_values <- base::as.integer(group_counts[base::names(color_values)])
  count_values[base::is.na(count_values)] <- 0L
  legend_labels <- stats::setNames(
    base::paste0(base::names(color_values), ' (n=', count_values, ')'),
    base::names(color_values))

  select_scatter_label_rows <- function(label_tbl, label_n) {
    if (label_n <= 0L || base::nrow(label_tbl) == 0L) {
      return(label_tbl[0L, , drop = FALSE])
    }
    label_tbl <- label_tbl[base::order(-label_tbl$label_score, label_tbl$padj_x, label_tbl$padj_y, label_tbl$go_description), , drop = FALSE]
    if (quadrant != 'all') {
      return(utils::head(label_tbl, label_n))
    }
    label_tbl$label_group_rank <- ave(
      -label_tbl$label_score,
      label_tbl$scatter_quadrant,
      FUN = rank,
      ties.method = 'first')
    label_tbl <- label_tbl[base::order(label_tbl$label_group_rank, label_tbl$scatter_quadrant), , drop = FALSE]
    label_tbl$label_group_rank <- NULL
    utils::head(label_tbl, label_n)
  }

  label_requests <- combine_label_requests(label_terms = label_terms)
  if (!base::is.null(label_requests) && base::length(label_requests) > 0L) {
    label_tbl <- select_gsea_labels(plot_tbl, label_terms = label_requests)
  } else {
    label_x_tbl <- select_scatter_label_rows(
      plot_tbl[label_n_x_only > 0L & plot_tbl$significant_x & !plot_tbl$significant_y, , drop = FALSE],
      label_n_x_only)
    label_y_tbl <- select_scatter_label_rows(
      plot_tbl[label_n_y_only > 0L & !plot_tbl$significant_x & plot_tbl$significant_y, , drop = FALSE],
      label_n_y_only)

    label_both_tbl <- plot_tbl[
      label_n_both > 0L &
        plot_tbl$significant_x &
        plot_tbl$significant_y,
      ,
      drop = FALSE]
    label_both_tbl <- select_scatter_label_rows(label_both_tbl, label_n_both)

    label_tbl <- base::rbind(label_x_tbl, label_y_tbl, label_both_tbl)
  }
  label_tbl <- label_tbl[!base::duplicated(label_tbl$go_term_id), , drop = FALSE]
  label_tbl$label_text <- wrap_gsea_label(label_tbl$go_description, words_per_line = label_words_per_line)
  label_tbl$label_nudge_x <- base::ifelse(label_tbl$NES_x >= 0, 0.08, -0.08)
  label_tbl$label_nudge_y <- base::ifelse(label_tbl$NES_y >= 0, 0.08, -0.08)

  axis_x_limits <- NULL
  axis_y_limits <- NULL
  if (equal_axis_limits) {
    axis_limit <- base::max(base::abs(base::c(plot_tbl$NES_x, plot_tbl$NES_y)), na.rm = TRUE)
    if (!base::is.finite(axis_limit) || axis_limit <= 0) {
      axis_limit <- 1
    }
    axis_limit <- axis_limit * 1.08
    quadrant_min_abs_nes <- base::max(quadrant_min_abs_nes, 0, na.rm = TRUE)
    if (axis_limit <= quadrant_min_abs_nes) {
      axis_limit <- quadrant_min_abs_nes + 0.5
    }
    if (quadrant == 'q1') {
      axis_x_limits <- base::c(quadrant_min_abs_nes, axis_limit)
      axis_y_limits <- base::c(quadrant_min_abs_nes, axis_limit)
    } else if (quadrant == 'q3') {
      axis_x_limits <- base::c(-axis_limit, -quadrant_min_abs_nes)
      axis_y_limits <- base::c(-axis_limit, -quadrant_min_abs_nes)
    } else {
      axis_x_limits <- base::c(-axis_limit, axis_limit)
      axis_y_limits <- base::c(-axis_limit, axis_limit)
    }
  }
  axis_breaks <- NULL
  if (!base::is.null(axis_x_limits)) {
    axis_breaks <- scales::breaks_pretty(n = 5)(axis_x_limits)
    if (axis_x_limits[[1L]] <= 0 && axis_x_limits[[2L]] >= 0) {
      axis_breaks <- base::sort(base::unique(base::c(axis_breaks, 0)))
    }
  }
  if (base::is.null(axis_x_limits) && (!base::is.null(x_min) || !base::is.null(x_max))) {
    axis_x_limits <- calculate_axis_limits(plot_tbl$NES_x, buffer_fraction = 0.06)
  }
  if (base::is.null(axis_y_limits) && (!base::is.null(y_min) || !base::is.null(y_max))) {
    axis_y_limits <- calculate_axis_limits(plot_tbl$NES_y, buffer_fraction = 0.06)
  }
  if (!base::is.null(x_min)) {
    axis_x_limits[[1L]] <- x_min
  }
  if (!base::is.null(x_max)) {
    axis_x_limits[[2L]] <- x_max
  }
  if (!base::is.null(y_min)) {
    axis_y_limits[[1L]] <- y_min
  }
  if (!base::is.null(y_max)) {
    axis_y_limits[[2L]] <- y_max
  }
  if (base::is.null(axis_breaks) && !base::is.null(axis_x_limits)) {
    axis_breaks <- scales::breaks_pretty(n = 5)(axis_x_limits)
  }

  fit_line_layer <- if (show_fit_line && base::nrow(plot_tbl) > 1L) {
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
  background_groups <- base::c('Not significant in either', 'Not significant', 'Other')
  background_rows <- base::is.na(plot_tbl$point_group) | plot_tbl$point_group %in% background_groups
  background_tbl <- plot_tbl[background_rows, , drop = FALSE]
  foreground_tbl <- plot_tbl[!background_rows, , drop = FALSE]
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
      data = background_tbl,
      ggplot2::aes(color = .data$point_group),
      size = point_size,
      alpha = 1,
      stroke = 0) +
    ggplot2::geom_point(
      data = foreground_tbl,
      ggplot2::aes(color = .data$point_group),
      size = point_size,
      alpha = 1,
      stroke = 0) +
    label_layer +
    ggplot2::scale_color_manual(values = color_values, labels = legend_labels, drop = FALSE) +
    ggplot2::scale_x_continuous(
      breaks = axis_breaks,
      position = if (quadrant == 'q3') 'top' else 'bottom',
      expand = if (base::is.null(axis_x_limits)) ggplot2::expansion(mult = 0.06) else ggplot2::expansion(mult = 0)) +
    ggplot2::scale_y_continuous(
      breaks = axis_breaks,
      position = if (quadrant == 'q3') 'right' else 'left',
      expand = if (base::is.null(axis_y_limits)) ggplot2::expansion(mult = 0.06) else ggplot2::expansion(mult = 0)) +
    coord_layer +
    ggplot2::labs(x = x_label, y = y_label, color = if (color_by == 'significance') 'GSEA enrichment' else 'GO term category') +
    ggplot2::guides(color = ggplot2::guide_legend(
      nrow = legend_nrow,
      byrow = TRUE,
      title.position = 'top',
      title.hjust = 0.5,
      override.aes = base::list(size = point_size + 0.5))) +
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
      legend.title = ggplot2::element_text(size = 8, color = 'black', face = 'bold', hjust = 0.5),
      legend.text = ggplot2::element_text(size = 7.5, color = 'black'),
      legend.margin = ggplot2::margin(0, 0, 8, 0),
      legend.key.size = grid::unit(0.25, 'cm'),
      plot.margin = ggplot2::margin(10, 10, 8, 10))
}

plot_gsea_nes_scatter <- make_gsea_scatterplot
