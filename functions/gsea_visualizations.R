# ----
# author:
# - Zoheb Khan
#
# script path:
# - functions/gsea_visualizations.R
#
# public functions:
# - read_gsea_result_csv()
# - plot_gsea_waterfall()
# - plot_gsea_volcano()
# - plot_gsea_half_volcano()
# - plot_gsea_nes_scatter()
#
# This consolidated file contains the shared validation utilities and all
# public plotting functions. Download and source this file once.
# ----


# ---- shared validation utilities ----

# 1.0 normalize GSEA inputs -----------------

# independently supplied GSEA tables require an explicit unique term key
.gsea_standardize_results <- function(gsea_results,
  term_col = 'go_description',
  nes_col = 'NES',
  pvalue_col = NULL,
  padj_col = 'padj',
  id_col = NULL) {
  if (is.null(id_col) && 'go_term_id' %in% names(gsea_results)) {
    id_col <- 'go_term_id'
  }
  if (is.null(pvalue_col) && 'pval' %in% names(gsea_results)) {
    pvalue_col <- 'pval'
  }
  missing_cols <- setdiff(
    c(term_col, nes_col, pvalue_col, padj_col, id_col),
    names(gsea_results))
  if (length(missing_cols) > 0L) {
    stop(
      'Missing required GSEA plotting column(s): ',
      paste(missing_cols, collapse = ', '),
      call. = FALSE)
  }

  go_description <- as.character(gsea_results[[term_col]])
  go_term_id <- if (is.null(id_col)) {
    go_description
  } else {
    as.character(gsea_results[[id_col]])
  }

  if (anyNA(go_description) || any(!nzchar(trimws(go_description)))) {
    stop('GSEA term descriptions must be complete and non-empty.', call. = FALSE)
  }
  if (anyNA(go_term_id) || any(!nzchar(trimws(go_term_id)))) {
    stop('GSEA term keys must be complete and non-empty.', call. = FALSE)
  }
  duplicate_ids <- unique(go_term_id[duplicated(go_term_id)])
  if (length(duplicate_ids) > 0L) {
    stop(
      'GSEA term keys must be unique; duplicated key(s): ',
      paste(utils::head(duplicate_ids, 5L), collapse = ', '),
      call. = FALSE)
  }

  numeric_columns <- c(nes_col, padj_col, pvalue_col)
  if (!all(vapply(
    gsea_results[numeric_columns],
    is.numeric,
    logical(1)))) {
    stop('GSEA NES, adjusted p-value, and optional nominal p-value columns must be numeric.', call. = FALSE)
  }

  pvalue <- if (is.null(pvalue_col)) {
    rep(NA_real_, nrow(gsea_results))
  } else {
    as.numeric(gsea_results[[pvalue_col]])
  }
  plot_data <- data.frame(
    go_term_id = go_term_id,
    go_description = go_description,
    NES = as.numeric(gsea_results[[nes_col]]),
    pvalue = pvalue,
    padj = as.numeric(gsea_results[[padj_col]]),
    stringsAsFactors = FALSE)

  if (anyNA(plot_data$NES) || anyNA(plot_data$padj) ||
      any(!is.finite(plot_data$NES)) || any(!is.finite(plot_data$padj))) {
    stop('GSEA term, NES, and adjusted p-value columns must be complete and finite.', call. = FALSE)
  }
  if (!all(is.na(plot_data$pvalue)) &&
      (anyNA(plot_data$pvalue) || any(!is.finite(plot_data$pvalue)))) {
    stop('The optional nominal p-value column must be complete and finite when supplied.', call. = FALSE)
  }
  invalid_probability <- plot_data$padj < 0 | plot_data$padj > 1
  if (!all(is.na(plot_data$pvalue))) {
    invalid_probability <- invalid_probability |
      plot_data$pvalue < 0 | plot_data$pvalue > 1
  }
  if (any(invalid_probability)) {
    stop('GSEA p-values must be between 0 and 1.', call. = FALSE)
  }
  plot_data$pval <- plot_data$pvalue

  plot_data
}

#' Read one precomputed GSEA result table
#'
#' The input must contain one row per unique term and numeric NES and adjusted
#' p-value columns. A nominal p-value column is optional. No filtering or GSEA
#' calculation is performed. Matching uses `id_col` when supplied, then
#' `go_term_id` when present, and otherwise `term_col`.
#'
#' @param path Path to one CSV file
#' @param term_col Column containing term descriptions
#' @param nes_col Column containing normalized enrichment scores
#' @param pvalue_col Optional column containing nominal GSEA p-values. When
#'   omitted, `pval` is used automatically when present.
#' @param padj_col Column containing adjusted p-values
#' @param id_col Optional column containing unique stable term identifiers;
#'   overrides an existing `go_term_id` column
#'
#' @return A data frame with `go_term_id`, `go_description`, `NES`, `pvalue`,
#'   `pval`, and `padj` columns. The nominal p-value columns contain `NA` when
#'   no nominal p-value column was supplied.
#'
read_gsea_result_csv <- function(path,
  term_col = 'go_description',
  nes_col = 'NES',
  pvalue_col = NULL,
  padj_col = 'padj',
  id_col = NULL) {
  .gsea_standardize_results(
    gsea_results = utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE),
    term_col = term_col,
    nes_col = nes_col,
    pvalue_col = pvalue_col,
    padj_col = padj_col,
    id_col = id_col)
}

# 2.0 define shared plot utilities -----------------

.gsea_neg_log10 <- function(x) {
  -log10(pmax(x, .Machine$double.xmin))
}

.gsea_neutral_color <- '#9E9E9E'
.gsea_nonsignificant_color <- '#D9D9D9'

.gsea_category_colors <- c(
  blue = '#0055AAFF',
  red = '#BE3428FF',
  purple = '#6F4FA3FF',
  green = '#2CB11BFF')

.gsea_direction_colors <- c(
  'Positive' = unname(.gsea_category_colors[['red']]),
  'Negative' = unname(.gsea_category_colors[['blue']]),
  'Not significant' = .gsea_nonsignificant_color)

.gsea_default_category_colors <- function(group_names) {
  named_defaults <- c(
    'Ion channel' = unname(.gsea_category_colors[['blue']]),
    'Metabolism' = unname(.gsea_category_colors[['red']]),
    'Muscle contraction' = unname(.gsea_category_colors[['purple']]),
    'Heart development' = unname(.gsea_category_colors[['green']]),
    'Ribosomal' = unname(.gsea_category_colors[['blue']]),
    'Mitosis' = unname(.gsea_category_colors[['red']]),
    'DNA replication' = unname(.gsea_category_colors[['purple']]))
  default_colors <- unname(named_defaults[group_names])
  missing_colors <- is.na(default_colors)
  if (any(missing_colors)) {
    fallback_pool <- unname(.gsea_category_colors)
    missing_count <- sum(missing_colors)
    fallback_colors <- if (missing_count <= length(fallback_pool)) {
      fallback_pool[seq_len(missing_count)]
    } else {
      grDevices::colorRampPalette(fallback_pool)(missing_count)
    }
    default_colors[missing_colors] <- fallback_colors
  }
  stats::setNames(default_colors, group_names)
}

.gsea_p_value_axis_title <- function(p_col) {
  p_value_label <- if (p_col == 'padj') 'adjusted p-value' else 'nominal p-value'
  bquote(-log[10](.(p_value_label)))
}

.gsea_validate_number <- function(value,
  parameter_name,
  minimum = -Inf,
  maximum = Inf) {
  if (length(value) != 1L || !is.numeric(value) || is.na(value) ||
      !is.finite(value) || value < minimum || value > maximum) {
    stop('`', parameter_name, '` has an invalid numeric value.', call. = FALSE)
  }
  as.numeric(value)
}

.gsea_filter_by_padj <- function(plot_tbl,
  padj_threshold,
  padj_cols = 'padj') {
  if (is.null(padj_threshold)) {
    return(plot_tbl)
  }
  padj_threshold <- .gsea_validate_number(
    padj_threshold,
    'padj_threshold',
    minimum = 0,
    maximum = 1)
  keep_rows <- Reduce(
    `|`,
    lapply(padj_cols, function(padj_col) plot_tbl[[padj_col]] < padj_threshold))
  plot_tbl[keep_rows, , drop = FALSE]
}

.gsea_validate_count <- function(value, parameter_name, minimum = 0L) {
  if (length(value) != 1L || !is.numeric(value) || is.na(value) ||
      !is.finite(value) || value != floor(value) || value < minimum ||
      value > .Machine$integer.max) {
    stop('`', parameter_name, '` must be an integer of at least ', minimum, '.', call. = FALSE)
  }
  as.integer(value)
}

.gsea_validate_flag <- function(value, parameter_name) {
  if (length(value) != 1L || !is.logical(value) || is.na(value)) {
    stop('`', parameter_name, '` must be `TRUE` or `FALSE`.', call. = FALSE)
  }
  value
}

.gsea_validate_limits <- function(lower, upper, lower_name, upper_name) {
  if (!is.null(lower)) {
    lower <- .gsea_validate_number(lower, lower_name)
  }
  if (!is.null(upper)) {
    upper <- .gsea_validate_number(upper, upper_name)
  }
  if (!is.null(lower) && !is.null(upper) && lower >= upper) {
    stop('`', lower_name, '` must be less than `', upper_name, '`.', call. = FALSE)
  }
  invisible(NULL)
}

.gsea_validate_term_groups <- function(term_groups,
  parameter_name = 'term_groups') {
  valid_values <- is.list(term_groups) && length(term_groups) > 0L &&
    all(vapply(term_groups, function(group) {
      length(group) > 0L && !anyNA(group) && all(nzchar(trimws(as.character(group))))
    }, logical(1)))
  valid_names <- !is.null(names(term_groups)) &&
    !anyNA(names(term_groups)) &&
    all(nzchar(trimws(names(term_groups)))) &&
    !anyDuplicated(names(term_groups)) &&
    !'Other' %in% names(term_groups)
  if (!valid_values || !valid_names) {
    stop(
      '`', parameter_name,
      '` must be a non-empty uniquely named list of non-empty seeds; ',
      '`Other` is reserved for unmatched terms.',
      call. = FALSE)
  }
  invisible(NULL)
}

.gsea_resolve_named_colors <- function(color_values,
  required_names,
  parameter_name) {
  if (is.null(color_values)) {
    return(NULL)
  }
  if (is.null(names(color_values)) || anyDuplicated(names(color_values))) {
    stop('`', parameter_name, '` must have unique names.', call. = FALSE)
  }
  missing_colors <- setdiff(required_names, names(color_values))
  if (length(missing_colors) > 0L) {
    stop(
      '`',
      parameter_name,
      '` is missing color(s) for: ',
      paste(missing_colors, collapse = ', '),
      call. = FALSE)
  }
  color_values[required_names]
}

.gsea_resolve_term_group_colors <- function(term_groups = NULL,
  term_group_colors = NULL) {
  if (is.null(term_groups)) {
    return(c('Other' = .gsea_neutral_color))
  }
  group_names <- names(term_groups)
  if (is.null(term_group_colors)) {
    term_group_colors <- .gsea_default_category_colors(group_names)
  }

  term_group_colors <- .gsea_resolve_named_colors(
    color_values = term_group_colors,
    required_names = group_names,
    parameter_name = 'term_group_colors')
  c(term_group_colors[group_names], 'Other' = .gsea_neutral_color)
}

.gsea_wrap_label <- function(labels,
  words_per_line = 2L) {
  words_per_line <- .gsea_validate_count(words_per_line, 'label_words_per_line', minimum = 1L)

  vapply(as.character(labels), function(label) {
    words <- strsplit(trimws(label), '\\s+', perl = TRUE)[[1L]]
    if (length(words) <= words_per_line) {
      return(paste(words, collapse = ' '))
    }
    word_groups <- split(words, ceiling(seq_along(words) / words_per_line))
    paste(
      vapply(word_groups, paste, character(1L), collapse = ' '),
      collapse = '\n')
  }, character(1L), USE.NAMES = FALSE)
}

.gsea_split_point_layers <- function(plot_tbl,
  group_col,
  background_groups) {
  group_values <- as.character(plot_tbl[[group_col]])
  background_rows <- is.na(group_values) | group_values %in% background_groups
  list(
    background = plot_tbl[background_rows, , drop = FALSE],
    foreground = plot_tbl[!background_rows, , drop = FALSE])
}

# keep terms neutral when callers do not supply category seeds
.gsea_classify_terms <- function(go_description) {
  factor(rep('Other', length(go_description)), levels = 'Other')
}

.gsea_match_term_groups <- function(gsea_results,
  term_groups,
  parameter_name = 'term_groups') {
  .gsea_validate_term_groups(term_groups, parameter_name = parameter_name)

  match_text_seed <- function(seed, descriptions) {
    escaped_seed <- gsub('([][{}()+*^$|\\\\?.])', '\\\\\\1', seed)
    grepl(paste0('\\b', escaped_seed, '\\b'), descriptions, ignore.case = TRUE, perl = TRUE)
  }

  term_group <- rep('Other', nrow(gsea_results))

  for (group_name in names(term_groups)) {
    seeds <- as.character(term_groups[[group_name]])
    matched <- rep(FALSE, nrow(gsea_results))
    for (seed in seeds) {
      id_match <- tolower(gsea_results$go_term_id) == tolower(seed)
      term_match <- match_text_seed(seed, gsea_results$go_description)
      matched <- matched | id_match | term_match
    }
    replace_idx <- term_group == 'Other' & matched
    term_group[replace_idx] <- group_name
  }

  term_group
}

.gsea_assign_term_groups <- function(gsea_results,
  term_groups = NULL,
  parameter_name = 'term_groups') {
  if (is.null(term_groups)) {
    return(.gsea_classify_terms(gsea_results$go_description))
  }
  matched_group <- .gsea_match_term_groups(
    gsea_results = gsea_results,
    term_groups = term_groups,
    parameter_name = parameter_name)
  factor(matched_group, levels = c(names(term_groups), 'Other'))
}

.gsea_select_labels <- function(plot_tbl,
  label_terms = NULL,
  label_n = 0L,
  rank_col = NULL) {
  if (nrow(plot_tbl) == 0L) {
    return(plot_tbl[0L, , drop = FALSE])
  }

  if (!is.null(label_terms) && length(label_terms) > 0L) {
    label_rows <- integer()
    for (requested_label in label_terms) {
      if (is.numeric(requested_label)) {
        if (is.null(rank_col) || !rank_col %in% names(plot_tbl)) {
          warning('Rank label requested, but no rank column is available; skipping: ', requested_label, call. = FALSE)
          next
        }
        matched_row <- which(plot_tbl[[rank_col]] == requested_label)
      } else {
        requested_label <- as.character(requested_label)
        matched_row <- which(tolower(plot_tbl$go_term_id) == tolower(requested_label))
        if (length(matched_row) == 0L) {
          matched_row <- which(tolower(plot_tbl$go_description) == tolower(requested_label))
        }
        if (length(matched_row) > 1L) {
          stop(
            'Requested GO description matches multiple term keys: ',
            requested_label,
            '. Supply a unique term ID instead.',
            call. = FALSE)
        }
      }
      if (length(matched_row) == 0L) {
        warning('Requested GO label was not found in the plotted terms; skipping: ', requested_label, call. = FALSE)
        next
      }
      label_rows <- c(label_rows, matched_row[[1L]])
    }
    label_tbl <- plot_tbl[unique(label_rows), , drop = FALSE]
    return(label_tbl)
  }

  if (label_n <= 0L) {
    return(plot_tbl[0L, , drop = FALSE])
  }
  if (is.null(rank_col) || !rank_col %in% names(plot_tbl)) {
    plot_tbl <- plot_tbl[order(-plot_tbl$label_score, plot_tbl$padj, plot_tbl$go_description), , drop = FALSE]
    return(utils::head(plot_tbl, label_n))
  }

  plot_tbl$label_bin <- cut(
    plot_tbl[[rank_col]],
    breaks = seq(1, max(plot_tbl[[rank_col]]) + 1, length.out = label_n + 1L),
    include.lowest = TRUE,
    labels = FALSE)
  label_tbl <- do.call(rbind, lapply(sort(unique(plot_tbl$label_bin)), function(label_bin) {
    bin_tbl <- plot_tbl[plot_tbl$label_bin == label_bin, , drop = FALSE]
    bin_tbl <- bin_tbl[order(-bin_tbl$label_score, bin_tbl$padj, bin_tbl$go_description), , drop = FALSE]
    utils::head(bin_tbl, 1L)
  }))
  label_tbl$label_bin <- NULL
  label_tbl[order(label_tbl[[rank_col]]), , drop = FALSE]
}

.gsea_select_volcano_labels <- function(plot_tbl,
  label_n = 12L) {
  if (label_n <= 0L || nrow(plot_tbl) == 0L) {
    return(plot_tbl[0L, , drop = FALSE])
  }

  positive_label_n <- ceiling(label_n / 2L)
  negative_label_n <- floor(label_n / 2L)
  positive_tbl <- plot_tbl[plot_tbl$significant & plot_tbl$NES > 0, , drop = FALSE]
  negative_tbl <- plot_tbl[plot_tbl$significant & plot_tbl$NES < 0, , drop = FALSE]
  positive_tbl <- positive_tbl[order(positive_tbl$NES, positive_tbl$padj), , drop = FALSE]
  negative_tbl <- negative_tbl[order(negative_tbl$NES, negative_tbl$padj), , drop = FALSE]
  positive_tbl$label_rank <- seq_len(nrow(positive_tbl))
  negative_tbl$label_rank <- seq_len(nrow(negative_tbl))

  label_tbl <- rbind(
    .gsea_select_labels(positive_tbl, label_n = positive_label_n, rank_col = 'label_rank'),
    .gsea_select_labels(negative_tbl, label_n = negative_label_n, rank_col = 'label_rank'))
  if (nrow(label_tbl) < label_n) {
    remaining_tbl <- plot_tbl[plot_tbl$significant & !plot_tbl$go_term_id %in% label_tbl$go_term_id, , drop = FALSE]
    remaining_tbl <- remaining_tbl[order(-remaining_tbl$label_score, remaining_tbl$padj), , drop = FALSE]
    label_tbl <- rbind(label_tbl, utils::head(remaining_tbl, label_n - nrow(label_tbl)))
  }
  label_tbl$label_rank <- NULL
  label_tbl[order(label_tbl$NES), , drop = FALSE]
}

.gsea_geom_text_repel <- function(data,
  mapping,
  label_color = NULL,
  ...) {
  layer_args <- list(data = data, mapping = mapping, ...)
  if (!is.null(label_color)) {
    layer_args$color <- label_color
  }
  do.call(ggrepel::geom_text_repel, layer_args)
}

.gsea_normalize_rank_by <- function(rank_by) {
  rank_by <- match.arg(rank_by, c('NES', 'padj', 'pvalue', 'pval'))
  if (rank_by == 'pval') {
    return('pvalue')
  }
  rank_by
}

.gsea_rank_direction_terms <- function(plot_tbl,
  direction,
  rank_by = 'NES') {
  if (rank_by == 'NES') {
    if (direction == 'positive') {
      return(plot_tbl[order(-plot_tbl$NES, plot_tbl$padj, plot_tbl$go_description), , drop = FALSE])
    }
    return(plot_tbl[order(plot_tbl$NES, plot_tbl$padj, plot_tbl$go_description), , drop = FALSE])
  }

  rank_values <- plot_tbl[[rank_by]]
  nes_tiebreak <- if (direction == 'positive') -plot_tbl$NES else plot_tbl$NES
  plot_tbl[order(rank_values, nes_tiebreak, plot_tbl$go_description, na.last = TRUE), , drop = FALSE]
}

.gsea_axis_limits <- function(values,
  buffer_fraction = 0.045) {
  value_range <- range(values)
  value_span <- diff(value_range)
  if (value_span == 0) {
    value_span <- max(abs(value_range), 1)
  }
  value_range + c(-1, 1) * value_span * buffer_fraction
}

.gsea_resolve_axis_limits <- function(values,
  lower = NULL,
  upper = NULL,
  buffer_fraction = 0.18) {
  axis_limits <- .gsea_axis_limits(values, buffer_fraction = buffer_fraction)
  if (!is.null(lower)) {
    axis_limits[[1L]] <- lower
  }
  if (!is.null(upper)) {
    axis_limits[[2L]] <- upper
  }
  axis_limits
}

.gsea_theme <- function(font_family) {
  ggplot2::theme_classic(base_family = font_family, base_size = 9) +
    ggplot2::theme(
      axis.title = ggplot2::element_text(size = 9, color = 'black'),
      axis.text = ggplot2::element_text(size = 8, color = 'black'),
      axis.line = ggplot2::element_line(color = 'black', linewidth = 0.24),
      axis.ticks = ggplot2::element_line(color = 'black', linewidth = 0.20),
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      legend.title = ggplot2::element_text(size = 8, color = 'black', face = 'bold'),
      legend.text = ggplot2::element_text(size = 7.5, color = 'black'),
      legend.key.size = grid::unit(0.25, 'cm'))
}

# ---- waterfall plot function ----

# 1.0 create waterfall plot -----------------

#' Make a GSEA waterfall plot
#'
#' Positive and negative plots exclude zero NES values. Terms are ranked before
#' `top_n` is applied. Key precedence is explicit `id_col`, then `go_term_id`
#' when present, otherwise `term_col`.
#'
#' @param gsea_results GSEA result table.
#' @param term_col Column containing term descriptions and the fallback key.
#' @param nes_col Column containing normalized enrichment scores.
#' @param pvalue_col Optional column containing nominal GSEA p-values. It is
#'   needed only when `rank_by = 'pvalue'` and the input does not already use
#'   the standardized `pval` column.
#' @param padj_col Column containing adjusted p-values.
#' @param padj_threshold Optional adjusted p-value threshold applied before
#'   direction filtering, ranking, and plotting. Terms must satisfy
#'   `padj < padj_threshold`. Use `NULL` to disable this filter.
#' @param id_col Optional column containing unique stable term identifiers;
#'   overrides an existing `go_term_id` column.
#' @param NES_direction NES direction to plot.
#' @param top_n Number of terms to show.
#' @param rank_by Column used to rank terms before taking `top_n`. Use `NES`,
#'   `padj`, `pvalue`, or `pval`.
#' @param x_axis_title Optional x-axis title. Defaults to
#'   `Ranked GO terms (by <rank_by>)`.
#' @param title Optional plot title. Defaults to a direction-specific positive
#'   or negative NES waterfall title.
#' @param label_by_groups Optional exact term descriptions or GO IDs to label.
#' @param label_ranks Optional plotted ranks to label after ranking and
#'   filtering.
#' @param label_n Number of labels to automatically select across the ranked
#'   terms when neither `label_by_groups` nor `label_ranks` is supplied.
#' @param color_by_groups Optional named GO category seeds for coloring. When
#'   omitted, all terms are assigned to the neutral `Other` category.
#' @param term_group_colors Optional named colors for `color_by_groups`.
#' @param point_size Point size.
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
plot_gsea_waterfall <- function(gsea_results,
  term_col = 'go_description',
  nes_col = 'NES',
  pvalue_col = NULL,
  padj_col = 'padj',
  padj_threshold = NULL,
  id_col = NULL,
  NES_direction = c('positive', 'negative'), # nolint: object_name_linter.
  top_n = 100,
  rank_by = c('NES', 'padj', 'pvalue', 'pval'),
  x_axis_title = NULL,
  title = NULL,
  label_by_groups = NULL,
  label_ranks = NULL,
  label_n = 12L,
  color_by_groups = NULL,
  term_group_colors = NULL,
  point_size = 1.80,
  label_size = 8.0,
  label_fontface = 'plain',
  label_words_per_line = 2L,
  label_repel_direction = 'y',
  label_x_nudge = NULL,
  label_y_nudge_fraction = 0.11,
  y_min = NULL,
  y_max = NULL,
  legend_position = NULL,
  font_family = 'Nimbus Sans') {
  NES_direction <- match.arg(NES_direction) # nolint: object_name_linter.
  rank_by <- .gsea_normalize_rank_by(rank_by)
  top_n <- .gsea_validate_count(top_n, 'top_n', minimum = 1L)
  point_size <- .gsea_validate_number(point_size, 'point_size', minimum = 0)
  if (length(label_by_groups) == 0L && length(label_ranks) == 0L) {
    label_n <- .gsea_validate_count(label_n, 'label_n')
  }
  .gsea_validate_limits(y_min, y_max, 'y_min', 'y_max')
  if (!is.null(label_ranks) &&
      (!is.numeric(label_ranks) || anyNA(label_ranks) ||
          any(!is.finite(label_ranks)) || any(label_ranks != floor(label_ranks)) ||
          any(label_ranks < 1) || any(label_ranks > .Machine$integer.max))) {
    stop('`label_ranks` must contain positive integers.', call. = FALSE)
  }
  gsea_results <- .gsea_standardize_results(
    gsea_results = gsea_results,
    term_col = term_col,
    nes_col = nes_col,
    pvalue_col = pvalue_col,
    padj_col = padj_col,
    id_col = id_col)
  gsea_results <- .gsea_filter_by_padj(
    plot_tbl = gsea_results,
    padj_threshold = padj_threshold)
  if (rank_by == 'pvalue' && all(is.na(gsea_results$pvalue))) {
    stop("`rank_by = 'pvalue'` requires an optional nominal p-value column.", call. = FALSE)
  }

  if (NES_direction == 'positive') {
    plot_tbl <- gsea_results[gsea_results$NES > 0, , drop = FALSE]
  } else {
    plot_tbl <- gsea_results[gsea_results$NES < 0, , drop = FALSE]
  }

  if (nrow(plot_tbl) == 0L) {
    stop(
      'No GSEA terms remain after applying the adjusted p-value threshold and requested direction.',
      call. = FALSE)
  }

  plot_tbl <- .gsea_rank_direction_terms(
    plot_tbl,
    direction = NES_direction,
    rank_by = rank_by)
  plot_tbl <- utils::head(plot_tbl, top_n)
  plot_tbl$waterfall_rank <- seq_len(nrow(plot_tbl))
  plot_tbl$label_score <- abs(plot_tbl$NES) * .gsea_neg_log10(plot_tbl$padj)
  plot_tbl$term_group <- .gsea_assign_term_groups(
    gsea_results = plot_tbl,
    term_groups = color_by_groups,
    parameter_name = 'color_by_groups')
  color_values <- .gsea_resolve_term_group_colors(
    term_groups = color_by_groups,
    term_group_colors = term_group_colors)
  label_tbl <- .gsea_select_labels(
    plot_tbl = plot_tbl,
    label_terms = c(as.list(label_by_groups), as.list(as.integer(label_ranks))),
    label_n = label_n,
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
    if (NES_direction == 'negative') {
      label_tbl$label_nudge_y <- -label_tbl$label_nudge_y
    }
  }
  x_breaks <- unique(c(0, scales::breaks_pretty(n = 4)(c(1, nrow(plot_tbl)))))
  if (is.null(x_axis_title)) {
    x_axis_title <- paste0('Ranked GO terms (by ', rank_by, ')')
  }
  if (is.null(title)) {
    title <- if (NES_direction == 'positive') {
      'Positive NES waterfall'
    } else {
      'Negative NES waterfall'
    }
  }
  if (is.null(legend_position)) {
    legend_position <- if (NES_direction == 'positive') c(0.98, 0.98) else c(0.03, 0.98)
  }
  legend_anchor <- if (NES_direction == 'positive') c(1, 1) else c(0, 1)
  point_layers <- .gsea_split_point_layers(
    plot_tbl = plot_tbl,
    group_col = 'term_group',
    background_groups = 'Other')

  ggplot2::ggplot(plot_tbl, ggplot2::aes(x = .data$waterfall_rank, y = .data$NES)) +
    ggplot2::geom_point(
      data = point_layers$background,
      ggplot2::aes(color = .data$term_group),
      size = point_size,
      alpha = 1,
      stroke = 0) +
    ggplot2::geom_point(
      data = point_layers$foreground,
      ggplot2::aes(color = .data$term_group),
      size = point_size,
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
    ggplot2::labs(
      title = title,
      x = x_axis_title,
      y = 'Normalized enrichment score (NES)',
      color = 'GO:BP category') +
    ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(size = point_size))) +
    .gsea_theme(font_family) +
    ggplot2::theme(
      axis.ticks.length = grid::unit(1.4, 'pt'),
      legend.position = legend_position,
      legend.justification = legend_anchor,
      legend.background = ggplot2::element_blank(),
      legend.box.background = ggplot2::element_blank(),
      legend.title = ggplot2::element_text(size = 8, hjust = 0),
      legend.text = ggplot2::element_text(size = 7.5, color = 'black'),
      legend.key.size = grid::unit(0.25, 'cm'),
      legend.box.just = 'left',
      plot.title = ggplot2::element_text(size = 10, face = 'bold', color = 'black', hjust = 0),
      plot.margin = ggplot2::margin(8, 10, 8, 8))
}

# ---- symmetric volcano plot function ----

# 1.0 create symmetric volcano plot -----------------

# retain the established SVG panel footprint with the portable y-axis label
.gsea_volcano_y_title_margin_pt <- 5.42

#' Make a symmetric GSEA volcano plot
#'
#' All validated terms are retained. Significance is defined as
#' `padj < padj_cutoff`, and nonsignificant terms are shown in grey. Zero
#' adjusted p-values are displayed at the finite cap controlled by `y_max`;
#' supplied values are not modified in the input table.
#'
#' @param gsea_results GSEA result table.
#' @param term_col Column containing term descriptions and the fallback key.
#' @param nes_col Column containing normalized enrichment scores.
#' @param padj_col Column containing adjusted p-values.
#' @param id_col Optional column containing unique stable term identifiers;
#'   overrides an existing `go_term_id` column.
#' @param padj_cutoff Adjusted p-value cutoff for coloring and shaded regions.
#' @param label_n Number of significant terms to label.
#' @param label_terms Optional exact term descriptions or GO IDs to label.
#' @param y_max Finite display cap for -log10 adjusted p-values. The coordinate
#'   range includes additional annotation space above this value.
#' @param y_min Lower y-axis limit for capped -log10 adjusted p-values.
#' @param point_size Point size.
#' @param point_colors Optional named colors for `Positive`, `Negative`, and
#'   `Not significant`.
#' @param label_size Label font size in points.
#' @param label_color Optional text color for all GO term labels. When `NULL`,
#'   labels use the color of their corresponding point group.
#' @param label_nudge_x Horizontal starting offset for labeled terms.
#' @param label_nudge_y Vertical starting offset for labeled terms.
#' @param count_label_size Font size in points for significant-count labels.
#' @param contrast_label Optional concise contrast label for the NES axis.
#' @param label_words_per_line Number of GO term words shown on each label line.
#' @param legend_position Legend position passed to `ggplot2::theme()`.
#' @param font_family Figure font family.
#'
#' @return A ggplot object without file-system side effects
#'
plot_gsea_volcano <- function(gsea_results,
  term_col = 'go_description',
  nes_col = 'NES',
  padj_col = 'padj',
  id_col = NULL,
  padj_cutoff = 0.05,
  label_n = 14L,
  label_terms = NULL,
  y_max = 25,
  y_min = 0,
  point_size = 1.3,
  point_colors = NULL,
  label_size = 8.0,
  label_color = NULL,
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
  point_size <- .gsea_validate_number(point_size, 'point_size', minimum = 0)
  if (length(label_terms) == 0L) {
    label_n <- .gsea_validate_count(label_n, 'label_n')
  }
  .gsea_validate_limits(y_min, y_max, 'y_min', 'y_max')
  plot_tbl <- .gsea_standardize_results(
    gsea_results = gsea_results,
    term_col = term_col,
    nes_col = nes_col,
    padj_col = padj_col,
    id_col = id_col)
  plot_tbl$neg_log10_padj <- .gsea_neg_log10(plot_tbl$padj)
  plot_tbl$plot_neg_log10_padj <- pmin(plot_tbl$neg_log10_padj, y_max)
  plot_tbl$significant <- plot_tbl$padj < padj_cutoff
  plot_tbl$point_group <- ifelse(
    plot_tbl$significant & plot_tbl$NES >= 0,
    'Positive',
    ifelse(plot_tbl$significant & plot_tbl$NES < 0, 'Negative', 'Not significant'))
  plot_tbl$point_group <- factor(
    plot_tbl$point_group,
    levels = c('Positive', 'Negative', 'Not significant'))
  color_values <- .gsea_direction_colors[levels(plot_tbl$point_group)]
  custom_colors <- .gsea_resolve_named_colors(
    color_values = point_colors,
    required_names = levels(plot_tbl$point_group),
    parameter_name = 'point_colors')
  if (!is.null(custom_colors)) {
    color_values <- custom_colors
  }
  plot_tbl$label_score <- abs(plot_tbl$NES) * plot_tbl$plot_neg_log10_padj
  if (length(label_terms) > 0L) {
    label_tbl <- .gsea_select_labels(plot_tbl, label_terms = as.list(label_terms))
  } else {
    label_tbl <- .gsea_select_volcano_labels(plot_tbl, label_n = label_n)
  }
  x_limit <- max(abs(plot_tbl$NES), na.rm = TRUE)
  x_limit <- max(1, x_limit) * 1.10
  x_limits <- c(-x_limit, x_limit)
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
  repel_tbl$label_nudge_x[label_match] <- ifelse(
    label_tbl$NES < 0,
    label_nudge_x,
    -label_nudge_x)
  repel_tbl$label_nudge_y[label_match] <- label_nudge_y +
    rep(c(0.18, -0.12, 0.08, -0.18), length.out = length(label_match))
  negative_count <- sum(plot_tbl$significant & plot_tbl$NES < 0, na.rm = TRUE)
  positive_count <- sum(plot_tbl$significant & plot_tbl$NES > 0, na.rm = TRUE)
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

  label_mapping <- if (is.null(label_color)) {
    ggplot2::aes(label = .data$label_text, color = .data$point_group)
  } else {
    ggplot2::aes(label = .data$label_text)
  }

  ggplot2::ggplot(plot_tbl, ggplot2::aes(x = .data$NES, y = .data$plot_neg_log10_padj)) +
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
    .gsea_geom_text_repel(
      data = repel_tbl,
      mapping = label_mapping,
      label_color = label_color,
      nudge_x = repel_tbl$label_nudge_x,
      nudge_y = repel_tbl$label_nudge_y,
      family = font_family,
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
      label = paste0('Negative: ', negative_count),
      hjust = 0,
      vjust = 0,
      family = font_family,
      size = count_label_size / ggplot2::.pt,
      fontface = 'bold',
      color = 'white',
      fill = color_values[['Negative']],
      linewidth = 0,
      label.padding = grid::unit(0.12, 'lines'),
      label.r = grid::unit(0.08, 'lines')) +
    ggplot2::annotate(
      'label',
      x = x_limits[[2L]],
      y = count_y,
      label = paste0('Positive: ', positive_count),
      hjust = 1,
      vjust = 0,
      family = font_family,
      size = count_label_size / ggplot2::.pt,
      fontface = 'bold',
      color = 'white',
      fill = color_values[['Positive']],
      linewidth = 0,
      label.padding = grid::unit(0.12, 'lines'),
      label.r = grid::unit(0.08, 'lines')) +
    ggplot2::scale_color_manual(values = color_values, drop = TRUE) +
    ggplot2::scale_x_continuous(
      trans = scales::pseudo_log_trans(sigma = 1),
      breaks = scales::breaks_pretty(n = 7)(x_limits),
      expand = ggplot2::expansion(mult = 0.02)) +
    ggplot2::scale_y_continuous(
      trans = scales::pseudo_log_trans(sigma = 1),
      breaks = c(0, 1, 2, 5, 10, 25),
      expand = c(0, 0)) +
    ggplot2::coord_cartesian(xlim = x_limits, ylim = y_limits, clip = 'off') +
    ggplot2::labs(x = x_axis_label, y = .gsea_p_value_axis_title('padj')) +
    ggplot2::guides(color = ggplot2::guide_legend(
      title = 'fgsea NES direction',
      nrow = 1,
      byrow = TRUE,
      title.position = 'top',
      title.hjust = 0.5,
      override.aes = list(size = point_size + 0.8))) +
    .gsea_theme(font_family) +
    ggplot2::theme(
      axis.title.y = ggplot2::element_text(margin = ggplot2::margin(r = .gsea_volcano_y_title_margin_pt)),
      legend.position = legend_position,
      legend.direction = 'horizontal',
      legend.box = 'horizontal',
      legend.justification = 'center',
      plot.margin = ggplot2::margin(10, 18, 8, 12))
}

# ---- directional half-volcano plot function ----

# 1.0 create directional volcano plot -----------------

#' Make a directional half-volcano plot for GSEA results
#'
#' The selected direction and inclusive `inner_nes_limit` are applied before
#' plotting. Significant and nonsignificant terms within that boundary are
#' retained. Significance is defined as `padj < padj_cutoff`, and a threshold
#' line is drawn only when adjusted p-values are shown on the y-axis.
#'
#' @param gsea_results GSEA result table.
#' @param term_col Column containing term descriptions and the fallback key.
#' @param nes_col Column containing normalized enrichment scores.
#' @param pvalue_col Optional column containing nominal GSEA p-values. It is
#'   needed only when `p_col = 'pvalue'` and the input does not already use the
#'   standardized `pval` column.
#' @param padj_col Column containing adjusted p-values.
#' @param id_col Optional column containing unique stable term identifiers;
#'   overrides an existing `go_term_id` column.
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
#'   `Not significant` when `term_groups` is omitted. Direction labels are
#'   `Positive` and `Negative`.
#' @param label_size Label font size in points.
#' @param label_color Optional text color for all GO term labels. When `NULL`,
#'   labels use the color of their corresponding point group.
#' @param contrast_label Optional concise contrast label for the NES axis.
#' @param label_words_per_line Number of GO term words shown on each label line.
#' @param y_min Optional lower y-axis limit.
#' @param y_max Optional upper y-axis limit.
#' @param legend_position Legend position passed to `ggplot2::theme()`.
#' @param font_family Figure font family.
#'
#' @return A ggplot object without file-system side effects
#'
plot_gsea_half_volcano <- function(gsea_results,
  term_col = 'go_description',
  nes_col = 'NES',
  pvalue_col = NULL,
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
  point_size = 1.6,
  point_colors = NULL,
  label_size = 8.0,
  label_color = NULL,
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
  point_size <- .gsea_validate_number(point_size, 'point_size', minimum = 0)
  if (length(label_terms) == 0L) {
    label_n <- .gsea_validate_count(label_n, 'label_n')
  }
  inner_nes_limit <- .gsea_validate_number(
    inner_nes_limit,
    'inner_nes_limit',
    minimum = 0)
  .gsea_validate_limits(y_min, y_max, 'y_min', 'y_max')
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
  if (p_col == 'pvalue' && all(is.na(plot_tbl$pvalue))) {
    stop("`p_col = 'pvalue'` requires an optional nominal p-value column.", call. = FALSE)
  }

  direction_rows <- if (direction == 'positive') plot_tbl$NES > 0 else plot_tbl$NES < 0
  plot_tbl <- plot_tbl[direction_rows, , drop = FALSE]
  boundary_rows <- if (direction == 'positive') {
    plot_tbl$NES >= inner_nes_limit
  } else {
    plot_tbl$NES <= -inner_nes_limit
  }
  plot_tbl <- plot_tbl[boundary_rows, , drop = FALSE]
  if (nrow(plot_tbl) == 0L) {
    stop('No GSEA terms remain after applying the half-volcano NES boundary.', call. = FALSE)
  }
  plot_tbl$neg_log10_p <- .gsea_neg_log10(plot_tbl[[p_col]])
  plot_tbl$significant <- plot_tbl$padj < padj_cutoff
  if (is.null(term_groups)) {
    direction_label <- if (direction == 'positive') 'Positive' else 'Negative'
    plot_tbl$point_group <- ifelse(plot_tbl$significant, direction_label, 'Not significant')
    plot_tbl$point_group <- factor(plot_tbl$point_group, levels = c(direction_label, 'Not significant'))
    color_values <- .gsea_direction_colors[levels(plot_tbl$point_group)]
    custom_colors <- .gsea_resolve_named_colors(
      color_values = point_colors,
      required_names = levels(plot_tbl$point_group),
      parameter_name = 'point_colors')
    if (!is.null(custom_colors)) {
      color_values <- custom_colors
    }
    legend_title <- 'fgsea NES direction'
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
    color_values <- c(color_values, 'Not significant' = .gsea_nonsignificant_color)
    legend_title <- 'GO:BP category'
    background_groups <- c('Other', 'Not significant')
  }
  plot_tbl$label_score <- abs(plot_tbl$NES) * plot_tbl$neg_log10_p
  significant_tbl <- plot_tbl[plot_tbl$significant, , drop = FALSE]
  significant_tbl <- significant_tbl[order(abs(significant_tbl$NES), significant_tbl$padj), , drop = FALSE]
  significant_tbl$label_rank <- seq_len(nrow(significant_tbl))
  if (length(label_terms) > 0L) {
    label_tbl <- .gsea_select_labels(plot_tbl, label_terms = as.list(label_terms))
  } else {
    label_tbl <- .gsea_select_labels(significant_tbl, label_n = label_n, rank_col = 'label_rank')
  }
  label_tbl$label_text <- .gsea_wrap_label(label_tbl$go_description, words_per_line = label_words_per_line)
  repel_tbl <- plot_tbl
  repel_tbl$label_text <- ''
  repel_tbl$label_nudge_x <- 0
  repel_tbl$label_nudge_y <- 0
  label_match <- match(label_tbl$go_term_id, repel_tbl$go_term_id)
  repel_tbl$label_text[label_match] <- label_tbl$label_text
  label_direction <- if (direction == 'positive') 1 else -1
  repel_tbl$label_nudge_x[label_match] <- -label_direction * 0.10
  y_limits <- .gsea_resolve_axis_limits(
    plot_tbl$neg_log10_p,
    lower = y_min,
    upper = y_max,
    buffer_fraction = 0.035)
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
  label_y_limits <- y_limits
  if (p_col == 'padj' && padj_cutoff > 0) {
    label_y_limits[[1L]] <- max(label_y_limits[[1L]], -log10(padj_cutoff))
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
  significance_line <- if (p_col == 'padj' && padj_cutoff > 0) {
    ggplot2::geom_hline(
      yintercept = -log10(padj_cutoff),
      color = 'black',
      linewidth = 0.22,
      linetype = 'dashed')
  } else {
    NULL
  }

  label_mapping <- if (is.null(label_color)) {
    ggplot2::aes(label = .data$label_text, color = .data$point_group)
  } else {
    ggplot2::aes(label = .data$label_text)
  }

  ggplot2::ggplot(plot_tbl, ggplot2::aes(x = .data$NES, y = .data$neg_log10_p)) +
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
    significance_line +
    .gsea_geom_text_repel(
      data = repel_tbl,
      mapping = label_mapping,
      label_color = label_color,
      nudge_x = repel_tbl$label_nudge_x,
      nudge_y = repel_tbl$label_nudge_y,
      family = font_family,
      size = label_size / ggplot2::.pt,
      lineheight = 0.9,
      hjust = 0.5,
      min.segment.length = 0,
      segment.color = 'black',
      segment.size = 0.16,
      box.padding = 0.45,
      point.padding = 1.20,
      xlim = x_limits,
      ylim = label_y_limits,
      direction = 'both',
      force = 3.0,
      force_pull = 8.0,
      max.time = 12,
      max.iter = 30000,
      max.overlaps = Inf,
      seed = 42,
      show.legend = FALSE) +
    ggplot2::scale_color_manual(values = color_values, drop = TRUE) +
    ggplot2::scale_x_continuous(breaks = x_breaks, expand = ggplot2::expansion(mult = 0.015)) +
    ggplot2::scale_y_continuous(expand = c(0, 0), position = if (direction == 'negative') 'right' else 'left') +
    ggplot2::coord_cartesian(xlim = x_limits, ylim = y_limits, clip = 'off') +
    ggplot2::labs(x = x_axis_label, y = .gsea_p_value_axis_title(p_col), color = legend_title) +
    ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(size = point_size + 0.5))) +
    .gsea_theme(font_family) +
    ggplot2::theme(
      legend.position = legend_position,
      legend.justification = legend_anchor,
      legend.background = ggplot2::element_rect(fill = 'white', color = NA),
      plot.margin = ggplot2::margin(8, 10, 8, 10))
}

# ---- comparative NES scatterplot function ----

# 1.0 create comparative NES scatterplot -----------------

#' Make a comparative GSEA NES scatterplot
#'
#' Both inputs must contain complete, unique term keys. The plot uses their
#' exact inner intersection, retains descriptions from `gsea_x`, and does not
#' impute unmatched terms. All shared terms are retained by default;
#' significance is defined as `padj < padj_cutoff`.
#'
#' @param gsea_x First GSEA result table.
#' @param gsea_y Second GSEA result table.
#' @param x_label X-axis label.
#' @param y_label Y-axis label.
#' @param term_col Column containing term descriptions and the fallback key in
#'   both tables.
#' @param nes_col Column containing normalized enrichment scores in both tables.
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
#' @param equal_axis_limits Use matching x/y limits. For `q1` and `q3` this
#'   zooms to the selected quadrant while keeping the x and y ranges equal.
#'   Cannot be combined with explicit axis limits.
#' @param quadrant_min_abs_nes Inner absolute NES boundary for Q1 and Q3 plots
#'   when `equal_axis_limits = TRUE`.
#' @param x_min Optional lower x-axis limit.
#' @param x_max Optional upper x-axis limit.
#' @param y_min Optional lower y-axis limit.
#' @param y_max Optional upper y-axis limit.
#' @param show_fit_line Draw a descriptive linear best-fit line through plotted
#'   terms. The line is fitted after quadrant selection and is disabled by
#'   default.
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
  x_label = 'NES in x contrast',
  y_label = 'NES in y contrast',
  term_col = 'go_description',
  nes_col = 'NES',
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
  label_size = 8.0,
  label_color = 'black',
  label_words_per_line = 2L,
  padj_cutoff = 0.05,
  equal_axis_limits = FALSE,
  quadrant_min_abs_nes = 0,
  x_min = NULL,
  x_max = NULL,
  y_min = NULL,
  y_max = NULL,
  show_fit_line = FALSE,
  legend_nrow = NULL,
  legend_position = 'top',
  font_family = 'Nimbus Sans') {
  color_by <- match.arg(color_by)
  quadrant <- match.arg(quadrant)
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
    padj_col = padj_col,
    id_col = id_col)
  y_tbl <- .gsea_standardize_results(
    gsea_results = gsea_y,
    term_col = term_col,
    nes_col = nes_col,
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
  plot_tbl$label_score <- pmax(abs(plot_tbl$NES_x), abs(plot_tbl$NES_y)) *
    pmax(.gsea_neg_log10(plot_tbl$padj_x), .gsea_neg_log10(plot_tbl$padj_y))
  plot_tbl$scatter_quadrant <- ifelse(
    plot_tbl$NES_x >= 0 & plot_tbl$NES_y >= 0,
    'Q1',
    ifelse(
      plot_tbl$NES_x < 0 & plot_tbl$NES_y >= 0,
      'Q2',
      ifelse(plot_tbl$NES_x < 0 & plot_tbl$NES_y < 0, 'Q3', 'Q4')))

  if (color_by == 'significance') {
    plot_tbl$point_group <- ifelse(
      plot_tbl$significant_x & plot_tbl$significant_y,
      'Significant in both',
      ifelse(
        plot_tbl$significant_x & !plot_tbl$significant_y,
        paste0('Significant in ', x_name, ' only'),
        ifelse(
          !plot_tbl$significant_x & plot_tbl$significant_y,
          paste0('Significant in ', y_name, ' only'),
          'Not significant in either')))
    group_names <- c(
      'Significant in both',
      paste0('Significant in ', x_name, ' only'),
      paste0('Significant in ', y_name, ' only'),
      'Not significant in either')
    if (is.null(point_colors)) {
      color_values <- stats::setNames(
        c(.gsea_category_colors[['blue']],
          .gsea_category_colors[['red']],
          .gsea_category_colors[['purple']],
          .gsea_nonsignificant_color),
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
    label_order <- order(
      -label_tbl$label_score,
      label_tbl$padj_x,
      label_tbl$padj_y,
      label_tbl$go_description)
    label_tbl <- label_tbl[label_order, , drop = FALSE]
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
  repel_tbl <- plot_tbl
  repel_tbl$label_text <- ''
  repel_tbl$label_nudge_x <- 0
  repel_tbl$label_nudge_y <- 0
  label_match <- match(label_tbl$go_term_id, repel_tbl$go_term_id)
  repel_tbl$label_text[label_match] <- label_tbl$label_text
  repel_tbl$label_nudge_x[label_match] <- ifelse(label_tbl$NES_x >= 0, 0.05, -0.05)
  repel_tbl$label_nudge_y[label_match] <- ifelse(label_tbl$NES_y >= 0, 0.05, -0.05)

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
  label_x_limits <- if (is.null(axis_x_limits)) c(NA, NA) else axis_x_limits
  label_y_limits <- if (is.null(axis_y_limits)) c(NA, NA) else axis_y_limits
  label_layer <- ggrepel::geom_text_repel(
    data = repel_tbl,
    ggplot2::aes(label = .data$label_text),
    nudge_x = repel_tbl$label_nudge_x,
    nudge_y = repel_tbl$label_nudge_y,
    family = font_family,
    color = label_color,
    size = label_size / ggplot2::.pt,
    lineheight = 0.95,
    hjust = 0.5,
    min.segment.length = 0,
    segment.color = 'black',
    segment.size = 0.16,
    box.padding = 0.55,
    point.padding = 0.30,
    xlim = label_x_limits,
    ylim = label_y_limits,
    direction = 'both',
    force = 2.2,
    force_pull = 4.0,
    max.time = 12,
    max.iter = 30000,
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
    ggplot2::labs(
      x = x_label,
      y = y_label,
      color = if (color_by == 'significance') 'GSEA enrichment' else 'GO:BP category') +
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
