# ----
# author:
# - Zoheb Khan
#
# script path:
# - functions/gsea_plot_utils.R
# ----

# 1.0 normalize GSEA inputs -----------------

# independently supplied GSEA tables require an explicit unique term key
.gsea_standardize_results <- function(gsea_results,
                                      term_col = 'go_description',
                                      nes_col = 'NES',
                                      pvalue_col = 'pval',
                                      padj_col = 'padj',
                                      id_col = NULL) {
  if (is.null(id_col) && 'go_term_id' %in% names(gsea_results)) {
    id_col <- 'go_term_id'
  }
  required_cols <- c(term_col, nes_col, pvalue_col, padj_col, id_col)
  missing_cols <- setdiff(required_cols, names(gsea_results))
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

  numeric_columns <- c(nes_col, pvalue_col, padj_col)
  if (!all(vapply(gsea_results[numeric_columns], is.numeric, logical(1)))) {
    stop('GSEA NES and p-value columns must be numeric.', call. = FALSE)
  }

  plot_data <- data.frame(
    go_term_id = go_term_id,
    go_description = go_description,
    NES = as.numeric(gsea_results[[nes_col]]),
    pvalue = as.numeric(gsea_results[[pvalue_col]]),
    padj = as.numeric(gsea_results[[padj_col]]),
    stringsAsFactors = FALSE)

  if (anyNA(plot_data) || any(!is.finite(plot_data$NES))) {
    stop('GSEA plotting columns must contain complete finite values.', call. = FALSE)
  }
  invalid_probability <-
    plot_data$pvalue < 0 | plot_data$pvalue > 1 |
      plot_data$padj < 0 | plot_data$padj > 1
  if (any(invalid_probability)) {
    stop('GSEA p-values and adjusted p-values must be between 0 and 1.', call. = FALSE)
  }
  plot_data$pval <- plot_data$pvalue

  plot_data
}

#' Read one precomputed GSEA result table
#'
#' The input must contain one row per unique term and numeric NES, nominal
#' p-value, and adjusted p-value columns. No filtering or GSEA calculation is
#' performed. When `id_col` is omitted, `term_col` is the unique matching key.
#'
#' @param path Path to one CSV file
#' @param term_col Column containing term descriptions
#' @param nes_col Column containing normalized enrichment scores
#' @param pvalue_col Column containing nominal GSEA p-values
#' @param padj_col Column containing adjusted p-values
#' @param id_col Optional column containing unique stable term identifiers
#'
#' @return A data frame with `go_term_id`, `go_description`, `NES`, `pvalue`,
#'   `pval`, and `padj` columns
#'
#' @export
read_gsea_result_csv <- function(path,
                                 term_col = 'go_description',
                                 nes_col = 'NES',
                                 pvalue_col = 'pval',
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

.gsea_neutral_gray <- function() {
  '#9E9E9E'
}

.gsea_light_gray <- function() {
  '#D9D9D9'
}

.gsea_default_colors <- function() {
  c(
    'Significantly up' = '#CC79A7',
    'Significantly down' = '#0072B2',
    'Not significant' = .gsea_light_gray(),
    'Significant in both' = '#2A9D8F',
    'Significant in x only' = '#CC79A7',
    'Significant in y only' = '#0072B2',
    'Not significant in either' = .gsea_light_gray(),
    'Other' = .gsea_neutral_gray())
}

.gsea_category_colors <- function() {
  c(
    'Development/morphogenesis' = '#1B9E77',
    'Neuronal/signaling' = '#4F61BD',
    'Cell cycle/genome' = '#D73027',
    'Metabolism/translation' = '#7A5195',
    'Other' = .gsea_neutral_gray())
}

.gsea_validate_number <- function(value,
                                  parameter_name,
                                  minimum = -Inf,
                                  maximum = Inf,
                                  minimum_inclusive = TRUE,
                                  maximum_inclusive = TRUE) {
  if (length(value) != 1L || !is.numeric(value) || is.na(value) || !is.finite(value)) {
    stop('`', parameter_name, '` has an invalid numeric value.', call. = FALSE)
  }
  valid_minimum <- if (minimum_inclusive) value >= minimum else value > minimum
  valid_maximum <- if (maximum_inclusive) value <= maximum else value < maximum
  if (!valid_minimum || !valid_maximum) {
    stop('`', parameter_name, '` has an invalid numeric value.', call. = FALSE)
  }
  as.numeric(value)
}

.gsea_validate_count <- function(value, parameter_name, minimum = 0L) {
  if (length(value) != 1L || !is.numeric(value) || is.na(value) ||
    !is.finite(value) || value != floor(value) || value < minimum ||
    value > .Machine$integer.max) {
    stop('`', parameter_name, '` must be an integer of at least ', minimum, '.', call. = FALSE)
  }
  as.integer(value)
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

.gsea_validate_positive_integer_vector <- function(value, parameter_name) {
  if (is.null(value)) {
    return(invisible(NULL))
  }
  if (!is.numeric(value) || anyNA(value) || any(!is.finite(value)) ||
    any(value != floor(value)) || any(value < 1) ||
    any(value > .Machine$integer.max)) {
    stop('`', parameter_name, '` must contain positive integers.', call. = FALSE)
  }
  invisible(NULL)
}

.gsea_validate_term_groups <- function(term_groups = NULL) {
  if (is.null(term_groups)) {
    return(invisible(NULL))
  }
  valid_values <- is.list(term_groups) && length(term_groups) > 0L &&
    all(vapply(term_groups, function(group) {
      length(group) > 0L && !anyNA(group) && all(nzchar(trimws(as.character(group))))
    }, logical(1)))
  valid_names <- !is.null(names(term_groups)) &&
    !anyNA(names(term_groups)) &&
    all(nzchar(trimws(names(term_groups)))) &&
    !anyDuplicated(names(term_groups))
  if (!valid_values || !valid_names) {
    stop('`term_groups` must be a non-empty uniquely named list of non-empty seeds.', call. = FALSE)
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
                                            term_group_colors = NULL,
                                            default_group = 'Other') {
  if (is.null(term_groups)) {
    return(.gsea_category_colors())
  }
  group_names <- names(term_groups)
  if (is.null(term_group_colors)) {
    if (!requireNamespace('RColorBrewer', quietly = TRUE)) {
      stop('The RColorBrewer package is required for default term-group colors.', call. = FALSE)
    }
    brewer_n <- max(3L, length(group_names))
    term_group_colors <- stats::setNames(
      RColorBrewer::brewer.pal(brewer_n, 'Set1')[seq_along(group_names)],
      group_names)
  }

  term_group_colors <- .gsea_resolve_named_colors(
    color_values = term_group_colors,
    required_names = group_names,
    parameter_name = 'term_group_colors')
  c(term_group_colors[group_names], stats::setNames(.gsea_neutral_gray(), default_group))
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

# fallback categories keep ungrouped scatterplots interpretable
.gsea_classify_terms <- function(go_description) {
  term <- tolower(as.character(go_description))
  category <- rep('Other', length(term))

  category[grepl(
    'cell cycle|mitotic|meiotic|chromosom|spindle|centromer|kinetochore|dna replication|dna repair|recombination',
    term)] <- 'Cell cycle/genome'
  category[category == 'Other' & grepl(
    'development|morphogenesis|differentiation|mesoderm|mesenchyme|embryo|organogenesis|migration|pattern specification',
    term)] <- 'Development/morphogenesis'
  category[category == 'Other' & grepl(
    'synap|neuro|axon|dendrit|transmitter|ion channel|receptor|signaling|signal transduction|membrane potential|secretion',
    term)] <- 'Neuronal/signaling'
  category[category == 'Other' & grepl(
    'ribosom|translation|rrna|metabolic|biosynthetic|mitochond|oxidative|protein folding|rna processing',
    term)] <- 'Metabolism/translation'

  factor(category, levels = names(.gsea_category_colors()))
}

.gsea_match_term_groups <- function(gsea_results,
                                    term_groups,
                                    default_group = 'Other') {
  .gsea_validate_term_groups(term_groups)

  match_text_seed <- function(seed, descriptions) {
    escaped_seed <- gsub('([][{}()+*^$|\\\\?.])', '\\\\\\1', seed)
    grepl(paste0('\\b', escaped_seed, '\\b'), descriptions, ignore.case = TRUE, perl = TRUE)
  }

  mapping <- data.frame(
    go_term_id = gsea_results$go_term_id,
    term_group = default_group,
    stringsAsFactors = FALSE)

  for (group_name in names(term_groups)) {
    seeds <- as.character(term_groups[[group_name]])
    matched <- rep(FALSE, nrow(gsea_results))
    for (seed in seeds) {
      id_match <- tolower(gsea_results$go_term_id) == tolower(seed)
      term_match <- match_text_seed(seed, gsea_results$go_description)
      matched <- matched | id_match | term_match
    }
    replace_idx <- mapping$term_group == default_group & matched
    mapping$term_group[replace_idx] <- group_name
  }

  mapping
}

.gsea_assign_term_groups <- function(gsea_results,
                                     term_groups = NULL,
                                     default_group = 'Other') {
  if (is.null(term_groups)) {
    return(.gsea_classify_terms(gsea_results$go_description))
  }
  term_mapping <- .gsea_match_term_groups(
    gsea_results = gsea_results,
    term_groups = term_groups,
    default_group = default_group)
  matched_group <- term_mapping$term_group[match(gsea_results$go_term_id, term_mapping$go_term_id)]
  factor(matched_group, levels = c(names(term_groups), default_group))
}

.gsea_combine_label_requests <- function(label_terms = NULL,
                                         label_ranks = NULL) {
  requests <- list()
  if (!is.null(label_terms) && length(label_terms) > 0L) {
    requests <- c(requests, as.list(label_terms))
  }
  if (!is.null(label_ranks) && length(label_ranks) > 0L) {
    requests <- c(requests, as.list(as.integer(label_ranks)))
  }
  if (length(requests) == 0L) {
    return(NULL)
  }
  requests
}

.gsea_select_labels <- function(plot_tbl,
                                label_terms = NULL,
                                label_n = 0L,
                                score_col = 'label_score',
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
        matched_row <- which(tolower(plot_tbl$go_description) == tolower(requested_label))
        if (length(matched_row) == 0L) {
          matched_row <- which(tolower(plot_tbl$go_term_id) == tolower(requested_label))
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
    plot_tbl <- plot_tbl[order(-plot_tbl[[score_col]], plot_tbl$padj, plot_tbl$go_description), , drop = FALSE]
    return(utils::head(plot_tbl, label_n))
  }

  plot_tbl$label_bin <- cut(
    plot_tbl[[rank_col]],
    breaks = seq(1, max(plot_tbl[[rank_col]]) + 1, length.out = label_n + 1L),
    include.lowest = TRUE,
    labels = FALSE)
  label_tbl <- do.call(rbind, lapply(sort(unique(plot_tbl$label_bin)), function(label_bin) {
    bin_tbl <- plot_tbl[plot_tbl$label_bin == label_bin, , drop = FALSE]
    bin_tbl <- bin_tbl[order(-bin_tbl[[score_col]], bin_tbl$padj, bin_tbl$go_description), , drop = FALSE]
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

  per_direction <- max(1L, floor(label_n / 2L))
  positive_tbl <- plot_tbl[plot_tbl$significant & plot_tbl$NES > 0, , drop = FALSE]
  negative_tbl <- plot_tbl[plot_tbl$significant & plot_tbl$NES < 0, , drop = FALSE]
  positive_tbl <- positive_tbl[order(positive_tbl$NES, positive_tbl$padj), , drop = FALSE]
  negative_tbl <- negative_tbl[order(negative_tbl$NES, negative_tbl$padj), , drop = FALSE]
  positive_tbl$label_rank <- seq_len(nrow(positive_tbl))
  negative_tbl$label_rank <- seq_len(nrow(negative_tbl))

  label_tbl <- rbind(
    .gsea_select_labels(positive_tbl, label_n = per_direction, score_col = 'label_score', rank_col = 'label_rank'),
    .gsea_select_labels(negative_tbl, label_n = per_direction, score_col = 'label_score', rank_col = 'label_rank'))
  if (nrow(label_tbl) < label_n) {
    remaining_tbl <- plot_tbl[plot_tbl$significant & !plot_tbl$go_term_id %in% label_tbl$go_term_id, , drop = FALSE]
    remaining_tbl <- remaining_tbl[order(-remaining_tbl$label_score, remaining_tbl$padj), , drop = FALSE]
    label_tbl <- rbind(label_tbl, utils::head(remaining_tbl, label_n - nrow(label_tbl)))
  }
  label_tbl$label_rank <- NULL
  label_tbl[order(label_tbl$NES), , drop = FALSE]
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
  rank_by <- .gsea_normalize_rank_by(rank_by)
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
  values <- values[is.finite(values)]
  if (length(values) == 0L) {
    return(c(NA_real_, NA_real_))
  }

  value_range <- range(values)
  value_span <- diff(value_range)
  if (value_span == 0) {
    value_span <- max(abs(value_range), 1)
  }
  value_range + c(-1, 1) * value_span * buffer_fraction
}

.gsea_waterfall_y_limits <- function(values,
                                     y_min = NULL,
                                     y_max = NULL,
                                     buffer_fraction = 0.18) {
  y_limits <- .gsea_axis_limits(values, buffer_fraction = buffer_fraction)
  if (!is.null(y_min)) {
    y_limits[[1L]] <- y_min
  }
  if (!is.null(y_max)) {
    y_limits[[2L]] <- y_max
  }
  y_limits
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
