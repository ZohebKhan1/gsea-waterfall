# shared utilities for gsea-waterfall plots

# keep input requirements explicit instead of guessing many column schemas
standardize_gsea_results <- function(gsea_results,
                                     term_col = 'go_description',
                                     nes_col = 'NES',
                                     pvalue_col = 'pval',
                                     padj_col = 'padj',
                                     id_col = NULL,
                                     ontology_col = NULL) {
  required_cols <- base::c(term_col, nes_col, pvalue_col, padj_col)
  missing_cols <- base::setdiff(required_cols, base::names(gsea_results))
  if (base::length(missing_cols) > 0L) {
    base::stop(
      'Missing required GSEA plotting column(s): ',
      base::paste(missing_cols, collapse = ', '),
      call. = FALSE)
  }

  n_terms <- base::nrow(gsea_results)
  go_description <- base::as.character(gsea_results[[term_col]])
  go_term_id <- if (!base::is.null(id_col) && id_col %in% base::names(gsea_results)) {
    base::as.character(gsea_results[[id_col]])
  } else {
    base::make.unique(go_description, sep = '_')
  }
  go_term_class <- if (!base::is.null(ontology_col) && ontology_col %in% base::names(gsea_results)) {
    base::as.character(gsea_results[[ontology_col]])
  } else {
    base::rep(NA_character_, n_terms)
  }

  out <- base::data.frame(
    go_term_id = go_term_id,
    go_term_class = go_term_class,
    go_description = go_description,
    NES = base::as.numeric(gsea_results[[nes_col]]),
    pvalue = base::as.numeric(gsea_results[[pvalue_col]]),
    padj = base::as.numeric(gsea_results[[padj_col]]),
    stringsAsFactors = FALSE)
  out$pval <- out$pvalue

  out[stats::complete.cases(out[, base::c('go_description', 'NES', 'pvalue', 'padj')]), , drop = FALSE]
}

read_gsea_result_csvs <- function(...,
                                  term_col = 'go_description',
                                  nes_col = 'NES',
                                  pvalue_col = 'pval',
                                  padj_col = 'padj',
                                  id_col = NULL,
                                  ontology_col = NULL) {
  paths <- base::unlist(base::list(...), use.names = FALSE)
  paths <- paths[!base::is.na(paths) & base::nzchar(paths)]
  if (base::length(paths) == 0L) {
    base::stop('Provide at least one GSEA CSV path.', call. = FALSE)
  }

  tables <- base::lapply(paths, function(path) {
    standardize_gsea_results(
      gsea_results = utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE),
      term_col = term_col,
      nes_col = nes_col,
      pvalue_col = pvalue_col,
      padj_col = padj_col,
      id_col = id_col,
      ontology_col = ontology_col)
  })
  base::do.call(base::rbind, tables)
}

safe_neg_log10 <- function(x) {
  -base::log10(base::pmax(x, .Machine$double.xmin))
}

gsea_neutral_gray <- function() {
  '#9E9E9E'
}

gsea_light_gray <- function() {
  '#D9D9D9'
}

gsea_default_colors <- function() {
  base::c(
    'Significantly up' = '#CC79A7',
    'Significantly down' = '#0072B2',
    'Not significant' = gsea_light_gray(),
    'Significant in both' = '#2A9D8F',
    'Significant in x only' = '#CC79A7',
    'Significant in y only' = '#0072B2',
    'Not significant in either' = gsea_light_gray(),
    'Other' = gsea_neutral_gray())
}

gsea_category_colors <- function() {
  base::c(
    'Development/morphogenesis' = '#1B9E77',
    'Neuronal/signaling' = '#4F61BD',
    'Cell cycle/genome' = '#D73027',
    'Metabolism/translation' = '#7A5195',
    'Other' = gsea_neutral_gray())
}

validate_max_categories <- function(term_groups = NULL,
                                    max_categories = 6L) {
  if (base::is.null(term_groups)) {
    return(base::invisible(NULL))
  }
  max_categories <- base::as.integer(max_categories)
  if (base::is.na(max_categories) || max_categories < 1L) {
    base::stop('`max_categories` must be a positive integer.', call. = FALSE)
  }
  if (base::length(term_groups) > max_categories) {
    base::stop('`term_groups` supports up to ', max_categories, ' named categories.', call. = FALSE)
  }
  base::invisible(NULL)
}

validate_term_groups <- function(term_groups = NULL,
                                 category_n = NULL) {
  if (base::is.null(term_groups)) {
    return(base::invisible(NULL))
  }
  if (base::is.null(base::names(term_groups)) || base::any(!base::nzchar(base::names(term_groups)))) {
    base::stop('`term_groups` must be a named list.', call. = FALSE)
  }
  if (!base::is.null(category_n) && category_n != base::length(term_groups)) {
    base::stop('`category_n` must equal the number of groups in `term_groups`.', call. = FALSE)
  }
  base::invisible(NULL)
}

resolve_named_colors <- function(color_values,
                                 required_names,
                                 parameter_name) {
  if (base::is.null(color_values)) {
    return(NULL)
  }
  missing_colors <- base::setdiff(required_names, base::names(color_values))
  if (base::length(missing_colors) > 0L) {
    base::stop(
      '`',
      parameter_name,
      '` is missing color(s) for: ',
      base::paste(missing_colors, collapse = ', '),
      call. = FALSE)
  }
  color_values[required_names]
}

resolve_term_group_colors <- function(term_groups = NULL,
                                      term_group_colors = NULL,
                                      default_group = 'Other') {
  if (base::is.null(term_groups)) {
    return(gsea_category_colors())
  }
  group_names <- base::names(term_groups)
  if (base::is.null(term_group_colors)) {
    if (!base::requireNamespace('RColorBrewer', quietly = TRUE)) {
      base::stop('The RColorBrewer package is required for default term-group colors.', call. = FALSE)
    }
    brewer_n <- base::max(3L, base::length(group_names))
    term_group_colors <- stats::setNames(
      RColorBrewer::brewer.pal(brewer_n, 'Set1')[base::seq_along(group_names)],
      group_names)
  }

  missing_colors <- base::setdiff(group_names, base::names(term_group_colors))
  if (base::length(missing_colors) > 0L) {
    base::stop(
      '`term_group_colors` is missing color(s) for: ',
      base::paste(missing_colors, collapse = ', '),
      call. = FALSE)
  }
  base::c(term_group_colors[group_names], stats::setNames(gsea_neutral_gray(), default_group))
}

wrap_gsea_label <- function(labels,
                            words_per_line = 2L) {
  words_per_line <- base::as.integer(words_per_line)
  if (base::is.na(words_per_line) || words_per_line < 1L) {
    base::stop('`words_per_line` must be a positive integer.', call. = FALSE)
  }

  base::vapply(base::as.character(labels), function(label) {
    words <- base::strsplit(base::trimws(label), '\\s+', perl = TRUE)[[1L]]
    if (base::length(words) <= words_per_line) {
      return(base::paste(words, collapse = ' '))
    }
    word_groups <- base::split(words, base::ceiling(base::seq_along(words) / words_per_line))
    base::paste(
      base::vapply(word_groups, base::paste, character(1L), collapse = ' '),
      collapse = '\n')
  }, character(1L), USE.NAMES = FALSE)
}

split_point_layers <- function(plot_tbl,
                               group_col,
                               background_groups) {
  group_values <- base::as.character(plot_tbl[[group_col]])
  background_rows <- base::is.na(group_values) | group_values %in% background_groups
  base::list(
    background = plot_tbl[background_rows, , drop = FALSE],
    foreground = plot_tbl[!background_rows, , drop = FALSE])
}

# fallback categories keep ungrouped scatterplots interpretable
classify_go_terms <- function(go_description) {
  term <- base::tolower(base::as.character(go_description))
  category <- base::rep('Other', base::length(term))

  category[base::grepl(
    'cell cycle|mitotic|meiotic|chromosom|spindle|centromer|kinetochore|dna replication|dna repair|recombination',
    term)] <- 'Cell cycle/genome'
  category[category == 'Other' & base::grepl(
    'development|morphogenesis|differentiation|mesoderm|mesenchyme|embryo|organogenesis|migration|pattern specification',
    term)] <- 'Development/morphogenesis'
  category[category == 'Other' & base::grepl(
    'synap|neuro|axon|dendrit|transmitter|ion channel|receptor|signaling|signal transduction|membrane potential|secretion',
    term)] <- 'Neuronal/signaling'
  category[category == 'Other' & base::grepl(
    'ribosom|translation|rrna|metabolic|biosynthetic|mitochond|oxidative|protein folding|rna processing',
    term)] <- 'Metabolism/translation'

  base::factor(category, levels = base::names(gsea_category_colors()))
}

match_go_term_groups <- function(gsea_results,
                                 term_groups,
                                 default_group = 'Other') {
  validate_term_groups(term_groups)

  match_text_seed <- function(seed, descriptions) {
    escaped_seed <- base::gsub('([][{}()+*^$|\\\\?.])', '\\\\\\1', seed)
    base::grepl(base::paste0('\\b', escaped_seed, '\\b'), descriptions, ignore.case = TRUE, perl = TRUE)
  }

  mapping <- base::data.frame(
    go_term_id = gsea_results$go_term_id,
    term_group = default_group,
    stringsAsFactors = FALSE)

  for (group_name in base::names(term_groups)) {
    seeds <- base::as.character(term_groups[[group_name]])
    matched <- base::rep(FALSE, base::nrow(gsea_results))
    for (seed in seeds) {
      id_match <- base::tolower(gsea_results$go_term_id) == base::tolower(seed)
      term_match <- match_text_seed(seed, gsea_results$go_description)
      matched <- matched | id_match | term_match
    }
    replace_idx <- mapping$term_group == default_group & matched
    mapping$term_group[replace_idx] <- group_name
  }

  mapping
}

assign_plot_term_groups <- function(gsea_results,
                                    term_groups = NULL,
                                    default_group = 'Other') {
  if (base::is.null(term_groups)) {
    return(classify_go_terms(gsea_results$go_description))
  }
  term_mapping <- match_go_term_groups(
    gsea_results = gsea_results,
    term_groups = term_groups,
    default_group = default_group)
  matched_group <- term_mapping$term_group[base::match(gsea_results$go_term_id, term_mapping$go_term_id)]
  base::factor(matched_group, levels = base::c(base::names(term_groups), default_group))
}

combine_label_requests <- function(label_terms = NULL,
                                   label_ranks = NULL) {
  requests <- base::list()
  if (!base::is.null(label_terms) && base::length(label_terms) > 0L) {
    requests <- base::c(requests, base::as.list(label_terms))
  }
  if (!base::is.null(label_ranks) && base::length(label_ranks) > 0L) {
    requests <- base::c(requests, base::as.list(base::as.integer(label_ranks)))
  }
  if (base::length(requests) == 0L) {
    return(NULL)
  }
  requests
}

select_gsea_labels <- function(plot_tbl,
                               label_terms = NULL,
                               label_n = 0L,
                               score_col = 'label_score',
                               rank_col = NULL) {
  if (base::nrow(plot_tbl) == 0L) {
    return(plot_tbl[0L, , drop = FALSE])
  }

  if (!base::is.null(label_terms) && base::length(label_terms) > 0L) {
    label_rows <- base::integer()
    for (requested_label in label_terms) {
      if (base::is.numeric(requested_label)) {
        if (base::is.null(rank_col) || !rank_col %in% base::names(plot_tbl)) {
          base::warning('Rank label requested, but no rank column is available; skipping: ', requested_label, call. = FALSE)
          next
        }
        matched_row <- base::which(plot_tbl[[rank_col]] == requested_label)
      } else {
        requested_label <- base::as.character(requested_label)
        matched_row <- base::which(base::tolower(plot_tbl$go_description) == base::tolower(requested_label))
        if (base::length(matched_row) == 0L) {
          matched_row <- base::which(base::tolower(plot_tbl$go_term_id) == base::tolower(requested_label))
        }
      }
      if (base::length(matched_row) == 0L) {
        base::warning('Requested GO label was not found in the plotted terms; skipping: ', requested_label, call. = FALSE)
        next
      }
      label_rows <- base::c(label_rows, matched_row[[1L]])
    }
    label_tbl <- plot_tbl[base::unique(label_rows), , drop = FALSE]
    return(label_tbl)
  }

  if (label_n <= 0L) {
    return(plot_tbl[0L, , drop = FALSE])
  }
  if (base::is.null(rank_col) || !rank_col %in% base::names(plot_tbl)) {
    plot_tbl <- plot_tbl[base::order(-plot_tbl[[score_col]], plot_tbl$padj, plot_tbl$go_description), , drop = FALSE]
    return(utils::head(plot_tbl, label_n))
  }

  plot_tbl$label_bin <- base::cut(
    plot_tbl[[rank_col]],
    breaks = base::seq(1, base::max(plot_tbl[[rank_col]]) + 1, length.out = label_n + 1L),
    include.lowest = TRUE,
    labels = FALSE)
  label_tbl <- do.call(base::rbind, base::lapply(base::sort(base::unique(plot_tbl$label_bin)), function(label_bin) {
    bin_tbl <- plot_tbl[plot_tbl$label_bin == label_bin, , drop = FALSE]
    bin_tbl <- bin_tbl[base::order(-bin_tbl[[score_col]], bin_tbl$padj, bin_tbl$go_description), , drop = FALSE]
    utils::head(bin_tbl, 1L)
  }))
  label_tbl$label_bin <- NULL
  label_tbl[base::order(label_tbl[[rank_col]]), , drop = FALSE]
}

select_volcano_labels <- function(plot_tbl,
                                  label_n = 12L) {
  if (label_n <= 0L || base::nrow(plot_tbl) == 0L) {
    return(plot_tbl[0L, , drop = FALSE])
  }

  per_direction <- base::max(1L, base::floor(label_n / 2L))
  positive_tbl <- plot_tbl[plot_tbl$significant & plot_tbl$NES > 0, , drop = FALSE]
  negative_tbl <- plot_tbl[plot_tbl$significant & plot_tbl$NES < 0, , drop = FALSE]
  positive_tbl <- positive_tbl[base::order(positive_tbl$NES, positive_tbl$padj), , drop = FALSE]
  negative_tbl <- negative_tbl[base::order(negative_tbl$NES, negative_tbl$padj), , drop = FALSE]
  positive_tbl$label_rank <- base::seq_len(base::nrow(positive_tbl))
  negative_tbl$label_rank <- base::seq_len(base::nrow(negative_tbl))

  label_tbl <- base::rbind(
    select_gsea_labels(positive_tbl, label_n = per_direction, score_col = 'label_score', rank_col = 'label_rank'),
    select_gsea_labels(negative_tbl, label_n = per_direction, score_col = 'label_score', rank_col = 'label_rank'))
  if (base::nrow(label_tbl) < label_n) {
    remaining_tbl <- plot_tbl[plot_tbl$significant & !plot_tbl$go_term_id %in% label_tbl$go_term_id, , drop = FALSE]
    remaining_tbl <- remaining_tbl[base::order(-remaining_tbl$label_score, remaining_tbl$padj), , drop = FALSE]
    label_tbl <- base::rbind(label_tbl, utils::head(remaining_tbl, label_n - base::nrow(label_tbl)))
  }
  label_tbl$label_rank <- NULL
  label_tbl[base::order(label_tbl$NES), , drop = FALSE]
}

normalize_rank_by <- function(rank_by) {
  rank_by <- base::match.arg(rank_by, base::c('NES', 'padj', 'pvalue', 'pval'))
  if (rank_by == 'pval') {
    return('pvalue')
  }
  rank_by
}

rank_gsea_direction_terms <- function(plot_tbl,
                                      direction,
                                      rank_by = 'NES') {
  rank_by <- normalize_rank_by(rank_by)
  if (rank_by == 'NES') {
    if (direction == 'positive') {
      return(plot_tbl[base::order(-plot_tbl$NES, plot_tbl$padj, plot_tbl$go_description), , drop = FALSE])
    }
    return(plot_tbl[base::order(plot_tbl$NES, plot_tbl$padj, plot_tbl$go_description), , drop = FALSE])
  }

  rank_values <- plot_tbl[[rank_by]]
  nes_tiebreak <- if (direction == 'positive') -plot_tbl$NES else plot_tbl$NES
  plot_tbl[base::order(rank_values, nes_tiebreak, plot_tbl$go_description, na.last = TRUE), , drop = FALSE]
}

calculate_axis_limits <- function(values,
                                  buffer_fraction = 0.045) {
  values <- values[base::is.finite(values)]
  if (base::length(values) == 0L) {
    return(base::c(NA_real_, NA_real_))
  }

  value_range <- base::range(values)
  value_span <- base::diff(value_range)
  if (value_span == 0) {
    value_span <- base::max(base::abs(value_range), 1)
  }
  value_range + base::c(-1, 1) * value_span * buffer_fraction
}

calculate_waterfall_y_limits <- function(values,
                                         y_min = NULL,
                                         y_max = NULL,
                                         buffer_fraction = 0.18) {
  y_limits <- calculate_axis_limits(values, buffer_fraction = buffer_fraction)
  if (!base::is.null(y_min)) {
    y_limits[[1L]] <- y_min
  }
  if (!base::is.null(y_max)) {
    y_limits[[2L]] <- y_max
  }
  y_limits
}
