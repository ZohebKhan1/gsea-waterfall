# standalone helpers for gsea plotting functions
if (!base::exists('standardize_gsea_results', mode = 'function')) {
  # comparative ranked gsea plotting utilities

  # plotting requires term names, nes, nominal p-values, and adjusted p-values
  standardize_gsea_results <- function(gsea_results) {
    resolve_col <- function(candidates) {
      found <- candidates[candidates %in% base::names(gsea_results)]
      if (base::length(found) == 0L) {
        return(NA_character_)
      }
      found[[1L]]
    }

    col_map <- base::list(
      go_term_id = resolve_col(base::c('go_term_id', 'term_id', 'pathway', 'ID')),
      go_term_class = resolve_col(base::c('go_term_class', 'gene_set_collection', 'go_ontology', 'ontology')),
      go_description = resolve_col(base::c('go_description', 'term_name', 'description', 'pathway_name', 'Description')),
      genes_in_go = resolve_col(base::c('genes_in_go', 'go_term_size', 'size')),
      genes_scored_in_go = resolve_col(base::c('genes_scored_in_go', 'total_ranked_genes_in_go')),
      leading_edge_in_go = resolve_col(base::c('leading_edge_in_go', 'leading_edge_gene_count')),
      ES = resolve_col(base::c('ES', 'enrichment_score')),
      NES = resolve_col(base::c('NES', 'normalized_enrichment_score')),
      pvalue = resolve_col(base::c('pvalue', 'pval', 'p_value', 'p.value', 'fgsea_nominal_p_value')),
      padj = resolve_col(base::c('padj', 'qvalue', 'q_value', 'FDR', 'fdr', 'gsea_adjusted_p_value')),
      leading_edge_genes = resolve_col(base::c('leading_edge_genes', 'leading_edge_gene_symbols', 'leadingEdge')))
    required_cols <- base::c('go_description', 'NES', 'pvalue', 'padj')
    missing_cols <- required_cols[base::is.na(base::unlist(col_map[required_cols]))]
    if (base::length(missing_cols) > 0L) {
      base::stop(
        'Missing required GSEA plotting column(s): ',
        base::paste(missing_cols, collapse = ', '),
        call. = FALSE)
    }

    go_description <- base::as.character(gsea_results[[col_map$go_description]])
    go_term_id <- if (base::is.na(col_map$go_term_id)) {
      base::make.unique(go_description, sep = '_')
    } else {
      base::as.character(gsea_results[[col_map$go_term_id]])
    }
    go_term_class <- if (base::is.na(col_map$go_term_class)) {
      base::rep(NA_character_, base::length(go_description))
    } else {
      base::as.character(gsea_results[[col_map$go_term_class]])
    }

    optional_numeric <- function(col_name, coerce_fn) {
      if (base::is.na(col_name)) {
        return(base::rep(NA_real_, base::length(go_description)))
      }
      coerce_fn(gsea_results[[col_name]])
    }
    optional_character <- function(col_name) {
      if (base::is.na(col_name)) {
        return(base::rep(NA_character_, base::length(go_description)))
      }
      base::as.character(gsea_results[[col_name]])
    }

    out <- base::data.frame(
      go_term_id = go_term_id,
      go_term_class = go_term_class,
      go_description = go_description,
      genes_in_go = optional_numeric(col_map$genes_in_go, base::as.integer),
      genes_scored_in_go = optional_numeric(col_map$genes_scored_in_go, base::as.integer),
      leading_edge_in_go = optional_numeric(col_map$leading_edge_in_go, base::as.integer),
      ES = optional_numeric(col_map$ES, base::as.numeric),
      NES = base::as.numeric(gsea_results[[col_map$NES]]),
      pvalue = base::as.numeric(gsea_results[[col_map$pvalue]]),
      padj = base::as.numeric(gsea_results[[col_map$padj]]),
      leading_edge_genes = optional_character(col_map$leading_edge_genes),
      stringsAsFactors = FALSE)

    optional_cols <- base::intersect(
      base::c('contrast_id', 'deg_contrast_full_name', 'comparison_label'),
      base::names(gsea_results))
    if (base::length(optional_cols) > 0L) {
      out$contrast_id <- base::as.character(gsea_results[[optional_cols[[1L]]]])
    }

    out[stats::complete.cases(out[, base::c('go_description', 'NES', 'pvalue', 'padj')]), , drop = FALSE]
  }

  read_gsea_result_csvs <- function(...) {
    paths <- base::unlist(base::list(...), use.names = FALSE)
    paths <- paths[!base::is.na(paths) & base::nzchar(paths)]
    if (base::length(paths) == 0L) {
      base::stop('Provide at least one GSEA CSV path.', call. = FALSE)
    }

    tables <- base::lapply(paths, function(path) {
      standardize_gsea_results(utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE))
    })
    base::do.call(base::rbind, tables)
  }

  safe_neg_log10 <- function(x) {
    -base::log10(base::pmax(x, .Machine$double.xmin))
  }

  resolve_font_path <- function(candidates) {
    found <- candidates[base::file.exists(candidates)]
    if (base::length(found) == 0L) {
      return(NA_character_)
    }
    found[[1L]]
  }

  plot_font_paths <- function(font_family = 'Nimbus Sans') {
    if (font_family == 'Nimbus Sans') {
      font_dir <- 'tutorial/assets/fonts'
      return(base::c(
        regular = resolve_font_path(base::c(
          base::file.path(font_dir, 'NimbusSans-Regular.otf'),
          '/usr/share/fonts/opentype/urw-base35/NimbusSans-Regular.otf')),
        bold = resolve_font_path(base::c(
          base::file.path(font_dir, 'NimbusSans-Bold.otf'),
          '/usr/share/fonts/opentype/urw-base35/NimbusSans-Bold.otf')),
        italic = resolve_font_path(base::c(
          base::file.path(font_dir, 'NimbusSans-Italic.otf'),
          '/usr/share/fonts/opentype/urw-base35/NimbusSans-Italic.otf')),
        bolditalic = resolve_font_path(base::c(
          base::file.path(font_dir, 'NimbusSans-BoldItalic.otf'),
          '/usr/share/fonts/opentype/urw-base35/NimbusSans-BoldItalic.otf'))))
    }
    if (font_family == 'Latin Modern Sans') {
      return(base::c(
        regular = base::path.expand('~/.local/share/fonts/latin-modern-sans/LatinModernSans-Regular.otf'),
        bold = base::path.expand('~/.local/share/fonts/latin-modern-sans/LatinModernSans-Bold.otf'),
        italic = base::path.expand('~/.local/share/fonts/latin-modern-sans/LatinModernSans-Oblique.otf'),
        bolditalic = base::path.expand('~/.local/share/fonts/latin-modern-sans/LatinModernSans-BoldOblique.otf')))
    }
    base::c(regular = NA_character_, bold = NA_character_, italic = NA_character_, bolditalic = NA_character_)
  }

  register_plot_font <- function(font_family = 'Nimbus Sans') {
    font_paths <- plot_font_paths(font_family = font_family)
    if (base::is.na(font_paths[['regular']]) || !base::file.exists(font_paths[['regular']])) {
      return(base::invisible(FALSE))
    }

    registered <- tryCatch(
      systemfonts::register_font(
        name = font_family,
        plain = font_paths[['regular']],
        bold = if (base::file.exists(font_paths[['bold']])) font_paths[['bold']] else font_paths[['regular']],
        italic = if (base::file.exists(font_paths[['italic']])) font_paths[['italic']] else font_paths[['regular']],
        bolditalic = if (base::file.exists(font_paths[['bolditalic']])) {
          font_paths[['bolditalic']]
        } else {
          font_paths[['regular']]
        }),
      error = function(e) FALSE)
    base::invisible(registered)
  }

  plot_user_fonts <- function(font_family = 'Nimbus Sans') {
    font_paths <- plot_font_paths(font_family = font_family)
    if (base::is.na(font_paths[['regular']]) || !base::file.exists(font_paths[['regular']])) {
      return(base::list())
    }

    stats::setNames(base::list(font_paths[['regular']]), font_family)
  }

  nimbus_svg_web_fonts <- function() {
    base::list(
      svglite::font_face(
        family = 'Nimbus Sans',
        otf = '../fonts/NimbusSans-Regular.otf',
        weight = 400,
        style = 'normal'),
      svglite::font_face(
        family = 'Nimbus Sans',
        otf = '../fonts/NimbusSans-Bold.otf',
        weight = 700,
        style = 'normal'),
      svglite::font_face(
        family = 'Nimbus Sans',
        otf = '../fonts/NimbusSans-Italic.otf',
        weight = 400,
        style = 'italic'),
      svglite::font_face(
        family = 'Nimbus Sans',
        otf = '../fonts/NimbusSans-BoldItalic.otf',
        weight = 700,
        style = 'italic'))
  }

  font_data_uri <- function(path) {
    if (!base::requireNamespace('base64enc', quietly = TRUE)) {
      base::stop('The base64enc package is required to embed SVG fonts.', call. = FALSE)
    }
    if (!base::file.exists(path)) {
      base::stop('Font file not found: ', path, call. = FALSE)
    }
    extension <- base::tolower(tools::file_ext(path))
    mime_type <- if (extension == 'ttf') 'font/ttf' else 'font/otf'
    base::paste0('data:', mime_type, ';base64,', base64enc::base64encode(path))
  }

  embed_svg_font <- function(path, svg_font_url, font_file, font_name) {
    svg_text <- base::readLines(path, warn = FALSE)
    font_format <- if (base::tolower(tools::file_ext(font_file)) == 'ttf') 'truetype' else 'opentype'
    font_uri <- font_data_uri(font_file)
    svg_text <- base::gsub(
      pattern = base::paste0('url\\("', svg_font_url, '"\\) format\\("', font_format, '"\\)'),
      replacement = base::paste0('url("', font_uri, '") format("', font_format, '")'),
      x = svg_text,
      fixed = FALSE)
    if (!base::any(base::grepl(font_uri, svg_text, fixed = TRUE))) {
      base::stop('failed to embed ', font_name, ' in ', path, call. = FALSE)
    }
    base::writeLines(svg_text, path, useBytes = TRUE)
  }

  embed_nimbus_svg_fonts <- function(path) {
    font_dir <- base::paste0(base::dirname(base::dirname(path)), '/fonts')
    embed_svg_font(
      path,
      '../fonts/NimbusSans-Regular.otf',
      base::paste0(font_dir, '/NimbusSans-Regular.otf'),
      'Nimbus Sans Regular')
    embed_svg_font(
      path,
      '../fonts/NimbusSans-Bold.otf',
      base::paste0(font_dir, '/NimbusSans-Bold.otf'),
      'Nimbus Sans Bold')
    embed_svg_font(
      path,
      '../fonts/NimbusSans-Italic.otf',
      base::paste0(font_dir, '/NimbusSans-Italic.otf'),
      'Nimbus Sans Italic')
    embed_svg_font(
      path,
      '../fonts/NimbusSans-BoldItalic.otf',
      base::paste0(font_dir, '/NimbusSans-BoldItalic.otf'),
      'Nimbus Sans Bold Italic')
    base::invisible(path)
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
      base::stop(
        '`term_groups` supports up to ',
        max_categories,
        ' named categories.',
        call. = FALSE)
    }
    base::invisible(NULL)
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

  normalize_rank_by <- function(rank_by) {
    rank_by <- base::match.arg(rank_by, base::c('NES', 'padj', 'pvalue', 'pval'))
    if (rank_by == 'pval') {
      rank_by <- 'pvalue'
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
      base::paste(base::vapply(word_groups, base::paste, character(1L), collapse = ' '), collapse = '\n')
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

  get_go_offspring <- function(go_ids,
                               ontology) {
    go_ids <- base::unique(base::as.character(go_ids))
    go_ids <- go_ids[base::grepl('^GO:[0-9]+$', go_ids)]
    ontology <- base::toupper(stats::na.omit(base::unique(ontology)))
    if (base::length(go_ids) == 0L || base::length(ontology) != 1L || !ontology %in% base::c('BP', 'CC', 'MF')) {
      return(go_ids)
    }
    if (!base::requireNamespace('GO.db', quietly = TRUE)) {
      base::warning('GO.db is not installed; GO descendant matching was skipped.', call. = FALSE)
      return(go_ids)
    }

    offspring_map <- base::switch(
      ontology,
      BP = base::as.list(GO.db::GOBPOFFSPRING),
      CC = base::as.list(GO.db::GOCCOFFSPRING),
      MF = base::as.list(GO.db::GOMFOFFSPRING))
    base::unique(base::c(go_ids, base::unlist(offspring_map[go_ids], use.names = FALSE)))
  }

  match_go_term_groups <- function(gsea_results,
                                   term_groups,
                                   include_descendants = FALSE,
                                   default_group = 'Other') {
    gsea_results <- standardize_gsea_results(gsea_results)
    if (base::is.null(base::names(term_groups)) || base::any(!base::nzchar(base::names(term_groups)))) {
      base::stop('`term_groups` must be a named list.', call. = FALSE)
    }

    match_text_seed <- function(seed, descriptions) {
      escaped_seed <- base::gsub('([][{}()+*^$|\\\\?.])', '\\\\\\1', seed)
      base::grepl(
        base::paste0('\\b', escaped_seed, '\\b'),
        descriptions,
        ignore.case = TRUE,
        perl = TRUE)
    }

    mapping <- base::data.frame(go_term_id = gsea_results$go_term_id, term_group = default_group, stringsAsFactors = FALSE)
    direct_matches <- stats::setNames(base::vector('list', base::length(term_groups)), base::names(term_groups))
    for (group_name in base::names(term_groups)) {
      seeds <- base::as.character(term_groups[[group_name]])
      explicit_go_ids <- seeds[base::grepl('^GO:[0-9]+$', seeds)]
      text_seeds <- seeds[!base::grepl('^GO:[0-9]+$', seeds)]
      text_go_ids <- base::character()

      for (seed in text_seeds) {
        seed_match <- match_text_seed(seed, gsea_results$go_description)
        text_go_ids <- base::unique(base::c(text_go_ids, gsea_results$go_term_id[seed_match]))
      }

      direct_matches[[group_name]] <- base::unique(base::c(explicit_go_ids, text_go_ids))
      replace_idx <- mapping$term_group == default_group & mapping$go_term_id %in% direct_matches[[group_name]]
      mapping$term_group[replace_idx] <- group_name
    }

    if (!base::isTRUE(include_descendants)) {
      return(mapping)
    }

    for (group_name in base::names(term_groups)) {
      matched_ids <- direct_matches[[group_name]]
      if (base::length(matched_ids) > 0L) {
        expanded_ids <- base::character()
        for (ontology in base::intersect(base::unique(gsea_results$go_term_class), base::c('BP', 'CC', 'MF'))) {
          ontology_ids <- base::intersect(matched_ids, gsea_results$go_term_id[gsea_results$go_term_class == ontology])
          expanded_ids <- base::unique(base::c(expanded_ids, get_go_offspring(ontology_ids, ontology = ontology)))
        }
        matched_ids <- base::unique(base::c(matched_ids, expanded_ids))
      }

      replace_idx <- mapping$term_group == default_group & mapping$go_term_id %in% matched_ids
      mapping$term_group[replace_idx] <- group_name
    }

    mapping
  }

  assign_plot_term_groups <- function(gsea_results,
                                      term_groups = NULL,
                                      include_descendants = FALSE,
                                      default_group = 'Other') {
    if (base::is.null(term_groups)) {
      return(classify_go_terms(gsea_results$go_description))
    }
    term_mapping <- match_go_term_groups(
      gsea_results = gsea_results,
      term_groups = term_groups,
      include_descendants = include_descendants,
      default_group = default_group)
    matched_group <- term_mapping$term_group[base::match(gsea_results$go_term_id, term_mapping$go_term_id)]
    base::factor(matched_group, levels = base::c(base::names(term_groups), default_group))
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
      for (label_idx in base::seq_along(label_terms)) {
        requested_label <- label_terms[[label_idx]]
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
      label_tbl$label_order <- base::seq_len(base::nrow(label_tbl))
      return(label_tbl[base::order(label_tbl$label_order), , drop = FALSE])
    }

    if (label_n > 0L) {
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
      label_tbl <- label_tbl[base::order(label_tbl[[rank_col]]), , drop = FALSE]
      label_tbl$label_bin <- NULL
      return(utils::head(label_tbl, label_n))
    }

    plot_tbl[0L, , drop = FALSE]
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
}

if (!base::exists('gsea_light_gray', mode = 'function')) {
  gsea_light_gray <- function() {
    '#D9D9D9'
  }
}


#' Make a symmetric GSEA volcano plot
#'
#' @param gsea_results Standardized or legacy GSEA result table.
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
  register_plot_font(font_family = font_family)
  plot_tbl <- standardize_gsea_results(gsea_results)
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
