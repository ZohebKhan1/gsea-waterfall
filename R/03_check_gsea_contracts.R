#!/usr/bin/env Rscript

# ----
# author:
# - Zoheb Khan
#
# script path:
# - R/03_check_gsea_contracts.R
#
# functions:
# - functions/gsea_plot_utils.R
# - functions/gsea_waterfall_plot.R
# - functions/gsea_volcano_plot.R
# - functions/gsea_half_volcano_plot.R
# - functions/gsea_nes_scatter_plot.R
#
# input data:
# - tutorial/data/GSE122380_gsea_cardiomyocyte_vs_mesoderm.csv
# - tutorial/data/GSE122380_gsea_day9_vs_day6.csv
# - tutorial/data/GSE122380_gsea_day3_vs_day1.csv
# ----

# 0.0 verify project root and source functions -----------------

if (!file.exists('gsea.Rproj')) {
  stop('Run R/03_check_gsea_contracts.R from the gsea repository root.', call. = FALSE)
}

# source GSEA plotting functions
source('functions/gsea_plot_utils.R')
source('functions/gsea_waterfall_plot.R')
source('functions/gsea_volcano_plot.R')
source('functions/gsea_half_volcano_plot.R')
source('functions/gsea_nes_scatter_plot.R')

# 1.0 verify fixed input and scientific contracts -----------------

gsea_cardiomyocyte_vs_mesoderm <- read_gsea_result_csv(
  'tutorial/data/GSE122380_gsea_cardiomyocyte_vs_mesoderm.csv')
gsea_day9_vs_day6 <- read_gsea_result_csv(
  'tutorial/data/GSE122380_gsea_day9_vs_day6.csv')
gsea_day3_vs_day1 <- read_gsea_result_csv(
  'tutorial/data/GSE122380_gsea_day3_vs_day1.csv')

stopifnot(
  nrow(gsea_cardiomyocyte_vs_mesoderm) == 2482L,
  sum(gsea_cardiomyocyte_vs_mesoderm$NES < 0) == 744L,
  nrow(gsea_day9_vs_day6) == 2477L,
  sum(gsea_day9_vs_day6$NES < 0) == 828L,
  nrow(gsea_day3_vs_day1) == 2477L,
  sum(gsea_day3_vs_day1$NES < 0) == 896L)

waterfall_plot <- plot_gsea_waterfall(
  gsea_results = gsea_cardiomyocyte_vs_mesoderm,
  direction = 'positive',
  top_n = 100L)
volcano_plot <- plot_gsea_volcano(
  gsea_results = gsea_cardiomyocyte_vs_mesoderm)
nes_scatter_plot <- plot_gsea_nes_scatter(
  gsea_x = gsea_day9_vs_day6,
  gsea_y = gsea_day3_vs_day1,
  x_label = 'Day 9 vs. Day 6 NES',
  y_label = 'Day 3 vs. Day 1 NES',
  x_name = 'Day 9 vs. Day 6',
  y_name = 'Day 3 vs. Day 1',
  include_nonsignificant = TRUE,
  equal_axis_limits = TRUE)

stopifnot(
  nrow(waterfall_plot$data) == 100L,
  sum(volcano_plot$data$point_group == 'Significantly down') == 272L,
  sum(volcano_plot$data$point_group == 'Significantly up') == 445L)

expected_nes_scatter_counts <- c(
  'Not significant in either' = 1327L,
  'Significant in Day 3 vs. Day 1 only' = 601L,
  'Significant in both' = 280L,
  'Significant in Day 9 vs. Day 6 only' = 269L)
nes_scatter_counts <- table(nes_scatter_plot$data$point_group)
stopifnot(all(
  as.integer(nes_scatter_counts[names(expected_nes_scatter_counts)]) ==
    expected_nes_scatter_counts))

# 2.0 verify public input and identifier boundaries -----------------

.expect_error <- function(expression) {
  inherits(try(force(expression), silent = TRUE), 'try-error')
}

invalid_duplicate <- gsea_cardiomyocyte_vs_mesoderm[1:2, ]
invalid_duplicate$go_term_id[[2L]] <- invalid_duplicate$go_term_id[[1L]]
stopifnot(.expect_error(plot_gsea_waterfall(
  invalid_duplicate,
  id_col = 'go_term_id')))

invalid_probability <- gsea_cardiomyocyte_vs_mesoderm[1:2, ]
invalid_probability$padj[[1L]] <- 1.1
stopifnot(.expect_error(plot_gsea_waterfall(invalid_probability)))

invalid_numeric <- gsea_cardiomyocyte_vs_mesoderm[1:2, ]
invalid_numeric$NES <- as.character(invalid_numeric$NES)
stopifnot(.expect_error(plot_gsea_waterfall(invalid_numeric)))

invalid_description <- gsea_cardiomyocyte_vs_mesoderm[1:2, ]
invalid_description$go_description[[1L]] <- ' '
stopifnot(.expect_error(plot_gsea_waterfall(invalid_description)))

ambiguous_description <- gsea_cardiomyocyte_vs_mesoderm[1:2, ]
ambiguous_description$go_description[[2L]] <- ambiguous_description$go_description[[1L]]
stopifnot(.expect_error(plot_gsea_waterfall(
  ambiguous_description,
  label_terms = ambiguous_description$go_description[[1L]])))

stopifnot(
  .expect_error(plot_gsea_waterfall(
    gsea_cardiomyocyte_vs_mesoderm,
    top_n = 0L)),
  .expect_error(plot_gsea_volcano(
    gsea_cardiomyocyte_vs_mesoderm,
    padj_cutoff = 1.1)),
  .expect_error(plot_gsea_waterfall(
    gsea_cardiomyocyte_vs_mesoderm,
    term_groups = list('duplicate' = 'one', 'duplicate' = 'two'))),
  .expect_error(plot_gsea_waterfall(
    gsea_cardiomyocyte_vs_mesoderm,
    term_groups = list('Other' = 'development'))),
  .expect_error(plot_gsea_nes_scatter(
    gsea_x = gsea_day9_vs_day6,
    gsea_y = gsea_day3_vs_day1,
    x_label = 'x',
    y_label = 'y',
    include_nonsignificant = NA)))

x_results <- data.frame(
  term = c('A', 'B'),
  id = c('one', 'two'),
  NES = c(1, 2),
  p = c(0.1, 0.2),
  q = c(0.2, 0.3))
y_results <- data.frame(
  term = c('renamed A', 'renamed B'),
  id = c('one', 'two'),
  NES = c(1.5, 2.5),
  p = c(0.1, 0.2),
  q = c(0.2, 0.3))

stable_id_nes_scatter_plot <- plot_gsea_nes_scatter(
  x_results,
  y_results,
  x_label = 'x',
  y_label = 'y',
  term_col = 'term',
  nes_col = 'NES',
  pvalue_col = 'p',
  padj_col = 'q',
  id_col = 'id',
  include_nonsignificant = TRUE)
stopifnot(nrow(stable_id_nes_scatter_plot$data) == 2L)

y_results$id <- c('three', 'four')
stopifnot(.expect_error(plot_gsea_nes_scatter(
  x_results,
  y_results,
  x_label = 'x',
  y_label = 'y',
  term_col = 'term',
  nes_col = 'NES',
  pvalue_col = 'p',
  padj_col = 'q',
  id_col = 'id')))

message('GSEA input, scientific, plotting, identifier, and failure-boundary contracts passed.')
