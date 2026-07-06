#!/usr/bin/env Rscript

# Created:
# 2026-06-19
#
# Inputs:
# - functions/gsea_waterfall_plot.R: plots ranked NES waterfalls
# - functions/gsea_volcano_plot.R: plots symmetric NES volcanoes
# - functions/gsea_half_volcano_plot.R: plots directional half-volcanoes
# - functions/gsea_nes_scatter_plot.R: plots cross-contrast NES scatterplots
# - tutorial/data/GSE122380_gsea_cardiomyocyte_vs_mesoderm.csv: minimal example GSEA table
# - tutorial/data/GSE122380_gsea_day9_vs_day6.csv: minimal scatter x-axis GSEA table
# - tutorial/data/GSE122380_gsea_day3_vs_day1.csv: minimal scatter y-axis GSEA table
#
# Outputs:
# - tutorial/assets/figures/GSE122380_gsea_waterfall_cardiomyocyte_vs_mesoderm_positive.svg
# - tutorial/assets/figures/GSE122380_gsea_waterfall_cardiomyocyte_vs_mesoderm_negative.svg
# - tutorial/assets/figures/GSE122380_gsea_volcano_cardiomyocyte_vs_mesoderm.svg
# - tutorial/assets/figures/GSE122380_gsea_half_volcano_cardiomyocyte_vs_mesoderm_positive.svg
# - tutorial/assets/figures/GSE122380_gsea_half_volcano_cardiomyocyte_vs_mesoderm_negative.svg
# - tutorial/assets/figures/GSE122380_gsea_scatter_day9_vs_day6_x_day3_vs_day1_all_quadrants.svg
#
# Purpose:
# Generate the waterfall, volcano, half-volcano, and comparative NES scatterplot
# figures used by the tutorial site.
#
# Notes:
# The script starts from precomputed GSEA result tables. It does not run DESeq2,
# build ranked gene vectors, or calculate GSEA enrichment.

# 1.1 source files -----------------
base::source('functions/gsea_waterfall_plot.R')
base::source('functions/gsea_volcano_plot.R')
base::source('functions/gsea_half_volcano_plot.R')
base::source('functions/gsea_nes_scatter_plot.R')

# 1.2 load canonical inputs -----------------
gsea_cardiomyocyte_vs_mesoderm <- read_gsea_result_csvs(
  'tutorial/data/GSE122380_gsea_cardiomyocyte_vs_mesoderm.csv')
gsea_day9_vs_day6 <- read_gsea_result_csvs(
  'tutorial/data/GSE122380_gsea_day9_vs_day6.csv')
gsea_day3_vs_day1 <- read_gsea_result_csvs(
  'tutorial/data/GSE122380_gsea_day3_vs_day1.csv')

# 1.3 define parameters -----------------
figure_family = 'Nimbus Sans'
display_contrast_label = 'Cardiomyocyte vs. Mesoderm'
figure_width = 9
figure_width_waterfall = 7.4
figure_width_volcano = 4.8
figure_width_directional = 7.4
figure_height_volcano = 4.9
figure_height_waterfall = 4.3
figure_height_directional = 4.3
figure_height_scatter = 6.375

# 1.4 define local helper functions -----------------
save_svg <- function(plot, path, width, height) {
  ggplot2::ggsave(
    filename = path,
    plot = plot,
    width = width,
    height = height,
    device = function(filename, width, height, ...) {
      svglite::svglite(
        filename = filename,
        width = width,
        height = height,
        web_fonts = nimbus_svg_web_fonts(),
        user_fonts = plot_user_fonts(font_family = figure_family),
        ...)
    },
    fix_text_size = FALSE,
    bg = 'white',
    limitsize = FALSE)
  embed_nimbus_svg_fonts(path)
}

save_figure_svg <- function(plot, path_stub, width, height) {
  save_svg(
    plot = plot,
    path = base::paste0(path_stub, '.svg'),
    width = width,
    height = height)
}

# 1.5 create directories -----------------
base::dir.create('tutorial/assets/figures', recursive = TRUE, showWarnings = FALSE)

# 2.0 define reusable term groups and labels -----------------
positive_term_groups <- base::list(
  'Ion channel' = base::c(
    'Calcium', 'Sodium', 'Potassium', 'Ion', 'Voltage', 'Channel',
    'electrical', 'Action potential', 'transmembrane', 'conduction'),
  'Metabolism' = base::c(
    'Metabolic', 'Metabolism', 'Electron transport chain', 'Mitochondria',
    'Mitochondrial', 'Oxidative', 'Phosphorylation', 'ATP', 'lipid',
    'biosynthetic', 'respiration'),
  'Muscle contraction' = base::c(
    'Muscle contraction', 'Muscle cell differentiation', 'Striated',
    'Sarcomere', 'Myofibril', 'Actin', 'Muscle'),
  'Heart development' = base::c(
    'atrial', 'ventricular', 'heart', 'cardiac'))
positive_label_terms <- base::c(
  'myofibril assembly',
  'cardiac muscle contraction',
  'cardiac muscle cell action potential',
  'calcium-mediated signaling',
  'electron transport chain',
  'muscle cell differentiation',
  'ATP metabolic process',
  'cardiac cell development',
  'cardiac chamber morphogenesis',
  'muscle contraction',
  'cardiac conduction',
  'cardiac muscle cell contraction',
  'heart growth',
  'cardiac atrium development',
  'heart trabecula morphogenesis')
negative_term_groups <- base::list(
  'Ribosomal' = base::c(
    'ribosome', 'ribosomal', 'spliceosomal', 'rRNA', 'RNA'),
  'Mitosis' = base::c(
    'Mitosis', 'Meiosis', 'Nuclear division', 'cell cycle', 'telomere',
    'spindle', 'Mitotic', 'Meiotic'),
  'DNA replication' = base::c(
    'chromosome', 'DNA', 'base', 'repair', 'DNA-'))
negative_term_group_colors <- base::c(
  'Ribosomal' = '#377EB8',
  'Mitosis' = '#984EA3',
  'DNA replication' = '#E41A1C')
positive_volcano_label_terms <- base::c(
  'cardiac muscle contraction',
  'muscle cell differentiation',
  'cardiac chamber morphogenesis',
  'heart growth',
  'heart valve development')
negative_volcano_label_terms <- base::c(
  'mitotic nuclear division',
  'DNA replication initiation',
  'mitotic G2/M transition checkpoint',
  'centriole assembly',
  'ribosome assembly')
positive_half_volcano_label_terms <- base::c(
  'myofibril assembly',
  'cardiac muscle contraction',
  'cardiac muscle cell action potential',
  'electron transport chain',
  'muscle cell differentiation',
  'cardiac cell development',
  'cardiac chamber morphogenesis',
  'muscle contraction',
  'cardiac conduction',
  'heart growth',
  'heart trabecula morphogenesis')
negative_label_terms <- base::c(
  'DNA-templated DNA replication',
  'DNA replication',
  'cell cycle checkpoint signaling',
  'mitotic cell cycle checkpoint signaling',
  'ribosome biogenesis',
  'recombinational repair',
  'DNA replication initiation',
  'cell cycle DNA replication',
  'spindle organization')
all_quadrants_scatter_label_terms <- base::c(
  'mesenchyme development',
  'mesoderm development',
  'mesodermal cell differentiation',
  'myofibril assembly',
  'heart process',
  'oxidative phosphorylation',
  'regulation of nuclear division')

# 3.0 save ranked waterfall figures -----------------
save_figure_svg(
  plot_gsea_waterfall(
    gsea_results = gsea_cardiomyocyte_vs_mesoderm,
    direction = 'positive',
    top_n = 100L,
    rank_by = 'NES',
    x_label = 'Ranked GO terms by NES for Cardiomyocyte vs. Mesoderm',
    label_n = 10L,
    label_terms = positive_label_terms,
    term_groups = positive_term_groups,
    category_n = 4L),
  'tutorial/assets/figures/GSE122380_gsea_waterfall_cardiomyocyte_vs_mesoderm_positive',
  width = figure_width_waterfall,
  height = figure_height_waterfall)

save_figure_svg(
  plot_gsea_waterfall(
    gsea_results = gsea_cardiomyocyte_vs_mesoderm,
    direction = 'negative',
    top_n = 100L,
    rank_by = 'NES',
    x_label = 'Ranked GO terms by NES for Cardiomyocyte vs. Mesoderm',
    label_n = 10L,
    term_groups = negative_term_groups,
    term_group_colors = negative_term_group_colors),
  'tutorial/assets/figures/GSE122380_gsea_waterfall_cardiomyocyte_vs_mesoderm_negative',
  width = figure_width_waterfall,
  height = figure_height_waterfall)

# 4.0 save volcano figures -----------------
save_figure_svg(
  plot_gsea_volcano(
    gsea_results = gsea_cardiomyocyte_vs_mesoderm,
    label_n = 10L,
    label_terms = base::c(positive_volcano_label_terms, negative_volcano_label_terms),
    contrast_label = display_contrast_label,
    label_words_per_line = 3L),
  'tutorial/assets/figures/GSE122380_gsea_volcano_cardiomyocyte_vs_mesoderm',
  width = figure_width_volcano,
  height = figure_height_volcano)

save_figure_svg(
  plot_gsea_half_volcano(
    gsea_results = gsea_cardiomyocyte_vs_mesoderm,
    direction = 'positive',
    label_n = 10L,
    label_terms = positive_half_volcano_label_terms,
    term_groups = positive_term_groups,
    category_n = 4L,
    contrast_label = display_contrast_label),
  'tutorial/assets/figures/GSE122380_gsea_half_volcano_cardiomyocyte_vs_mesoderm_positive',
  width = figure_width_directional,
  height = figure_height_directional)

save_figure_svg(
  plot_gsea_half_volcano(
    gsea_results = gsea_cardiomyocyte_vs_mesoderm,
    direction = 'negative',
    label_n = 10L,
    label_terms = negative_label_terms,
    term_groups = negative_term_groups,
    term_group_colors = negative_term_group_colors,
    contrast_label = display_contrast_label),
  'tutorial/assets/figures/GSE122380_gsea_half_volcano_cardiomyocyte_vs_mesoderm_negative',
  width = figure_width_directional,
  height = figure_height_directional)

# 5.0 save comparative NES scatterplots -----------------
save_figure_svg(
  plot_gsea_nes_scatter(
    gsea_x = gsea_day9_vs_day6,
    gsea_y = gsea_day3_vs_day1,
    x_label = 'Day 9 vs. Day 6 NES',
    y_label = 'Day 3 vs. Day 1 NES',
    color_by = 'significance',
    quadrant = 'all',
    x_name = 'Day 9 vs. Day 6',
    y_name = 'Day 3 vs. Day 1',
    label_terms = all_quadrants_scatter_label_terms,
    include_nonsignificant = TRUE,
    equal_axis_limits = TRUE),
  'tutorial/assets/figures/GSE122380_gsea_scatter_day9_vs_day6_x_day3_vs_day1_all_quadrants',
  width = figure_width,
  height = figure_height_scatter)

# 6.0 remove script scratch objects -----------------
base::rm(
  list = base::c(
    'gsea_cardiomyocyte_vs_mesoderm',
    'gsea_day9_vs_day6',
    'gsea_day3_vs_day1',
    'figure_family',
    'display_contrast_label',
    'figure_width',
    'figure_width_waterfall',
    'figure_width_volcano',
    'figure_width_directional',
    'figure_height_volcano',
    'figure_height_waterfall',
    'figure_height_directional',
    'figure_height_scatter',
    'save_svg',
    'save_figure_svg',
    'positive_term_groups',
    'positive_label_terms',
    'negative_term_groups',
    'negative_term_group_colors',
    'positive_volcano_label_terms',
    'negative_volcano_label_terms',
    'positive_half_volcano_label_terms',
    'negative_label_terms',
    'all_quadrants_scatter_label_terms'),
  envir = base::.GlobalEnv)
