#!/usr/bin/env Rscript

# ----
# author:
# - Zoheb Khan
#
# script path:
# - R/01_generate_gsea_figures.R
#
# functions:
# - functions/gsea_visualizations.R
#
# input data:
# - tutorial/data/GSE122380_gsea_cardiomyocyte_vs_mesoderm.csv
# - tutorial/data/GSE122380_gsea_day9_vs_day6.csv
# - tutorial/data/GSE122380_gsea_day3_vs_day1.csv
# - tutorial/assets/fonts/NimbusSans-Regular.otf
# - tutorial/assets/fonts/NimbusSans-Bold.otf
# - tutorial/assets/fonts/NimbusSans-Italic.otf
# - tutorial/assets/fonts/NimbusSans-BoldItalic.otf
#
# outputs:
# - tutorial/assets/figures/GSE122380_gsea_waterfall_cardiomyocyte_vs_mesoderm_positive.svg
# - tutorial/assets/figures/GSE122380_gsea_waterfall_cardiomyocyte_vs_mesoderm_negative.svg
# - tutorial/assets/figures/GSE122380_gsea_volcano_cardiomyocyte_vs_mesoderm.svg
# - tutorial/assets/figures/GSE122380_gsea_half_volcano_cardiomyocyte_vs_mesoderm_positive.svg
# - tutorial/assets/figures/GSE122380_gsea_half_volcano_cardiomyocyte_vs_mesoderm_negative.svg
# - tutorial/assets/figures/GSE122380_gsea_scatter_day9_vs_day6_x_day3_vs_day1_all_quadrants.svg
# ----

# 0.0 verify project root and source functions -----------------

if (!file.exists('gsea.Rproj')) {
  stop('Run R/01_generate_gsea_figures.R from the gsea repository root.', call. = FALSE)
}

# source GSEA plotting functions
source('functions/gsea_visualizations.R')

# 1.0 define figure parameters and paths -----------------

gsea_figure_dir = 'tutorial/assets/figures'

gsea_figure_font_family = 'GSEA Nimbus Sans'
gsea_svg_font_family = 'Nimbus Sans'
gsea_figure_font_paths <- c(
  plain = 'tutorial/assets/fonts/NimbusSans-Regular.otf',
  bold = 'tutorial/assets/fonts/NimbusSans-Bold.otf',
  italic = 'tutorial/assets/fonts/NimbusSans-Italic.otf',
  bold_italic = 'tutorial/assets/fonts/NimbusSans-BoldItalic.otf')

waterfall_figure_width = 7.40
waterfall_figure_height = 4.30
volcano_figure_width = 4.80
volcano_figure_height = 4.90
half_volcano_figure_width = 7.40
half_volcano_figure_height = 4.30
nes_scatter_figure_width = 9.00
nes_scatter_figure_height = 6.38

gsea_padj_cutoff = 0.05
waterfall_top_n = 100L
half_volcano_inner_nes_limit = 1
waterfall_point_size = 1.30
volcano_point_size = 1.30
half_volcano_point_size = 1.60
nes_scatter_point_size = 1.30

gsea_contrast_label = 'Cardiomyocyte vs. Mesoderm'

positive_term_groups <- list(
  'Ion channel' = c(
    'Calcium', 'Sodium', 'Potassium', 'Ion', 'Voltage', 'Channel',
    'electrical', 'Action potential', 'transmembrane', 'conduction'),
  'Metabolism' = c(
    'Metabolic', 'Metabolism', 'Electron transport chain', 'Mitochondria',
    'Mitochondrial', 'Oxidative', 'Phosphorylation', 'ATP', 'lipid',
    'biosynthetic', 'respiration'),
  'Muscle contraction' = c(
    'Muscle contraction', 'Muscle cell differentiation', 'Striated',
    'Sarcomere', 'Myofibril', 'Actin', 'Muscle'),
  'Heart development' = c('atrial', 'ventricular', 'heart', 'cardiac'))

negative_term_groups <- list(
  'Ribosomal' = c('ribosome', 'ribosomal', 'spliceosomal', 'rRNA', 'RNA'),
  'Mitosis' = c(
    'Mitosis', 'Meiosis', 'Nuclear division', 'cell cycle', 'telomere',
    'spindle', 'Mitotic', 'Meiotic'),
  'DNA replication' = c('chromosome', 'DNA', 'base', 'repair', 'DNA-'))

negative_term_group_colors <- c(
  'Ribosomal' = '#377EB8',
  'Mitosis' = '#984EA3',
  'DNA replication' = '#E41A1C')

positive_waterfall_label_terms <- c(
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

volcano_label_terms <- c(
  'cardiac muscle contraction',
  'muscle cell differentiation',
  'cardiac chamber morphogenesis',
  'heart growth',
  'heart valve development',
  'mitotic nuclear division',
  'DNA replication initiation',
  'mitotic G2/M transition checkpoint',
  'centriole assembly',
  'ribosome assembly')

positive_half_volcano_label_terms <- c(
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

negative_half_volcano_label_terms <- c(
  'DNA-templated DNA replication',
  'DNA replication',
  'cell cycle checkpoint signaling',
  'mitotic cell cycle checkpoint signaling',
  'ribosome biogenesis',
  'recombinational repair',
  'DNA replication initiation',
  'cell cycle DNA replication',
  'spindle organization')

nes_scatter_label_terms <- c(
  'mesenchyme development',
  'mesoderm development',
  'mesodermal cell differentiation',
  'myofibril assembly',
  'heart process',
  'oxidative phosphorylation',
  'regulation of nuclear division')

# 1.1 read GSEA inputs -----------------

gsea_cardiomyocyte_vs_mesoderm <- read_gsea_result_csv(
  'tutorial/data/GSE122380_gsea_cardiomyocyte_vs_mesoderm.csv')
gsea_day9_vs_day6 <- read_gsea_result_csv(
  'tutorial/data/GSE122380_gsea_day9_vs_day6.csv')
gsea_day3_vs_day1 <- read_gsea_result_csv(
  'tutorial/data/GSE122380_gsea_day3_vs_day1.csv')

# 1.2 define SVG writer -----------------

# create inline web-font CSS before opening the SVG device
make_embedded_font_face <- function(font_path, weight, style) {
  font_uri <- base64enc::dataURI(file = font_path, mime = 'font/otf')
  structure(
    paste0(
      '    @font-face {\n',
      '      font-family: "', gsea_svg_font_family, '";\n',
      '      src: url("', font_uri, '") format("opentype");\n',
      '      font-weight: ', weight, ';\n',
      '      font-style: ', style, ';\n',
      '    }'),
    class = c('font_face', 'character'))
}

systemfonts::register_font(
  name = gsea_figure_font_family,
  plain = gsea_figure_font_paths[['plain']],
  bold = gsea_figure_font_paths[['bold']],
  italic = gsea_figure_font_paths[['italic']],
  bolditalic = gsea_figure_font_paths[['bold_italic']])

gsea_svg_font_faces <- list(
  make_embedded_font_face(
    font_path = gsea_figure_font_paths[['plain']],
    weight = '400',
    style = 'normal'),
  make_embedded_font_face(
    font_path = gsea_figure_font_paths[['bold']],
    weight = '700',
    style = 'normal'),
  make_embedded_font_face(
    font_path = gsea_figure_font_paths[['italic']],
    weight = '400',
    style = 'italic'),
  make_embedded_font_face(
    font_path = gsea_figure_font_paths[['bold_italic']],
    weight = '700',
    style = 'italic'))

save_gsea_figure <- function(plot, figure_path, figure_width, figure_height) {
  dir.create(dirname(figure_path), recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(
    filename = figure_path,
    plot = plot,
    device = function(filename, width, height, bg, ...) {
      svglite::svglite(
        filename = filename,
        width = width,
        height = height,
        bg = bg,
        web_fonts = gsea_svg_font_faces,
        fix_text_size = FALSE,
        ...)
    },
    width = figure_width,
    height = figure_height,
    units = 'in',
    bg = 'white')
  invisible(figure_path)
}

# 2.0 create and save tutorial figures -----------------

save_gsea_figure(
  plot = plot_gsea_waterfall(
    gsea_results = gsea_cardiomyocyte_vs_mesoderm,
    direction = 'positive',
    top_n = waterfall_top_n,
    rank_by = 'NES',
    x_label = 'Ranked GO terms by NES for Cardiomyocyte vs. Mesoderm',
    label_terms = positive_waterfall_label_terms,
    term_groups = positive_term_groups,
    point_size = waterfall_point_size,
    font_family = gsea_figure_font_family),
  figure_path = file.path(
    gsea_figure_dir,
    'GSE122380_gsea_waterfall_cardiomyocyte_vs_mesoderm_positive.svg'),
  figure_width = waterfall_figure_width,
  figure_height = waterfall_figure_height)

save_gsea_figure(
  plot = plot_gsea_waterfall(
    gsea_results = gsea_cardiomyocyte_vs_mesoderm,
    direction = 'negative',
    top_n = waterfall_top_n,
    rank_by = 'NES',
    x_label = 'Ranked GO terms by NES for Cardiomyocyte vs. Mesoderm',
    label_n = 10L,
    term_groups = negative_term_groups,
    term_group_colors = negative_term_group_colors,
    point_size = waterfall_point_size,
    font_family = gsea_figure_font_family),
  figure_path = file.path(
    gsea_figure_dir,
    'GSE122380_gsea_waterfall_cardiomyocyte_vs_mesoderm_negative.svg'),
  figure_width = waterfall_figure_width,
  figure_height = waterfall_figure_height)

save_gsea_figure(
  plot = plot_gsea_volcano(
    gsea_results = gsea_cardiomyocyte_vs_mesoderm,
    padj_cutoff = gsea_padj_cutoff,
    label_terms = volcano_label_terms,
    contrast_label = gsea_contrast_label,
    label_words_per_line = 3L,
    point_size = volcano_point_size,
    font_family = gsea_figure_font_family),
  figure_path = file.path(
    gsea_figure_dir,
    'GSE122380_gsea_volcano_cardiomyocyte_vs_mesoderm.svg'),
  figure_width = volcano_figure_width,
  figure_height = volcano_figure_height)

save_gsea_figure(
  plot = plot_gsea_half_volcano(
    gsea_results = gsea_cardiomyocyte_vs_mesoderm,
    direction = 'positive',
    p_col = 'padj',
    padj_cutoff = gsea_padj_cutoff,
    label_terms = positive_half_volcano_label_terms,
    term_groups = positive_term_groups,
    inner_nes_limit = half_volcano_inner_nes_limit,
    contrast_label = gsea_contrast_label,
    point_size = half_volcano_point_size,
    font_family = gsea_figure_font_family),
  figure_path = file.path(
    gsea_figure_dir,
    'GSE122380_gsea_half_volcano_cardiomyocyte_vs_mesoderm_positive.svg'),
  figure_width = half_volcano_figure_width,
  figure_height = half_volcano_figure_height)

save_gsea_figure(
  plot = plot_gsea_half_volcano(
    gsea_results = gsea_cardiomyocyte_vs_mesoderm,
    direction = 'negative',
    p_col = 'padj',
    padj_cutoff = gsea_padj_cutoff,
    label_terms = negative_half_volcano_label_terms,
    term_groups = negative_term_groups,
    term_group_colors = negative_term_group_colors,
    inner_nes_limit = half_volcano_inner_nes_limit,
    contrast_label = gsea_contrast_label,
    point_size = half_volcano_point_size,
    font_family = gsea_figure_font_family),
  figure_path = file.path(
    gsea_figure_dir,
    'GSE122380_gsea_half_volcano_cardiomyocyte_vs_mesoderm_negative.svg'),
  figure_width = half_volcano_figure_width,
  figure_height = half_volcano_figure_height)

save_gsea_figure(
  plot = plot_gsea_nes_scatter(
    gsea_x = gsea_day9_vs_day6,
    gsea_y = gsea_day3_vs_day1,
    x_label = 'Day 9 vs. Day 6 NES',
    y_label = 'Day 3 vs. Day 1 NES',
    color_by = 'significance',
    quadrant = 'all',
    x_name = 'Day 9 vs. Day 6',
    y_name = 'Day 3 vs. Day 1',
    label_terms = nes_scatter_label_terms,
    padj_cutoff = gsea_padj_cutoff,
    include_nonsignificant = TRUE,
    equal_axis_limits = TRUE,
    show_fit_line = TRUE,
    point_size = nes_scatter_point_size,
    font_family = gsea_figure_font_family),
  figure_path = file.path(
    gsea_figure_dir,
    'GSE122380_gsea_scatter_day9_vs_day6_x_day3_vs_day1_all_quadrants.svg'),
  figure_width = nes_scatter_figure_width,
  figure_height = nes_scatter_figure_height)
