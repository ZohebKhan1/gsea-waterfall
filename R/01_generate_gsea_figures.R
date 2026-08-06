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

# source GSEA plotting functions
source('functions/gsea_visualizations.R')

# 1.0 define relevant directory paths -----------------

gsea_figure_dir = 'tutorial/assets/figures'

gsea_figure_font_family = 'GSEA Nimbus Sans'
gsea_svg_font_family = 'Nimbus Sans'
gsea_figure_font_paths <- c(
  plain = 'tutorial/assets/fonts/NimbusSans-Regular.otf',
  bold = 'tutorial/assets/fonts/NimbusSans-Bold.otf',
  italic = 'tutorial/assets/fonts/NimbusSans-Italic.otf',
  bold_italic = 'tutorial/assets/fonts/NimbusSans-BoldItalic.otf'
)

# 1.1 define relevant figure parameter parameters -----------------

waterfall_figure_width = 7.40
waterfall_figure_height = 4.30
volcano_figure_width = 4.80
volcano_figure_height = 4.90
half_volcano_figure_width = 7.40
half_volcano_figure_height = 4.30
nes_scatter_figure_width = 9.00
nes_scatter_figure_height = 6.38

gsea_padj_cutoff = 0.05
waterfall_top_n = 100
half_volcano_inner_nes_limit = 1
volcano_point_size = 1.30
half_volcano_point_size = 1.60
nes_scatter_point_size = 1.30

# 1.1 define labels, groupings for all plots -----------------

## waterfall plot labeling/categorization

gsea_contrast_label = 'Cardiomyocyte vs. Mesoderm'

color_by_group_categories <- list(
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
  'Heart development' = c(
    'atrial', 'ventricular', 'heart', 'cardiac')
)

positive_waterfall_label_terms <- c(
  'myofibril assembly',
  'cardiac muscle cell action potential',
  'calcium-mediated signaling',
  'electron transport chain',
  'muscle cell differentiation',
  'cardiac cell development',
  'cardiac chamber morphogenesis',
  'muscle contraction',
  'cardiac conduction',
  'cardiac muscle cell contraction',
  'heart growth',
  'cardiac atrium development',
  'heart trabecula morphogenesis'
)

negative_color_by_group_categories <- list(
  'Ribosomal' = c(
    'ribosome', 'ribosomal', 'spliceosomal', 'rRNA', 'RNA'),
  'Mitosis' = c(
    'Mitosis', 'Meiosis', 'Nuclear division', 'cell cycle', 'telomere',
    'spindle', 'Mitotic', 'Meiotic'),
  'DNA replication' = c(
    'chromosome', 'DNA', 'base', 'repair', 'DNA-')
)

negative_waterfall_label_terms <- c(
  'DNA replication',
  'chromosome segregation',
  'double-strand break repair',
  'cell cycle G2/M phase transition',
  'maturation of 5.8S rRNA',
  'mitotic nuclear division',
  'regulation of cell cycle phase transition',
  'regulation of mitotic cell cycle phase transition',
  'spliceosomal complex assembly'
)

## symmetrical volcano plot GO term labels

positive_volcano_label_terms <- c(
  'muscle tissue development',
  'heart process',
  'muscle cell development',
  'myofibril assembly'
)

negative_volcano_label_terms <- c(
  'DNA-templated DNA replication',
  'double-strand break repair',
  'interstrand cross-link repair',
  'positive regulation of mitotic cell cycle'
)
volcano_label_terms <- c(
  negative_volcano_label_terms,
  positive_volcano_label_terms
)

## half volcano plot GO term labels

positive_half_volcano_label_terms <- c(
  'response to oxidative stress',
  'lipid localization',
  'ribose phosphate biosynthetic process',
  'smooth muscle cell proliferation',
  'cardiac ventricle development',
  'ATP biosynthetic process',
  'heart process',
  'myofibril assembly'
)

negative_half_volcano_label_terms <- c(
  'DNA-templated DNA replication',
  'mitotic cell cycle checkpoint signaling',
  'spliceosomal tri-snRNP complex assembly',
  'base-excision repair',
  'positive regulation of mitotic cell cycle',
  'negative regulation of telomere maintenance',
  'regulation of RNA splicing',
  'DNA-templated transcription elongation'
)

## comparative NES scatterplot GO term labels

nes_scatter_label_terms <- c(
  'mesenchyme development',
  'mesoderm development',
  'mesodermal cell differentiation',
  'myofibril assembly',
  'heart process',
  'oxidative phosphorylation',
  'regulation of nuclear division'
)

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
    padj_threshold = gsea_padj_cutoff,
    NES_direction = 'positive',
    top_n = waterfall_top_n,
    rank_by = 'NES',
    label_by_groups = positive_waterfall_label_terms,
    color_by_groups = color_by_group_categories,
    font_family = gsea_figure_font_family),
  figure_path = file.path(
    gsea_figure_dir,
    'GSE122380_gsea_waterfall_cardiomyocyte_vs_mesoderm_positive.svg'),
  figure_width = waterfall_figure_width,
  figure_height = waterfall_figure_height)

save_gsea_figure(
  plot = plot_gsea_waterfall(
    gsea_results = gsea_cardiomyocyte_vs_mesoderm,
    padj_threshold = gsea_padj_cutoff,
    NES_direction = 'negative',
    top_n = waterfall_top_n,
    rank_by = 'NES',
    label_by_groups = negative_waterfall_label_terms,
    color_by_groups = negative_color_by_group_categories,
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
    term_groups = color_by_group_categories,
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
    term_groups = negative_color_by_group_categories,
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
    equal_axis_limits = TRUE,
    show_fit_line = TRUE,
    point_size = nes_scatter_point_size,
    font_family = gsea_figure_font_family),
  figure_path = file.path(
    gsea_figure_dir,
    'GSE122380_gsea_scatter_day9_vs_day6_x_day3_vs_day1_all_quadrants.svg'),
  figure_width = nes_scatter_figure_width,
  figure_height = nes_scatter_figure_height)
