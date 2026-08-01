#!/usr/bin/env Rscript

# ----
# author:
# - Zoheb Khan
#
# script path:
# - R/02_render_github_pages.R
#
# input data:
# - tutorial/tutorial.Rmd
# - tutorial/_bookdown.yml
# - tutorial/_output.yml
# - tutorial/style.css
# - tutorial/code-collapse.html
# - tutorial/assets/figures/*.svg
# - tutorial/assets/fonts/*.otf
#
# outputs:
# - docs/index.html
# - docs/search_index.json
# - docs/style.css
# - docs/.nojekyll
# - docs/libs/**
# - docs/assets/**
# ----

# 0.0 verify project root -----------------

if (!file.exists('gsea.Rproj')) {
  stop('Run R/02_render_github_pages.R from the gsea repository root.', call. = FALSE)
}

# 1.0 replace the generated GitHub Pages site -----------------

# remove only bookdown-owned site artifacts before regeneration
unlink(c('docs/index.html', 'docs/search_index.json', 'docs/style.css'))
unlink('docs/libs', recursive = TRUE)
unlink('docs/assets', recursive = TRUE)

bookdown::render_book(
  input = 'tutorial',
  output_format = 'bookdown::gitbook',
  clean = TRUE)

unlink(c('docs/tutorial.md', 'docs/reference-keys.txt', 'tutorial/index.rds'))
unlink(list.files(
  'docs',
  pattern = '^\\.DS_Store$',
  all.files = TRUE,
  full.names = TRUE,
  recursive = TRUE))
file.create('docs/.nojekyll')
