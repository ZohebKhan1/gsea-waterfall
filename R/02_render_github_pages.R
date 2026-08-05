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
# - tutorial/assets/fonts/LatinModernSans-Regular.otf
# - tutorial/assets/fonts/LatinModernSans-Bold.otf
# - tutorial/assets/fonts/LatinModernSans-Oblique.otf
# - tutorial/assets/fonts/LatinModernSans-BoldOblique.otf
#
# outputs:
# - docs/index.html
# - docs/search_index.json
# - docs/style.css
# - docs/.nojekyll
# - docs/libs/**
# - docs/assets/**
# ----

# 1.0 replace the generated GitHub Pages site -----------------

# remove stale copied dependencies and assets before regeneration
unlink(c('docs/libs', 'docs/assets'), recursive = TRUE)

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
