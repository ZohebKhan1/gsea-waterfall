#!/usr/bin/env Rscript

strip_html_text <- function(html) {
  text <- base::gsub('<script[^>]*>.*?</script>', ' ', html)
  text <- base::gsub('<style[^>]*>.*?</style>', ' ', text)
  text <- base::gsub('<[^>]+>', ' ', text)
  text <- base::gsub('&quot;', '"', text, fixed = TRUE)
  text <- base::gsub('&lt;', '<', text, fixed = TRUE)
  text <- base::gsub('&gt;', '>', text, fixed = TRUE)
  text <- base::gsub('&amp;', '&', text, fixed = TRUE)
  text <- base::gsub('[[:space:]]+', ' ', text)
  base::trimws(text)
}

write_search_index <- function(html_path, output_path) {
  html <- base::paste(base::readLines(html_path, warn = FALSE), collapse = ' ')
  title_match <- base::regmatches(html, base::regexpr('<title>[^<]+</title>', html))
  title <- base::sub('^<title>', '', base::sub('</title>$', '', title_match))
  search_index <- base::list(base::list('index.html', title, strip_html_text(html)))
  base::writeLines(jsonlite::toJSON(search_index, auto_unbox = TRUE, pretty = TRUE), output_path)
}

strip_trailing_whitespace <- function(path) {
  lines <- base::readLines(path, warn = FALSE)
  base::writeLines(base::sub('[[:space:]]+$', '', lines), path, useBytes = TRUE)
}

sync_asset_directory <- function(source_dir, output_dir) {
  if (!base::dir.exists(output_dir)) {
    base::dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  }
  source_files <- base::list.files(source_dir, recursive = TRUE, all.files = TRUE, no.. = TRUE)
  output_files <- base::list.files(output_dir, recursive = TRUE, all.files = TRUE, no.. = TRUE)
  stale_files <- base::setdiff(output_files, source_files)
  if (base::length(stale_files) > 0L) {
    base::unlink(base::file.path(output_dir, stale_files), recursive = TRUE)
  }
  base::file.copy(source_dir, base::dirname(output_dir), recursive = TRUE, overwrite = TRUE)
}

base::dir.create('docs', recursive = TRUE, showWarnings = FALSE)
base::dir.create('docs/assets', recursive = TRUE, showWarnings = FALSE)
base::dir.create('docs/libs', recursive = TRUE, showWarnings = FALSE)

# render the single-page tutorial before syncing assets
bookdown::render_book('tutorial')

rendered_html <- if (base::file.exists('tutorial/index.html')) 'tutorial/index.html' else 'tutorial/_site/index.html'
bookdown_lib_dir <- if (base::dir.exists('tutorial/_site/libs')) 'tutorial/_site/libs' else 'tutorial/libs'

base::file.copy(rendered_html, 'docs/index.html', overwrite = TRUE)
base::file.copy('tutorial/style.css', 'docs/style.css', overwrite = TRUE)
strip_trailing_whitespace('docs/index.html')

# keep github pages assets aligned with the tutorial source
sync_asset_directory('tutorial/assets', 'docs/assets')
if (base::dir.exists(bookdown_lib_dir)) {
  sync_asset_directory(bookdown_lib_dir, 'docs/libs')
}
write_search_index('docs/index.html', 'docs/search_index.json')

# remove bookdown intermediates after docs has been synced
base::unlink('tutorial/index.html')
base::unlink('tutorial/libs', recursive = TRUE)
base::unlink('tutorial/_site', recursive = TRUE)

base::rm(
  list = base::c(
    'strip_html_text',
    'strip_trailing_whitespace',
    'sync_asset_directory',
    'write_search_index',
    'rendered_html',
    'bookdown_lib_dir'),
  envir = base::.GlobalEnv)
