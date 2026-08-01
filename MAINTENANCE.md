# Development and maintenance

## Project scope

Open `gsea.Rproj` or run commands from the repository root. This is a
standalone plotting toolkit rather than an R package. Reusable functions live
under `functions/`; scripts under `R/` generate the example figures, render the
GitHub Pages site, and run regression checks.

The toolkit accepts precomputed GSEA results. Experimental design, contrasts,
gene-set universes, ranking statistics, identifier mapping, and multiple-testing
procedures are determined upstream. The plotting functions preserve the supplied
NES values, p-values, adjusted p-values, term descriptions, and row membership.

The maintained examples use these decisions:

- one GO biological-process term per input row;
- key precedence of explicit `id_col`, then `go_term_id`, then a complete and
  unique term description;
- significance defined as `padj < 0.05`;
- waterfall ranking by NES before selecting 100 terms per direction; and
- exact shared-key intersections for comparative plots, without imputation.

Changes to these decisions should include updated regression checks and a clear
explanation in the commit message or release notes.

## Public interface

Source `gsea_plot_utils.R` before the plot modules. The supported functions are:

- `read_gsea_result_csv()`
- `plot_gsea_waterfall()`
- `plot_gsea_volcano()`
- `plot_gsea_half_volcano()`
- `plot_gsea_nes_scatter()`

Plot functions return ggplot objects and do not write files. The reader returns
normalized `go_term_id`, `go_description`, `NES`, `pvalue`, `pval`, and `padj`
columns. The `pval` alias allows reader output to pass directly to plot functions
that use the common GSEA column name.

## Local checks

Install the runtime and development dependencies once:

```r
install.packages(c(
  "ggplot2", "ggrepel", "lintr", "RColorBrewer", "scales"
))
```

Parse and lint the maintained R files:

```sh
Rscript -e 'files <- list.files(c("R", "functions"), "[.]R$", full.names = TRUE); invisible(lapply(files, parse))'
Rscript -e 'lints <- c(lintr::lint_dir("R"), lintr::lint_dir("functions")); print(lints); quit(status = as.integer(length(lints) > 0L))'
```

Run the fixed-input and public-interface checks:

```sh
Rscript R/03_check_gsea_contracts.R
```

The same checks run in `.github/workflows/r-checks.yaml` for pushes and pull
requests targeting `main`.

## Figures and website

The maintained inputs are the three CSV files under `tutorial/data/`. Generate
the six SVG figures, then render the tutorial to `docs/`:

```sh
Rscript R/01_generate_gsea_figures.R
Rscript R/02_render_github_pages.R
```

Edit `tutorial/tutorial.Rmd`, the tutorial configuration, and files under
`tutorial/assets/`; do not hand-edit generated files under `docs/`. After
rendering, inspect the figures and confirm that local links, fonts, the search
index, and responsive layout work as expected.

## Regression reference

`R/03_check_gsea_contracts.R` asserts the maintained row counts, NES direction
counts, significance groups, identifier matching, and invalid-input behavior.
Update those expectations only when the example inputs or intended behavior
changes, and explain the reason alongside the code change.

## Release checklist

1. Run the parse, lint, and regression commands above.
2. Regenerate figures and the site when their source changes.
3. Review source and generated diffs separately.
4. Confirm the README example and rendered tutorial still work.
5. Summarize user-visible or scientific changes in the release notes.
