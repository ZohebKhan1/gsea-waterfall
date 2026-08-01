# Maintenance contract

## Scope and scientific decisions

Open `gsea.Rproj` or run commands from the repository root. This repository is
a standalone plotting toolkit, not an R package. Files under `functions/` are
the copyable public implementation; scripts under `R/` own the example figures,
contract checks, and generated GitHub Pages site.

The toolkit receives precomputed GSEA tables. Experimental units, study blocks,
models, contrasts, filters, enrichment universes, ranking metrics, identifier
mapping, and multiple-testing procedures belong to the upstream analysis and
are outside this repository's authority. The maintained example preserves the
supplied NES, nominal p-values, adjusted p-values, term descriptions, and row
membership.

The maintained presentation decisions are:

- one GO biological-process term per input row;
- key precedence of explicit `id_col`, then `go_term_id` when present, otherwise
  a complete, unique term description;
- adjusted-p-value significance defined as `padj < 0.05`;
- waterfall ranking by NES before selecting 100 terms per direction;
- the established term groups and label selections; and
- the exact shared-key intersection for comparative plots.

Changing any item above is an analytical or visual-contract change and requires
appropriate before-and-after evidence.

## Public API

Users source `gsea_plot_utils.R` first, followed by the plot modules they need.
The supported public functions are:

- `read_gsea_result_csv()`
- `plot_gsea_waterfall()`
- `plot_gsea_volcano()`
- `plot_gsea_half_volcano()`
- `plot_gsea_nes_scatter()`

Plot functions return ggplot objects and do not write files. The reader returns
normalized `go_term_id`, `go_description`, `NES`, `pvalue`, `pval`, and `padj`
columns. `pval` is retained as an input-compatible alias so reader output can be
passed directly to plot functions that default to the common GSEA column name.
No compatibility function aliases or multi-file reader are supported.

## Output ownership

| Artifact | Unit and format | Writer | Consumer | Empty-result behavior |
| :-- | :-- | :-- | :-- | :-- |
| `tutorial/data/*.csv` | one GO term per row, CSV | bundled fixed input | reader, contract checks, figure script | not applicable |
| `tutorial/assets/figures/*.svg` | one editable plot, SVG | `R/01_generate_gsea_figures.R` | README and site renderer | plotting fails before a file is written |
| `docs/index.html`, `docs/search_index.json`, `docs/style.css`, `docs/.nojekyll`, `docs/libs/**`, `docs/assets/**` | one single-page bookdown site | `R/02_render_github_pages.R` | GitHub Pages | rendering fails; no empty site is admitted |

`tutorial/tutorial.Rmd` owns site content. The tutorial YAML, CSS, HTML include,
figures, and fonts are maintained build and presentation sources. Do not edit
generated files under `docs/`; the site writer removes and recreates only its
declared generated artifacts.

## Canonical verification

R Codex Utils is installed outside the repository and invoked from `PATH`.
Select the repository guideline explicitly for every canonical run:

```sh
export R_CODEX_UTILS_GUIDELINES=coding_guidelines_v8.0.md
r_codex_utils check parse --format text R/*.R functions/*.R
```

Run the scientific and public-boundary regression check with its dependencies
declared:

```sh
r_codex_utils script \
  --input functions/gsea_plot_utils.R \
  --input functions/gsea_waterfall_plot.R \
  --input functions/gsea_volcano_plot.R \
  --input functions/gsea_half_volcano_plot.R \
  --input functions/gsea_nes_scatter_plot.R \
  --input tutorial/data/GSE122380_gsea_cardiomyocyte_vs_mesoderm.csv \
  --input tutorial/data/GSE122380_gsea_day9_vs_day6.csv \
  --input tutorial/data/GSE122380_gsea_day3_vs_day1.csv \
  R/03_check_gsea_contracts.R
```

Run `R/01_generate_gsea_figures.R` with the source modules, CSV inputs, font
files, and six exact SVG outputs declared. Use `--replace-output` only for those
six SVG files. Run `R/02_render_github_pages.R` with the tutorial sources and
key site outputs declared; the script owns replacement of generated site
families.

After rendering, verify editable SVG text, four embedded Nimbus Sans faces, no
raster images or empty font sources, exact source-to-site figure copies, valid
search JSON, resolved local links, responsive layout, and absence of machine
paths or stale generated assets.

## Regression baseline

| Contract | Expected value |
| :-- | :-- |
| Cardiomyocyte vs. Mesoderm | 2,482 terms; 744 negative NES |
| Day 9 vs. Day 6 | 2,477 terms; 828 negative NES |
| Day 3 vs. Day 1 | 2,477 terms; 896 negative NES |
| Symmetric volcano | 272 significantly down; 445 significantly up |
| Comparative scatter | 1,327 neither; 601 Day 3-only; 280 both; 269 Day 9-only |

The six maintained SVG dimensions are stable. The scatter height is 459.36
points because its 6.38-inch source constant follows the two-decimal physical-
dimension rule; this presentation-only value does not alter plotted rows,
classifications, or coordinates.
