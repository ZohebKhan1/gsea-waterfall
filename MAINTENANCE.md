# Maintenance contract

## Project and scientific boundaries

Open `gsea.Rproj` or run commands from the repository root. The repository is a
standalone plotting toolkit, not an R package. The files in `functions/` are the
copyable public implementation; the scripts in `R/` own the example figures and
generated GitHub Pages site.

The three example CSV files contain one row per unique GO biological-process
term and one precomputed GSEA contrast. The plotting workflow does not select
samples, fit a model, calculate enrichment, adjust p-values, or reinterpret
missing values. It preserves the supplied NES, nominal p-value, and adjusted
p-value values. The example uses the existing adjusted p-value cutoff of 0.05,
term groups, label selections, direction filters, NES ranking, and cross-contrast
term intersection.

## Public API

Users source `gsea_plot_utils.R` first, followed by any plotting module they
need. The supported public functions are:

- `read_gsea_result_csv()`
- `plot_gsea_waterfall()`
- `plot_gsea_volcano()`
- `plot_gsea_half_volcano()`
- `plot_gsea_nes_scatter()`

The August 2026 cleanup removed undocumented `make_gsea_*()` aliases, the
duplicate `category_n` argument, the speculative `max_categories` limit, the
scatter alias `label_n_one_contrast`, and the multi-file CSV reader behavior.
All repository callers and tutorial examples use the API above.

## Output ownership

| Artifact | Row or visual unit | Writer | Consumer |
| :-- | :-- | :-- | :-- |
| `tutorial/data/*.csv` | one GO term per contrast | fixed example input | plotting functions and tutorial |
| `tutorial/assets/figures/*.svg` | one editable example figure | `R/01_generate_gsea_figs.R` | README and tutorial |
| `docs/index.html`, `docs/search_index.json`, `docs/style.css`, `docs/libs/**`, `docs/assets/**` | one single-page site and its linked assets | `R/02_render_github_pgs.R` via bookdown | GitHub Pages |

`tutorial/tutorial.Rmd` is the site source of truth. Do not edit generated files
under `docs/` to correct source content.

## Canonical verification

R Codex Utils is installed outside the repository and invoked from `PATH`; the
repository does not vendor its launcher. During the August 2026 modernization,
the local `coding_guidelines_v8.0.md` file was selected explicitly with
`R_CODEX_UTILS_GUIDELINES=coding_guidelines_v8.0.md`.

Set the reviewed local guideline, then parse all maintained R sources before
execution:

```sh
export R_CODEX_UTILS_GUIDELINES=coding_guidelines_v8.0.md
r_codex_utils check parse --format text R/*.R functions/*.R
```

Run the durable scientific/API regression check with its dependencies declared:

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

Run the figure and site scripts with the sources, inputs, and output families
declared. `R/01_generate_gsea_figs.R` lists its six exact SVG outputs;
`R/02_render_github_pgs.R` lists the generated site families. Use
`--replace-output` for the six SVG files only. The site script itself removes
and recreates only its owned `docs/` artifacts.

After rendering, verify that `docs/` contains no machine-specific paths or
stale assets, inspect every SVG for editable text and four embedded Nimbus Sans
faces, and visually compare the six figures and site against the checkpoint
baseline.

## August 2026 verification record

The modernization was closed against the checkpoint baseline with these
results:

- all eight maintained R files parsed successfully, and the durable contract
  script completed without an error;
- the fixed contrasts retained 2,482, 2,477, and 2,477 GO-term rows, with 744,
  828, and 896 negative NES values, respectively;
- the volcano retained 272 significantly down and 445 significantly up terms;
  the scatter retained 1,327 nonsignificant, 601 Day 3-only, 280 shared, and 269
  Day 9-only terms;
- all six SVGs retained editable text, contained four embedded Nimbus Sans font
  faces, and contained no raster image or empty font source;
- the two waterfall rasters matched the checkpoint exactly, while the other
  normalized raster comparisons had RMSE below 0.04. The volcano panel geometry
  matched exactly;
- the scatter height is 459.36 points rather than 459.00 points because its
  6.38-inch maintained constant follows the two-decimal dimension rule. This
  0.36-point presentation-only change does not affect rows, classifications, or
  plotted coordinates; and
- the generated site contained 34 owned files, all 28 local file references
  resolved, its source and generated styles and figures matched, and desktop
  and 430-pixel responsive layouts were visually inspected. No machine paths,
  stale copied tools, or generated `.DS_Store` files remained at closure.
