# gsea-waterfall

`gsea-waterfall` is a small, copyable R toolkit for visualizing precomputed Gene
Set Enrichment Analysis (GSEA) results. It provides ranked NES waterfall plots,
symmetric and directional volcano plots, and cross-contrast NES scatterplots.
The functions return ggplot objects and do not run GSEA or write files.

[Rendered tutorial](https://zohebkhan1.github.io/gsea-waterfall/)

## Quick start

Install the plotting dependencies:

```r
utils::install.packages(c('ggplot2', 'ggrepel', 'RColorBrewer', 'scales'))
```

Download the shared utilities first, followed by the plot modules:

```r
base_url <- 'https://raw.githubusercontent.com/ZohebKhan1/gsea-waterfall/main/functions'

utils::download.file(paste0(base_url, '/gsea_plot_utils.R'), 'gsea_plot_utils.R')
utils::download.file(paste0(base_url, '/gsea_waterfall_plot.R'), 'gsea_waterfall_plot.R')
utils::download.file(paste0(base_url, '/gsea_volcano_plot.R'), 'gsea_volcano_plot.R')
utils::download.file(paste0(base_url, '/gsea_half_volcano_plot.R'), 'gsea_half_volcano_plot.R')
utils::download.file(paste0(base_url, '/gsea_nes_scatter_plot.R'), 'gsea_nes_scatter_plot.R')

source('gsea_plot_utils.R')
source('gsea_waterfall_plot.R')
source('gsea_volcano_plot.R')
source('gsea_half_volcano_plot.R')
source('gsea_nes_scatter_plot.R')
```

Read one precomputed contrast and build a plot:

```r
gsea_results <- read_gsea_result_csv('my_gsea_results.csv')

waterfall_plot <- plot_gsea_waterfall(
  gsea_results = gsea_results,
  direction = 'positive',
  top_n = 100L,
  rank_by = 'NES'
)

ggplot2::ggsave(
  'gsea_waterfall.svg', waterfall_plot,
  width = 7.4, height = 4.3, units = 'in'
)
```

## Input contract

Each input table represents one GSEA contrast with one row per unique term:

| Column | Required | Meaning |
| :-- | :-- | :-- |
| `go_description` | yes | term description and fallback matching key |
| `go_term_id` | no | stable matching key used automatically when present |
| `NES` | yes | normalized enrichment score |
| `pval` | yes | nominal GSEA p-value in `[0, 1]` |
| `padj` | yes | adjusted p-value in `[0, 1]` |

Required columns must be complete, and NES values must be finite. Key precedence
is explicit `id_col`, then `go_term_id` when present, otherwise
`go_description`. Comparative plots use the exact shared-key intersection and
do not impute absent terms.

If your GSEA table uses different column names, map them when reading it:

```r
fgsea_results <- read_gsea_result_csv(
  'fgsea_results.csv',
  term_col = 'pathway',
  nes_col = 'NES',
  pvalue_col = 'pval',
  padj_col = 'padj'
)
```

For comparative plots, use `id_col` with a stable term identifier whenever one
is available. Without stable IDs, term descriptions must match exactly across
tables.

## Most-used parameters

The following controls cover most user workflows. Styling and layout arguments
remain available in the function definitions for users who need them.

| Function | Parameter | What it changes |
| :-- | :-- | :-- |
| `plot_gsea_waterfall()` | `direction` | Plots positive or negative NES terms. |
|  | `top_n` | Keeps this many terms after direction filtering and ranking. |
|  | `rank_by` | Ranks by NES, adjusted p-value, or nominal p-value. |
|  | `label_terms`, `label_n` | Labels selected terms or automatically selects labels. |
|  | `term_groups` | Applies caller-defined keyword/ID groups; omitted groups remain neutral. |
| `plot_gsea_volcano()` | `padj_cutoff` | Defines the significance colors and count labels. |
|  | `label_terms`, `label_n` | Selects exact labels or the strongest significant labels. |
| `plot_gsea_half_volcano()` | `direction` | Shows one NES direction. |
|  | `p_col` | Chooses nominal or adjusted p-value for the y-axis. |
|  | `inner_nes_limit` | Removes terms whose absolute NES is below the boundary. |
| `plot_gsea_nes_scatter()` | `include_nonsignificant` | Keeps terms that are not significant in either contrast. |
|  | `quadrant` | Shows all terms, Q1, or Q3. |
|  | `equal_axis_limits` | Uses matched x/y limits for direct comparison. |
|  | `show_fit_line` | Adds an optional descriptive linear fit to the plotted terms. |

All plot functions return `ggplot` objects. They only select rows for display;
they do not rewrite the supplied GSEA statistics or write files.

## Plot functions

### `plot_gsea_waterfall()`

Ranks one NES direction. Terms are neutral by default; supply `term_groups` to
color caller-defined categories.

The displayed assets use explicit category seeds and selected labels from
`R/01_generate_gsea_figures.R`; the parameter examples above are intentionally
compact for reuse with other GSEA tables.

<img src="tutorial/assets/figures/GSE122380_gsea_waterfall_cardiomyocyte_vs_mesoderm_positive.svg" alt="Positive NES waterfall" width="760">

**Figure 1.** Positive NES GO biological-process waterfall for Cardiomyocyte vs.
Mesoderm.

<img src="tutorial/assets/figures/GSE122380_gsea_waterfall_cardiomyocyte_vs_mesoderm_negative.svg" alt="Negative NES waterfall" width="760">

**Figure 2.** Negative NES waterfall highlighting ribosomal, mitotic, and
DNA-replication term groups.

### `plot_gsea_volcano()`

Shows positive, negative, and non-significant terms in one symmetric NES
volcano plot.

<img src="tutorial/assets/figures/GSE122380_gsea_volcano_cardiomyocyte_vs_mesoderm.svg" alt="Symmetric GSEA volcano" width="560">

**Figure 3.** Symmetric GO biological-process GSEA volcano for Cardiomyocyte vs.
Mesoderm.

### `plot_gsea_half_volcano()`

Expands one NES direction and optionally applies caller-supplied term groups.

<img src="tutorial/assets/figures/GSE122380_gsea_half_volcano_cardiomyocyte_vs_mesoderm_positive.svg" alt="Positive NES half-volcano" width="760">

**Figure 4.** Positive NES half-volcano.

<img src="tutorial/assets/figures/GSE122380_gsea_half_volcano_cardiomyocyte_vs_mesoderm_negative.svg" alt="Negative NES half-volcano" width="760">

**Figure 5.** Negative NES half-volcano.

### `plot_gsea_nes_scatter()`

Compares matched term NES values across two contrasts and can color terms by
significance pattern or a supplied biological grouping.

<img src="tutorial/assets/figures/GSE122380_gsea_scatter_day9_vs_day6_x_day3_vs_day1_all_quadrants.svg" alt="Comparative NES scatterplot" width="900">

**Figure 6.** Comparative NES scatterplot for Day 9 vs. Day 6 and Day 3 vs. Day
1.

The [tutorial](https://zohebkhan1.github.io/gsea-waterfall/) walks through the
same workflow with copyable calls and figure-specific parameter explanations.

## Reproduce the figures and site

Open `gsea.Rproj` or run from the repository root:

```r
utils::install.packages(c(
  'base64enc', 'bookdown', 'knitr', 'rmarkdown', 'svglite', 'systemfonts'
))
```

```sh
Rscript R/01_generate_gsea_figures.R
Rscript R/02_render_github_pages.R
Rscript R/03_check_gsea_contracts.R
```

The figure script writes six editable, font-embedded SVGs under
`tutorial/assets/figures/`. The site script renders `tutorial/tutorial.Rmd`
directly to `docs/`; generated files under `docs/` should not be edited by
hand. The contract script verifies fixed scientific counts, stable-key matching,
and invalid-input boundaries.

## Repository layout

```text
functions/          downloadable plotting functions and shared utilities
R/                  figure generation, site rendering, and contract checks
tutorial/data/      fixed precomputed example GSEA tables
tutorial/assets/    generated tutorial figures and bundled fonts
tutorial/tutorial.Rmd
docs/               generated GitHub Pages site
```

## Scientific context

The waterfall presentation was inspired by pathway-level figures in
[Ciceri et al., *Nature*, 2024](https://doi.org/10.1038/s41586-023-06984-8),
[Xu et al., *Nature Cell Biology*, 2025](https://doi.org/10.1038/s41556-025-01751-5),
[Vuong et al., *Science*, 2026](https://doi.org/10.1126/science.aea1259), and
[Risgaard et al., *Science*, 2026](https://doi.org/10.1126/science.aea1549).

The underlying GSEA method is described by
[Subramanian et al.](https://www.pnas.org/doi/10.1073/pnas.0506580102). The
[fgsea](https://bioconductor.org/packages/fgsea/) package is one source of
precomputed tables compatible with this plotting contract.
