# gsea-waterfall

`gsea-waterfall` visualizes precomputed Gene Set Enrichment Analysis (GSEA)
results as ranked NES waterfall plots, volcano plots, and cross-contrast
comparative NES scatterplots. The functions return `ggplot` objects. They are a
visualization layer applied after GSEA has been performed; they do not run GSEA
or write files.

[GSEA visualization examples shown here](https://zohebkhan1.github.io/gsea-waterfall/)

## Example data and preprocessing

The bundled tables are derived from the public [NCBI GEO accession
GSE122380](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE122380), a
bulk RNA-seq time course of human iPSC differentiation into cardiomyocytes.
Basic sequencing and sample-level QC was performed before differential
expression analysis. A DESeq2 comparison of day 9 (cardiomyocyte) versus day 3
(mesoderm) generated the cardiomyocyte-versus-mesoderm contrast used for the
bundled GSEA results.

## Use case

GSEA results are commonly shown as running enrichment plots for individual
pathways or as dotplot/barplot-style summaries of a short list of terms. Those
views work well for a focused pathway or compact top-term summary, but they
provide limited room to inspect the broader ranked result.

The waterfall plots are intended for that complementary use case. They allow
roughly 100–200 ranked GO terms to be viewed together, including terms beyond
the usual top 10–20. This can make related parent/child terms and broader
functional patterns visible in the same ranked context. Caller-defined
`term_groups` can organize displayed terms by biological theme; they are
annotations only and do not perform redundancy reduction or additional
statistics.

The visual design was informed by pathway-level figures from Ciceri et al.
(2024), Xu et al. (2025), Vuong et al. (2026), and Risgaard et al. (2026). Full
citations are provided at the end of this README.

## Repository layout

```text
functions/          downloadable plotting functions and shared utilities
R/                  figure generation, site rendering, and contract checks
tutorial/data/      fixed precomputed example GSEA tables
tutorial/assets/    generated tutorial figures and bundled fonts
tutorial/tutorial.Rmd
docs/               generated GitHub Pages site
```

## Quick start

Install the plotting dependencies:

```r
install.packages(c('ggplot2', 'ggrepel', 'RColorBrewer', 'scales'))
```

Download the shared utilities first, followed by the plot modules:

```r
base_url <- 'https://raw.githubusercontent.com/ZohebKhan1/gsea-waterfall/main/functions'

download.file(paste0(base_url, '/gsea_plot_utils.R'), 'gsea_plot_utils.R')
download.file(paste0(base_url, '/gsea_waterfall_plot.R'), 'gsea_waterfall_plot.R')
download.file(paste0(base_url, '/gsea_volcano_plot.R'), 'gsea_volcano_plot.R')
download.file(paste0(base_url, '/gsea_half_volcano_plot.R'), 'gsea_half_volcano_plot.R')
download.file(paste0(base_url, '/gsea_nes_scatter_plot.R'), 'gsea_nes_scatter_plot.R')

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
```

## Input data

Each input table represents one GSEA contrast with one row per unique term:

The completed result table can come from any GSEA workflow, including `fgsea`,
`clusterProfiler`, or another package. If the source uses different column
names, map them when reading the table as shown below.

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

```r
fgsea_results <- read_gsea_result_csv(
  'fgsea_results.csv',
  term_col = 'pathway',
  nes_col = 'NES',
  pvalue_col = 'pval',
  padj_col = 'padj'
)
```

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

## Public API

Each plot function returns a `ggplot` object and does not write files. The
source files contain the complete argument documentation; the table below is
the practical overview of the public functions.

| Function | Purpose |
| :-- | :-- |
| `read_gsea_result_csv()` | Reads and validates one completed GSEA result table, with optional source-column mapping. |
| `plot_gsea_waterfall()` | Ranks positive or negative NES terms and optionally groups or labels them. |
| `plot_gsea_volcano()` | Shows both NES directions against adjusted-p-value significance in one panel. |
| `plot_gsea_half_volcano()` | Shows one NES direction with nominal or adjusted p-value on the y-axis. |
| `plot_gsea_nes_scatter()` | Compares matching term NES values across two GSEA contrasts. |

The [GSEA visualization examples shown here](https://zohebkhan1.github.io/gsea-waterfall/)
walk through the same workflow with copyable calls and figure-specific
parameter explanations.

## Citations

1. Subramanian A, Tamayo P, Mootha VK, et al. Gene set enrichment analysis: a knowledge-based approach for interpreting genome-wide expression profiles. *PNAS*. 2005. [Link](https://www.pnas.org/doi/10.1073/pnas.0506580102)
2. Ciceri G, Baggiolini A, Cho HS, et al. An epigenetic barrier sets the timing of human neuronal maturation. *Nature*. 2024;626:881–890. [Link](https://doi.org/10.1038/s41586-023-06984-8)
3. Xu N, Cho HS, Hackland JOS, et al. Genome-wide CRISPR screen identifies Menin and SUZ12 as regulators of human developmental timing. *Nature Cell Biology*. 2025;27:1411–1421. [Link](https://doi.org/10.1038/s41556-025-01751-5)
4. Vuong CK, Weber A, Seong P, et al. A single-cell multiomic analysis identifies molecular and gene-regulatory mechanisms dysregulated in developing Down syndrome neocortex. *Science*. 2026;392:eaea1259. [Link](https://doi.org/10.1126/science.aea1259)
5. Risgaard RD, et al. Molecular and cellular processes disrupted in the early postnatal Down syndrome prefrontal cortex. *Science*. 2026;392:eaea1549. [Link](https://doi.org/10.1126/science.aea1549)

The [`fgsea`](https://bioconductor.org/packages/fgsea/) package is one source
of precomputed tables compatible with this plotting contract.
