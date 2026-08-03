# gsea-waterfall

`gsea-waterfall` visualizes precomputed Gene Set Enrichment Analysis (GSEA)
results as ranked NES waterfall plots, volcano plots, and cross-contrast
comparative NES scatterplots. The functions return `ggplot` objects. They are a
visualization layer applied after GSEA has been performed; they do not run GSEA
or write files.

[GSEA visualization examples shown here](https://zohebkhan1.github.io/gsea-waterfall/)

## What is GSEA?

Gene Set Enrichment Analysis evaluates whether the genes in a predefined gene
set are concentrated toward the top or bottom of a ranked gene list. Rather
than testing genes one at a time, GSEA summarizes coordinated changes across
biological processes or pathways. A typical result table contains one row per
gene set with a term description, normalized enrichment score (NES), nominal
p-value, and adjusted p-value.

GSEA can be performed with any appropriate tool, including `fgsea`,
`clusterProfiler`, or another package. This repository accepts completed GSEA
result tables and does not select gene sets, calculate enrichment statistics,
or perform multiple-testing correction. At minimum, a waterfall needs a term
description and NES. The shared reader in this repository also validates
nominal and adjusted p-values so the same standardized table can support
volcano and comparative plots. If your source uses different column names,
map them when reading the table as shown below.

NES is the primary direction and ranking measure: positive and negative values
indicate enrichment on opposite sides of the ranked gene list. Adjusted
p-values describe statistical significance and are used by the volcano and
comparative plots for coloring and display decisions. NES and adjusted
p-values answer different questions and should not be interpreted
interchangeably. The plots use these fields differently: waterfalls rank by
NES or a p-value, volcanoes show NES against adjusted p-value, half-volcanoes
can show nominal or adjusted p-value, and comparative scatterplots use NES with
adjusted-p-value significance classes.

## Repository layout

```text
functions/          downloadable plotting functions and shared utilities
R/                  figure generation, site rendering, and contract checks
tutorial/data/      fixed precomputed example GSEA tables
tutorial/assets/    generated tutorial figures and bundled fonts
tutorial/tutorial.Rmd
docs/               generated GitHub Pages site
```

## Why use GSEA waterfall plots?

GSEA results are often presented as a running-score enrichment plot when one
pathway is the focus, or as a dotplot/barplot-style figure in which NES or
adjusted p-value is encoded on the x-axis and other GSEA metrics are shown by
point size or color. These compact plots are useful for a short list of terms,
but they become crowded when many results are displayed and commonly focus on
only the top 10–20 terms.

The top-ranked terms can also be functionally redundant: related GO parent and
child terms may appear next to one another because they represent overlapping
gene sets. A waterfall plot does not remove this redundancy automatically, but
it provides a broader view of the ranked enrichment landscape, often allowing
100–200 terms to be inspected rather than only the first few rows.

The waterfall functions also allow displayed GO terms to be grouped using
caller-supplied keyword or identifier categories. In the cardiomyocyte iPSC
time-course bulk RNA-seq example, these annotations highlight processes such
as cardiac development, metabolism, ribosomal activity, mitosis, and DNA
replication. They are visualization annotations, not additional statistical
tests; terms remain neutral by default unless `term_groups` is supplied.

The initial idea for the waterfall presentation was informed by pathway-level
figures from Ciceri et al. (2024), Xu et al. (2025), Vuong et al. (2026), and
Risgaard et al. (2026). Full citations are provided at the end of this README.

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
