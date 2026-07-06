# gsea-waterfall

Sourceable R functions for visualizing precomputed Gene Set Enrichment Analysis
(GSEA) results. The repository is intentionally not an R package; users can
download and source the plotting functions directly.

GitHub Pages: `https://zohebkhan1.github.io/gsea-waterfall/`

## Repository Layout

```text
.
├── README.md
├── R/
│   ├── 01_generate_tutorial_figures.R
│   └── render_tutorial_site.R
├── docs/
│   └── index.html
├── functions/
│   ├── gsea_half_volcano_plot.R
│   ├── gsea_nes_scatter_plot.R
│   ├── gsea_volcano_plot.R
│   └── gsea_waterfall_plot.R
└── tutorial/
    ├── data/
    ├── assets/
    │   ├── figures/
    │   └── fonts/
    ├── style.css
    └── tutorial.Rmd
```

## Get The Functions

Download only the plotting functions from GitHub:

```r
# set the raw github function path
base_url <- "https://raw.githubusercontent.com/ZohebKhan1/gsea-waterfall/main/functions"

download.file(paste0(base_url, "/gsea_waterfall_plot.R"), "gsea_waterfall_plot.R")
download.file(paste0(base_url, "/gsea_volcano_plot.R"), "gsea_volcano_plot.R")
download.file(paste0(base_url, "/gsea_half_volcano_plot.R"), "gsea_half_volcano_plot.R")
download.file(paste0(base_url, "/gsea_nes_scatter_plot.R"), "gsea_nes_scatter_plot.R")

source("gsea_waterfall_plot.R")
source("gsea_volcano_plot.R")
source("gsea_half_volcano_plot.R")
source("gsea_nes_scatter_plot.R")
```

## Function Overview

`plot_gsea_waterfall()` ranks enriched GO terms in one NES direction and supports
manual biological category coloring.

`plot_gsea_volcano()` shows positive, negative, and non-significant GSEA terms in
one symmetric NES volcano plot.

`plot_gsea_half_volcano()` expands one NES direction into a larger directional
volcano panel.

`plot_gsea_nes_scatter()` compares NES values for matched GO terms across two
GSEA contrasts and can color terms by significance class or biological category.

The helper `read_gsea_result_csvs()` reads one or more GSEA CSV files and
standardizes the columns used by the plotting functions.

## Required Input

The public plotting workflow requires a precomputed GSEA result table with these
columns:

| Column           | Required | Meaning                     |
| :--------------- | :------- | :-------------------------- |
| `go_description` | yes      | GO term name                |
| `NES`            | yes      | normalized enrichment score |
| `pval`           | yes      | nominal GSEA p-value        |
| `padj`           | yes      | adjusted p-value            |

GO IDs, ontology labels, leading-edge genes, and term sizes can be present, but
they are optional for the plotting workflow.

## Reproduce The Tutorial

```r
source("R/01_generate_tutorial_figures.R")
source("R/render_tutorial_site.R")
```

The rendered site is written to `docs/index.html` for GitHub Pages.

## Figures

### Figure 1. Positive NES waterfall

<img src="tutorial/assets/figures/GSE122380_gsea_waterfall_cardiomyocyte_vs_mesoderm_positive.svg" alt="Positive NES waterfall" width="760">

**Figure 1.** Positive NES GO biological-process waterfall for Cardiomyocyte vs. Mesoderm.
Terms are ranked by NES, colored by manually defined cardiac-relevant GO term
groups, and labeled with selected high-priority terms.

### Figure 2. Negative NES waterfall

<img src="tutorial/assets/figures/GSE122380_gsea_waterfall_cardiomyocyte_vs_mesoderm_negative.svg" alt="Negative NES waterfall" width="760">

**Figure 2.** Negative NES GO biological-process waterfall for the same contrast. The selected
groups highlight ribosomal, mitotic, and DNA-replication programs.

### Figure 3. Symmetric GSEA volcano

<img src="tutorial/assets/figures/GSE122380_gsea_volcano_cardiomyocyte_vs_mesoderm.svg" alt="Symmetric GSEA volcano" width="560">

**Figure 3.** Symmetric GO biological-process GSEA volcano for Cardiomyocyte vs. Mesoderm. NES
is shown on the x-axis, adjusted p-value on the y-axis, and point color marks
significant positive, significant negative, or non-significant enrichment.

### Figure 4. Positive NES half-volcano

<img src="tutorial/assets/figures/GSE122380_gsea_half_volcano_cardiomyocyte_vs_mesoderm_positive.svg" alt="Positive NES half-volcano" width="760">

**Figure 4.** Positive NES half-volcano for Cardiomyocyte vs. Mesoderm. The enlarged
directional view uses the same labeled positive terms as the positive NES
waterfall plot.

### Figure 5. Negative NES half-volcano

<img src="tutorial/assets/figures/GSE122380_gsea_half_volcano_cardiomyocyte_vs_mesoderm_negative.svg" alt="Negative NES half-volcano" width="760">

**Figure 5.** Negative NES half-volcano for Cardiomyocyte vs. Mesoderm. The selected terms
emphasize DNA replication, checkpoint, and ribosome-associated programs.

### Figure 6. Comparative NES scatterplot

<img src="tutorial/assets/figures/GSE122380_gsea_scatter_day9_vs_day6_x_day3_vs_day1_all_quadrants.svg" alt="Comparative NES scatterplot" width="900">

**Figure 6.** Comparative GO biological-process NES scatterplot for Day 9 vs. Day 6 and Day 3
vs. Day 1. The fitted line summarizes the NES relationship across plotted GO
terms, while point color indicates whether each GO term is significant in one,
both, or neither contrast.

## Notes On GSEA And fgsea

These functions plot GSEA results; they do not run DESeq2, build ranked gene
lists, or calculate enrichment. The example tables can come from any GSEA
workflow that reports GO term names, NES, nominal p-values, and adjusted
p-values.

The original GSEA method is described by Subramanian et al. The `fgsea`
Bioconductor package provides a fast implementation for preranked GSEA and can
produce input tables for these plots.

## References

- Subramanian A, Tamayo P, Mootha VK, et al. Gene set enrichment analysis: a
  knowledge-based approach for interpreting genome-wide expression profiles.
  PNAS. 2005. https://www.pnas.org/doi/10.1073/pnas.0506580102
- Korotkevich G, Sukhov V, Budin N, Shpak B, Artyomov MN, Sergushichev A. Fast
  gene set enrichment analysis. https://www.biorxiv.org/content/10.1101/060012v3
- Bioconductor `fgsea`: Fast Gene Set Enrichment Analysis.
  https://bioconductor.org/packages/fgsea/
- The Gene Ontology Consortium. The Gene Ontology resource.
  https://geneontology.org/
