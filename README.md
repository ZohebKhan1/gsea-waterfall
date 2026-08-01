# gsea-waterfall

`gsea-waterfall` is a visualization package for Gene Set Enrichment Analysis
(GSEA) results. It provides R functions for ranked NES waterfall plots,
symmetric volcano plots, directional half-volcano plots, and cross-contrast NES
scatterplots.

GitHub Pages: `https://zohebkhan1.github.io/gsea-waterfall/`

## Download The Functions

Download the plotting functions and source them in R:

```r
# set the raw github function path
base_url <- "https://raw.githubusercontent.com/ZohebKhan1/gsea-waterfall/main/functions"

download.file(paste0(base_url, "/gsea_plot_utils.R"), "gsea_plot_utils.R")
download.file(paste0(base_url, "/gsea_waterfall_plot.R"), "gsea_waterfall_plot.R")
download.file(paste0(base_url, "/gsea_volcano_plot.R"), "gsea_volcano_plot.R")
download.file(paste0(base_url, "/gsea_half_volcano_plot.R"), "gsea_half_volcano_plot.R")
download.file(paste0(base_url, "/gsea_nes_scatter_plot.R"), "gsea_nes_scatter_plot.R")

source("gsea_plot_utils.R")
source("gsea_waterfall_plot.R")
source("gsea_volcano_plot.R")
source("gsea_half_volcano_plot.R")
source("gsea_nes_scatter_plot.R")
```

## Function Overview

### `plot_gsea_waterfall()`

Ranks enriched GO terms in one NES direction and colors manually defined
biological categories. This is the main ranked-pathway view for showing broad
positive or negative GSEA structure without reducing the result to only a small
top-term list.

<img src="tutorial/assets/figures/GSE122380_gsea_waterfall_cardiomyocyte_vs_mesoderm_positive.svg" alt="Positive NES waterfall" width="760">

**Figure 1.** Positive NES GO biological-process waterfall for Cardiomyocyte vs. Mesoderm.
Terms are ranked by NES, colored by manually defined cardiac-relevant GO term
groups, and labeled with selected high-priority terms.

<img src="tutorial/assets/figures/GSE122380_gsea_waterfall_cardiomyocyte_vs_mesoderm_negative.svg" alt="Negative NES waterfall" width="760">

**Figure 2.** Negative NES GO biological-process waterfall for the same contrast. The selected
groups highlight ribosomal, mitotic, and DNA-replication programs.

### `plot_gsea_volcano()`

Shows positive, negative, and non-significant GSEA terms in one symmetric NES
volcano plot. This view is useful when the user needs the full bidirectional
contrast in a single panel.

<img src="tutorial/assets/figures/GSE122380_gsea_volcano_cardiomyocyte_vs_mesoderm.svg" alt="Symmetric GSEA volcano" width="560">

**Figure 3.** Symmetric GO biological-process GSEA volcano for Cardiomyocyte vs. Mesoderm. NES
is shown on the x-axis, adjusted p-value on the y-axis, and point color marks
significant positive, significant negative, or non-significant enrichment.

### `plot_gsea_half_volcano()`

Expands one NES direction into a larger directional volcano panel. This is meant
for closer inspection of either enriched or depleted programs while retaining
the same biological grouping logic used by the waterfall plots.

<img src="tutorial/assets/figures/GSE122380_gsea_half_volcano_cardiomyocyte_vs_mesoderm_positive.svg" alt="Positive NES half-volcano" width="760">

**Figure 4.** Positive NES half-volcano for Cardiomyocyte vs. Mesoderm. The enlarged
directional view uses the same labeled positive terms as the positive NES
waterfall plot.

<img src="tutorial/assets/figures/GSE122380_gsea_half_volcano_cardiomyocyte_vs_mesoderm_negative.svg" alt="Negative NES half-volcano" width="760">

**Figure 5.** Negative NES half-volcano for Cardiomyocyte vs. Mesoderm. The selected terms
emphasize DNA replication, checkpoint, and ribosome-associated programs.

### `plot_gsea_nes_scatter()`

Compares NES values for matched GO terms across two GSEA contrasts. This is the
cross-contrast view for identifying GO terms that are shared, contrast-specific,
or directionally divergent.

<img src="tutorial/assets/figures/GSE122380_gsea_scatter_day9_vs_day6_x_day3_vs_day1_all_quadrants.svg" alt="Comparative NES scatterplot" width="900">

**Figure 6.** Comparative GO biological-process NES scatterplot for Day 9 vs. Day 6 and Day 3
vs. Day 1. The fitted line summarizes the NES relationship across plotted GO
terms, while point color indicates whether each GO term is significant in one,
both, or neither contrast.

### `read_gsea_result_csvs()`

Reads one or more GSEA CSV files and standardizes the columns used by the
plotting functions. It is a small input helper, not a plotting function.

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

## Waterfall Plot Rationale

The ranked GSEA waterfall visualization was heavily inspired by pathway-level
waterfall plots used in recent developmental and stem-cell genomics studies.
This project turns that visual style into a parameterizable, reusable R function
for precomputed GSEA results.

Examples include [Ciceri, _Nature_, 2024](https://doi.org/10.1038/s41586-023-06984-8),
[Xu, _Nature Cell Biology_, 2025](https://doi.org/10.1038/s41556-025-01751-5),
[Vuong, _Science_, 2026](https://doi.org/10.1126/science.aea1259), and
[Risgaard, _Science_, 2026](https://doi.org/10.1126/science.aea1549).

The advantage of these plots is breadth. Instead of showing only a small subset
of enriched terms, such as the top 10 GO terms, the waterfall view can rank up
to 200 GO terms by NES and color them by biological class. This gives a broader
overview of enriched and depleted pathways while still allowing selected terms
to be labeled directly.

## Reproduce The Tutorial

```r
source("R/01_generate_gsea_figs.R")
source("R/02_render_github_pgs.R")
```

The rendered site is written to `docs/index.html` for GitHub Pages.

## Function Reference

### `plot_gsea_waterfall()`

Ranks GO terms in one NES direction and colors selected biological categories.

| Parameter           | Meaning                                                            |
| :------------------ | :----------------------------------------------------------------- |
| `gsea_results`      | GSEA result table.                                                 |
| `term_col`          | Column containing GO term names.                                   |
| `nes_col`           | Column containing NES values.                                      |
| `pvalue_col`        | Column containing nominal GSEA p-values.                           |
| `padj_col`          | Column containing adjusted p-values.                               |
| `id_col`            | Optional column containing GO IDs.                                 |
| `direction`         | `"positive"` or `"negative"`.                                      |
| `top_n`             | Number of ranked terms to plot.                                    |
| `rank_by`           | Column used for ordering terms, usually `"NES"`.                   |
| `x_label`           | X-axis title.                                                      |
| `label_n`           | Number of automatic labels when `label_terms` is not supplied.     |
| `label_terms`       | Exact GO descriptions or GO IDs to label.                          |
| `term_groups`       | Named list of biological groups and matching text seeds or GO IDs. |
| `category_n`        | Expected number of groups in `term_groups`.                        |
| `term_group_colors` | Optional named colors for biological groups.                       |

### `plot_gsea_volcano()`

Plots positive, negative, and non-significant GSEA terms in one symmetric NES
volcano.

| Parameter              | Meaning                                             |
| :--------------------- | :-------------------------------------------------- |
| `gsea_results`         | GSEA result table.                                  |
| `term_col`             | Column containing GO term names.                    |
| `nes_col`              | Column containing NES values.                       |
| `pvalue_col`           | Column containing nominal GSEA p-values.            |
| `padj_col`             | Column containing adjusted p-values.                |
| `id_col`               | Optional column containing GO IDs.                  |
| `padj_cutoff`          | Adjusted p-value cutoff for significance.           |
| `label_n`              | Number of significant terms to label automatically. |
| `label_terms`          | Exact GO descriptions or GO IDs to label.           |
| `contrast_label`       | Short contrast name used in the x-axis title.       |
| `label_words_per_line` | Number of words per line in wrapped GO labels.      |
| `point_size`           | Point size.                                         |
| `label_size`           | GO label text size.                                 |

### `plot_gsea_half_volcano()`

Expands one NES direction into a directional volcano panel.

| Parameter           | Meaning                                                            |
| :------------------ | :----------------------------------------------------------------- |
| `gsea_results`      | GSEA result table.                                                 |
| `term_col`          | Column containing GO term names.                                   |
| `nes_col`           | Column containing NES values.                                      |
| `pvalue_col`        | Column containing nominal GSEA p-values.                           |
| `padj_col`          | Column containing adjusted p-values.                               |
| `id_col`            | Optional column containing GO IDs.                                 |
| `direction`         | `"positive"` or `"negative"`.                                      |
| `p_col`             | P-value column shown on the y-axis, usually `"padj"`.              |
| `padj_cutoff`       | Adjusted p-value cutoff for significance.                          |
| `label_terms`       | Exact GO descriptions or GO IDs to label.                          |
| `term_groups`       | Named list of biological groups and matching text seeds or GO IDs. |
| `term_group_colors` | Optional named colors for biological groups.                       |
| `inner_nes_limit`   | Inner NES boundary shown on the x-axis.                            |
| `contrast_label`    | Short contrast name used in the x-axis title.                      |

### `plot_gsea_nes_scatter()`

Compares NES values for matched GO terms across two GSEA contrasts.

| Parameter                | Meaning                                                 |
| :----------------------- | :------------------------------------------------------ |
| `gsea_x`, `gsea_y`       | GSEA result tables for the x- and y-axis contrasts.     |
| `x_label`, `y_label`     | Axis titles.                                            |
| `term_col`               | Column containing GO term names in both tables.         |
| `nes_col`                | Column containing NES values in both tables.            |
| `pvalue_col`             | Column containing nominal GSEA p-values in both tables. |
| `padj_col`               | Column containing adjusted p-values in both tables.     |
| `id_col`                 | Optional column containing GO IDs in both tables.       |
| `color_by`               | `"significance"` or `"term_group"`.                     |
| `x_name`, `y_name`       | Short contrast names used in legend labels.             |
| `label_terms`            | Exact GO descriptions or GO IDs to label.               |
| `include_nonsignificant` | Keep terms not significant in either contrast.          |
| `equal_axis_limits`      | Use matched x/y limits.                                 |
| `show_fit_line`          | Draw a linear best-fit line.                            |

## Notes On GSEA And fgsea

The example tables can come from any GSEA workflow that reports GO term names,
NES, nominal p-values, and adjusted p-values.

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
- Ciceri G, Baggiolini A, Cho HS, et al. An epigenetic barrier sets the timing
  of human neuronal maturation. Nature. 2024;626:881-890.
  https://doi.org/10.1038/s41586-023-06984-8
- Xu N, Cho HS, Hackland JOS, et al. Genome-wide CRISPR screen identifies Menin
  and SUZ12 as regulators of human developmental timing. Nature Cell Biology.
  2025;27:1411-1421. https://doi.org/10.1038/s41556-025-01751-5
- Vuong CK, Weber A, Seong P, et al. A single-cell multiomic analysis identifies
  molecular and gene-regulatory mechanisms dysregulated in developing Down
  syndrome neocortex. Science. 2026;392:eaea1259.
  https://doi.org/10.1126/science.aea1259
- Risgaard RD, et al. Molecular and cellular processes disrupted in the early
  postnatal Down syndrome prefrontal cortex. Science. 2026;392:eaea1549.
  https://doi.org/10.1126/science.aea1549
