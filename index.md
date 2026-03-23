# proteoDA

proteoDA is a streamlined, user-friendly R package for quantitative
proteomics data preprocessing and statistical analysis. The package uses
a custom S3 class `DAList` structure allowing for provenance tracking
and reproducibility.

proteoDA version 2.0 extends the previous version with additional
features such as:

1.  [`align_data_and_metadata()`](https://byrumlab.github.io/proteoDA/reference/align_data_and_metadata.md)
    to check sample name, group, and data alignment
2.  [`filter_proteins_per_contrast()`](https://byrumlab.github.io/proteoDA/reference/filter_proteins_per_contrast.md)
    to apply filtering proteins separately for each comparison
3.  [`perseus_impute()`](https://byrumlab.github.io/proteoDA/reference/perseus_impute.md)
    to apply Perseus style missing value imputation
4.  [`run_filtered_limma_analysis()`](https://byrumlab.github.io/proteoDA/reference/run_filtered_limma_analysis.md)
    runs limma and also calculates the movingSD and logFC_zscore
5.  [`write_norm_report()`](https://byrumlab.github.io/proteoDA/reference/write_norm_report.md)
    has been updated to provide additional diagnostic plots
6.  [`normalize_data()`](https://byrumlab.github.io/proteoDA/reference/normalize_data.md)
    now has an option for within group normalization with cyclic loess
7.  [`write_limma_tables()`](https://byrumlab.github.io/proteoDA/reference/write_limma_tables.md)
    has been updated to include an overview worksheet and other
    customization options for the final Excel file.
8.  [`write_limma_plots()`](https://byrumlab.github.io/proteoDA/reference/write_limma_plots.md)
    includes a new parameter to highlight proteins on the static volcano
    plots

## Installation

`proteoDA` is available for install from GitHub via the `remotes`
package.

### Install proteoDA version 1.0.0

``` r
install.packages("remotes")
remotes::install_github("ByrumLab/proteoDA@v1.0.0",
                        dependencies = TRUE)
```

### Install proteoDA version 2.0.0

``` r

install.packages("devtools")
devtools::install_github("ByrumLab/proteoDA",
                         dependencies = TRUE)
 # or
install.packages("remotes")
remotes::install_github("ByrumLab/proteoDA@v2.0.0",
                        dependencies = TRUE)
```

Once `proteoDA` is installed, load it into R:

``` r
library(proteoDA)
```

## Workflow

![proteoDA workflow
flowchart](./data-raw/proteoDA_flowchart.png?raw=true)

proteoDA workflow flowchart

## Example pipeline

For a detailed explanation of the pipeline, check out the
`proteoDA_v2.0_workflow` vignette.

1.  Setup parameters in the `proteoDA_params.R` file and prepare data to
    load into the DAList object
2.  Run `proteoDA_analysis_version2.0_2026.03.19.R` analysis.

## Getting help

Documentation and tutorials are available on the package website:
<https://byrumlab.github.io/proteoDA> Additional information can be
found in the documentation for each function. If you need further
assistance, [file an issue on
GitHub](https://github.com/ByrumLab/proteoDA/issues).

## Reporting issues

If you find any bugs or unexpected behaviors, [file an issue on
GitHub](https://github.com/ByrumLab/proteoDA/issues). It is helpful if
you can include a minimal reproducible example (reprex) that triggers
the issue, check out the `reprex` [R
package](https://reprex.tidyverse.org/) for more information and tools
on creating reproducible examples.

## Contributing

We welcome code contributions from users. To contribute, [open a pull
request](https://github.com/ByrumLab/proteoDA/pulls) against the main
branch. Please note that the proteoDA project is released with a
[Contributor Code of
Conduct](https://contributor-covenant.org/version/2/1/CODE_OF_CONDUCT.html).
By contributing to this project, you agree to abide by its terms.
