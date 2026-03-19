
<!-- README.md is generated from README.Rmd. Please edit README.Rmd file -->

# proteoDA

<!-- badges: start -->

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.7962306.svg)](https://doi.org/10.5281/zenodo.7962306)
[![DOI](https://joss.theoj.org/papers/10.21105/joss.05184/status.svg)](https://doi.org/10.21105/joss.05184)
[![R-CMD-check](https://github.com/ByrumLab/proteoDA/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/ByrumLab/proteoDA/actions/workflows/R-CMD-check.yaml)
[![Codecov test
coverage](https://codecov.io/github/ByrumLab/proteoDA/branch/main/graph/badge.svg)](https://app.codecov.io/github/ByrumLab/proteoDA?branch=main)
<!-- badges: end -->

proteoDA is a streamlined, user-friendly R package for quantitative
proteomics data preprocessing and statistical analysis. The package uses
a custom S3 class `DAList` structure allowing for provenance tracking
and reproducibility.

proteoDA version 2.0 extends the previous version with additional
features such as:

1.  `align_data_and_metadata()` to check sample name, group, and data
    alignment
2.  `filter_proteins_per_contrast()` to apply filtering proteins
    separately for each comparison
3.  `perseus_impute()` to apply Perseus style missing value imputation
4.  `run_filtered_limma_analysis()` runs limma and also calculates the
    movingSD and logFC_zscore
5.  `write_norm_report()` has been updated to provide additional
    diagnostic plots
6.  `normalize_data()` now has an option for within group normalization
    with cyclic loess
7.  `write_limma_tables()` has been updated to include an overview
    worksheet and other customization options for the final Excel file.
8.  `write_limma_plots()` includes a new parameter to highlight proteins
    on the static volcano plots

## Installation

`proteoDA` is available for install from GitHub via the `remotes`
package.

### Install proteoDA version 1.0.0

``` r
install.packages("remotes")
remotes::install_github("ByrumLab/proteoDA@v1.0.0",
                        dependencies = TRUE,
                        build_vignettes = TRUE,
                        force = TRUE)
```

### Install proteoDA version 2.0.0

``` r

install.packages("devtools")
devtools::install_github("ByrumLab/proteoDA",
                         dependencies = TRUE,
                         build_vignettes = TRUE,
                         force = TRUE)
 # or
install.packages("remotes")
remotes::install_github("ByrumLab/proteoDA@v2.0.0",
                        dependencies = TRUE,
                        build_vignettes = TRUE,
                        force = TRUE)
```

Using the `build_vignettes = TRUE` argument will build the tutorial
vignettes when you install, which you can access by running
`browseVignettes(package = "proteoDA")`. However, building the vignettes
requires some additional software dependencies. If you run into issues
when installing the vignettes, you can set `build_vignettes = FALSE`.

Once `proteoDA` is installed, load it into R:

``` r
library(proteoDA)
```

## Workflow

<figure>
<img src="./data-raw/proteoDA_flowchart.png?raw=true"
alt="proteoDA workflow flowchart" />
<figcaption aria-hidden="true">proteoDA workflow flowchart</figcaption>
</figure>

## Example pipeline

For a detailed explanation of the pipeline, check out the
`proteoDA_v2.0_workflow` vignette.

1.  Setup parameters in the `proteoDA_params.R` file and prepare data to
    load into the DAList object
2.  Run `proteoDA_analysis_version2.0_2026.03.19.R` analysis.

## Getting help

For general help on using `proteoDA`, check out the vignettes by running
`browseVignettes(package = "proteoDA")`. If you did not build the
vignette upon install, you can find a pre-built `.html` version of the
vignette in the `vignettes` folder on
[GitHub](https://github.com/ByrumLab/proteoDA) in the `release/1.0`
branch. Additional information can be found in the documentation for
each function. If you need further assistance, [file an issue on
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
