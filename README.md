
<!-- README.md is generated from README.Rmd. Please edit README.Rmd file -->

# proteoDA

proteoDA is a streamlined, user-friendly R package for the analysis of
high resolution mass spectrometry protein data. The package uses a
custom S3 class that keeps the R objects consistent across the pipeline
and is easily pipe-able, so minimal R knowledge is required.

## Installation

`proteoDA` is not yet on CRAN, but it is available for install from
GitHub via the `devtools` package.

Install `devtools` if you haven’t already:

``` r
install.packages("devtools")
```

Then install and load `proteoDA`:

``` r
devtools::install_github(ByrumLab/proteoDA)
```

## Example pipeline

Add workflow figure

## Prepare data to import into `proteoDA`

`proteoDA` is built around a `DAList` object: a list which holds all the
information for a single quantitative proteomics experiment. Most
`proteoDA` functions take in a `DAList` as input and output a `DAList`.
A `DAList` is a list with 7 slots:

1.  **data**- MS intensity data where each row is a protein and each
    column is a sample. It must be numeric.
2.  **annotation**- Protein annotation information such as the accession
    id, description, gene symbol, etc. “uniprot_id” is a required
    column. It must have the same number of rows as data. Each column is
    a separate piece of annotation info.
3.  **metadata**- A data frame of sample information such as
    sample_name, group, batch, gender, paired, etc. The row names must
    match the column names in data where there is one row per sample.
4.  **design**- A list which holds information on the statistical
    design: the design formula and matrix, and contrasts.
5.  **eBayes_fit**- The model fit object from running `limma` models on
    the data.
6.  **results**- Statistical tables of differential abundance (DA)
    results for each of the statistical terms/contrasts analyzed.
7.  **tags**- A “miscellaneous” slot used to keep track of information
    on filtering, normalization, function parameters, etc.

The first three slots are required by the user to create the `DAList`
object. After creating a new object with the `DAList()` function,
further`proteoDA` functions are used to process the data and add the
rest of the slots.

``` r
# example data has 4 protein annotation columns and 17 sample columns
input <- read.csv("vignettes/DIA_data.csv")  # subset the data and annotation
data <- input[,5:21]                         # subset to extract sample columns
anno <- input[,1:4]                          # subset to extract annotation
meta <- read.csv("vignettes/metafile.csv")   # sample information
row.names(meta) <- meta$data_column_name     # "data_column_name" is first column of metafile.csv
```

All of the functions in the pipeline check for a proper structure of the
DAList object, both when it is input into the function and before it
returns the result.

## For more details, please check out the vignette.
