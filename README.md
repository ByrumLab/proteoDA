
<!-- README.md is generated from README.Rmd. Please edit that file -->

# proteoDA

<!-- badges: start -->
<!-- badges: end -->

proteoDA is a public-friendly R package for the analysis of high resolution
mass spectrometry protein data. The package utilizes a custom S3 class
that keeps the R objects consistent across the pipeline and is easily
chain/pipe-able.

## Installation

The Github repository where the package is stored is public. You will need the
`devtools` package installed and the name of the repository. 

``` r
install.packages("devtools")
library(devtools)
devtools::install_github(ByrumLab/proteoDA)
                         
```
## Example pipeline

Add workflow figure

## Prepare data to import into `proteoDA`

`DAList` holds all the information for a single quantitative proteomics
experiment and all our functions take it as input. It is a list with 7
slots:

1)  <b>data-</b> MS intensity data where each row is a protein and each column
    is a sample. It must be numeric. 
2)  <b>annotation-</b> Protein annotation information such as the accession id, 
    description, gene symbol, etc. "uniprot_id" is a required column. 
    It must have the same number of rows as data. 
    Each column is a separate piece of annotation info.
3)  <b>metadata-</b> A data frame of sample information such as sample_name,
    group, batch, gender, paired, etc.  
    The rownames must match the column names in data where 
    there is one row per sample.
4)  <b>design-</b> A list which holds information on the statistical design:
    the design formula and matrix, and contrasts.
5)  <b>eBayes_fit</b>- The model fit object from running `limma` models on
    the data.
6)  <b>results-</b> Statistical tables of differential abundance (DA) results 
    for each of the statistical terms/contrasts analyzed.
7)  <b>tags-</b> A “miscellaneous” slot used to keep track of information
    on filtering, normalization, function parameters, etc.

The first three slots are required by the user to create the `DAList` object.
After which, `proteoDA` functions are used to process the data and add
the rest of the slots. 

```r

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
