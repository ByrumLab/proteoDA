# Prepare limma contrasts matrix

Create the contrasts matrix, for use in a limma model. This function
utilizes the function
[limma::makeContrasts](https://rdrr.io/pkg/limma/man/makeContrasts.html)
with a user provided file or vector of a list of comparisons. Note: The
label on the plots is defined by what is written in the contrast
statement prior to the equal sign.

## Usage

``` r
add_contrasts(DAList, contrasts_vector = NULL, contrasts_file = NULL)
```

## Arguments

- DAList:

  A DAList. Must have a non-empty statistical design.

- contrasts_vector:

  A vector of contrasts.

- contrasts_file:

  The path to the contrasts file listing the desired contrasts. Must be
  a .csv, .tsv, or .txt file.

## Value

A DAList with added contrasts associated with the limma design

## Examples

``` r
if (FALSE) { # \dontrun{
# An example of a .csv file with three comparisons
data -> add_contrasts(data, contrasts_file = "path/to/file.csv")

# file info
Treatment1_vs_Control= Treatment1 - Control
Treatment2_vs_Control= Treatment2 - Control
Treatment2_vs_Treatment1= Treatment2 - Treatment1

# An example using a vector
data <- add_contrasts(contrasts_vector = c("Treat1_vs_Control=Treat1-Control",
     "Treat2_vs_Control=Treat2-Control"))
} # }
```
