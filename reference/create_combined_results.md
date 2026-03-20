# Combine statistical results

Internal function used to construct a data frame of combined results.

## Usage

``` r
create_combined_results(annotation, data, statlist)
```

## Arguments

- annotation:

  A data frame of annotation data for each gene/protein.

- data:

  A data frame of normalized intensity data for each sample.

- statlist:

  A list of per-contrast DE results.

## Value

A data frame of the combined results.
