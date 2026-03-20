# Summarize the number of DA proteins in a contrast

Internal function to summarize the number of DE genes/proteins for a
given contrast.

## Usage

``` r
summarize_contrast_DA(contrast_name, contrast_res_list)
```

## Arguments

- contrast_name:

  The name of the contrast to summarize.

- contrast_res_list:

  A list of per-contrast DA results.

## Value

A data frame summarizing differential abundance for the given contrast.
