# Scree Plot for PCA

Generates a scree plot showing the proportion of variance for each
principal component.

## Usage

``` r
qc_pca_scree_plot7(pca, max_pc = NULL, text.sizes = c(16, 14, 14))
```

## Arguments

- pca:

  A prcomp object resulting from PCA analysis.

- max_pc:

  The maximum number of PCs to display in the scree plot. Default: NULL
  (shows all PCs)

- text.sizes:

  A numeric vector specifying text sizes for different plot elements:
  title, axis titles, axis labels. Default: c(16, 14, 14)

## Value

A list containing the scree plot and the data used for plotting.

## Details

This function visualizes the variance explained by each principal
component in the PCA analysis.

## See also

[`cli_abort`](https://cli.r-lib.org/reference/cli_abort.html)

## Examples

``` r
if (FALSE) { # \dontrun{
if(interactive()){
 # Example usage
 qc_pca_scree_plot7(pca)
 }
} # }
```
