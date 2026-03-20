# Extract differential abundance results from a model fit

Extracts tables of statistical results from the model fit created by
[`fit_limma_model`](https://byrumlab.github.io/proteoDA/reference/fit_limma_model.md)
and add them to the results slot of the DAList. The results are a list
of tables, one for each contrast/term specific in the statistical
design. In models with intercepts these are ignored by default, but can
be output by setting extract_intercept to TRUE. See
[`limma::decideTests`](https://rdrr.io/pkg/limma/man/decideTests.html)
and [`limma::topTable`](https://rdrr.io/pkg/limma/man/toptable.html) for
information on the statistical results.

## Usage

``` r
extract_DA_results(
  DAList,
  pval_thresh = 0.05,
  lfc_thresh = 1,
  adj_method = "BH",
  extract_intercept = F
)
```

## Arguments

- DAList:

  A DAList, which must contain a model fit.

- pval_thresh:

  The p-value threshold used to determine significance (significant when
  p \< pval_thresh). Default is 0.05.

- lfc_thresh:

  The logFC threshold used to determine significance (significant when
  \|logFC\| \> lfc.tresh). Default is 1. LogFC are base 2.

- adj_method:

  The method used for adjusting P-values. Possible values are "none",
  "BH", "BY", and "holm". Default is "BH", for the Benjamini-Hochberg
  correction. See
  [`stats::p.adjust`](https://rdrr.io/r/stats/p.adjust.html) for
  details.

- extract_intercept:

  For models with an intercept term, should results for the intercept be
  extracted? Default if FALSE.

## Value

A DAList object, with differential abundance results added to the
results slot.

## Examples

``` r
if (FALSE) { # \dontrun{
# Using default thresholds and p-value adjustment
results <- extract_DA_results(DAList)

# Relax significance and log fold-change thresholds
results <- extract_DA_results(DAList,
                              pval_thresh = 0.1,
                              lfc_thresh = 0)

# Include intercept term in results
results <- extract_DA_results(DAList,
                              extract_intercept = T)
} # }
```
