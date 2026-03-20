# Fit the limma differential abundance model

Fits the limma differential abundance model to the intensity data,
following the specified design and (optional) contrast matrices. When a
random factor is included, uses
[`limma::duplicateCorrelation`](https://rdrr.io/pkg/limma/man/dupcor.html)
to estimate the intra-block correlation within groups. Uses
[`limma::lmFit`](https://rdrr.io/pkg/limma/man/lmFit.html) to fit the
initial model, optionally re-parameterizes the results in terms of
contrasts with
[`limma::contrasts.fit`](https://rdrr.io/pkg/limma/man/contrasts.fit.html),
and then recomputes moderated statistics following limma's empirical
Bayes model with
[`limma::eBayes`](https://rdrr.io/pkg/limma/man/ebayes.html).

## Usage

``` r
fit_limma_model(DAList)
```

## Arguments

- DAList:

  A DAList, which must contain a statistical design.

## Value

A DAList object, with the model fit added in the eBayes_fit slot.

## Examples

``` r
if (FALSE) { # \dontrun{
model_fit <- fit_limma_model(DAList)
} # }
```
