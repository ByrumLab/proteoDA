# Prepare per-contrast model data for plotting

Internal function used to prepare a results data frame for both static
and interactive plots in reports.

## Usage

``` r
prep_plot_model_data(model_results, contrast)
```

## Arguments

- model_results:

  The results slot of a DAList object.

- contrast:

  The name of the contrast for which to prep the model data.

## Value

A data frame of model results for the given contrast.
