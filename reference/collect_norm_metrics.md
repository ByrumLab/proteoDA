# Collect basic normalization metrics from a normList and groups

Computes summary metrics for each method in a `normList`:

- PCV: pooled coefficient of variation (computed on linear scale)

- PMAD: pooled median absolute deviation (log2 scale)

- PEV: pooled estimate of variance (log2 scale)

- COR: mean pairwise Pearson correlation within groups (log2 scale)

## Usage

``` r
collect_norm_metrics(normList, groups, contrast = NA_character_)
```

## Arguments

- normList:

  Named list of matrices (methods x samples). Names are method labels.

- groups:

  Character vector of group labels aligned to columns of matrices.

- contrast:

  Character scalar to label the contrast in the output.

## Value

data.frame with columns: contrast, method, metric, value
