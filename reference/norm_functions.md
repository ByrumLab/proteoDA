# Cyclic loess normalization (optionally within groups) Expects log2 data. If `groups` is given, normalizes each group independently.

Cyclic loess normalization (optionally within groups) Expects log2 data.
If `groups` is given, normalizes each group independently.

## Usage

``` r
log2Norm(dat)

medianNorm(logDat)

meanNorm(logDat)

vsnNorm(dat)

quantileNorm(logDat)

cycloessNorm(logDat, groups = NULL, method = "fast", ...)

rlrNorm(logDat)

giNorm(dat)
```
