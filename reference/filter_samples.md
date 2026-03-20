# Filter samples from a DAList

This function is used to remove samples from a DAList, filtering using
data in the metadata data frame of the DAList. Samples which do not
produce a value of TRUE for the supplied condition are removed from the
data, annotation, and metadata slots of the DAList. If condition
evaluates to NA, the function will return an error.

## Usage

``` r
filter_samples(DAList, condition)
```

## Arguments

- DAList:

  A DAList object to be filtered.

- condition:

  An expression that returns a logical value, defined in terms of
  variables present in the metadata data frame of the supplied DAList.
  Samples are kept if the condition is TRUE for that sample.

## Value

A DAList, with samples that do not meet the condition removed.

## Examples

``` r
if (FALSE) { # \dontrun{
# Suppose the DAList$metadata data frame contains three columns:
# sample_ID = An alpha-numeric ID uniquely identifying a sample
# treatment = A character listing the treatment group a sample belongs to
# year = A numeric value listing the year a sample was collected

# Remove a specific sample by ID
filtered <- filter_samples(DAList,
                           sample_ID != "abc123")

# Keep any sample which contains the string "control" in the treatment:
filtered <- filter_samples(DAList,
                           grepl(pattern = "control",
                                 x = treatment))

# Remove any sample from before 2010
filtered <- filter_samples(DAList,
                           year >= 2010)


# Filtering functions can be chained together
filtered <- DAList |>
  filter_samples(grepl(pattern = "control",
                       x = treatment)) |>
  filter_samples(year >= 2010)
} # }
```
