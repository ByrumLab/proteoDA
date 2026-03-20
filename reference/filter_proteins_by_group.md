# Filter protein data by number of quantified samples in a group

This function is used to remove proteins from a DAList, filtering out
proteins based on levels of missing values in the "data" data frame of
the DAList. The grouping_column must be a column in the metadata of the
DAList, which lists the group membership for each sample. The min_reps
and min_groups arguments that determine the number of replicates/samples
per group (min_reps) and number of groups (min_groups) in which a
protein must have non-missing intensity values in order to be retained.
This function assumes that all missing values are encoded as NA. See
[`zero_to_missing`](https://byrumlab.github.io/proteoDA/reference/missing_data.md)
and
[`missing_to_zero`](https://byrumlab.github.io/proteoDA/reference/missing_data.md)
for helper functions to convert missing values to and from 0.

## Usage

``` r
filter_proteins_by_group(
  DAList,
  min_reps = NULL,
  min_groups = NULL,
  grouping_column = "group"
)
```

## Arguments

- DAList:

  A DAList object to be filtered.

- min_reps:

  The minimum number of replicates/samples within a group that need to
  have a non-missing intensity value for a given protein in order for
  that protein to be considered as quantified within a group.

- min_groups:

  The minimum number of groups that must have at least min_reps non-zero
  samples for a given protein to be retained.

- grouping_column:

  The name of the column in the metadata which provides the group
  membership for each sample. Default is "group".

## Value

A DAList, with proteins that are not present in sufficient samples and
groups removed.

## Examples

``` r
if (FALSE) { # \dontrun{
  # Suppose the DAList contains data from 20 samples across 4
  # experimental groups (5 samples per group), with the group membership
  # listed in a column named "group"

  # Strict filtering:
  # no missing data
  # Proteins must be present in all samples in all groups
  filtered <- filter_proteins_by_group(DAList,
                                       min_reps = 5,
                                       min_groups = 4,
                                       grouping_column = "group")
  # Lax filtering:
  # protein must be present in at least one sample in each group
  filtered <- filter_proteins_by_group(DAList,
                                       min_reps = 1,
                                       min_groups = 4,
                                       grouping_column = "group")

 # Filtering functions can be chained together
 filtered <- DAList |>
   filter_proteins_by_annotation(!grepl(pattern = "keratin",
                                        x = protein_name)) |>
   filter_proteins_by_group(min_reps = 1,
                            min_groups = 4,
                            grouping_column = "group")
} # }
```
