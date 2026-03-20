# Align and validate sample correspondence between data and metadata

In many quantitative proteomics workflows, mass spectrometry sample
files are named using numeric suffixes (e.g. `Sample1`, `Sample2`, ...,
`Sample10`). When read into R or other languages, these can be
automatically sorted lexicographically (`Sample1`, `Sample10`,
`Sample11`, `Sample2`, ...), causing incorrect sample-to-group mapping
and misaligned metadata.

This helper ensures that the `data` and `metadata` objects are correctly
matched and ordered before creating a `DAList`. It automatically fixes
ordering issues, groups replicates of the same condition together, and
performs natural ("human") sorting of sample names containing numbers.

## Usage

``` r
align_data_and_metadata(
  data,
  metadata,
  sample_col = NULL,
  group_col = NULL,
  prefer_group_blocks = TRUE,
  strict = FALSE
)
```

## Arguments

- data:

  A numeric matrix or data frame of protein intensities, where rows are
  proteins and columns are samples.

- metadata:

  A data frame of sample metadata, one row per sample.

- sample_col:

  Optional name of the column in `metadata` that contains the sample
  identifiers. If NULL, row names are used.

- group_col:

  Optional column name in `metadata` used to group replicates together
  (e.g. "group", "condition").

- prefer_group_blocks:

  Logical; if TRUE (default) and `group_col` is supplied, both `data`
  and `metadata` will be reordered to keep replicates of the same group
  contiguous and sorted naturally. If FALSE, only the metadata will be
  reordered to match data columns.

- strict:

  Logical; if TRUE, stops with an error if sample orders differ. If
  FALSE (default), automatically reorders when sets match.

## Value

A list with three elements:

- data:

  The reordered numeric intensity data.

- metadata:

  The reordered sample metadata.

- changes:

  A data frame summarizing the old vs new order of samples.

## Details

The function performs three major checks:

1.  Ensures the same sample IDs appear in both `data` and `metadata`.

2.  Detects and removes any ordering mismatches.

3.  (Optionally) groups replicates of the same experimental condition
    together using natural sorting.

If sample names differ only by order, the function reorders them safely.
If samples are missing or duplicated, it stops with a descriptive error.

## Examples

``` r
if (FALSE) { # \dontrun{
# Example data with misordered samples
dat <- data.frame(
  Sample1  = c(10, 20, 30),
  Sample10 = c(11, 21, 31),
  Sample2  = c(12, 22, 32)
)
rownames(dat) <- paste0("Prot", 1:3)

meta <- data.frame(
  sample = c("Sample1", "Sample2", "Sample10"),
  group  = c("A", "A", "B")
)
rownames(meta) <- meta$sample

# Align automatically, grouping by condition
aligned <- align_data_and_metadata(
  data = dat,
  metadata = meta,
  sample_col = "sample",
  group_col = "group",
  prefer_group_blocks = TRUE
)

# Check results
aligned$changes
head(aligned$data)
head(aligned$metadata)
} # }
```
