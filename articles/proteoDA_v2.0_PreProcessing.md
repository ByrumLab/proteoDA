# proteoDA v2.0 Preprocessing

## Introduction

Preprocessing is a critical first step in any proteomics differential
abundance (DA) workflow. In `proteoDA`, we separate preprocessing
(cleaning and filtering) from modeling (limma-based statistical
analysis) so that you can:

- Explicitly document which samples and proteins were removed and why
- Apply different filtering rules for all groups vs contrast-specific
  analyses
- Track preprocessing decisions via the DAList\$tags slot for
  reproducibility

This vignette focuses on preprocessing functions that operate on a
`DAList`:

- Converting embedded zeros to missing values
- Filtering samples using metadata
- Filtering proteins using annotation
- Filtering proteins using missingness patterns across groups
- Filtering proteins per contrast using a contrasts file

Normalization, imputation, and limma designs/statistical thresholds are
covered in separate vignettes:

- **Overall workflow:**
  [`vignette("proteoDA_v2.0_workflow", package = "proteoDA")`](https://byrumlab.github.io/proteoDA/articles/proteoDA_v2.0_workflow.md)
- **Normalization and imputation:**
  `vignette("proteoDA_v2.0_Norm_Impute", package = "proteoDA")`
- **Limma designs and statistical thresholds:**
  [`vignette("proteoDA_v2.0_design", package = "proteoDA")`](https://byrumlab.github.io/proteoDA/articles/proteoDA_v2.0_design.md)

Typical `proteoDA` preprocessing flow:

1.  Import data → create a `DAList`
2.  Convert zeros →
    [`zero_to_missing()`](https://byrumlab.github.io/proteoDA/reference/missing_data.md)
3.  Optional sample-level and annotation-based filtering
    ([`filter_samples()`](https://byrumlab.github.io/proteoDA/reference/filter_samples.md),
    [`filter_proteins_by_annotation()`](https://byrumlab.github.io/proteoDA/reference/filter_proteins_by_annotation.md))
4.  Protein-level missingness filtering
    ([`filter_proteins_by_group()`](https://byrumlab.github.io/proteoDA/reference/filter_proteins_by_group.md)
    or
    [`filter_proteins_by_proportion()`](https://byrumlab.github.io/proteoDA/reference/filter_proteins_by_proportion.md))
5.  Optional contrast-specific filtering
    ([`filter_proteins_per_contrast()`](https://byrumlab.github.io/proteoDA/reference/filter_proteins_per_contrast.md))
6.  Normalize / impute → fit limma models See the normalization vignette
    for details on available normalization strategies and diagnostics:
    `vignette("proteoDA_v2.0_Norm_Impute", package = "proteoDA")`
7.  Write tables and plots →
    [`write_limma_tables()`](https://byrumlab.github.io/proteoDA/reference/write_limma_tables.md)
    and
    [`write_limma_plots()`](https://byrumlab.github.io/proteoDA/reference/write_limma_plots.md)

------------------------------------------------------------------------

## Overview of preprocessing functions

### Summary table

| Function                      | Purpose                                                            | Key_Args                                                       | When_to_use                                                                   | Output                                                               |
|:------------------------------|:-------------------------------------------------------------------|:---------------------------------------------------------------|:------------------------------------------------------------------------------|:---------------------------------------------------------------------|
| zero_to_missing               | Convert embedded zeros to NA for missing-aware downstream steps    | \-                                                             | Early: almost always, before filtering and modeling                           | DAList with NAs representing missingness                             |
| missing_to_zero               | Convert NA back to 0 for export or external tools                  | \-                                                             | Late: before exporting to tools that require zeros                            | DAList with zeros instead of NAs                                     |
| filter_samples                | Remove samples based on metadata conditions (e.g., QC failures)    | condition (expression on metadata columns)                     | When specific samples should be dropped (e.g., outliers, failed runs)         | DAList with reduced sample set                                       |
| filter_proteins_by_annotation | Remove proteins based on annotation (e.g., contaminants, keratins) | condition (expression on annotation columns)                   | When specific proteins should be dropped (e.g., keratins, decoys, spikes)     | DAList with reduced protein set                                      |
| filter_proteins_by_group      | Global missingness filtering across groups using replicate counts  | min_reps, min_groups, grouping_column                          | When all groups should share a common missingness rule                        | DAList (globally filtered)                                           |
| filter_proteins_by_proportion | Global missingness filtering across groups using proportions       | min_prop, grouping_column                                      | When group sizes are imbalanced but you want comparable missingness per group | DAList (globally filtered)                                           |
| filter_proteins_per_contrast  | Contrast-specific filtering using non-missing counts in groups     | min_reps, require_both_groups, grouping_column, contrasts_file | When different contrasts need different missingness rules                     | DAList with per-contrast filtered protein sets and data_per_contrast |

Preprocessing options in proteoDA v2.0

### General usage patterns

All preprocessing functions accept and return a `DAList`, allowing you
to chain them:

``` r
library(proteoDA)

DAList_clean <- DAList_raw |>
            zero_to_missing() |>
            filter_samples(batch != "bad_batch") |>
            filter_proteins_by_annotation(!grepl("KERATIN", protein_name)) |>
            filter_proteins_by_group(
                min_reps       = 2,
                min_groups     = 2,
                grouping_column = "group"
)
```

Most functions also update `DAList$tags` so decisions can be inspected
later.

------------------------------------------------------------------------

## Handling zeros and missing values

### zero_to_missing()

Embedded zeros often represent missing intensities rather than true zero
abundance. Statistical modeling in `proteoDA` assumes missingness is
encoded as `NA`.

Typical usage:

``` r
DAList <- zero_to_missing(DAList)
```

When to use

- Early in the pipeline, before filtering by missingness or running
  limma
- Any time upstream tools encode absence as 0 rather than NA

What it does

- Replaces 0 values in `DAList$data` with NA
- Leaves metadata and annotation unchanged
- Optionally records a tag in `DAList$tags$zero_to_missing`
  (implementation specific)

### missing_to_zero()

In some cases you may need to convert NA back to 0, e.g., when:

- Exporting to tools that cannot handle NA, such as PCA.
- Re-using data in pipelines that assume zeros

Typical usage:

``` r
DAList_export <- missing_to_zero(DAList)
```

You should not run limma on zero-filled data; use zeros only for export
or visualization where appropriate.

------------------------------------------------------------------------

## Sample-level filtering

### filter_samples(): drop samples using metadata

[`filter_samples()`](https://byrumlab.github.io/proteoDA/reference/filter_samples.md)
removes samples from a `DAList` based on logical conditions evaluated on
`DAList$metadata`.

Function recap

``` r
DAList <- filter_samples(DAList,
                  condition)
```

- condition is an expression using columns from DAList\$metadata
- Samples are kept where condition evaluates to TRUE

Example: remove a specific sample

``` r
DAList <- filter_samples(DAList,
                  sample_ID != "Sample_001")
```

Example: keep only “control” and “treated” groups

``` r
DAList <- filter_samples(DAList, 
                         group %in% c("Control", "Treated"))
```

Example: remove all samples in a group

``` r
DAList <- filter_samples(DAList,
                         group != "Pool")
```

Example: remove pre-2018 samples

``` r
DAList <- filter_samples(DAList,
                      year >= 2018)
```

Behavior and messages

- Updates `DAList$metadata` and `DAList$data` to keep only selected
  samples (`DAList$data` columns are re-subset to the remaining row
  names in `metadata`)

- Prints how many samples were removed and shows the removed metadata

- Errors if the filter produces any NA values (to avoid silent
  mis-filtering)

------------------------------------------------------------------------

## Protein-level filtering using annotation

### filter_proteins_by_annotation()

[`filter_proteins_by_annotation()`](https://byrumlab.github.io/proteoDA/reference/filter_proteins_by_annotation.md)
removes proteins using logical conditions on `DAList$annotation`. This
is useful for removing contaminants, such as keratins, decoys, spike-in
standards, etc.

Function recap

``` r
# general command
DAList <- filter_proteins_by_annotation(DAList,
            condition)
```

- `condition` is evaluated in the context of `DAList$annotation`
- Proteins are kept where condition is `TRUE`

Example: remove a single protein

``` r
DAList <- filter_proteins_by_annotation(DAList,
            uniprot_id != "P12345")
```

Example: remove keratins and contaminants

``` r
# Remove any protein which contains "keratin" in the name:
DAList <- DAList |>
      filter_proteins_by_annotation(!grepl("KERATIN", protein_name)) |>
      filter_proteins_by_annotation(contaminant_flag == FALSE)
```

Example: remove proteins with molecular weight \< 30 kDa

``` r
# Remove any protein with molecular weight < 30 kDa
DAList <- filter_proteins_by_annotation(DAList,
                                          molecular_weight > 30)
```

Behavior and tags

- Updates `DAList$annotation` and re-subsets `DAList$data` rows
- On success, `DAList$tags$filter_proteins_by_annotation` is updated
  with the condition call
- Errors if the condition evaluates to `NA` for any row

------------------------------------------------------------------------

## Protein-level filtering by missingness across groups

For DA analysis, you typically require some minimum level of non-missing
data per group. `proteoDA` provides two global strategies:

- By `counts` (number of replicates in number of groups) →
  [`filter_proteins_by_group()`](https://byrumlab.github.io/proteoDA/reference/filter_proteins_by_group.md)
- By `proportion` (fraction of samples per group) →
  [`filter_proteins_by_proportion()`](https://byrumlab.github.io/proteoDA/reference/filter_proteins_by_proportion.md)

Both depend on a grouping column in `DAList$metadata` (e.g. “group”).

### `filter_proteins_by_group()`: numeric replicate thresholds

**General Rule of Thumb:** Require at least 2/3 of biological replicates
to have values in at least one group.

Limma requires at least 2 replicates but power increases with increased
replicates due to better estimation of model variance.

For each protein:

- Count how many samples per group have non-missing values
- Require at least `min_reps` replicates within at least `min_groups`
  groups

Function recap

``` r
DAList_filtered <- filter_proteins_by_group(
                DAList,
                min_reps       = 2,
                min_groups     = 1,
                grouping_column = "group"
)
```

Arguments

- `min_reps` (required)  
  Minimum number of non-missing samples within a group for a protein to
  be considered present

- `min_groups` (required)  
  Minimum number of groups that must meet min_reps

- `grouping_column`  
  Column in `DAList$metadata` with group labels (default “group”)

Example: strict filtering (no missingness allowed)

- data contains 5 biological replicates and 4 conditions or “groups”

``` r
DAList_strict <- filter_proteins_by_group(
          DAList,
          min_reps       = 5,     # group size
          min_groups     = 4,     # all groups
          grouping_column = "group"
)
```

Example: at least one replicate in all groups

``` r
DAList_lax <- filter_proteins_by_group(
              DAList,
              min_reps       = 1,
              min_groups     = 4,
              grouping_column = "group"
)
```

Behavior and diagnostics

- Checks that `min_reps` is not larger than the smallest group size  
  (otherwise errors with an informative message)
- Checks that enough `groups` can meet `min_reps` to satisfy
  `min_groups`
- Builds an internal matrix of non-missing counts per protein × group
- Returns a filtered `DAList` and logs parameters in
  `DAList$tags$filter_proteins_by_group`

### filter_proteins_by_proportion(): proportional thresholds

When group sizes are imbalanced, you may prefer to specify a proportion
rather than an absolute number of replicates.

For each protein and group:

- Compute group-specific thresholds = ceiling(group_size \* min_prop)
- Require that the protein is present in at least threshold samples in
  each group

Function recap

``` r
DAList_prop <- filter_proteins_by_proportion(
            DAList,
            min_prop       = 0.66,
            grouping_column = "group"
)
```

Arguments

- `min_prop`  
  Minimum fraction of samples within each group that must be
  non-missing  
  (must be between 0 and 1, e.g. 0.5, 0.66, 1.0)

- `grouping_column`  
  Column in `DAList$metadata` with group labels

Example: strict (no missing)

``` r
DAList_all_present <- filter_proteins_by_proportion(
        DAList,
        min_prop       = 1,
        grouping_column = "group"
)
```

Example: moderate (≥50% of samples per group)

``` r
DAList_half <- filter_proteins_by_proportion(
      DAList,
      min_prop       = 0.5,
      grouping_column = "group"
)
```

Behavior and tags

- Computes group-specific thresholds via ceiling(n_group \* min_prop)
- Requires all groups to meet their thresholds for a protein to be kept
- Updates `DAList$data`, `DAList$annotation`, and logs parameters in
  `DAList$tags$filter_proteins_by_proportion`

------------------------------------------------------------------------

## Contrast-specific missingness filtering

Global filters may be too strict or too lenient for some contrasts,
especially in heterogeneous designs.

[`filter_proteins_per_contrast()`](https://byrumlab.github.io/proteoDA/reference/filter_proteins_per_contrast.md)
allows you to:

- Use the same raw `DAList` for all contrasts
- Apply `per-contrast` missingness rules based on the groups involved in
  each contrast
- Store results in `DAList$data_per_contrast` and
  `DAList$annotation_per_contrast` for downstream limma runs

### **filter_proteins_per_contrast()**

Concept

Given:

- A `contrasts_file` where each row looks like
  `GroupA_vs_GroupB=GroupA - GroupB`
- A `grouping_column` in metadata (e.g., “group”)

For each contrast:

- Identify the two groups in the contrast (e.g., “GroupA” and “GroupB”)
- Count non-missing values per protein in each group
- Apply `min_reps` and `require_both_groups` to decide which proteins to
  keep
- Store contrast-specific data and annotation

Function recap

``` r
DAList_filtered <- filter_proteins_per_contrast(
        DAList           = DAList,
        contrasts_file   = "data/contrasts.csv",
        min_reps         = 2,
        require_both_groups = TRUE,
        grouping_column  = "group"
)
```

Arguments

- `DAList`  
  Input `DAList` with raw or pre-filtered data, metadata, annotation

- `contrasts_file`  
  Path to CSV with lines like `GroupA_vs_GroupB=GroupA - GroupB` (no
  header)

- `min_reps`  
  Minimum replicates required per group in the contrast

- `require_both_groups`  
  **TRUE:** both groups must meet `min_reps`  
  **FALSE:** at least one group must meet `min_reps`

- `grouping_column`  
  Column in `metadata` with `group` labels

Outputs and bookkeeping

On success, the function:

- Adds `DAList$filtered_proteins_per_contrast` (named list of retained
  IDs per contrast)

- Adds `DAList$data_per_contrast` (named list of contrast-specific data
  matrices)

- Adds `DAList$annotation_per_contrast` (named list of contrast-specific
  annotation data frames)

- Adds `DAList$tags$retention_summary` (data frame with Contrast and
  Proteins_Retained)

These per-contrast slots are then used by downstream functions such as
[`run_filtered_limma_analysis()`](https://byrumlab.github.io/proteoDA/reference/run_filtered_limma_analysis.md)
or plotting/reporting functions that operate per contrast.

*NOTE:* If you use `filtered_proteins_per_contrast()` and then apply
[`normalize_data()`](https://byrumlab.github.io/proteoDA/reference/normalize_data.md),
the normalized values will be added to `DAList$data_per_contrast` and
`DAList$data` will remain untouched.

Example: heterogeneous design

``` r
DAList <- zero_to_missing(DAList)

DAList <- filter_proteins_per_contrast(
        DAList             = DAList,
        contrasts_file     = "design/contrasts.csv",
        min_reps           = 2,
        require_both_groups = FALSE,
        grouping_column    = "group"
)

DAList$tags$retention_summary
```

------------------------------------------------------------------------

## Example workflows

### Simple balanced design

``` r
DAList_pp <- DAList_raw |>
        zero_to_missing() |>
        filter_samples(group %in% c("Control", "Treatment")) |>
        filter_proteins_by_annotation(!grepl("KERATIN", protein_name)) |>
        filter_proteins_by_group(
            min_reps       = 2,
            min_groups     = 2,
            grouping_column = "group"
)
```

### Imbalanced groups with many samples

``` r
DAList_pp <- DAList_raw |>
    zero_to_missing() |>
    filter_samples(!qc_fail) |>
    filter_proteins_by_proportion(
        min_prop       = 0.5,
        grouping_column = "group"
)
```

### Mixed design with contrast-specific filters

``` r
DAList_pp <- DAList_raw |>
    zero_to_missing() |>
    filter_proteins_by_annotation(!decoy) |>
    filter_proteins_per_contrast(
        DAList             = _,
        contrasts_file     = "data/contrasts.csv",
        min_reps           = 2,
        require_both_groups = FALSE,
        grouping_column    = "group"
)
```

------------------------------------------------------------------------

## Summary

In `proteoDA`, preprocessing is a separate, explicit step designed to
make downstream limma modeling more robust and more transparent. The key
ideas from this vignette are:

- **Handle zeros early**  
  Use
  [`zero_to_missing()`](https://byrumlab.github.io/proteoDA/reference/missing_data.md)
  so that missingness is consistently encoded as `NA` before any
  filtering or modeling.

- **Filter at the sample level when needed**  
  Use
  [`filter_samples()`](https://byrumlab.github.io/proteoDA/reference/filter_samples.md)
  with metadata-based conditions to remove poor-quality, failed, or
  out-of-scope samples before touching protein-level filters.

- **Use annotation to remove unwanted proteins**  
  Use
  [`filter_proteins_by_annotation()`](https://byrumlab.github.io/proteoDA/reference/filter_proteins_by_annotation.md)
  to remove contaminants, keratins, decoys, spike-ins, or any proteins
  flagged by annotation.

- **Control missingness across groups**

  - [`filter_proteins_by_group()`](https://byrumlab.github.io/proteoDA/reference/filter_proteins_by_group.md)
    for replicate-count–based thresholds (`min_reps`, `min_groups`)  
  - [`filter_proteins_by_proportion()`](https://byrumlab.github.io/proteoDA/reference/filter_proteins_by_proportion.md)
    when group sizes are imbalanced and a proportion-based rule
    (`min_prop`) is easier to reason about.

- **Apply contrast-specific rules when designs are heterogeneous**  
  [`filter_proteins_per_contrast()`](https://byrumlab.github.io/proteoDA/reference/filter_proteins_per_contrast.md)
  allows you to tailor missingness thresholds to each contrast and
  stores per-contrast data and annotation in `DAList$data_per_contrast`
  and `DAList$annotation_per_contrast`.

- **Track decisions via tags**  
  Most preprocessing functions write into `DAList$tags`, making it
  easier to document which filters were applied and with what
  parameters.

Together, these tools let you build a reproducible preprocessing
pipeline that matches your experimental design, balances sensitivity and
robustness, and feeds clean, well-documented inputs into the `proteoDA`
modeling and reporting functions.

------------------------------------------------------------------------

## FAQ

### Should I always run `zero_to_missing()?`

In most DIA/DDA proteomics workflows, yes:

- Upstream tools often encode non-detected values as 0, NA, or NaN
- proteoDA has a function in the beginning of the workflow that finds
  blanks and NAs and replaces those with a numeric `0` so it is clear
  what input data values are included
- Downstream modeling in `proteoDA` assumes missing values are `NA`
- Running
  [`zero_to_missing()`](https://byrumlab.github.io/proteoDA/reference/missing_data.md)
  early ensures consistent missing value handling

The only time you might not need it is if your input data are already
NA-coded or you want to run limma using zero values specifically.

------------------------------------------------------------------------

### When should I filter samples vs proteins?

Use
[`filter_samples()`](https://byrumlab.github.io/proteoDA/reference/filter_samples.md)
for:

- Poor-quality or flagged samples
- Known batch failures
- Removing specific sample subgroups

Use *protein-level* filters (`filter_proteins_by_annotation`,
`filter_proteins_by_group`, `filter_proteins_by_proportion`,
`filter_proteins_per_contrast`) for:

- Contaminants, decoys, spike-ins
- Insufficient data density (missingness)

In practice, you typically:

1.  Filter samples first
2.  Then apply protein-level filters based on the cleaned sample set

------------------------------------------------------------------------

### How strict should my missingness filters be?

There is no single rule, but some common patterns:

- Exploratory analysis
  - filter_proteins_by_group(min_reps = 1, min_groups = n_groups)  
    or  
  - filter_proteins_by_proportion(min_prop = 0.5)  
- Confirmatory analysis / publication-ready
  - filter_proteins_by_group(min_reps = 2 or 3, min_groups \>= 2)
  - Possibly combine with imputation restricted to filtered proteins

We recommend trying a few thresholds and checking how many proteins
remain and how model diagnostics look (e.g., moving SD plots).

------------------------------------------------------------------------

### When should I use `filter_proteins_per_contrast()`?

Use contrast-specific filtering when: - Your design includes
heterogeneous groups or many unused samples for a given contrast - Some
contrasts have fewer replicates per group and need a slightly different
rule - Some groups have more missing values - You want to maximize
proteins per contrast without globally over-filtering

If your experiment is simple and balanced, global filters may be enough.
The QC missing value heatmap can help decide how much filtering is
required.

------------------------------------------------------------------------

### How do preprocessing choices interact with normalization and imputation?

A typical (recommended) order:

1.  Import → create DAList
2.  zero_to_missing()
3.  Sample-level and annotation-level filtering
4.  Protein-level missingness filtering (filter_proteins_by_group /
    filter_proteins_by_proportion / filter_proteins_per_contrast)
5.  Normalization (see normalization vignette)
6.  Imputation (if desired, after normalization)
7.  Fit limma models (see design/statistics vignettes)

See the normalization vignette for details on available normalization
strategies and diagnostics:
`vignette("proteoDA_v2.0_Norm_Impute", package = "proteoDA")`

The choice of limma design and contrasts is described in more detail in:
[`vignette("proteoDA_v2.0_design", package = "proteoDA")`](https://byrumlab.github.io/proteoDA/articles/proteoDA_v2.0_design.md)

Consistency is key: once you decide on a pipeline, keep the same
sequence for all datasets and document thresholds in the project
reports.

------------------------------------------------------------------------

### How can I check what filtering was applied?

Many preprocessing functions store parameters under `DAList$tags`. For
example:

``` r
str(DAList$tags, max.level = 1)

DAList$tags$filter_proteins_by_group

DAList$tags$filter_proteins_by_proportion

DAList$tags$retention_summary
```

You can also directly inspect the dimensions of the data and annotation:

``` r
dim(DAList$data)

dim(DAList$annotation)
```

For contrast-specific filtering:

``` r
names(DAList$data_per_contrast)

DAList$tags$retention_summary
```

These summaries help you reproduce and justify filtering decisions in
manuscripts and reports.
