# Create a DAList

Creates a DAList from existing raw data. DALists are S3 objects of class
DAList and contain 7 slots which hold data and results for a
quantitative proteomics project (see Details for a full description of
each slot). The first three slots (data, annotation, and metadata) are
required to create a DAList. The other slots can be supplied, but are
generally added via their respective analysis functions. *NB*- you must
ensure that the data, annotation, and metadata are in the proper order.
To help with this, `DAList` requires the row names of the metadata to
match the column names of the data. However, DAList cannot know if the
rows in the data and the annotation are in the proper order:
double-check that they are!

## Usage

``` r
DAList(
  data,
  annotation,
  metadata,
  design = NULL,
  eBayes_fit = NULL,
  results = NULL,
  tags = NULL
)
```

## Arguments

- data:

  A data frame or matrix containing protein intensity data for each
  sample. Rows are proteins, columns are samples.

- annotation:

  A data frame containing protein annotation data. Must contain a column
  named "uniprot_id" containing unique entries for each row.

- metadata:

  A data frame containing metadata on each sample.

- design:

  Optional: a list of design information, see Details.

- eBayes_fit:

  Optional: a model fit object, see Details.

- results:

  Optional: a list of data frame(s) summarizing differential abundance
  statistics, see Details.

- tags:

  Optional: a list containing flags, parameter values, and other
  settings, see Details.

## Value

A DAList object

## Details

A DAList contains 7 slots:

1.  data- A matrix or data frame containing protein intensity data. Rows
    are proteins, columns are samples. This array should contain numeric
    information only, and the number of rows (proteins) should match the
    number of rows in the annotation data frame. The columns (samples)
    in the data should be in the same order as the rows of samples in
    the metadata data frame, with column names that match the row names
    of the metadata.

2.  annotation- A data frame of protein information, containing the same
    number of rows as the data slot and in the same order. This data
    frame must contain a column named "uniprot_id", which contains a
    unique value for each row/protein (ideally, the UniProtID). The
    values in this column will be used to uniquely identify proteins in
    the DAList object. Additional columns containing other information
    (gene names, descriptions, molecular weights, etc) can be supplied.

3.  metadata- A data frame containing information on sample metadata,
    one row per sample. It should contain the same number of rows
    (samples) as there are columns in the data data frame, in the same
    order, with row names that match the column names of the data.

4.  design- A list containing information on the statistical design used
    for analyzing differential abundance. Set using
    [`add_design`](https://byrumlab.github.io/proteoDA/reference/add_design.md)
    and
    [`add_design`](https://byrumlab.github.io/proteoDA/reference/add_design.md).

5.  eBayes_fit- A list containing the model fit object. This is added
    using
    [`fit_limma_model`](https://byrumlab.github.io/proteoDA/reference/fit_limma_model.md),
    which uses functions from the limma package to perform differential
    abundance analysis. Specifically, this slot should be an MAarrayLM
    object output from
    [`limma::eBayes`](https://rdrr.io/pkg/limma/man/ebayes.html), which
    is used within
    [`fit_limma_model`](https://byrumlab.github.io/proteoDA/reference/fit_limma_model.md).

6.  results- A list of data frames, in which each data frame contains
    differential abundance results for a given statistical term or
    contrast. These tables are derived from the output of
    [`limma::decideTests`](https://rdrr.io/pkg/limma/man/decideTests.html)
    and
    [`limma::topTable`](https://rdrr.io/pkg/limma/man/toptable.html).

7.  tags- A list containing flags, parameter values, and other settings
    for the analysis. Generally, should not be set or modified by hand,
    but will be set by various pipeline functions.

## Examples

``` r
if (FALSE) { # \dontrun{

# Prepare data
# In most cases, these would be imported from external
# files and processed.
# Here, an example with 4 samples and 3 proteins.

# Raw numeric data
raw_data <- data.frame(sampleA = c(10000, 15000, 12500),
                       sampleB = c(20000, 30000, 25000),
                       sampleC = c(18000, 23000, 20500),
                       sampleD = c(36000, 46000, 41000))

# Annotation data
# Must be in same order as numeric values in data
# and requires one column named uniprot_id containing unique values
protein_info <- data.frame(uniprot_id = c("A0A023ZSD2", "M9NH73", "H9AZR4"),
                           gene = c("Mcf", "shh", "WntA"))

# sample metadata
# should be in same order as samples in the data,
# with row names that match the column names of the data
sample_info <- data.frame(sample_id = c("sampleA", "sampleB",
                                        "sampleC", "sampleD"),
                          group = c("control",   "control",
                                    "treatment", "treatment"))
rownames(sample_info) <- colnames(raw_data)


# Assemble the DAList from existing raw data.
raw <- DAList(data = raw_data,
              annotation = protein_info,
              metadata = sample_info)
} # }
```
