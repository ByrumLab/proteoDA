
<!-- README.md is generated from README.Rmd. Please edit that file -->

# proteoDA

<!-- badges: start -->
<!-- badges: end -->

proteoDA is a public-friendly R package for the analysis of high resolution
mass spectrometry protein data. The package utilizes a custom S3 class
that keeps the R objects consistent across the pipeline and is easily
chain/pipe-able.

## Installation

The Github repository where the package is stored is public. You will need the
`devtools` package installed and the name of the repository. 

``` r
install.packages("devtools")
library(devtools)
devtools::install_github(ByrumLab/proteoDA)
                         
```

## DAList() S3 object structure

`DAList` holds all the information for a single quantitative proteomics
experiment and all our functions take it as input. It is a list with 7
slots:

1)  <b>data-</b> MS intensity data where each row is a protein and each column
    is a sample. It must be numeric. 
2)  <b>annotation-</b> Protein annotation information such as the accession id, 
    description, gene symbol, etc. "uniprot_id" is a required column. 
    It must have the same number of rows as data. 
    Each column is a separate piece of annotation info.
3)  <b>metadata-</b> A data frame of sample information such as sample_name,
    group, batch, gender, paired, etc.  
    The rownames must match the column names in data where 
    there is one row per sample.
4)  design- A list which holds information on the statistical design:
    the design formula and matrix, contrasts, etc.
5)  eBayes_fit- The model fit object from running `limma`’s models on
    the data.
6)  results- Statistical tables of DE results for each of the
    statistical terms/contrasts analyzed.
7)  tags- Basically a “miscellaneous” slot, used to keep track of info
    on filtering, normalization, etc.

For the public R package, users will have to make the first three slots
all at once, and then our functions are used to process the data and add
the rest of the slots. For our internal usage, we have some functions to
import data and metadata and assemble a DAList object.

All of the functions in the pipeline check for a proper structure of the
DAList object, both when it is input into the function and before it
returns the result. We still need to improve these checks and think of
all the edge cases to test for.

## Example pipeline

With the switch to a new S3 object, the analysis pipeline has changed a
little bit. Some of these changes are to accommodate the new data
structure, and others are to make things a little more straightforward
for the public R package. So, here’s a new example pipeline, with some
explanation.

### Data Import

For our internal use, we have a `read_DIA_data()` function that will
read in our data, based on the formatting of our Samples Report files
from MaxQuant:

``` r
data <- read_DIA_data(input_file = "path/to/Samples Report of from Maxquant.csv") 
```

This function also filters out contaminants and decoys during import,
before converting everything into a DAList.

After importing the raw data, we import metadata and add it to the
DAList object:

``` r
data <- add_metadata(DAList = data, metadata_file = "path/to/metdata.csv") 
```

One thing to note: all the functions in the pipeline take the DAList
object as their first argument/input. So, most functions can be piped
together, using either the base R pipe `|>` or the tidyverse/magrittr
`%>%`. So, the above code can be condensed to:

``` r
data <- read_DIA_data(input_file = "path/to/Samples Report of from Maxquant.csv")  |> 
  add_metadata(metadata_file = "path/to/metdata.csv") 
```

When we split our package into the public and internal versions, these
functions will be in the internal one.

### Filtering

#### Sample filtering

The filtering functions have been rewritten to hopefully be easier to
use and combine. First, there’s functionality to filter out samples
based on their metadata, with the `filter_samples()` function.
Specifying filtering is easier, and now works like core R functions such
as `subset()` or `dplyr::filter()`: you just supply a logical expression
that gives the filtering criteria. So, we can remove all samples where
the group is equal to “pool”

``` r
data <- filter_samples(data, group != "Pool")
```

Or where a sample ID contains a certain string:

``` r
data <- filter_samples(data, !stringr::str_detect(sample, "sampleX"))
```

Again, these can be chained together:

``` r
data <- filter_samples(data, group != "Pool") |> 
  filter_samples(!stringr::str_detect(sample, "sampleX"))
```

The filtering functions that are applied are tracked in the `tags` slot
of the DAList object. Any samples that are filtered out are
automatically removed from the data as well.

#### Protein/gene filtering

We also have some new functions for filtering out proteins/genes. First,
we can filter proteins by information in their annotation data:

``` r
data <- filter_proteins_by_annotation(data, Protein.Name != "bad protein")
```

We can also filter proteins based on their degree of missing data across
samples. One function, `filter_proteins_by_group()`, replicates the
filtering we used to apply in the old `process_data()` function (e.g.,
requiring a non-zero intensity in a minimum number of samples across a
minimum number of groups).

``` r
data <- filter_proteins_by_group(data, min_reps = 5, min_groups = 3, grouping_column = "group")
```

We can also filter by proportion, which is useful when sample sizes are
uneven across groups:

``` r
data <- filter_proteins_by_proportion(data, 
                                      min_prop = 1, # No missing data allowed
                                      grouping_column = "group")
```

As with sample filtering, the protein filter functions that are applied
are tracked in the `tags` slot. With both types of filtering,
proteins/genes are removed form both the annotation and data slots.

This more modular way of filtering hopefully makes it easier to add
other filtering methods. If there’s one you think we should add, make a
GitHub issue!

### Normalization report

Next, we make the normalization report:

``` r
write_norm_report(data,
                  grouping_column = "group",
                  output_dir = "directory/to/save/",
                  file = "filename.pdf", # See note below
                  overwrite = T)
```

In the S3 version of the package, this function has no defaults for
`output_dir` and `file`: you need to specify those yourself. In the
example above, I’ve explicitly chosen generic ones: these are not the
structure that Stephanie wants. When we split our package into public
and internal packages, we can set the defaults we want for our internal
package. For now, you’ll need to specify them by hand.

This function invisibly returns the input data/DAList, so it can be
chained if desired (though you may not want to).

### Normalize data

After looking at the normalization report and deciding what
normalization method to use, you normalize the data:

``` r
norm_data <- data |>
  normalize_data("quantile")
```

This function replaces the unnormalized data in the `data` slot of the
`DAList` with normalized data, and updates the `tags` to record that the
dat are normalized and the method of normalization.

### QC report

Then, we make the QC report. The QC report will actually work on either
normalized or non-normalized data, and will use the appropriate y-axis
names (intensity vs. normalized intensity). But we generally want to run
this on normalized data.

``` r
write_qc_report(norm_data,
                color_column = "group",
                output_dir = "directory/to/save",
                file = "QC_report.pdf",
                overwrite = T)
```

This function is like the normalization report: it has no defaults for
`output_dir` and `file`: you need to specify those yourself. In the
example above, I’ve explicitly chosen generic ones: these are not the
structure that Stephanie wants. When we split our package into public
and internal packages, we can set the defaults we want for our internal
package. For now, you’ll need to specify them by hand.

This function invisibly returns the input data/DAList, so it can be
chained if desired (though you may not want to).

### Specifying the statistical model

This is an area where there are some big changes, which hopefully will
make it easier to specify more complicated models. The first step is to
use `add_design()` to add the statistical model. Besides the DAList, the
only argument this function takes in is a formula of the desired
statistical model. In the previous code of the package, this was done by
specifying column names, and there were some limitations around model
specification (intercept vs. no-intercept, inclusion of interactions,
etc). Now, there’s much more flexibility. The “default” model we often
use (a no intercept model with “group” as the only statistical factor)
would be:

``` r
norm_data <- add_design(norm_data,
                        design_formula = ~ 0 + group)
```

You could also specify a model with intercepts:

``` r
norm_data <- add_design(norm_data,
                        design_formula = ~ group)
```

More complicated things like interactions:

``` r
norm_data <- add_design(norm_data,
                        design_formula = ~ treatment*genotype)
```

Or even go crazy including random factors (though you can only include
one, and only for intercepts):

``` r
norm_data <- add_design(norm_data,
                        design_formula = ~ treatment*genotype + (1 | batch))
```

The model equation syntax is the same as functions like
`lm()`/`glm()`/etc., with a limited random effect syntax in the style of
`lme4`/`brms`. The function does checking to make sure you’re inputting
a valid equation, and if you don’t the error messages are hopefully
clear about what to fix.

Then, you can optionally add contrasts, with `add_contrasts()`. These
can be specified from a file, or just a string of contrasts:

``` r
# Specify file
norm_data <- add_contrasts(norm_data,
                           contrasts_file = "path/to/contrasts.txt")

# Specify by hand
norm_data <- norm_data |>
  add_contrasts(contrasts_vector = c("example=groupA-groupB"))
```

Specifying contrasts is optional. If you don’t add other contrasts, then
by default you’ll get estimates of the model terms directly from your
design matrix. In cases where that isn’t what you want (e.g., you want
all pairwise comparisons among the different levels in a group), you
need to add contrasts. Adding contrasts has to be done after adding a
design with `add_design()`, and the function will perform checks to make
sure your contrasts and statistical design make sense together.

### Fitting the model

After you add the statistical design, you fit the limma model with
`fit_limma_model()` and get tables of differential expression results
with `extract_DA_results()`. `fit_limma_model()` adds the `eBayes_fit`
slot to the DAList, and `extract_DA_results()` adds the `results` slot
as a list with one element for each statistical term/contrast.

``` r
final <- fit_limma_model(norm_data) |> 
  extract_DA_results(pval_thresh = 0.055, lfc_thresh = 1, adj_method = "BH")
```

Besides the report functions above, the whole analysis could in theory
be chained together into one pipeline:

``` r
full_chain <- read_DIA_data(input_file = "path/to/Samples Report of from Maxquant.csv")  |> 
  add_metadata(metadata_file = "path/to/metdata.csv") |> 
  filter_samples(group != "Pool") |> 
  filter_samples(!stringr::str_detect(sample, "sampleX")) |> 
  filter_proteins_by_group(min_reps = 5, 
                           min_groups = 3, 
                           grouping_column = "group") |> 
  normalize_data("quantile") |> 
  add_design(~ treatment*genotype + (1 | batch)) |> 
  add_contrasts(contrasts_vector = c("example=groupA-groupB")) |> 
  fit_limma_model() |> 
  extract_DA_results()
```

### Writing final reports

Then, you write out results of the analysis in tables
(`write_limma_tables()`) and interactive reports/plots
(`write_limma_plots`). We might want to change those names, substituting
“results” for “limma”.

``` r
write_limma_tables(full_chain,
                   output_dir = "path/to/output/dir",
                   overwrite = T)
```

``` r
write_limma_plots(full_chain,
                  grouping_column = "group")
```

Confusingly, at the moment these functions do have some defaults for
output directories and filenames (though I can’t guarantee that they’re
exactly what Stephanie wants. Again, can always specify by hand). But,
this will change once we split into the public and internal R package:
the defaults for the public package will not be what we want to use
internally.

## Experiments other than DIA protein data

At the moment, everything is written specifically with DIA protein data
in mind. But, most of the middle functions should work regardless. The
data import functions probably won’t work with other data. But, if you
can get your data into the right structure to be put into our S3 object
(see above), the rest of the functions should work (though some text
labels might be incorrect: e.g., “protein” instead of “phospho”).
