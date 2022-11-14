
<!-- README.md is generated from README.Rmd. Please edit that file -->

# proteomicsDIA

<!-- badges: start -->
<!-- badges: end -->

Our package is getting closer to being ready for publication. One major
task we’ve been working on for publication is streamlining/standardizing
the input and output of all functions to use a single list with a custom
S3 class. This has required a pretty comprehensive re-write of all our
functions, but hopefully does so in a way that makes the package easier
to use.

We’re at the point where it would be helpful for folks to test out the
new S3-based pipeline. It would be especially helpful if folks could
re-run some old projects using the new S3 pipeline and make sure the
results match the old package version.

We also need a new name!

## Contributing to package development

I would highly recommend checking out the [R
Packages](https://r-pkgs.org/) book from Hadley Wickham and Jenny Bryan.
It is a pretty comprehensive book on creating R packages, and provides a
nice framework that we can all follow (and the framework which the
package is currently following).

You’ll want to install a few packages for development:

``` r
install.packages(c("devtools", "roxygen2", "testthat", "knitr", "rmarkdown", "usethis"))
```

You may want to set up your `.Rprofile` to automatically load some of
all of the packages when you start R (see
[here](https://r-pkgs.org/setup.html#personal-startup-configuration) for
more info), though I don’t usually bother and just load them as I need
them.

## Installation

The Github repository where the package is stored is private, so
installation is a little more complicated that usual. You’ll need the
`devtools` package installed (which you’ll want to have installed anyway
for developing the package, see below). You’ll also need the personal
access token (PAT) that gives you access to the repository. Email Tim or
Stephanie to get it. Once you have it, you can install the development
version of `proteomicsDIA` from [GitHub](https://github.com/) with the
`devtools::install_github()` function. To install the old style of the
package:

``` r
# install.packages("devtools")
devtools::install_github("ByrumLab/proteomicsDIA",
                         auth_token = "COPY_PAT_HERE")
```

And, to install the new version based on the S3 object:

``` r
# install.packages("devtools")
devtools::install_github("ByrumLab/proteomicsDIA@s3_object",
                         auth_token = "COPY_PAT_HERE")
```

You can only have one version installed at a time. See the end of the
README for an example workflow for the old pipeline.

## New s3 object structure

The new version of the pipeline is based around a new S3 class,
`DIAlist` (probably need to change the name to something more general).
It holds all the info for a single quantitative proteomics experiment,
and all our functions take it as input. Its a list, with 7 slots:

1)  data- MS intensity data, where each row is a protein and each column
    is a sample.
2)  annotation- Protein/gene annotation information. Same number of rows
    as data, Each column is a separate piece of annotation info.
3)  metadata- A data frame of sample metadata, 1 row per sample (should
    match the number of columns in data).
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
import data and metadata and assemble a DIAlist object.

All of the functions in the pipeline check for a proper structure of the
DIAlist object, both when it is input into the function and before it
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

One thing to notes: this function just reads in the data and splits it
into the intensity data and annotation data. It doesn’t do any filtering
out of decoys or contaminants, as the old data import function did.

After importing the raw data, we import metadata and add it to the
DIAlist object:

``` r
data <- add_metadata(DIAlist = data, metadata_file = "path/to/metdata.csv") 
```

One thing to note: all the functions in the pipeline take the DIAlist
object as their first argument/input. So, most functions can be piped
together, using either the base R pipe `|>` or the tidyverse/magrittr
`%>%`. So, the above code can be condensed to:

``` r
data <- read_DIA_data(input_file = "path/to/Samples Report of from Maxquant.csv")  %>% 
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
data <- filter_samples(data, group != "Pool") %>% 
  filter_samples(!stringr::str_detect(sample, "sampleX"))
```

The filtering functions that are applied are tracked in the `tags` slot
of the DIAlist object. Any samples that are filtered out are
automatically removed from the data as well.

#### Protein/gene filtering

We also have some new functions for filtering out proteins/genes. First,
we can filter proteins by information in their annotation data:

``` r
data <- filter_proteins_by_annotation(data, Protein.Name != "bad protein")
```

As with `filter_samples()`, this function takes an expression that
evaluates to a logical condition that determines filtering. For our
internal use, we have a function that filters out contaminant proteins
(this used to happen in the data import step):

``` r
data <- filter_proteins_contaminants(data) 
# this is just a wrapper around filter_proteins_by_annotation()
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

Again, these can be chained:

``` r
data <- filter_proteins_contaminants(data) %>% 
  filter_proteins_by_group(min_reps = 5, min_groups = 3, grouping_column = "group")
```

As with sample filtering, the protein filter functions that are applied
are tracked in the `tags` slot. With both types of filtering,
proteins/genes are removed form both the annotation and data slots.

This more modular way of filtering hopefully makes it easier to add
other filtering methods. If there’s one you think we should add, make a
GitHub issue!

### Normalization report

Next, we make the proteinorm normalization report:

``` r
write_proteinorm_report(data,
                        grouping_column = "group",
                        out_dir = "directory/to/save/",
                        file = "filename.pdf", # See note below
                        overwrite = T)
```

In the S3 version of the package, this function has no defaults for
`out_dir` and `file`: you need to specify those yourself. In the example
above, I’ve explicitly chosen generic ones: these are not the structure
that Stephanie wants. When we split our package into public and internal
packages, we can set the defaults we want for our internal package. For
now, you’ll need to specify them by hand.

Also note: this output from this function can’t currently be chained (it
doesn’t output the DIAlist object). I’m still thinking about whether to
change that.

### Normalize data

After looking at the normalization report and deciding what
normalization method to use, you normalize the data:

``` r
norm_data <- data %>%
  normalize_data("quantile")
```

This function replaces the unnormalized data in the `data` slot of the
`DIAlist` with normalized data, and updates the `tags` to record that
the dat are normalized and the method of normalization.

### QC report

Then, we make the QC report. The QC report will actually work on either
normalized or non-normalized data, and will use the appropriate y-axis
names (intensity vs. normalized intensity). But we generally want to run
this on normalized data.

``` r
write_qc_report(full_higgs_chain,
                grouping_column = "group",
                out_dir = "directory/to/save",
                file = "QC_report.pdf",
                overwrite = T)
```

This function is like the normalization report: it has no defaults for
`out_dir` and `file`: you need to specify those yourself. In the example
above, I’ve explicitly chosen generic ones: these are not the structure
that Stephanie wants. When we split our package into public and internal
packages, we can set the defaults we want for our internal package. For
now, you’ll need to specify them by hand.

Also note: this output from this function can’t currently be chained (it
doesn’t output the DIAlist object). I’m still thinking about whether to
change that.

### Specifying the statistical model

This is an area where there are some big changes, which hopefully will
make it easier to specify more complicated models. The first step is to
use `add_design()` to add the statistical model. Besides the DIAlist,
the only argument this function takes in is a formula of the desired
statistical model. In the previous code of the package, this was done by
specifying column names, and there were some limitations around model
specification (intercept vs. no-intercept, inclusion of interactions,
etc). Now, there’s much more flexibility. The “default” model we often
use (a no intercept model with “group” as the only statistical factor)
would be:

``` r
norm_data <- add_design(data,
                        design_formula = ~ 0 + group)
```

You could also specify a model with intercepts:

``` r
norm_data <- add_design(data,
                        design_formula = ~ group)
```

More complicated things like interactions:

``` r
norm_data <- add_design(data,
                        design_formula = ~ treatment*genotype)
```

Or even go crazy including random factors (though you can only include
one, and only for intercepts):

``` r
norm_data <- add_design(data,
                        design_formula = ~ treatment*genotype + (1 | batch))
```

The model equation syntax is the same as functions like
`lm()`/`glm()`/etc., with a limited random effect syntax of in the style
of `lme4`/`brms`. The function does checking to make sure you’re
inputting a valid equation, and if you don’t the error messages are
hopefully clear about what to fix.

Then, you can optionally add contrasts, with `add_contrasts()`. These
can be specified from a file, or just a string of contrasts:

``` r
# Specify file
norm_data <- add_contrasts(norm_data,
                           contrasts_file = "path/to/contrasts.txt")

# Specify by hand
norm_data <- norm_data %>%
  add_contrasts(contrasts_vector = c("example=groupA-groupB"))
```

Specifying contrasts is optional. If you don’t add other contrasts, then
by default you’ll get estimates of the model terms directly from your
design matrix In cases where that isn’t what you want (e.g., you want
all pairwise comparisons among the different levels in a group), you
need to add contrasts. Adding contrasts has to be done after adding a
design with `add_design()`, and the function will perform checks to make
sure your contrasts and statistical design make sense together.

### Fitting the model

After you add the statistical design, you fit the limma model with
`fit_limma_model()` and get tables of differential expression results
`extract_DE_results()`. `fit_limma_model()` adds the `eBayes_fit` slot
to the DIAlist, and `extract_DE_results()` adds the `results` slot as a
list with one element for each statistical term/contrast.

``` r
final <- fit_limma_model(norm_data) %>% 
  extract_DE_results(pval_thresh = 0.5, lfc.thresh = 1, adj.method = "BH")
```

Besides the report functions above, the whole analysis could in theory
be chained together into one pipeline:

``` r
full_chain <- read_DIA_data(input_file = "path/to/Samples Report of from Maxquant.csv")  %>% 
  add_metadata(metadata_file = "path/to/metdata.csv") %>% 
  filter_samples(group != "Pool") %>% 
  filter_samples(!stringr::str_detect(sample, "sampleX")) %>% 
  filter_proteins_contaminants() %>% 
  filter_proteins_by_group(min_reps = 5, min_groups = 3, grouping_column = "group") %>% 
  normalize_data("quantile") %>% 
  add_design(design_formula = ~ treatment*genotype + (1 | batch)) %>% 
  add_contrasts(contrasts_vector = c("example=groupA-groupB")) %>% 
  fit_limma_model() %>% 
  extract_DE_results(pval_thresh = 0.5, lfc.thresh = 1, adj.method = "BH")
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
write_limma_plots(results_higgs,
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

### Example workflow- Old style

Here’s an example workflow for using the “old” style of the package,
without the S3 object. I’ve left it here for reference, but once we make
the switch to the new style it won’t be relevant anymore.

``` r
# Load the proteomicsDIA package
library(proteomicsDIA)

# Extract Maxquant data
extracted_data <- read_DIA_data("path/to/Samples Report of from Maxquant.csv")

# make targets
targets <- make_targets(input_file = "path/to/metdata.csv",
                        sample_IDs = colnames(extracted_data$data))

# If you want, can save the generated targets dataframe for internal checking
write.csv(targets, file = "targets.csv")


# Subset targets
sub <- subset_targets(targets = targets, 
                      filter_list = list(group = "pool",
                                         sample = c("sampleA", "sampleB")))

# Process data
norm <- process_data(data = extracted_data$data,
                     targets = sub$targets,
                     min.reps = 5,
                     min.grps = 3)
# Make the proteinorm report
write_proteinorm_report(processed_data = norm,
                        grouping_column = "group")

# Make the QC report
write_qc_report(processed_data = norm, 
                chosen_norm_method = "vsn",
                grouping_column = "group")

# Make the design matrix
design <- make_design(targets=norm$targets,
                      group_column = "group",
                      factor_columns = NULL,
                      paired_column = NULL)
# Make the contrasts matrix
contrasts <- make_contrasts(file = "path/to/contrasts/file.csv",
                            design = design$design)

# Fit the model
fit <- fit_limma_model(data = norm$normList[["vsn"]], # choose your normalization method
                       design_obj = design,
                       contrasts_obj = contrasts)

# Extract the differential expression results
results <- extract_limma_DE_results(limma_fit = fit)

# Save the tables of results
write_limma_tables(model_results = results,
                    norm.method = "vsn",
                    annotation = extracted_data$annot,
                    ilab = "example_1234")
# And save the plots and interactive report
write_limma_plots(model_results = results,
                   annotation =  extracted_data$annot,
                   groups = norm$targets$group,
                   output_dir = "example_1234")

# If you see any typos or things that aren't clear, let me know!
```
