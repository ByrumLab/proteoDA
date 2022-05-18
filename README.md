
<!-- README.md is generated from README.Rmd. Please edit that file -->

# proteomicsDIA

<!-- badges: start -->
<!-- badges: end -->

The `proteomicsDIA` package will (eventually) be our internal R package
that supplies all the functions we need for doing proteomics analysis of
DIA data.

## Non-DIA experiments

At the moment, these functions should work fine for non-DIA experiments
as well. Most of the functions are the same for the different experiment
types. Some functions, like `extract_data()` and `make_targets()`, have
slightly different behavior depending on the pipeline. At the moment,
only the DIA code within these functions has been updated: you might
notice differences in the text printed to the console and in the output
objects. But, the non-DIA code is more or less unchanged from version 26
of Charity’s functions, so it should still all work properly. If you run
into any problems, [file an
issue](https://github.com/ByrumLab/proteomicsDIA/issues).

## Installation

The Github repository where the package is stored is private, so
installation is a little more complicated that usual. You’ll need the
`devtools` package installed (which you’ll want to have installed anyway
for developing the package, see below). You’ll also need the personal
access token (PAT) that gives you access to the repository. Email Tim or
Stephanie to get it. Once you have it, you can install the development
version of `proteomicsDIA` from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("ByrumLab/proteomicsDIA",
                         auth_token = "COPY_PAT_HERE")
```

## Example

The package is getting developed in pipeline order, adding each of the
functions and subfunctions needed for the various steps in the example
workflows (see the progress tracker below). The goal is to get all the
current functionality into the R package before we start adding new
features or making changes.

``` r
# Load the proteomicsDIA package
library(proteomicsDIA)

# Extract Maxquant data
extracted_data <- extract_data("path/to/Samples Report of from Maxquant.csv",
                               pipe = "DIA",
                               enrich = "protein")

# make targets
targets <- make_targets(file = "path/to/metdata.csv",
                        sampleIDs = colnames(extracted_data$data),
                        pipe = "DIA",
                        enrich = "protein")

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
make_proteinorm_report(normList = norm$normList,
                       groups = norm$targets$group, 
                       enrich = "protein")

# Make the QC report
make_qc_report(normList = norm$normList, 
               norm.method = "vsn",
               groups = norm$targets$group,
               batch = norm$targets$group,
               enrich = "protein")

# Make the design matrix
design <- make_design(targets=norm$targets,
                      group_column = "group",
                      factor_columns = NULL,
                      paired_column = NULL)
# Make the contrasts matrix
contrasts <- make_contrasts(file = "path/to/contrasts/file.csv",
                            design = design$design)

# Functions are documented, check them out:
?extract_data
?make_targets
?subset_targets
?process_data
?make_proteinorm_report
?make_qc_report
?make_design
?make_contrasts

# If you see any typos or things that aren't clear, let me know!
```

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

### Collaboration with Git

We need to figure this out! Deciding on the best path forward with this
depends a little on everyone’s familiarity and level of comfort with
Git. I (Tim) have used Git for the last couple years for all my
projects, though I have rarely used it to collaborate with others and
thus have rarely needed to worry about the optimal collaboration
workflow. We can work together to figure that out more, but in the
meantime here are some resources that might be helpful:

-   [Happy Git with R](https://happygitwithr.com/index.html)- A guide to
    using Git and Github with RStudio.
-   [Chapter 18 of R packages](https://r-pkgs.org/git.html) is on Git.
-   Git GUIs- Personally, I mostly avoid using Git in the command line
    and instead use GUIs. RStudio has a minimal one, explained in the
    first link above. Initially, I mostly used
    [GitKraken](https://www.gitkraken.com/git-client), which is free
    with a student license and has some nice touches (better diff
    viewing, better branch visualization). Lately, I’ve just been using
    the Git functionality in [Visual Studio
    Code](https://code.visualstudio.com/), which also has some nice
    extensions for using Git.

Once we’re all set up and comfortable with Git, we’ll need to figure out
how to manage workflows. That might be a topic for another day.

### Misc ideas

See the GitHub issue tracker for other ideas

-   Add a changelog.md file, to textually track/explain changes (with
    version #s).

-   Check out code coverage and CI tools, to automate testing. Some
    options in usethis. Might want to wait to do much testing until
    after we decide on whether we want to do much re-writing and
    reformatting of the code.

### Progress tracker

I’m working my way through the “top-level” functions in the example DIA
workflow sent by Charity. I’m not doing much on the TMT, LF, or
phosphoTMT side of things, though I have done a little code restyling
for those sometimes.

On the code side, I’ve been doing relatively little. I’ve been doing
some restyling of the code for clarity, and reduced redundancy. I’ve
dropped TODOs through the code to mark ideas for re-doing of code. The
main thing I’ve done is to rework print statements and error messages
(using the `cli` R package) to make these better formatted and more
informative.

On the documentation side, I’m just trying to get some initial
documentation for each main function and all its subfunctions, linking
them together. The goal is to explain all the inputs and outputs to get
a better idea of how they fit together.

On the test side, I started writing some unit tests for the data
extraction functions, but I think I might hold off on that for now. If
we do any changing of the code structure, we may have to re-write the
tests as well. Maybe better to just do it once.

| Function                   | Code               | Document           | Test               |
|----------------------------|--------------------|--------------------|--------------------|
| `extract_data`             | :heavy_check_mark: | :heavy_check_mark: | :heavy_check_mark: |
| `make_targets`             | :heavy_check_mark: | :heavy_check_mark: | :x:                |
| `subset_targets`           | :heavy_check_mark: | :heavy_check_mark: | :x:                |
| `process_data`             | :heavy_check_mark: | :heavy_check_mark: | :x:                |
| `make_proteinorm_report`   | :heavy_check_mark: | :heavy_check_mark: | :x:                |
| `make_qc_report`           | :heavy_check_mark: | :heavy_check_mark: | :x:                |
| `make_design`              | :heavy_check_mark: | :heavy_check_mark: | :x:                |
| `make_contrasts`           | :heavy_check_mark: | :heavy_check_mark: | :x:                |
| `fit_limmma_model`         | :heavy_check_mark: | :x:                | :x:                |
| `extract_limma_DE_results` | :heavy_check_mark: | :x:                | :x:                |
| `write_limma_results`      | :x:                | :x:                | :x:                |
| `make_limma_report`        | :x:                | :x:                | :x:                |

The final steps in the pipeline, for the limma analysis, are being
reworked a little. Originally, there were two functions:
`run_limma_analysis` fit the model, extracted data, saved some .csv
files, and made the glimma plots. `add_limma_results` added limma
results to an excel spreadsheet. I’m working to separate out some of the
functionality into independent functions. Now, we have:

-   `fit_limma_model`- Which just does model fitting, and returns a list
    of various model fit objects.
-   `extract_limma_DE_results`- Which extracts statistical results for
    each contrast from the list of limma model fits.
-   TBD- `write_limma_results`- Which will output the various .csv and
    excel files of results.
-   TBD- `make_limma_report`- Which will output the GLIMMA plots and
    html report for end users.
